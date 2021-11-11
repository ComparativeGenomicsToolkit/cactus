#!/usr/bin/env python3

#Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Script strings together all the components to make the basic pipeline for reconstruction.
"""

import os
import sys
import xml.etree.ElementTree as ET
import math
from argparse import ArgumentParser
from toil.lib.bioio import system
from toil.lib.bioio import getTempFile
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.lib.bioio import getLogLevelString
from sonLib.bioio import catFiles

from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import runCactusConsolidated
from cactus.shared.common import findRequiredNode
from cactus.shared.common import RoundedJob

from cactus.paf.local_alignment import make_paf_alignments

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper

############################################################
############################################################
############################################################
##Shared functions
############################################################
############################################################
############################################################

def extractNode(node):
    """Make an XML node free of its parent subtree
    """
    return ET.fromstring(ET.tostring(node))

def getJobNode(phaseNode, jobClass):
    """Gets a job node for a given job.
    """
    className = jobClass.__name__
    assert className != ''
    assert className.isalnum()
    return phaseNode.find(className)

class CactusJob(RoundedJob):
    """Base job for all cactus workflow jobs.
    """
    def __init__(self, phaseNode, constantsNode, overlarge=False,
                 checkpoint=False, preemptable=True):
        self.phaseNode = phaseNode
        self.constantsNode = constantsNode
        self.overlarge = overlarge
        self.jobNode = getJobNode(self.phaseNode, self.__class__)
        if self.jobNode is not None:
            logger.info("JobNode = %s" % self.jobNode.attrib)

        memory = None
        cores = None
        if hasattr(self, 'memoryPoly'):
            # Memory should be determined by a polynomial fit on the
            # input size
            memory = self.evaluateResourcePoly(self.memoryPoly)
            if hasattr(self, 'memoryCap'):
                memory = int(min(memory, self.memoryCap))

        disk = None
        if memory is None and overlarge:
            memory = self.getOptionalJobAttrib("overlargeMemory", typeFn=int,
                                               default=getOptionalAttrib(self.constantsNode, "defaultOverlargeMemory", int, default=sys.maxsize))
            cores = self.getOptionalJobAttrib("overlargeCpu", typeFn=int,
                                              default=getOptionalAttrib(self.constantsNode, "defaultOverlargeCpu", int, default=None))
        elif memory is None:
            memory = self.getOptionalJobAttrib("memory", typeFn=int,
                                               default=getOptionalAttrib(self.constantsNode, "defaultMemory", int, default=sys.maxsize))
            cores = self.getOptionalJobAttrib("cpu", typeFn=int,
                                              default=getOptionalAttrib(self.constantsNode, "defaultCpu", int, default=sys.maxsize))
        RoundedJob.__init__(self, memory=memory, cores=cores, disk=disk,
                            checkpoint=checkpoint, preemptable=preemptable)

    def evaluateResourcePoly(self, poly):
        """Evaluate a polynomial based on the total sequence size."""
        features = {'totalSequenceSize': self.cactusWorkflowArguments.totalSequenceSize}
        if hasattr(self, 'featuresFn'):
            features.update(self.featuresFn())
        if hasattr(self, 'feature'):
            x = features[self.feature]
        else:
            x = features['totalSequenceSize']
        resource = 0
        for degree, coefficient in enumerate(reversed(poly)):
            resource += coefficient * (x**degree)
        return int(resource)

    def getOptionalJobAttrib(self, attribName, typeFn=None, default=None, errorIfNotPresent=False):
        """Gets an optional attribute of the job node.
        """
        return getOptionalAttrib(node=self.jobNode, attribName=attribName, typeFn=typeFn, default=default, errorIfNotPresent=errorIfNotPresent)

    def addService(self, job):
        """Works around toil issue #1695, returning something we can index for multiple return values."""
        rv = super(CactusJob, self).addService(job)
        if hasattr(rv, '__getitem__'):
            # New Toil. Promise we got is indexable for more sub-Promises.
            # Return the indexable root return value promise, from which
            # promises for different indexes can be obtained.
            return rv
        else:
            # Running on old Toil. Return the whole service host job so we can
            # get more Promises out of it. TODO: Remove this when everyone has
            # upgraded Toil.
            return self._services[-1]

class CactusPhasesJob(CactusJob):
    """Base job for each workflow phase job.
    """
    def __init__(self, cactusWorkflowArguments=None, phaseName=None, topFlowerName=0,
                 checkpoint=False, preemptable=True, halID=None,
                 fastaID=None):
        self.phaseName = phaseName
        phaseNode = findRequiredNode(cactusWorkflowArguments.configNode, phaseName)
        constantsNode = findRequiredNode(cactusWorkflowArguments.configNode, "constants")
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.topFlowerName = topFlowerName
        self.halID = halID
        self.fastaID = fastaID
        CactusJob.__init__(self, phaseNode=phaseNode, constantsNode=constantsNode, overlarge=False,
                           checkpoint=checkpoint, preemptable=preemptable)

    def makeFollowOnPhaseJob(self, job, phaseName):
        return self.addFollowOn(job(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName=phaseName,
                                    topFlowerName=self.topFlowerName, halID=self.halID, fastaID=self.fastaID)).rv()

    def getPhaseNumber(self):
        return len(self.cactusWorkflowArguments.configNode.findall(self.phaseNode.tag))

############################################################
############################################################
############################################################
##The blast phase that uses the trimming strategy.
############################################################
############################################################
############################################################

def inverseJukesCantor(d):
    """Takes a substitution distance and calculates the number of expected changes per site (inverse jukes cantor)
    d = -3/4 * log(1 - 4/3 * p)
    exp(-4/3 * d) = 1 - 4/3 * p
    4/3 * p = 1 - exp(-4/3 * d)
    p = 3/4 * (1 - exp(-4/3 * d))
    """
    assert d >= 0.0
    return 0.75 * (1 - math.exp(-d * 4.0/3.0))

def setupFilteringByIdentity(cactusWorkflowArguments):
    #Filter by identity
    blastNode = findRequiredNode(cactusWorkflowArguments.configNode, "blast")
    if getOptionalAttrib(blastNode, "filterByIdentity", bool, False): #Do the identity filtering
        adjustedPath = max(float(blastNode.attrib["identityRatio"]) * cactusWorkflowArguments.longestPath,
        float(blastNode.attrib["minimumDistance"]))
        identity = str(100 - math.ceil(100 * inverseJukesCantor(adjustedPath)))
        blastNode.attrib["lastzArguments"] = blastNode.attrib["lastzArguments"] + (" --identity=%s" % identity)

############################################################
############################################################
############################################################
##The consolidate phase, which runs the setup, caf,
## bar, reference and cactus to hal algorithms in one job
## on a multi-node machine
############################################################
############################################################
############################################################

class CactusPafAlign(CactusPhasesJob):
    def __init__(self, standAlone = False, *args, **kwargs):
        self.standAlone = standAlone
        super(CactusPafAlign, self).__init__(*args, **kwargs)

    def run(self, fileStore):
        exp = self.cactusWorkflowArguments.experimentWrapper

        # run the alignment process
        self.cactusWorkflowArguments.alignmentsID = self.addChildJobFn(make_paf_alignments,
                                            self.cactusWorkflowArguments.speciesTree,
                                            {event: exp.getSequenceID(event) for event in exp.getGenomesWithSequence()},
                                            exp.getRootGenome(),
                                            self.cactusWorkflowArguments.configNode).rv()

        # calculate the total size of the sequences in the alignment problem
        self.cactusWorkflowArguments.totalSequenceSize = sum(os.stat(fileStore.readGlobalFile(exp.getSequenceID(i))).st_size for i in exp.getGenomesWithSequence())

        # hack to get cactus-blast standalone tool working
        if self.standAlone:
            return self.cactusWorkflowArguments

        # now call consolidated to use the resulting alignments
        return self.makeFollowOnPhaseJob(CactusConsolidated, phaseName="consolidated")

class CactusConsolidated(CactusPhasesJob):
    def __init__(self, *args, **kwargs):
        super(CactusConsolidated, self).__init__(*args, **kwargs)
        if "cactusWorkflowArguments" in kwargs and kwargs["cactusWorkflowArguments"].consCores:
            self.cores = kwargs["cactusWorkflowArguments"].consCores

        exp = self.cactusWorkflowArguments.experimentWrapper
        x = self.cactusWorkflowArguments.totalSequenceSize #  total size of sequences in problem

        # Disk requirements
        self.disk = int(3 * x)

        # this is the old caf job's memory function
        memoryPoly = [1.80395944e+01, 7.96042247e+07]
        memoryCap = 120e09
        resource = 0
        for degree, coefficient in enumerate(reversed(memoryPoly)):
            resource += coefficient * (x**degree)
        resource = min(resource, memoryCap)
        self.memory = int(resource)
    
    """Run cactus_consolidated (this spans everythring from bar to hal-genenerator)."""
    def run(self, fileStore):
        # Get the experiment object
        exp = self.cactusWorkflowArguments.experimentWrapper

        # Build up a genome -> fasta map.
        seqMap = {event: fileStore.readGlobalFile(exp.getSequenceID(event)) for event in exp.getGenomesWithSequence()}

        # Get the alignments file
        alignments = fileStore.readGlobalFile(self.cactusWorkflowArguments.alignmentsID)
        logger.info("Alignments file: %s" % alignments)

        # Split the alignments file into primary and secondary
        primary_alignment_file = fileStore.getLocalTempFile()
        system("grep -v 'tp:A:S' {} > {} || true".format(alignments, primary_alignment_file))
        secondary_alignment_file = fileStore.getLocalTempFile()
        system("grep 'tp:A:S' {} > {} || true".format(alignments, secondary_alignment_file))

        # Temporary place to store the output c2h file
        tmpHal = fileStore.getLocalTempFile()
        tmpFasta = fileStore.getLocalTempFile()
        tmpRef = fileStore.getLocalTempFile()

        # Keep inputs in same directory for the docker interface
        tmpConfig = fileStore.getLocalTempFile()
        self.cactusWorkflowArguments.configWrapper.writeXML(tmpConfig)

        messages = runCactusConsolidated(seqMap=seqMap,
                                         newickTreeString=self.cactusWorkflowArguments.speciesTree,
                                         cactusParams=tmpConfig,
                                         alignmentsFile=primary_alignment_file,
                                         outputFile=tmpHal,
                                         outputHalFastaFile=tmpFasta,
                                         outputReferenceFile=tmpRef,
                                         secondaryAlignmentsFile=secondary_alignment_file,
                                         constraintAlignmentsFile=None,
                                         logLevel=getLogLevelString(),
                                         outgroupEvents=exp.getOutgroupGenomes(),
                                         referenceEvent=exp.getRootGenome(),
                                         cores=self.cores)

        # Log back any messages
        for message in messages:
            fileStore.logToMaster(message)

        # Write the temporary output files to the final output
        # At top level--have the final .c2h file
        halID = fileStore.writeGlobalFile(tmpHal)
        fastaID = fileStore.writeGlobalFile(tmpFasta)
        referenceID = fileStore.writeGlobalFile(tmpRef)

        intermediateResultsUrl = getattr(self.cactusWorkflowArguments, 'intermediateResultsUrl', None)
        if intermediateResultsUrl is not None:
            # The user requested to keep the c2h files in a separate place. Export it there.
            url = intermediateResultsUrl + ".c2h"
            fileStore.exportFile(halID, url)

            # The user requested to keep the hal fasta files in a separate place. Export it there.
            url = intermediateResultsUrl + ".hal.fa"
            fileStore.exportFile(fastaID, url)

            # The user requested to keep the reference fasta files in a separate place. Export it there.
            url = intermediateResultsUrl + ".reference.fa"
            fileStore.exportFile(referenceID, url)

        exp.setHalID(halID)
        exp.setHalFastaID(fastaID)
        exp.setReferenceID(referenceID)

        return exp

############################################################
############################################################
############################################################
#Main function
############################################################
############################################################
############################################################
class CactusWorkflowArguments:
    """Object for representing a cactus workflow's arguments
    """
    def __init__(self, options, experimentFile, configNode, seqIDMap, consCores=None):
        #Get a local copy of the experiment file
        self.experimentFile = experimentFile
        self.experimentNode = ET.parse(self.experimentFile).getroot()
        self.experimentWrapper = ExperimentWrapper(self.experimentNode)
        self.alignmentsID = None
        for genome, seqID in list(seqIDMap.items()):
            print(('setting this', genome, seqID))
            self.experimentWrapper.setSequenceID(genome, seqID)
        #Get the species tree
        self.speciesTree = self.experimentNode.attrib["species_tree"]
        #Get any list of 'required species' for the blocks of the cactus.
        self.outgroupEventNames = getOptionalAttrib(self.experimentNode, "outgroup_events", errorIfNotPresent=False)
        #Constraints
        self.constraintsID = getOptionalAttrib(self.experimentNode, "constraintsID", errorIfNotPresent=False)
        #Space to put the path to the directory containing beds of
        #outgroup coverage on ingroups, so that any sequence aligning
        #to an outgroup can be rescued after bar phase
        self.ingroupCoverageIDs = None
        # Same, but for the final bed file
        self.ingroupCoverageID = None
        # If not None, a url prefix to dump database files to
        # (i.e. file:///path/to/prefix). The dumps will be labeled
        # -caf, -avg, etc.
        self.intermediateResultsUrl = options.intermediateResultsUrl
        # Number of cores for cactus_consolidated
        self.consCores = consCores

        #The config node
        self.configNode = configNode
        self.configWrapper = ConfigWrapper(self.configNode)
        #Now deal with the constants that need to be added here
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()
        self.configWrapper.setBuildHal(options.buildHal)
        self.configWrapper.setBuildFasta(options.buildFasta)

def addCactusWorkflowOptions(parser):
    parser.add_argument("--experiment", dest="experimentFile",
                      help="The file containing a link to the experiment parameters")
    parser.add_argument("--buildHal", dest="buildHal", action="store_true",
                      help="Build a hal file", default=False)
    parser.add_argument("--buildFasta", dest="buildFasta", action="store_true",
                      help="Build a fasta file of the input sequences (and reference sequence, used with hal output)",
                      default=False)
    parser.add_argument("--intermediateResultsUrl",
                        help="URL prefix to save intermediate results like DB dumps to (e.g. "
                        "prefix-dump-caf, prefix-dump-avg, etc.)", default=None)

class RunCactusPreprocessorThenCactusSetup(RoundedJob):
    def __init__(self, options, cactusWorkflowArguments):
        RoundedJob.__init__(self)
        self.options = options
        self.cactusWorkflowArguments = cactusWorkflowArguments

    def run(self, fileStore):
        exp = self.cactusWorkflowArguments.experimentWrapper
        genomes = exp.getGenomesWithSequence()
        originalSequenceIDs = [exp.getSequenceID(genome) for genome in genomes]
        preprocessedSeqIDs = self.addChild(CactusPreprocessor(originalSequenceIDs, self.cactusWorkflowArguments.configNode)).rv()
        self.addFollowOn(AfterPreprocessing(self.cactusWorkflowArguments, exp, genomes, preprocessedSeqIDs))

class AfterPreprocessing(RoundedJob):
    def __init__(self, cactusWorkflowArguments, eW, genomes, preprocessedSeqIDs):
        RoundedJob.__init__(self)
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.eW = eW
        self.genomes = genomes
        self.preprocessedSeqIDs = preprocessedSeqIDs

    def run(self, fileStore):
        #Now replace the input sequences with the preprocessed sequences in the experiment wrapper
        for genome, preprocessedSeqID in zip(self.genomes, self.preprocessedSeqIDs):
            self.eW.setSequenceID(genome, preprocessedSeqID)
        fileStore.logToMaster("doTrimStrategy() = %s, outgroupEventNames = %s" % (self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy(), self.cactusWorkflowArguments.outgroupEventNames))
        # Use the trimming strategy to blast ingroups vs outgroups.
        self.addFollowOn(CactusPafAlign(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName="blast"))


def runCactusWorkflow(args):
    ##########################################
    #Construct the arguments.
    ##########################################

    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    options = parser.parse_args(args)
    options.disableCaching = True
    set_logging_from_options(options)

    experimentWrapper = ExperimentWrapper(ET.parse(options.experimentFile).getroot())
    with Toil(options) as toil:
        seqIDMap = dict()
        for name in experimentWrapper.getGenomesWithSequence():
            fullSeq = getTempFile()
            seq = experimentWrapper.getSequenceID(name)
            if os.path.isdir(seq):
                catFiles([os.path.join(seq, seqFile) for seqFile in os.listdir(seq)], fullSeq)
            else:
                fullSeq = seq
            experimentWrapper.setSequenceID(name, toil.importFile(makeURL(fullSeq)))
            print((name, experimentWrapper.getSequenceID(name)))

        experimentWrapper.writeXML(options.experimentFile)

        configNode = ET.parse(experimentWrapper.getConfigPath()).getroot()
        print(seqIDMap)
        cactusWorkflowArguments = CactusWorkflowArguments(options, experimentFile=options.experimentFile, configNode=configNode, seqIDMap=seqIDMap)

        toil.start(RunCactusPreprocessorThenCactusSetup(options, cactusWorkflowArguments))

if __name__ == '__main__':
    runCactusWorkflow(sys.argv)
