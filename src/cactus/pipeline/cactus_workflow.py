#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework to structure all the related wrappers.
"""

import os
import sys
import xml.etree.ElementTree as ET
import math
import time
import copy
from argparse import ArgumentParser
from operator import itemgetter

from sonLib.bioio import newickTreeParser

from toil.lib.bioio import getTempFile
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions
from toil.lib.bioio import system
from sonLib.bioio import catFiles
from sonLib.bioio import getLogLevelString

from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.common import RunAsFollowOn
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusCaf
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import runCactusSplitFlowersBySecondaryGrouping
from cactus.shared.common import encodeFlowerNames
from cactus.shared.common import decodeFirstFlowerName
from cactus.shared.common import runCactusConvertAlignmentToCactus
from cactus.shared.common import runCactusPhylogeny
from cactus.shared.common import runCactusBar
from cactus.shared.common import runCactusMakeNormal
from cactus.shared.common import runCactusReference
from cactus.shared.common import runCactusAddReferenceCoordinates
from cactus.shared.common import runCactusCheck
from cactus.shared.common import runCactusHalGenerator
from cactus.shared.common import runCactusFlowerStats
from cactus.shared.common import runCactusSecondaryDatabase
from cactus.shared.common import runCactusFastaGenerator
from cactus.shared.common import findRequiredNode
from cactus.shared.common import runConvertAlignmentsToInternalNames
from cactus.shared.common import runStripUniqueIDs
from cactus.shared.common import RoundedJob
from cactus.shared.common import readGlobalFileWithoutCache

from cactus.blast.blast import BlastIngroupsAndOutgroups
from cactus.blast.blast import BlastOptions
from cactus.blast.mappingQualityRescoringAndFiltering import mappingQualityRescoring

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.ktserverToil import KtServerService
from cactus.pipeline.ktserverControl import stopKtserver

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

    def getOptionalPhaseAttrib(self, attribName, typeFn=None, default=None):
        """Gets an optional attribute of the phase node.
        """
        return getOptionalAttrib(node=self.phaseNode, attribName=attribName, typeFn=typeFn, default=default)

    def getOptionalJobAttrib(self, attribName, typeFn=None, default=None):
        """Gets an optional attribute of the job node.
        """
        return getOptionalAttrib(node=self.jobNode, attribName=attribName, typeFn=typeFn, default=default)

    def addService(self, job):
        """Works around toil issue #1695, returning a Job rather than a Promise."""
        super(CactusJob, self).addService(job)
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

    def makeRecursiveChildJob(self, job, launchSecondaryKtForRecursiveJob=False):
        newChild = job(phaseNode=extractNode(self.phaseNode),
                       constantsNode=extractNode(self.constantsNode),
                       cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString,
                       flowerNames=encodeFlowerNames((self.topFlowerName,)),
                       flowerSizes=[self.cactusWorkflowArguments.totalSequenceSize],
                       overlarge=True,
                       cactusWorkflowArguments=self.cactusWorkflowArguments)

        if launchSecondaryKtForRecursiveJob and ExperimentWrapper(self.cactusWorkflowArguments.experimentNode).getDbType() == "kyoto_tycoon":
            cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)
            memory = max(2500000000, self.evaluateResourcePoly([4.10201882, 2.01324291e+08]))
            cpu = cw.getKtserverCpu(default=0.1)
            dbElem = ExperimentWrapper(self.cactusWorkflowArguments.scratchDbElemNode)
            dbString = self.addService(KtServerService(dbElem=dbElem, isSecondary=True, memory=memory, cores=cpu)).rv(0)
            newChild.phaseNode.attrib["secondaryDatabaseString"] = dbString
            return self.addChild(newChild).rv()
        else:
            return self.addChild(newChild).rv()

    def makeFollowOnPhaseJob(self, job, phaseName):
        return self.addFollowOn(job(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName=phaseName,
                                    topFlowerName=self.topFlowerName, halID=self.halID, fastaID=self.fastaID)).rv()

    def runPhase(self, recursiveJob, nextPhaseJob, nextPhaseName, doRecursion=True, launchSecondaryKtForRecursiveJob=False):
        """
        Adds a recursive child job and then a follow-on phase job. Returns the result of the follow-on
        phase job.
        """
        logger.info("Starting %s phase job at %s seconds (recursing = %i)" % (self.phaseNode.tag, time.time(), doRecursion))
        if doRecursion:
            self.makeRecursiveChildJob(recursiveJob, launchSecondaryKtForRecursiveJob)
        return self.makeFollowOnPhaseJob(job=nextPhaseJob, phaseName=nextPhaseName)

    def makeFollowOnCheckpointJob(self, checkpointConstructor, phaseName, ktServerDump=None):
        """Add a follow-on checkpoint phase."""
        return self.addFollowOn(checkpointConstructor(\
                   phaseName=phaseName, ktServerDump=ktServerDump,
                   cactusWorkflowArguments=self.cactusWorkflowArguments,
                   topFlowerName=self.topFlowerName,
                   halID=self.halID, fastaID=self.fastaID)).rv()

    def getPhaseNumber(self):
        return len(self.cactusWorkflowArguments.configNode.findall(self.phaseNode.tag))

    def setupSecondaryDatabase(self):
        """Setup the secondary database
        """
        confXML = ET.fromstring(self.cactusWorkflowArguments.secondaryDatabaseString)
        dbElem = DbElemWrapper(confXML)
        if dbElem.getDbType() != "kyoto_tycoon":
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=True)

    def cleanupSecondaryDatabase(self):
        """Cleanup the secondary database
        """
        confXML = ET.fromstring(self.cactusWorkflowArguments.secondaryDatabaseString)
        dbElem = DbElemWrapper(confXML)
        if dbElem.getDbType() != "kyoto_tycoon":
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=False)

class CactusCheckpointJob(CactusPhasesJob):
    """Special "checkpoint" phase job, launching and restoring the primary Cactus DB.

    Inherit from this and run `runPhaseWithPrimaryDB` to start a new
    primary DB before running the given phase.

    Meant to provide a restore point in case of pipeline failure. Note
    that the checkpointed job is technically this job's child, which
    starts the database.
    """
    def __init__(self, ktServerDump=None, *args, **kwargs):
        self.ktServerDump = ktServerDump
        super(CactusCheckpointJob, self).__init__(*args, **kwargs)

    def runPhaseWithPrimaryDB(self, jobConstructor):
        """Start and load a new primary DB before running the given phase.
        """
        job = jobConstructor(cactusWorkflowArguments=self.cactusWorkflowArguments,
                             phaseName=self.phaseName, topFlowerName=self.topFlowerName)
        startDBJob = StartPrimaryDB(job, ktServerDump=self.ktServerDump,
                                    cactusWorkflowArguments=self.cactusWorkflowArguments,
                                    phaseName=self.phaseName, topFlowerName=self.topFlowerName)
        promise = self.addChild(startDBJob)
        return promise

class StartPrimaryDB(CactusPhasesJob):
    """Launches a primary Cactus DB."""
    def __init__(self, nextJob, ktServerDump=None, *args, **kwargs):
        self.nextJob = nextJob
        self.ktServerDump = ktServerDump
        kwargs['checkpoint'] = True
        kwargs['preemptable'] = False
        super(StartPrimaryDB, self).__init__(*args, **kwargs)

    def run(self, fileStore):
        cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)

        if self.cactusWorkflowArguments.experimentWrapper.getDbType() == "kyoto_tycoon":
            memory = max(2500000000, self.evaluateResourcePoly([4.10201882, 2.01324291e+08]))
            cores = cw.getKtserverCpu(default=0.1)
            dbElem = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
            service = self.addService(KtServerService(dbElem=dbElem,
                                                      existingSnapshotID=self.ktServerDump,
                                                      isSecondary=False,
                                                      memory=memory, cores=cores))
            dbString = service.rv(0)
            snapshotID = service.rv(1)
            self.nextJob.cactusWorkflowArguments.cactusDiskDatabaseString = dbString
            # TODO: This part needs to be cleaned up
            self.nextJob.cactusWorkflowArguments.snapshotID = snapshotID
            return self.addChild(self.nextJob).rv()
        else:
            return self.addFollowOn(self.nextJob).rv()

class SavePrimaryDB(CactusPhasesJob):
    """Saves the DB to a file and clears the DB."""
    def __init__(self, *args, **kwargs):
        super(SavePrimaryDB, self).__init__(*args, **kwargs)

    def run(self, fileStore):
        stats = runCactusFlowerStats(cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString,
                                     flowerName=0)
        fileStore.logToMaster("At end of %s phase, got stats %s" % (self.phaseName, stats))
        dbElem = DbElemWrapper(ET.fromstring(self.cactusWorkflowArguments.cactusDiskDatabaseString))
        # Send the terminate message
        stopKtserver(dbElem)
        # Wait for the file to appear in the right place. This may take a while
        while True:
            with fileStore.readGlobalFileStream(self.cactusWorkflowArguments.snapshotID) as f:
                if f.read(1) != b'':
                    # The file is no longer empty
                    break
            time.sleep(10)
        # We have the file now
        intermediateResultsUrl = getattr(self.cactusWorkflowArguments, 'intermediateResultsUrl', None)
        if intermediateResultsUrl is not None:
            # The user requested to keep the DB dumps in a separate place. Export it there.
            url = intermediateResultsUrl + "-dump-" + self.phaseName
            fileStore.exportFile(self.cactusWorkflowArguments.snapshotID, url)
        return self.cactusWorkflowArguments.snapshotID

class CactusRecursionJob(CactusJob):
    """Base recursive job for traversals up and down the cactus tree.
    """
    flowerFeatures = lambda self: {'flowerGroupSize': sum(self.flowerSizes),
                                   'maxFlowerSize': max(self.flowerSizes),
                                   'numFlowers': len(self.flowerSizes)}
    featuresFn = flowerFeatures
    feature = 'flowerGroupSize'
    maxSequenceSizeOfFlowerGroupingDefault = 1000000
    def __init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, flowerSizes, overlarge=False, precomputedAlignmentIDs=None, checkpoint = False, cactusWorkflowArguments=None, preemptable=True, memPoly=None):
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.flowerNames = flowerNames
        self.flowerSizes = flowerSizes
        self.cactusWorkflowArguments = cactusWorkflowArguments

        #need to do this because the alignment IDs are jobstore promises, and can't
        #be stored in the config XML until they are respolved into actual IDs, which doesn't
        #happen until the follow-on job after CactusBarWrapperLarge
        self.precomputedAlignmentIDs = precomputedAlignmentIDs

        CactusJob.__init__(self, phaseNode=phaseNode, constantsNode=constantsNode, overlarge=overlarge,
                           checkpoint=checkpoint, preemptable=preemptable)

    def makeFollowOnRecursiveJob(self, job, phaseNode=None):
        """Sets the followon to the given recursive job
        """
        if phaseNode == None:
            phaseNode = self.phaseNode
        return self.addFollowOn(job(phaseNode=phaseNode, constantsNode=self.constantsNode,
                                    cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                    flowerNames=self.flowerNames, flowerSizes=self.flowerSizes,
                                    overlarge=self.overlarge,
                                    precomputedAlignmentIDs=self.precomputedAlignmentIDs,
                                    cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()

    def makeFollowOnRecursiveJobWithPromisedRequirements(self, job, phaseNode=None):
        """
        Toil's PromisedRequirements don't actually work with real job
        classes, only functions. So this is a hacky way of working
        around that by introducing our own level of indirection.
        """
        if phaseNode == None:
            phaseNode = self.phaseNode
        return self.addFollowOn(RunAsFollowOn(job, phaseNode=phaseNode, constantsNode=self.constantsNode,
                                              cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                              flowerNames=self.flowerNames, flowerSizes=self.flowerSizes,
                                              overlarge=self.overlarge,
                                              precomputedAlignmentIDs=self.precomputedAlignmentIDs,
                                              cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()

    def makeChildJobs(self, flowersAndSizes, job, overlargeJob=None,
                      phaseNode=None):
        """Make a set of child jobs for a given set of flowers and chosen child job
        """
        if overlargeJob == None:
            overlargeJob = job
        if phaseNode == None:
            phaseNode = self.phaseNode

        logger.info("Make wrapper jobs: There are %i flowers" % len(flowersAndSizes))
        for overlarge, flowerNames, flowerSizes in flowersAndSizes:
            if overlarge: #Make sure large flowers are on their own, in their own job
                flowerStatsString = runCactusFlowerStats(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                                         flowerName=decodeFirstFlowerName(flowerNames))
                self._fileStore.logToMaster("Adding an oversize flower for job class %s and stats %s" \
                                 % (overlargeJob, flowerStatsString))
                self.addChild(overlargeJob(cactusDiskDatabaseString=
                                           self.cactusDiskDatabaseString,
                                           phaseNode=phaseNode,
                                           constantsNode=self.constantsNode,
                                           flowerNames=flowerNames,
                                           flowerSizes=flowerSizes,
                                           overlarge=True,
                                           cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
            else:
                logger.info("Adding recursive flower job")
                self.addChild(job(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                  phaseNode=phaseNode, constantsNode=self.constantsNode,
                                  flowerNames=flowerNames,
                                  flowerSizes=flowerSizes,
                                  overlarge=False,
                                  cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()

    def makeRecursiveJobs(self, fileStore=None, job=None, phaseNode=None):
        """Make a set of child jobs for a given set of parent flowers.
        """
        if job == None:
            job = self.__class__
        jobNode = getJobNode(self.phaseNode, job)
        flowersAndSizes=runCactusGetFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                            features=self.featuresFn(),
                                            jobName=job.__name__,
                                            fileStore=fileStore,
                                            flowerNames=self.flowerNames,
                                            minSequenceSizeOfFlower=getOptionalAttrib(jobNode, "minFlowerSize", int, 0),
                                            maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(jobNode, "maxFlowerGroupSize", int,
                                            default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault),
                                            maxSequenceSizeOfSecondaryFlowerGrouping=getOptionalAttrib(jobNode, "maxFlowerWrapperGroupSize", int,
                                            default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault))
        return self.makeChildJobs(flowersAndSizes=flowersAndSizes,
                              job=job, phaseNode=phaseNode)

    def makeExtendingJobs(self, job, fileStore=None, overlargeJob=None, phaseNode=None):
        """Make set of child jobs that extend the current cactus tree.
        """

        jobNode = getJobNode(self.phaseNode, job)
        flowersAndSizes=runCactusExtendFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                              features=self.featuresFn(),
                                              jobName=job.__name__,
                                              fileStore=fileStore,
                                              flowerNames=self.flowerNames,
                                              minSequenceSizeOfFlower=getOptionalAttrib(jobNode, "minFlowerSize", int, 0),
                                              maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(jobNode, "maxFlowerGroupSize", int,
                                              default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault))
        return self.makeChildJobs(flowersAndSizes=flowersAndSizes,
                                  job=job, overlargeJob=overlargeJob,
                                  phaseNode=phaseNode)

    def makeWrapperJobs(self, job, overlargeJob=None, phaseNode=None):
        """Takes the list of flowers for a recursive job and splits them up to fit the given wrapper job(s).
        """
        splitFlowerNames = runCactusSplitFlowersBySecondaryGrouping(self.flowerNames)
        # We've split the flower names up into groups, but now we need
        # to put the flower sizes in so that they correspond with
        # their flower.
        flowersAndSizes = []
        flowersSoFar = 0
        for overlarge, flowerNames in splitFlowerNames:
            # Number of flowers in this grouping.
            numFlowers = int(flowerNames.split()[0])
            flowersAndSizes += [(overlarge, flowerNames, self.flowerSizes[flowersSoFar:flowersSoFar + numFlowers])]
            flowersSoFar += numFlowers
        totalFlowers = int(self.flowerNames.split()[0])
        assert flowersSoFar == totalFlowers, \
               "Didn't process all flowers while going through a secondary grouping."
        return self.makeChildJobs(flowersAndSizes=flowersAndSizes,
                                  job=job, overlargeJob=overlargeJob,
                                  phaseNode=phaseNode)

############################################################
############################################################
############################################################
##The blast phase that uses the trimming strategy.
############################################################
############################################################
############################################################

def prependUniqueIDs(fas, outputDir):
    """Prepend unique ints to fasta headers.

    (prepend rather than append since trimmed outgroups have a start
    token appended, which complicates removal slightly)
    """
    uniqueID = 0
    ret = []
    for fa in fas:
        outPath = os.path.join(outputDir, os.path.basename(fa))
        out = open(outPath, 'w')
        for line in open(fa):
            if len(line) > 0 and line[0] == '>':
                header = "id=%d|%s" % (uniqueID, line[1:-1])
                out.write(">%s\n" % header)
            else:
                out.write(line)
        ret.append(outPath)
        uniqueID += 1
    return ret

def setupDivergenceArgs(cactusWorkflowArguments):
    #Adapt the config file to use arguments for the appropriate divergence distance
    cactusWorkflowArguments.longestPath = getLongestPath(newickTreeParser(cactusWorkflowArguments.speciesTree))
    if cactusWorkflowArguments.outgroupEventNames == None:
        distanceToAddToRootAlignment = getOptionalAttrib(cactusWorkflowArguments.configNode, "distanceToAddToRootAlignment", float, 0.0)
        cactusWorkflowArguments.longestPath += distanceToAddToRootAlignment
    cw = ConfigWrapper(cactusWorkflowArguments.configNode)
    cw.substituteAllDivergenceContolledParametersWithLiterals(cactusWorkflowArguments.longestPath)

def setupFilteringByIdentity(cactusWorkflowArguments):
    #Filter by identity
    cafNode = findRequiredNode(cactusWorkflowArguments.configNode, "caf")
    if getOptionalAttrib(cafNode, "filterByIdentity", bool, False): #Do the identity filtering
        adjustedPath = max(float(cafNode.attrib["identityRatio"]) * cactusWorkflowArguments.longestPath,
        float(cafNode.attrib["minimumDistance"]))
        identity = str(100 - math.ceil(100 * inverseJukesCantor(adjustedPath)))
        cafNode.attrib["lastzArguments"] = cafNode.attrib["lastzArguments"] + (" --identity=%s" % identity)


class CactusTrimmingBlastPhase(CactusPhasesJob):
    """Blast ingroups vs outgroups using the trimming strategy before
    running cactus setup.
    """
    def run(self, fileStore):
        fileStore.logToMaster("Running blast using the trimming strategy")

        exp = self.cactusWorkflowArguments.experimentWrapper
        ingroupsAndOriginalIDs = [(g, exp.getSequenceID(g)) for g in exp.getGenomesWithSequence() if g not in exp.getOutgroupGenomes()]
        outgroupsAndOriginalIDs = [(g, exp.getSequenceID(g)) for g in exp.getOutgroupGenomes()]
        from sonLib.nxnewick import NXNewick
        print((NXNewick().writeString(exp.getTree())))
        print((exp.getRootGenome()))
        print(ingroupsAndOriginalIDs)
        print(outgroupsAndOriginalIDs)
        sequences = [fileStore.readGlobalFile(id) for id in map(itemgetter(1), ingroupsAndOriginalIDs + outgroupsAndOriginalIDs)]
        self.cactusWorkflowArguments.totalSequenceSize = sum(os.stat(x).st_size for x in sequences)

        renamedInputSeqDir = fileStore.getLocalTempDir()
        uniqueFas = prependUniqueIDs(sequences, renamedInputSeqDir)
        uniqueFaIDs = [fileStore.writeGlobalFile(seq, cleanup=True) for seq in uniqueFas]

        # Set the uniquified IDs for the ingroups and outgroups
        ingroupsAndNewIDs = list(zip(list(map(itemgetter(0), ingroupsAndOriginalIDs)), uniqueFaIDs[:len(ingroupsAndOriginalIDs)]))
        outgroupsAndNewIDs = list(zip(list(map(itemgetter(0), outgroupsAndOriginalIDs)), uniqueFaIDs[len(ingroupsAndOriginalIDs):]))

        # Change the blast arguments depending on the divergence
        setupDivergenceArgs(self.cactusWorkflowArguments)
        setupFilteringByIdentity(self.cactusWorkflowArguments)

        cafNode = findRequiredNode(self.cactusWorkflowArguments.configNode, "caf")

        # FIXME: this is really ugly and steals the options from the caf tag
        blastJob = self.addChild(BlastIngroupsAndOutgroups(
            BlastOptions(chunkSize=getOptionalAttrib(cafNode, "chunkSize", int),
                         overlapSize=getOptionalAttrib(cafNode, "overlapSize", int),
                         lastzArguments=getOptionalAttrib(cafNode, "lastzArguments"),
                         compressFiles=getOptionalAttrib(cafNode, "compressFiles", bool),
                         realign=getOptionalAttrib(cafNode, "realign", bool),
                         realignArguments=getOptionalAttrib(cafNode, "realignArguments"),
                         memory=getOptionalAttrib(cafNode, "lastzMemory", int, sys.maxsize),
                         smallDisk=getOptionalAttrib(cafNode, "lastzSmallDisk", int, sys.maxsize),
                         largeDisk=getOptionalAttrib(cafNode, "lastzLargeDisk", int, sys.maxsize),
                         minimumSequenceLength=getOptionalAttrib(cafNode, "minimumSequenceLengthForBlast", int, 1),
                         trimFlanking=self.getOptionalPhaseAttrib("trimFlanking", int, 10),
                         trimMinSize=self.getOptionalPhaseAttrib("trimMinSize", int, 0),
                         trimThreshold=self.getOptionalPhaseAttrib("trimThreshold", float, 0.8),
                         trimWindowSize=self.getOptionalPhaseAttrib("trimWindowSize", int, 10),
                         trimOutgroupFlanking=self.getOptionalPhaseAttrib("trimOutgroupFlanking", int, 100),
                         trimOutgroupDepth=self.getOptionalPhaseAttrib("trimOutgroupDepth", int, 1),
                         keepParalogs=self.getOptionalPhaseAttrib("keepParalogs", bool, False)),
            list(map(itemgetter(0), ingroupsAndNewIDs)), list(map(itemgetter(1), ingroupsAndNewIDs)),
            list(map(itemgetter(0), outgroupsAndNewIDs)), list(map(itemgetter(1), outgroupsAndNewIDs))))

        # Alignment post processing to filter alignments
        if getOptionalAttrib(cafNode, "runMapQFiltering", bool, False):
            minimumMapQValue=getOptionalAttrib(cafNode, "minimumMapQValue", float, 0.0)
            maxAlignmentsPerSite=getOptionalAttrib(cafNode, "maxAlignmentsPerSite", int, 1)
            alpha=getOptionalAttrib(cafNode, "alpha", float, 1.0)
            fileStore.logToMaster("Running mapQ uniquifying with parameters, minimumMapQValue: %s, maxAlignmentsPerSite %s, alpha: %s" %
                                  (minimumMapQValue, maxAlignmentsPerSite, alpha))
            blastJob = blastJob.encapsulate() # Encapsulate to ensure that blast Job and all its successors
            # run before mapQ
            mapQJob = blastJob.addFollowOnJobFn(mappingQualityRescoring, blastJob.rv(0),
                                                minimumMapQValue=minimumMapQValue,
                                                maxAlignmentsPerSite=maxAlignmentsPerSite,
                                                alpha=alpha,
                                                logLevel=getLogLevelString(),
                                                preemptable=True)
            self.cactusWorkflowArguments.alignmentsID = mapQJob.rv(0)
            self.cactusWorkflowArguments.secondaryAlignmentsID = mapQJob.rv(1)
        else:
            fileStore.logToMaster("Not running mapQ filtering")
            self.cactusWorkflowArguments.alignmentsID = blastJob.rv(0)
            self.cactusWorkflowArguments.secondaryAlignmentsID = None

        self.cactusWorkflowArguments.outgroupFragmentIDs = blastJob.rv(1)
        self.cactusWorkflowArguments.ingroupCoverageIDs = blastJob.rv(2)

        # Update the experiment wrapper to point to the uniquified IDs
        # for the ingroups, and the trimmed results for the outgroups.
        for genome, seqID in ingroupsAndNewIDs:
            self.cactusWorkflowArguments.experimentWrapper.setSequenceID(genome, seqID)

        updateJob = blastJob.addFollowOnJobFn(updateExpWrapperForOutgroups, self.cactusWorkflowArguments.experimentWrapper,
                                              list(map(itemgetter(0), outgroupsAndNewIDs)), self.cactusWorkflowArguments.outgroupFragmentIDs)
        self.cactusWorkflowArguments.experimentWrapper = updateJob.rv()

        return self.makeFollowOnCheckpointJob(CactusSetupCheckpoint, "setup")

def updateExpWrapperForOutgroups(job, expWrapper, outgroupGenomes, outgroupFragmentIDs):
    for genome, outgroupFragmentID in zip(outgroupGenomes, outgroupFragmentIDs):
        expWrapper.setSequenceID(genome, outgroupFragmentID)
    return expWrapper

############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################

def getLongestPath(node, distance=0.0):
    """Identify the longest path from the mrca of the leaves of the species tree.
    """
    i, j = distance, distance
    if node.left != None:
        i = getLongestPath(node.left, abs(node.left.distance)) + distance
    if node.right != None:
        j = getLongestPath(node.right, abs(node.right.distance)) + distance
    return max(i, j)

class CactusSetupCheckpoint(CactusCheckpointJob):
    """Start a new DB, run the setup and CAF phases, save the DB, then launch the BAR checkpoint."""
    def run(self, fileStore):
        ktServerDump = self.runPhaseWithPrimaryDB(CactusSetupPhase).rv()
        return self.makeFollowOnCheckpointJob(CactusBarCheckpoint, "bar", ktServerDump=ktServerDump)

class CactusSetupPhase(CactusPhasesJob):
    """Initialises the cactus database and adapts the config file for the run."""
    def run(self, fileStore):
        experiment = self.cactusWorkflowArguments.experimentWrapper
        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupDivergenceArgs(self.cactusWorkflowArguments)
        # Build up a genome -> fasta map.
        seqIDMap = dict((genome, experiment.getSequenceID(genome)) for genome in experiment.getGenomesWithSequence())
        seqMap = dict((genome, fileStore.readGlobalFile(id)) for genome, id in list(seqIDMap.items()))

        messages = runCactusSetup(cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString,
                                  seqMap=seqMap,
                                  newickTreeString=self.cactusWorkflowArguments.speciesTree,
                                  outgroupEvents=experiment.getOutgroupGenomes(),
                                  makeEventHeadersAlphaNumeric=self.getOptionalPhaseAttrib("makeEventHeadersAlphaNumeric", bool, False))
        for message in messages:
            logger.info(message)
        return self.makeFollowOnPhaseJob(CactusCafPhase, "caf")

############################################################
############################################################
############################################################
#The CAF phase.
#
#Creates the reconstruction structure with blocks
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

class CactusCafPhase(CactusPhasesJob):
    memoryPoly = [2.51087392e+00, 4.49616219e+08]

    def __init__(self, **kwargs):
        # Ensure non-preemptability, since this job takes a long time
        # when there are a ton of alignments.
        kwargs['preemptable'] = False
        super(CactusCafPhase, self).__init__(**kwargs)

    def run(self, fileStore):
        if len(self.cactusWorkflowArguments.ingroupCoverageIDs) > 0:
            # Convert the bed files to use 64-bit cactus Names instead
            # of the headers. Ideally this should belong in the bar
            # phase but we run stripUniqueIDs before then.
            bedFiles = [fileStore.readGlobalFile(path) for path in self.cactusWorkflowArguments.ingroupCoverageIDs]
            tempFile = fileStore.getLocalTempFile()
            system("cat %s > %s" % (" ".join(bedFiles), tempFile))
            ingroupCoverageFile = fileStore.getLocalTempFile()
            runConvertAlignmentsToInternalNames(self.cactusWorkflowArguments.cactusDiskDatabaseString, tempFile, ingroupCoverageFile, self.topFlowerName, isBedFile=True)
            self.cactusWorkflowArguments.ingroupCoverageID = fileStore.writeGlobalFile(ingroupCoverageFile)

        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupFilteringByIdentity(self.cactusWorkflowArguments)
        #Setup any constraints
        if self.cactusWorkflowArguments.constraintsID != None: #Setup the constraints arg
            constraintsFile = fileStore.readGlobalFile(self.cactusWorkflowArguments.constraintsID)
            newConstraintsFile = fileStore.getLocalTempFile()
            runCactusConvertAlignmentToCactus(self.cactusWorkflowArguments.cactusDiskDatabaseString,
                                              constraintsFile, newConstraintsFile)
            self.cactusWorkflowArguments.constraintsID = fileStore.writeGlobalFile(newConstraintsFile, cleanup=True)

        assert self.getPhaseNumber() == 1

        # Convert the cigar files to use 64-bit cactus Names instead of the headers.

        # Primary alignments first
        alignmentsFile = fileStore.readGlobalFile(self.cactusWorkflowArguments.alignmentsID)
        convertedAlignmentsFile = fileStore.getLocalTempFile()
        runConvertAlignmentsToInternalNames(cactusDiskString=self.cactusWorkflowArguments.cactusDiskDatabaseString, alignmentsFile=alignmentsFile, outputFile=convertedAlignmentsFile, flowerName=self.topFlowerName)
        fileStore.logToMaster("Converted headers of cigar file %s to internal names, new file %s" % (self.cactusWorkflowArguments.alignmentsID, convertedAlignmentsFile))
        self.cactusWorkflowArguments.alignmentsID = fileStore.writeGlobalFile(convertedAlignmentsFile, cleanup=True)

        # Now any secondary alignments
        if self.cactusWorkflowArguments.secondaryAlignmentsID != None:
            secondaryAlignmentsFile = fileStore.readGlobalFile(self.cactusWorkflowArguments.secondaryAlignmentsID)
            convertedAlignmentsFile = fileStore.getLocalTempFile()
            runConvertAlignmentsToInternalNames(cactusDiskString=self.cactusWorkflowArguments.cactusDiskDatabaseString, alignmentsFile=secondaryAlignmentsFile, outputFile=convertedAlignmentsFile, flowerName=self.topFlowerName)
            fileStore.logToMaster("Converted headers of secondary cigar file %s to internal names, new file %s" % (self.cactusWorkflowArguments.secondaryAlignmentsID, convertedAlignmentsFile))
            self.cactusWorkflowArguments.secondaryAlignmentsID = fileStore.writeGlobalFile(convertedAlignmentsFile, cleanup=True)

        # While we're at it, remove the unique IDs prepended to
        # the headers inside the cactus DB.
        runStripUniqueIDs(self.cactusWorkflowArguments.cactusDiskDatabaseString)

        return self.runPhase(CactusCafWrapper, SavePrimaryDB, "caf")

class CactusCafWrapper(CactusRecursionJob):
    """Runs cactus_caf on one flower and one alignment file.
    """
    featuresFn = lambda self: {'alignmentsSize': self.cactusWorkflowArguments.alignmentsID.size}
    memoryPoly = [1.80395944e+01, 7.96042247e+07]
    memoryCap = 120e09
    feature = 'totalSequenceSize'

    def __init__(self, **kwargs):
        # We want to make sure this job is *not* preempted, because it's
        # the longest-running job in the entire process, taking up to
        # 12 hours.
        kwargs['preemptable'] = False
        super(CactusCafWrapper, self).__init__(**kwargs)

    def runCactusCafInWorkflow(self, alignmentFile, secondaryAlignmentFile, fileStore, constraints=None):
        debugFilePath = self.getOptionalPhaseAttrib("phylogenyDebugPrefix")
        exp = self.cactusWorkflowArguments.experimentWrapper
        if debugFilePath != None:
            debugFilePath += exp.getRootGenome()
        messages = runCactusCaf(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                          features=self.featuresFn(),
                          fileStore=fileStore,
                          jobName=self.__class__.__name__,
                          alignments=alignmentFile,
                          secondaryAlignments=secondaryAlignmentFile,
                          flowerNames=self.flowerNames,
                          constraints=constraints,
                          annealingRounds=self.getOptionalPhaseAttrib("annealingRounds"),
                          deannealingRounds=self.getOptionalPhaseAttrib("deannealingRounds"),
                          trim=self.getOptionalPhaseAttrib("trim"),
                          minimumTreeCoverage=self.getOptionalPhaseAttrib("minimumTreeCoverage", float),
                          blockTrim=self.getOptionalPhaseAttrib("blockTrim", float),
                          minimumBlockDegree=self.getOptionalPhaseAttrib("minimumBlockDegree", int),
                          minimumIngroupDegree=self.getOptionalPhaseAttrib("minimumIngroupDegree", int),
                          minimumOutgroupDegree=self.getOptionalPhaseAttrib("minimumOutgroupDegree", int),
                          alignmentFilter=self.getOptionalPhaseAttrib("alignmentFilter"),
                          lastzArguments=self.getOptionalPhaseAttrib("lastzArguments"),
                          minimumSequenceLengthForBlast=self.getOptionalPhaseAttrib("minimumSequenceLengthForBlast", int, 1),
                          maxAdjacencyComponentSizeRatio=self.getOptionalPhaseAttrib("maxAdjacencyComponentSizeRatio", float),
                          minLengthForChromosome=self.getOptionalPhaseAttrib("minLengthForChromosome", int),
                          proportionOfUnalignedBasesForNewChromosome=self.getOptionalPhaseAttrib("proportionOfUnalignedBasesForNewChromosome", float),
                          maximumMedianSequenceLengthBetweenLinkedEnds=self.getOptionalPhaseAttrib("maximumMedianSequenceLengthBetweenLinkedEnds", int),
                          realign=self.getOptionalPhaseAttrib("realign", bool),
                          realignArguments=self.getOptionalPhaseAttrib("realignArguments"),
                          phylogenyNumTrees=self.getOptionalPhaseAttrib("phylogenyNumTrees", int, 1),
                          phylogenyRootingMethod=self.getOptionalPhaseAttrib("phylogenyRootingMethod"),
                          phylogenyScoringMethod=self.getOptionalPhaseAttrib("phylogenyScoringMethod"),
                          phylogenyBreakpointScalingFactor=self.getOptionalPhaseAttrib("phylogenyBreakpointScalingFactor"),
                          phylogenySkipSingleCopyBlocks=self.getOptionalPhaseAttrib("phylogenySkipSingleCopyBlocks", bool),
                          phylogenyMaxBaseDistance=self.getOptionalPhaseAttrib("phylogenyMaxBaseDistance"),
                          phylogenyMaxBlockDistance=self.getOptionalPhaseAttrib("phylogenyMaxBlockDistance"),
                          phylogenyDebugFile=debugFilePath,
                          phylogenyKeepSingleDegreeBlocks=self.getOptionalPhaseAttrib("phylogenyKeepSingleDegreeBlocks", bool),
                          phylogenyTreeBuildingMethod=self.getOptionalPhaseAttrib("phylogenyTreeBuildingMethod"),
                          phylogenyCostPerDupPerBase=self.getOptionalPhaseAttrib("phylogenyCostPerDupPerBase"),
                          phylogenyCostPerLossPerBase=self.getOptionalPhaseAttrib("phylogenyCostPerLossPerBase"),
                          referenceEventHeader=exp.getRootGenome(),
                          phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce=self.getOptionalPhaseAttrib("phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce"),
                          numTreeBuildingThreads=self.getOptionalPhaseAttrib("numTreeBuildingThreads"),
                          doPhylogeny=self.getOptionalPhaseAttrib("doPhylogeny", bool, False),
                          minimumBlockHomologySupport=self.getOptionalPhaseAttrib("minimumBlockHomologySupport"),
                          minimumBlockDegreeToCheckSupport=self.getOptionalPhaseAttrib("minimumBlockDegreeToCheckSupport"),
                          phylogenyNucleotideScalingFactor=self.getOptionalPhaseAttrib("phylogenyNucleotideScalingFactor"),
                          removeRecoverableChains=self.getOptionalPhaseAttrib("removeRecoverableChains"),
                          minimumNumberOfSpecies=self.getOptionalPhaseAttrib("minimumNumberOfSpecies", int),
                          phylogenyHomologyUnitType=self.getOptionalPhaseAttrib("phylogenyHomologyUnitType"),
                          phylogenyDistanceCorrectionMethod=self.getOptionalPhaseAttrib("phylogenyDistanceCorrectionMethod"),
                          maxRecoverableChainsIterations=self.getOptionalPhaseAttrib("maxRecoverableChainsIterations", int),
                          maxRecoverableChainLength=self.getOptionalPhaseAttrib("maxRecoverableChainLength", int))
        for message in messages:
            logger.info(message)

    def run(self, fileStore):
        alignments = fileStore.readGlobalFile(self.cactusWorkflowArguments.alignmentsID)
        logger.info("Alignments file: %s" % alignments)

        secondaryAlignments = None
        if self.cactusWorkflowArguments.secondaryAlignmentsID != None:
            secondaryAlignments = fileStore.readGlobalFile(self.cactusWorkflowArguments.secondaryAlignmentsID)

        constraints = None
        if self.cactusWorkflowArguments.constraintsID is not None:
            constraints = fileStore.readGlobalFile(self.cactusWorkflowArguments.constraintsID)


        self.runCactusCafInWorkflow(alignmentFile=alignments, secondaryAlignmentFile=secondaryAlignments, fileStore=fileStore,
                                    constraints=constraints)


############################################################
############################################################
############################################################
#The BAR phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

class CactusBarCheckpoint(CactusCheckpointJob):
    """Load the DB, run the BAR, AVG, and normalization phases, save the DB, then run the reference checkpoint."""
    def run(self, fileStore):
        ktServerDump = self.runPhaseWithPrimaryDB(CactusBarPhase).rv()
        return self.makeFollowOnCheckpointJob(CactusReferenceCheckpoint, "reference", ktServerDump=ktServerDump)

class CactusBarPhase(CactusPhasesJob):
    """Runs bar algorithm."""
    def run(self, fileStore):
        return self.runPhase(CactusBarRecursion, CactusNormalPhase, "normal", doRecursion=self.getOptionalPhaseAttrib("runBar", bool, False))

class CactusBarRecursion(CactusRecursionJob):
    """This job does the get flowers down pass for the BAR alignment phase."""
    memoryPoly = [2e+09]

    def run(self, fileStore):
        self.makeRecursiveJobs(fileStore=fileStore)
        self.makeExtendingJobs(fileStore=fileStore,
                               job=CactusBarWrapper, overlargeJob=CactusBarWrapperLarge)

def runBarForJob(self, fileStore=None, features=None, calculateWhichEndsToComputeSeparately=False, endAlignmentsToPrecomputeOutputFile=None, precomputedAlignments=None):
    return runCactusBar(jobName=self.__class__.__name__,
                 fileStore=fileStore,
                 features=features,
                 cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                 flowerNames=self.flowerNames,
                 maximumLength=self.getOptionalPhaseAttrib("bandingLimit", float),
                 spanningTrees=self.getOptionalPhaseAttrib("spanningTrees", int),
                 gapGamma=self.getOptionalPhaseAttrib( "gapGamma", float),
                 matchGamma=self.getOptionalPhaseAttrib( "matchGamma", float),
                 splitMatrixBiggerThanThis=self.getOptionalPhaseAttrib("splitMatrixBiggerThanThis", int),
                 anchorMatrixBiggerThanThis=self.getOptionalPhaseAttrib("anchorMatrixBiggerThanThis", int),
                 repeatMaskMatrixBiggerThanThis=self.getOptionalPhaseAttrib("repeatMaskMatrixBiggerThanThis", int),
                 diagonalExpansion=self.getOptionalPhaseAttrib("diagonalExpansion"),
                 constraintDiagonalTrim=self.getOptionalPhaseAttrib("constraintDiagonalTrim", int),
                 minimumBlockDegree=self.getOptionalPhaseAttrib("minimumBlockDegree", int),
                 minimumIngroupDegree=self.getOptionalPhaseAttrib("minimumIngroupDegree", int),
                 minimumOutgroupDegree=self.getOptionalPhaseAttrib("minimumOutgroupDegree", int),
                 alignAmbiguityCharacters=self.getOptionalPhaseAttrib("alignAmbiguityCharacters", bool),
                 pruneOutStubAlignments=self.getOptionalPhaseAttrib("pruneOutStubAlignments", bool),
                 useProgressiveMerging=self.getOptionalPhaseAttrib("useProgressiveMerging", bool),
                 calculateWhichEndsToComputeSeparately=calculateWhichEndsToComputeSeparately,
                 endAlignmentsToPrecomputeOutputFile=endAlignmentsToPrecomputeOutputFile,
                 largeEndSize=self.getOptionalPhaseAttrib("largeEndSize", int),
                 precomputedAlignments=precomputedAlignments,
                 ingroupCoverageFile=self.cactusWorkflowArguments.ingroupCoverageID if self.getOptionalPhaseAttrib("rescue", bool) else None,
                 minimumSizeToRescue=self.getOptionalPhaseAttrib("minimumSizeToRescue"),
                 minimumCoverageToRescue=self.getOptionalPhaseAttrib("minimumCoverageToRescue"),
                 minimumNumberOfSpecies=self.getOptionalPhaseAttrib("minimumNumberOfSpecies", int))

class CactusBarWrapper(CactusRecursionJob):
    """Runs the BAR algorithm implementation.
    """
    memoryPoly = [2.81473430e-01, 2.96245523e+09]

    def run(self, fileStore):
        messages = runBarForJob(self, features=self.featuresFn(), fileStore=fileStore)
        for message in messages:
            fileStore.logToMaster(message)

class CactusBarWrapperLarge(CactusRecursionJob):
    """Breaks up the bar into a series of smaller bars, then runs them.
    """
    memoryPoly = [3e+09]

    def run(self, fileStore):
        logger.info("Starting the cactus bar preprocessor job to breakup the bar alignment")
        veryLargeEndSize=self.getOptionalPhaseAttrib("veryLargeEndSize", int, default=1000000)
        maxFlowerGroupSize = self.getOptionalJobAttrib("maxFlowerGroupSize", int,
                                            default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault)
        endsToAlign = []
        endSizes = []
        precomputedAlignmentIDs = []
        for line in runBarForJob(self, features=self.featuresFn(),
                                 fileStore=fileStore, calculateWhichEndsToComputeSeparately=True):
            endToAlign, sequencesInEndAlignment, basesInEndAlignment = line.split()
            sequencesInEndAlignment = int(sequencesInEndAlignment)
            basesInEndAlignment = int(basesInEndAlignment)

            # Sanity check
            assert len(endsToAlign) == len(endSizes), "End alignments and size lists out of sync"

            #If we have a really big end align separately
            if basesInEndAlignment >= veryLargeEndSize:
                alignmentID = self.addChild(CactusBarEndAlignerWrapper(self.phaseNode, self.constantsNode,
                                                        self.cactusDiskDatabaseString, self.flowerNames,
                                                        self.flowerSizes, True, [ endToAlign ], [ basesInEndAlignment ],
                                                        cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
                precomputedAlignmentIDs.append(alignmentID)
                logger.info("Precomputing very large end alignment for %s with %i caps and %i bases" % \
                             (endToAlign, sequencesInEndAlignment, basesInEndAlignment))
            else:
                endsToAlign.append(endToAlign)
                endSizes.append(basesInEndAlignment)
                if sum(endSizes) >= maxFlowerGroupSize:
                    alignmentID = self.addChild(CactusBarEndAlignerWrapper(self.phaseNode,
                                                       self.constantsNode,
                                                       self.cactusDiskDatabaseString,
                                                       self.flowerNames, self.flowerSizes, False,
                                                       endsToAlign, endSizes,
                                                       cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
                    precomputedAlignmentIDs.append(alignmentID)
                    endsToAlign = []
                    endSizes = []
        if len(endsToAlign) > 0:
            precomputedAlignmentIDs.append(self.addChild(CactusBarEndAlignerWrapper(
                self.phaseNode, self.constantsNode, self.cactusDiskDatabaseString, self.flowerNames,
                self.flowerSizes, False, endsToAlign, endSizes,
                cactusWorkflowArguments=self.cactusWorkflowArguments)).rv())
        self.precomputedAlignmentIDs = precomputedAlignmentIDs
        self.makeFollowOnRecursiveJobWithPromisedRequirements(CactusBarWrapperWithPrecomputedEndAlignments)
        logger.info("Breaking bar job into %i separate jobs" % \
                             (len(precomputedAlignmentIDs)))

class CactusBarEndAlignerWrapper(CactusRecursionJob):
    """Computes an end alignment."""
    def featuresFn(self):
        """Merges both end size features and flower features--they will both
        have an impact on resource usage."""
        d = {'endGroupSize': sum(self.endSizes),
             'maxEndSize': max(self.endSizes),
             'numEnds': len(self.endSizes)}
        d.update(self.flowerFeatures())
        return d

    memoryPoly = [1.495e-03, 4.87e+09]
    feature = 'maxEndSize'
    memoryCap = 40e09

    def __init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, flowerSizes,
                 overlarge, endsToAlign, endSizes, cactusWorkflowArguments):
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.endsToAlign = endsToAlign
        self.endSizes = endSizes
        CactusRecursionJob.__init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, flowerSizes, overlarge, cactusWorkflowArguments=self.cactusWorkflowArguments, preemptable=True)

    def run(self, fileStore):
        self.endsToAlign = [ int(i) for i in self.endsToAlign ]
        self.endsToAlign.sort()
        self.flowerNames = encodeFlowerNames((decodeFirstFlowerName(self.flowerNames),) + tuple(self.endsToAlign)) #The ends to align become like extra flower names
        alignmentFile = fileStore.getLocalTempFile()
        messages = runBarForJob(self, features=self.featuresFn(),
                                fileStore=fileStore,
                                endAlignmentsToPrecomputeOutputFile=alignmentFile)
        for message in messages:
            fileStore.logToMaster(message)
        return fileStore.writeGlobalFile(alignmentFile, cleanup=False)

class CactusBarWrapperWithPrecomputedEndAlignments(CactusRecursionJob):
    """Runs the BAR algorithm implementation with some precomputed end alignments."""
    featuresFn = lambda self: {'alignmentsSize': sum([fileID.size for fileID in self.precomputedAlignmentIDs])}
    feature = 'alignmentsSize'
    memoryPoly = [1.99700749e+00, 3.29659639e+08]

    def run(self, fileStore):
        if self.precomputedAlignmentIDs:
            precomputedAlignments = [readGlobalFileWithoutCache(fileStore, fileID) for fileID in self.precomputedAlignmentIDs]
            messages = runBarForJob(self, features=self.featuresFn(),
                                    fileStore=fileStore,
                                    precomputedAlignments=precomputedAlignments)
        else:
            messages = runBarForJob(self)
        for message in messages:
            fileStore.logToMaster(message)

############################################################
############################################################
############################################################
#Normalisation pass
############################################################
############################################################
############################################################

class CactusNormalPhase(CactusPhasesJob):
    """Phase to normalise the graph, ensuring all chains are maximal
    """
    def run(self, fileStore):
        normalisationIterations = self.getOptionalPhaseAttrib("iterations", int, default=0)
        if normalisationIterations > 0:
            self.phaseNode.attrib["normalised"] = "1"
            self.phaseNode.attrib["iterations"] = str(normalisationIterations-1)
            return self.runPhase(CactusNormalRecursion, CactusNormalPhase, "normal")
        else:
            return self.makeFollowOnPhaseJob(CactusAVGPhase, "avg")

class CactusNormalRecursion(CactusRecursionJob):
    """This job does the down pass for the normal phase.
    """
    def run(self, fileStore):
        self.makeRecursiveJobs(fileStore=fileStore)
        return self.makeFollowOnRecursiveJob(CactusNormalRecursion2)

class CactusNormalRecursion2(CactusRecursionJob):
    """This job sets up the normal wrapper in an up traversal of the tree.
    """
    def run(self, fileStore):
        return self.makeWrapperJobs(CactusNormalWrapper)

class CactusNormalWrapper(CactusRecursionJob):
    """This jobs run the normalisation script.
    """
    def run(self, fileStore):
        runCactusMakeNormal(self.cactusDiskDatabaseString, flowerNames=self.flowerNames,
                            maxNumberOfChains=self.getOptionalPhaseAttrib("maxNumberOfChains", int, default=30))

############################################################
############################################################
############################################################
#Phylogeny pass
############################################################
############################################################
############################################################

class CactusAVGPhase(CactusPhasesJob):
    """Phase to build avgs for each flower.
    """
    def run(self, fileStore):
        return self.runPhase(CactusAVGRecursion, SavePrimaryDB, 'avg', doRecursion=self.getOptionalPhaseAttrib("buildAvgs", bool, False))

class CactusAVGRecursion(CactusRecursionJob):
    """This job does the recursive pass for the AVG phase.
    """
    def run(self, fileStore):
        self.makeWrapperJobs(CactusAVGWrapper)
        return self.makeFollowOnRecursiveJob(CactusAVGRecursion2)

class CactusAVGRecursion2(CactusRecursionJob):
    """This job does the recursive pass for the AVG phase.
    """
    def run(self, fileStore):
        self.makeRecursiveJobs(job=CactusAVGRecursion, fileStore=fileStore)

class CactusAVGWrapper(CactusRecursionJob):
    """This job runs tree building
    """
    def run(self, fileStore):
        runCactusPhylogeny(self.cactusDiskDatabaseString, flowerNames=self.flowerNames)

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################


class CactusReferenceCheckpoint(CactusCheckpointJob):
    """Run reference phase, then save the DB and run the HAL phase."""
    def run(self, fileStore):
        child = self.runPhaseWithPrimaryDB(CactusReferencePhase)
        experiment = child.rv(0)
        ktServerDump = child.rv(1)
        self.cactusWorkflowArguments = copy.deepcopy(self.cactusWorkflowArguments)
        self.cactusWorkflowArguments.experimentWrapper = experiment
        return self.makeFollowOnCheckpointJob(CactusHalCheckpoint, "hal", ktServerDump=ktServerDump)

class CactusReferencePhase(CactusPhasesJob):
    """Runs the reference problem algorithm"""
    def run(self, fileStore):
        self.setupSecondaryDatabase()
        self.phaseNode.attrib["experimentPath"] = self.cactusWorkflowArguments.experimentFile
        self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
        exp = self.cactusWorkflowArguments.experimentWrapper
        return self.runPhase(CactusReferenceRecursion, CactusSetReferenceCoordinatesDownPhase, "reference",
                             doRecursion=exp.isRootReconstructed(),
                             launchSecondaryKtForRecursiveJob=True)

class CactusReferenceRecursion(CactusRecursionJob):
    """This job creates the wrappers to run the reference problem algorithm, the follow on job then recurses down.
    """
    memoryPoly = [2e+09]

    def run(self, fileStore):
        self.makeWrapperJobs(CactusReferenceWrapper)
        return self.makeFollowOnRecursiveJob(CactusReferenceRecursion2)

class CactusReferenceWrapper(CactusRecursionJob):
    """Actually run the reference code.
    """
    memoryPoly = [0.71709110685129696, 141266641]
    feature = 'maxFlowerSize'

    def run(self, fileStore):
        exp = self.cactusWorkflowArguments.experimentWrapper
        runCactusReference(fileStore=fileStore,
                       jobName=self.__class__.__name__,
                       features=self.featuresFn(),
                       cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                       flowerNames=self.flowerNames,
                       matchingAlgorithm=self.getOptionalPhaseAttrib("matchingAlgorithm"),
                       permutations=self.getOptionalPhaseAttrib("permutations", int),
                       referenceEventString=exp.getRootGenome(),
                       useSimulatedAnnealing=self.getOptionalPhaseAttrib("useSimulatedAnnealing", bool),
                       theta=self.getOptionalPhaseAttrib("theta", float),
                       phi=self.getOptionalPhaseAttrib("phi", float),
                       maxWalkForCalculatingZ=self.getOptionalPhaseAttrib("maxWalkForCalculatingZ", int),
                       ignoreUnalignedGaps=self.getOptionalPhaseAttrib("ignoreUnalignedGaps", bool),
                       wiggle=self.getOptionalPhaseAttrib("wiggle", float),
                       numberOfNs=self.getOptionalPhaseAttrib("numberOfNs", int),
                       minNumberOfSequencesToSupportAdjacency=self.getOptionalPhaseAttrib("minNumberOfSequencesToSupportAdjacency", int),
                       makeScaffolds=self.getOptionalPhaseAttrib("makeScaffolds", bool))

class CactusReferenceRecursion2(CactusRecursionJob):
    memoryPoly = [2e+09]
    def run(self, fileStore):
        self.makeRecursiveJobs(fileStore=fileStore,
                               job=CactusReferenceRecursion)
        return self.makeFollowOnRecursiveJob(CactusReferenceRecursion3)

class CactusReferenceRecursion3(CactusRecursionJob):
    """After completing the recursion for the reference algorithm, the up pass of adding in the reference coordinates is performed.
    """
    def run(self, fileStore):
        return self.makeWrapperJobs(CactusSetReferenceCoordinatesUpWrapper)

class CactusSetReferenceCoordinatesUpWrapper(CactusRecursionJob):
    """Does the up pass for filling in the reference sequence coordinates, once a reference has been established.
    """
    memoryPoly = [1.3030742924744299, 180741939.947]
    feature = 'maxFlowerSize'

    def run(self, fileStore):
        exp = self.cactusWorkflowArguments.experimentWrapper
        runCactusAddReferenceCoordinates(fileStore=fileStore, jobName=self.__class__.__name__,
                                         features=self.featuresFn(),
                                         cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                         secondaryDatabaseString=self.getOptionalPhaseAttrib("secondaryDatabaseString"),
                                         flowerNames=self.flowerNames,
                                         referenceEventString=exp.getRootGenome(),
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"),
                                         bottomUpPhase=True)

class CactusSetReferenceCoordinatesDownPhase(CactusPhasesJob):
    """This is the second part of the reference coordinate setting, the down pass.
    """
    def run(self, fileStore):
        self.cleanupSecondaryDatabase()
        exp = self.cactusWorkflowArguments.experimentWrapper
        return self.runPhase(CactusSetReferenceCoordinatesDownRecursion, CactusExtractReferencePhase, "reference", doRecursion=exp.isRootReconstructed())

class CactusSetReferenceCoordinatesDownRecursion(CactusRecursionJob):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """
    memoryPoly = [2e+09]

    def run(self, fileStore):
        self.makeWrapperJobs(CactusSetReferenceCoordinatesDownWrapper)
        return self.makeFollowOnRecursiveJob(CactusSetReferenceCoordinatesDownRecursion2)

class CactusSetReferenceCoordinatesDownRecursion2(CactusRecursionJob):
    memoryPoly = [2e+09]

    def run(self, fileStore):
        return self.makeRecursiveJobs(fileStore=fileStore,
                                      job=CactusSetReferenceCoordinatesDownRecursion)

class CactusSetReferenceCoordinatesDownWrapper(CactusRecursionJob):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """
    memoryPoly = [0.52844015396914878, 116287385]
    feature = 'maxFlowerSize'

    def run(self, fileStore):
        exp = self.cactusWorkflowArguments.experimentWrapper
        runCactusAddReferenceCoordinates(fileStore=fileStore, features=self.featuresFn(),
                                         jobName=self.__class__.__name__,
                                         cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                         flowerNames=self.flowerNames,
                                         referenceEventString=exp.getRootGenome(),
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"),
                                         bottomUpPhase=False)

class CactusExtractReferencePhase(CactusPhasesJob):
    memoryPoly = [2.24519561e+00, 4.70479486e+08]

    def run(self, fileStore):
        experiment = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
        if experiment.isRootReconstructed():
            fileStore.logToMaster("Starting Reference Extract Phase")
            eventName = experiment.getRootGenome()
            referencePath = fileStore.getLocalTempFile()
            cactus_call(parameters=["cactus_getReferenceSeq", "--cactusDisk",
                                    self.cactusWorkflowArguments.cactusDiskDatabaseString, "--flowerName", "0",
                                    "--referenceEventString", eventName, "--outputFile",
                                    os.path.basename(referencePath), "--logLevel", getLogLevelString()])
            referenceID = fileStore.writeGlobalFile(referencePath)
            experiment.setReferenceID(referenceID)
            intermediateResultsUrl = getattr(self.cactusWorkflowArguments, 'intermediateResultsUrl', None)
            if intermediateResultsUrl is not None:
                # The user requested to keep the hal fasta files in a separate place. Export it there.
                url = intermediateResultsUrl + ".reference.fa"
                fileStore.exportFile(referenceID, url)
        self.cactusWorkflowArguments.experimentWrapper = experiment
        return experiment, self.makeFollowOnPhaseJob(CactusCheckPhase, "check")

############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################

class CactusCheckPhase(CactusPhasesJob):
    """The check phase, where we verify everything is as it should be
    """
    def run(self, fileStore):
        normalNode = findRequiredNode(self.cactusWorkflowArguments.configNode, "normal")
        self.phaseNode.attrib["checkNormalised"] = getOptionalAttrib(normalNode, "normalised", default="0")
        return self.runPhase(CactusCheckRecursion, SavePrimaryDB, "check", doRecursion=self.getOptionalPhaseAttrib("runCheck", bool, False))

class CactusCheckRecursion(CactusRecursionJob):
    """This job does the recursive pass for the check phase.
    """
    def run(self, fileStore):
        self.makeRecursiveJobs(fileStore=fileStore)
        self.makeWrapperJobs(CactusCheckWrapper)

class CactusCheckWrapper(CactusRecursionJob):
    """Runs the actual check wrapper
    """
    def run(self, fileStore):
        runCactusCheck(self.cactusDiskDatabaseString, self.flowerNames, checkNormalised=self.getOptionalPhaseAttrib("checkNormalised", bool, False))

############################################################
############################################################
############################################################
#Hal generation
############################################################
############################################################
############################################################

class CactusHalCheckpoint(CactusCheckpointJob):
    """Load the DB and run the final reference and HAL phases."""
    def run(self, fileStore):
        return self.runPhaseWithPrimaryDB(CactusHalGeneratorPhase).rv()

class CactusHalGeneratorPhase(CactusPhasesJob):
    def run(self, fileStore):
        if self.getOptionalPhaseAttrib("buildFasta", bool, default=False):
            self.fastaID = self.makeRecursiveChildJob(CactusFastaGenerator)
        return self.makeFollowOnPhaseJob(CactusHalGeneratorPhase2, "hal")

class CactusHalGeneratorPhase2(CactusPhasesJob):
    def run(self, fileStore):
        if self.getOptionalPhaseAttrib("buildHal", bool, default=False):
            self.setupSecondaryDatabase()
            self.phaseNode.attrib["experimentPath"] = self.cactusWorkflowArguments.experimentFile
            self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
            self.phaseNode.attrib["outputFile"] = "1"

            self.halID = self.makeRecursiveChildJob(CactusHalGeneratorRecursion, launchSecondaryKtForRecursiveJob=True)

        return self.makeFollowOnPhaseJob(CactusHalGeneratorPhase3, "hal")

class CactusHalGeneratorPhase3(CactusPhasesJob):
    def run(self, fileStore):
        self.cactusWorkflowArguments.experimentWrapper.setHalID(self.halID)
        self.cactusWorkflowArguments.experimentWrapper.setHalFastaID(self.fastaID)
        return self.cactusWorkflowArguments.experimentWrapper

class CactusFastaGenerator(CactusRecursionJob):
    memoryPoly = [2.99160856e+00, 4.48507512e+08]
    feature = 'totalSequenceSize'

    def run(self, fileStore):
        experiment = self.cactusWorkflowArguments.experimentWrapper
        tmpFasta = fileStore.getLocalTempFile()
        runCactusFastaGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                flowerName=decodeFirstFlowerName(self.flowerNames),
                                outputFile=tmpFasta,
                                referenceEventString=experiment.getRootGenome())
        intermediateResultsUrl = getattr(self.cactusWorkflowArguments, 'intermediateResultsUrl', None)
        fastaID = fileStore.writeGlobalFile(tmpFasta)
        if intermediateResultsUrl is not None:
            # The user requested to keep the hal fasta files in a separate place. Export it there.
            url = intermediateResultsUrl + ".hal.fa"
            fileStore.exportFile(fastaID, url)
        return fastaID

class CactusHalGeneratorRecursion(CactusRecursionJob):
    """Generate the hal file by merging indexed hal files from the children.
    """
    memoryPoly = [2e+09]
    def run(self, fileStore):
        i = extractNode(self.phaseNode)
        if "outputFile" in i.attrib:
            i.attrib.pop("outputFile")

        self.makeRecursiveJobs(fileStore=fileStore, phaseNode=i)
        return self.makeFollowOnRecursiveJob(CactusHalGeneratorUpWrapper)

class CactusHalGeneratorUpWrapper(CactusRecursionJob):
    """Generate the .c2h strings for this flower, storing them in the secondary database."""
    memoryPoly = [4e+09]

    def run(self, fileStore):
        if self.getOptionalPhaseAttrib("outputFile"):
            tmpHal = fileStore.getLocalTempFile()
        else:
            tmpHal = None
        runCactusHalGenerator(jobName=self.__class__.__name__,
                              features=self.featuresFn(),
                              fileStore=fileStore,
                              cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                              secondaryDatabaseString=self.getOptionalPhaseAttrib("secondaryDatabaseString"),
                              flowerNames=self.flowerNames,
                              referenceEventString=self.cactusWorkflowArguments.experimentWrapper.getRootGenome(),
                              outputFile=tmpHal,
                              showOnlySubstitutionsWithRespectToReference=\
                              self.getOptionalPhaseAttrib("showOnlySubstitutionsWithRespectToReference", bool))
        if tmpHal:
            # At top level--have the final .c2h file
            intermediateResultsUrl = getattr(self.cactusWorkflowArguments, 'intermediateResultsUrl', None)
            halID = fileStore.writeGlobalFile(tmpHal)
            if intermediateResultsUrl is not None:
                # The user requested to keep the c2h files in a separate place. Export it there.
                url = intermediateResultsUrl + ".c2h"
                fileStore.exportFile(halID, url)
            return halID
        else:
            return None

class CactusHalGeneratorPhaseCleanup(CactusPhasesJob):
    """Cleanup the database used to build the hal
    """
    def run(self, fileStore):
        self.cleanupSecondaryDatabase()

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
    def __init__(self, options, experimentFile, configNode, seqIDMap):
        #Get a local copy of the experiment file
        self.experimentFile = experimentFile
        self.experimentNode = ET.parse(self.experimentFile).getroot()
        self.scratchDbElemNode = ET.parse(self.experimentFile).getroot()
        self.experimentWrapper = ExperimentWrapper(self.experimentNode)
        self.alignmentsID = None
        for genome, seqID in list(seqIDMap.items()):
            print(('setting this', genome, seqID))
            self.experimentWrapper.setSequenceID(genome, seqID)
        #Get the database string
        self.cactusDiskDatabaseString = ET.tostring(self.experimentNode.find("cactus_disk").find("st_kv_database_conf"), encoding='unicode').replace('\n', '')
        #Get the species tree
        self.speciesTree = self.experimentNode.attrib["species_tree"]
        #Get any list of 'required species' for the blocks of the cactus.
        self.outgroupEventNames = getOptionalAttrib(self.experimentNode, "outgroup_events")
        #Constraints
        self.constraintsID = getOptionalAttrib(self.experimentNode, "constraintsID")
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
        self.ktServerDump = None

        #Secondary, scratch DB
        secondaryConf = copy.deepcopy(self.experimentNode.find("cactus_disk").find("st_kv_database_conf"))
        secondaryElem = DbElemWrapper(secondaryConf)
        self.secondaryDatabaseString = secondaryElem.getConfString()

        #The config node
        self.configNode = configNode
        self.configWrapper = ConfigWrapper(self.configNode)
        #Now deal with the constants that need to be added here
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()
        self.configWrapper.setBuildHal(options.buildHal)
        self.configWrapper.setBuildFasta(options.buildFasta)

        #Now build the remaining options from the arguments
        findRequiredNode(self.configNode, "avg").attrib["buildAvgs"] = "1" if options.buildAvgs else "0"

def addCactusWorkflowOptions(parser):
    parser.add_argument("--experiment", dest="experimentFile",
                      help="The file containing a link to the experiment parameters")
    parser.add_argument("--buildAvgs", dest="buildAvgs", action="store_true",
                      help="Build trees", default=False)
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
        eW = self.cactusWorkflowArguments.experimentWrapper
        genomes = eW.getGenomesWithSequence()
        originalSequenceIDs = [eW.getSequenceID(genome) for genome in genomes]
        preprocessedSeqIDs = self.addChild(CactusPreprocessor(originalSequenceIDs, self.cactusWorkflowArguments.configNode)).rv()
        self.addFollowOn(AfterPreprocessing(self.cactusWorkflowArguments, eW, genomes, preprocessedSeqIDs))

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
        self.addFollowOn(CactusTrimmingBlastPhase(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName="trimBlast"))

def runCactusWorkflow(args):
    ##########################################
    #Construct the arguments.
    ##########################################

    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    options = parser.parse_args(args)
    options.disableCaching = True
    setLoggingFromOptions(options)

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
