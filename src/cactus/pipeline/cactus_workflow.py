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
from sonLib.bioio import catFiles
from sonLib.bioio import getLogLevelString

from toil.job import Job
from toil.common import Toil
from toil.realtimeLogger import RealtimeLogger

from cactus.shared.common import makeURL
from cactus.shared.common import RunAsFollowOn
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import runCactusConsolidated
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import runCactusSplitFlowersBySecondaryGrouping
from cactus.shared.common import encodeFlowerNames
from cactus.shared.common import decodeFirstFlowerName
from cactus.shared.common import runCactusFlowerStats
from cactus.shared.common import runCactusSecondaryDatabase
from cactus.shared.common import findRequiredNode
from cactus.shared.common import RoundedJob

from cactus.blast.blast import BlastIngroupsAndOutgroups
from cactus.blast.blast import BlastOptions
from cactus.blast.mappingQualityRescoringAndFiltering import mappingQualityRescoring

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.dbServerToil import DbServerService, getDbServer

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

    def makeRecursiveChildJob(self, job, launchSecondaryDbForRecursiveJob=False):
        newChild = job(phaseNode=extractNode(self.phaseNode),
                       constantsNode=extractNode(self.constantsNode),
                       cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString,
                       flowerNames=encodeFlowerNames((self.topFlowerName,)),
                       flowerSizes=[self.cactusWorkflowArguments.totalSequenceSize],
                       overlarge=True,
                       cactusWorkflowArguments=self.cactusWorkflowArguments)

        if launchSecondaryDbForRecursiveJob and ExperimentWrapper(self.cactusWorkflowArguments.experimentNode).getDbType() in ["kyoto_tycoon", "redis"]:
            cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)
            memory = max(2500000000, self.evaluateResourcePoly([4.10201882, 2.01324291e+08]))
            cpu = cw.getKtserverCpu(default=0.1)
            dbElem = ExperimentWrapper(self.cactusWorkflowArguments.scratchDbElemNode)
            dbString = self.addService(DbServerService(dbElem=dbElem, isSecondary=True, memory=memory, cores=cpu)).rv(0)
            newChild.phaseNode.attrib["secondaryDatabaseString"] = dbString
            return self.addChild(newChild).rv()
        else:
            return self.addChild(newChild).rv()

    def makeFollowOnPhaseJob(self, job, phaseName):
        return self.addFollowOn(job(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName=phaseName,
                                    topFlowerName=self.topFlowerName, halID=self.halID, fastaID=self.fastaID)).rv()

    def runPhase(self, recursiveJob, nextPhaseJob, nextPhaseName, doRecursion=True, launchSecondaryDbForRecursiveJob=False):
        """
        Adds a recursive child job and then a follow-on phase job. Returns the result of the follow-on
        phase job.
        """
        RealtimeLogger.info("Starting %s phase job at %s seconds (recursing = %i)" % (self.phaseNode.tag, time.time(), doRecursion))
        if doRecursion:
            self.makeRecursiveChildJob(recursiveJob, launchSecondaryDbForRecursiveJob)
        return self.makeFollowOnPhaseJob(job=nextPhaseJob, phaseName=nextPhaseName)

    def makeFollowOnCheckpointJob(self, checkpointConstructor, phaseName, dbServerDump=None):
        """Add a follow-on checkpoint phase."""
        return self.addFollowOn(checkpointConstructor(\
                   phaseName=phaseName, dbServerDump=dbServerDump,
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
        # TODO: It seems that it does nothing for the known databases. Maybe we can remove it.
        if dbElem.getDbType() not in ["kyoto_tycoon", "redis"]:
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=True)

    def cleanupSecondaryDatabase(self):
        """Cleanup the secondary database
        """
        confXML = ET.fromstring(self.cactusWorkflowArguments.secondaryDatabaseString)
        dbElem = DbElemWrapper(confXML)
        # TODO: It seems that it does nothing for the known databases. Maybe we can remove it.
        if dbElem.getDbType() not in ["kyoto_tycoon", "redis"]:
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=False)

class CactusCheckpointJob(CactusPhasesJob):
    """Special "checkpoint" phase job, launching and restoring the primary Cactus DB.

    Inherit from this and run `runPhaseWithPrimaryDB` to start a new
    primary DB before running the given phase.

    Meant to provide a restore point in case of pipeline failure. Note
    that the checkpointed job is technically this job's child, which
    starts the database.
    """
    def __init__(self, dbServerDump=None, *args, **kwargs):
        self.dbServerDump = dbServerDump
        super(CactusCheckpointJob, self).__init__(*args, **kwargs)

    def runPhaseWithPrimaryDB(self, jobConstructor):
        """Start and load a new primary DB before running the given phase.
        """
        job = jobConstructor(cactusWorkflowArguments=self.cactusWorkflowArguments,
                             phaseName=self.phaseName, topFlowerName=self.topFlowerName)
        startDBJob = StartPrimaryDB(job, dbServerDump=self.dbServerDump,
                                    cactusWorkflowArguments=self.cactusWorkflowArguments,
                                    phaseName=self.phaseName, topFlowerName=self.topFlowerName)
        promise = self.addChild(startDBJob)
        return promise

class StartPrimaryDB(CactusPhasesJob):
    """Launches a primary Cactus DB."""
    def __init__(self, nextJob, dbServerDump=None, *args, **kwargs):
        self.nextJob = nextJob
        self.dbServerDump = dbServerDump
        kwargs['checkpoint'] = True
        kwargs['preemptable'] = False
        super(StartPrimaryDB, self).__init__(*args, **kwargs)

    def run(self, fileStore):
        cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)

        if self.cactusWorkflowArguments.experimentWrapper.getDbType() in ["kyoto_tycoon", "redis"]:
            memory = max(2500000000, self.evaluateResourcePoly([4.10201882, 2.01324291e+08]))
            cores = cw.getKtserverCpu(default=0.1)
            dbElem = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
            service = self.addService(DbServerService(dbElem=dbElem,
                                                      existingSnapshotID=self.dbServerDump,
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
        getDbServer(dbElem, fileStore).stopServer()
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

def prependUniqueIDs(fas, outputDir, idMap=None, firstID=0):
    """Prepend unique ints to fasta headers.

    (prepend rather than append since trimmed outgroups have a start
    token appended, which complicates removal slightly)
    """
    uniqueID = firstID
    ret = []
    for fa in fas:
        outPath = os.path.join(outputDir, os.path.basename(fa))
        out = open(outPath, 'w')
        for line in open(fa):
            if len(line) > 0 and line[0] == '>':
                header = "id=%d|%s" % (uniqueID, line[1:-1])
                out.write(">%s\n" % header)
                if idMap is not None:
                    idMap[line[1:-1].rstrip()] = header.rstrip()
            else:
                out.write(line)
        ret.append(outPath)
        uniqueID += 1
    return ret

def getLongestPath(node, distance=0.0):
    """Identify the longest path from the mrca of the leaves of the species tree.
    """
    i, j = distance, distance
    if node.left != None:
        i = getLongestPath(node.left, abs(node.left.distance)) + distance
    if node.right != None:
        j = getLongestPath(node.right, abs(node.right.distance)) + distance
    return max(i, j)

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
    def __init__(self, standAlone = False, *args, **kwargs):
        self.standAlone = standAlone
        super(CactusTrimmingBlastPhase, self).__init__(*args, **kwargs)
        
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
                         keepParalogs=self.getOptionalPhaseAttrib("keepParalogs", bool, False),
                         gpuLastz=getOptionalAttrib(cafNode, "gpuLastz", bool, False)),
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

        # hack to get cactus-blast standalone tool working
        if self.standAlone:
            return self.cactusWorkflowArguments

        return self.makeFollowOnPhaseJob(CactusConsolidated1, phaseName="consolidated")

        #return self.makeFollowOnCheckpointJob(CactusSetupCheckpoint, "setup")

def updateExpWrapperForOutgroups(job, expWrapper, outgroupGenomes, outgroupFragmentIDs):
    for genome, outgroupFragmentID in zip(outgroupGenomes, outgroupFragmentIDs):
        expWrapper.setSequenceID(genome, outgroupFragmentID)
    return expWrapper

############################################################
############################################################
############################################################
##The optional consolidate phase, which runs the setup, caf,
## bar, reference and cactus to hal algorithms in one job
## on a multi-node machine
############################################################
############################################################
############################################################

class CactusConsolidated1(CactusPhasesJob):
    """Start the secondary DB."""
    def run(self, fileStore):
        # Get the experiment object
        experiment = self.cactusWorkflowArguments.experimentWrapper
        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupDivergenceArgs(self.cactusWorkflowArguments)

        # Build up a genome -> fasta map.
        seqIDMap = dict((genome, experiment.getSequenceID(genome)) for genome in experiment.getGenomesWithSequence())
        seqMap = dict((genome, fileStore.readGlobalFile(id)) for genome, id in list(seqIDMap.items()))

        # Get the alignment files
        alignments = fileStore.readGlobalFile(self.cactusWorkflowArguments.alignmentsID)
        logger.info("Alignments file: %s" % alignments)

        secondaryAlignments = None
        if self.cactusWorkflowArguments.secondaryAlignmentsID != None:
            secondaryAlignments = fileStore.readGlobalFile(self.cactusWorkflowArguments.secondaryAlignmentsID)

        constraints = None
        if self.cactusWorkflowArguments.constraintsID is not None:
            constraints = fileStore.readGlobalFile(self.cactusWorkflowArguments.constraintsID)

        # Temporary place to store the output c2h file
        tmpHal = fileStore.getLocalTempFile()
        tmpFasta = fileStore.getLocalTempFile()
        tmpRef = fileStore.getLocalTempFile()

        messages = runCactusConsolidated(cactusParams=experiment.getConfigPath(),
                                         seqMap=seqMap,
                                         newickTreeString=self.cactusWorkflowArguments.speciesTree,
                                         alignmentsFile=alignments,
                                         #secondaryDatabaseString=self.cactusWorkflowArguments.secondaryDatabaseString,
                                         outputFile=tmpHal,
                                         outputHalFastaFile=tmpFasta,
                                         outputReferenceFile=tmpRef,
                                         secondaryAlignmentsFile=secondaryAlignments,
                                         constraintAlignmentsFile=constraints,
                                         logLevel=None,
                                         outgroupEvents=experiment.getOutgroupGenomes(),
                                         referenceEvent=experiment.getRootGenome())

        # Log back any messages
        for message in messages:
            logger.info(message)

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

        self.cactusWorkflowArguments.experimentWrapper.setHalID(halID)
        self.cactusWorkflowArguments.experimentWrapper.setHalFastaID(fastaID)
        self.cactusWorkflowArguments.experimentWrapper.setReferenceID(referenceID)

        return self.cactusWorkflowArguments.experimentWrapper

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
        self.dbServerDump = None

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
