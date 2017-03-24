#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework to structure all the related wrappers.
"""

import os
import sys
import xml.etree.ElementTree as ET
import math
import time
import random
import copy
from argparse import ArgumentParser

from sonLib.bioio import newickTreeParser

from toil.lib.bioio import getTempFile
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions
from toil.lib.bioio import system
from sonLib.bioio import catFiles
from sonLib.bioio import getLogLevelString


from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.common import RunAsFollowOn

from toil.job import Job
from toil.common import Toil


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

from cactus.blast.cactus_blast import BlastIngroupsAndOutgroups
from cactus.blast.cactus_blast import BlastOptions

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.ktserverToil import KtServerService
from cactus.pipeline.ktserverControl import dumpKtServer, clearKtServer, restoreKtServer

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
    def __init__(self, phaseNode, constantsNode, overlarge=False, cactusSequencesID=None,
                 checkpoint=False, preemptable=True):
        self.phaseNode = phaseNode
        self.constantsNode = constantsNode
        self.overlarge = overlarge
        self.jobNode = getJobNode(self.phaseNode, self.__class__)
        self.cactusSequencesID = cactusSequencesID
        if self.jobNode:
            logger.info("JobNode = %s" % self.jobNode.attrib)

        memory = None
        cores = None
        if hasattr(self, 'memoryPoly'):
            # Memory should be determined by a polynomial fit on the
            # input size
            memory = self.evaluateResourcePoly(self.memoryPoly)

        if cactusSequencesID and hasattr(cactusSequencesID, 'size'):
            disk = 2*cactusSequencesID.size
        else:
            disk = None
        if memory is None and overlarge:
            memory = self.getOptionalJobAttrib("overlargeMemory", typeFn=int,
                                               default=getOptionalAttrib(self.constantsNode, "defaultOverlargeMemory", int, default=sys.maxint))
            cores = self.getOptionalJobAttrib("overlargeCpu", typeFn=int,
                                              default=getOptionalAttrib(self.constantsNode, "defaultOverlargeCpu", int, default=None))
        elif memory is None:
            memory = self.getOptionalJobAttrib("memory", typeFn=int,
                                               default=getOptionalAttrib(self.constantsNode, "defaultMemory", int, default=sys.maxint))
            cores = self.getOptionalJobAttrib("cpu", typeFn=int,
                                              default=getOptionalAttrib(self.constantsNode, "defaultCpu", int, default=sys.maxint))
        RoundedJob.__init__(self, memory=memory, cores=cores, disk=disk,
                            checkpoint=checkpoint, preemptable=preemptable)

    def evaluateResourcePoly(self, poly):
        """Evaluate a polynomial based on the total sequence size."""
        if hasattr(self, 'inputSizeFn'):
            size = self.inputSizeFn()
        else:
            size = self.cactusWorkflowArguments.totalSequenceSize
        resource = 0
        for degree, coefficient in enumerate(reversed(poly)):
            resource += coefficient * (size**degree)
        return int(resource)

    def getOptionalPhaseAttrib(self, attribName, typeFn=None, default=None):
        """Gets an optional attribute of the phase node.
        """
        return getOptionalAttrib(node=self.phaseNode, attribName=attribName, typeFn=typeFn, default=default)

    def getOptionalJobAttrib(self, attribName, typeFn=None, default=None):
        """Gets an optional attribute of the job node.
        """
        return getOptionalAttrib(node=self.jobNode, attribName=attribName, typeFn=typeFn, default=default)

class CactusPhasesJob(CactusJob):
    """Base job for each workflow phase job.
    """
    def __init__(self, cactusWorkflowArguments, phaseName, topFlowerName=0, index=0,
                 cactusSequencesID=None, checkpoint=False, preemptable=True, halID=None,
                 fastaID=None, ktServerDump=None):
        phaseNode = findRequiredNode(cactusWorkflowArguments.configNode, phaseName, index)
        constantsNode = findRequiredNode(cactusWorkflowArguments.configNode, "constants")
        self.ktServerDump = ktServerDump
        self.index = index
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.topFlowerName = topFlowerName
        self.halID = halID
        self.fastaID = fastaID
        CactusJob.__init__(self, phaseNode=phaseNode, constantsNode=constantsNode, overlarge=False,
                           cactusSequencesID=cactusSequencesID, checkpoint=checkpoint,
                           preemptable=preemptable)

    def makeRecursiveChildJob(self, job, launchSecondaryKtForRecursiveJob=False):
        newChild = job(phaseNode=extractNode(self.phaseNode), 
                       constantsNode=extractNode(self.constantsNode),
                       cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString, 
                       flowerNames=encodeFlowerNames((self.topFlowerName,)),
                       flowerSizes=[self.cactusWorkflowArguments.totalSequenceSize],
                       overlarge=True,
                       cactusSequencesID=self.cactusSequencesID,
                       cactusWorkflowArguments=self.cactusWorkflowArguments)

        if launchSecondaryKtForRecursiveJob and ExperimentWrapper(self.cactusWorkflowArguments.experimentNode).getDbType() == "kyoto_tycoon":
            cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)
            memory = self.evaluateResourcePoly([4.10201882, 2.01324291e+08])
            cpu = cw.getKtserverCpu(default=getOptionalAttrib(
                    self.constantsNode, "defaultCpu", int, default=0))
            dbElem = ExperimentWrapper(self.cactusWorkflowArguments.scratchDbElemNode)
            dbString = self.addService(KtServerService(dbElem = dbElem, isSecondary = True, memory=memory, cores=cpu))
            newChild.phaseNode.attrib["secondaryDatabaseString"] = dbString
            return self.addChild(newChild).rv()
        else:
            return self.addChild(newChild).rv()

    def makeFollowOnPhaseJob(self, job, phaseName, index=0, ktServerDump = None, checkpoint = False):
        return self.addFollowOn(job(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName=phaseName, 
                                    topFlowerName=self.topFlowerName, index=index, cactusSequencesID=self.cactusSequencesID, checkpoint=checkpoint, halID=self.halID, fastaID=self.fastaID, ktServerDump=ktServerDump)).rv()

    def runPhase(self, recursiveJob, nextPhaseJob, nextPhaseName, doRecursion=True, index=0, launchSecondaryKtForRecursiveJob=False, updateDatabase=False):
        """
        Adds a recursive child job and then a follow-on phase job. Returns the result of the follow-on
        phase job.
        """
        logger.info("Starting %s phase job with index %i at %s seconds (recursing = %i)" % (self.phaseNode.tag, self.getPhaseIndex(), time.time(), doRecursion))
        if doRecursion:
            cactusSequencesID = self.makeRecursiveChildJob(recursiveJob, launchSecondaryKtForRecursiveJob)
        if updateDatabase and doRecursion:
            self.cactusSequencesID = cactusSequencesID
        return self.makeFollowOnPhaseJob(job=nextPhaseJob, phaseName=nextPhaseName, index=index)

    def getPhaseIndex(self):
        return self.index

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

class CactusRecursionJob(CactusJob):
    """Base recursive job for traversals up and down the cactus tree.
    """
    maxSequenceSizeOfFlowerGroupingDefault = 1000000
    def __init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, flowerSizes, overlarge=False, precomputedAlignmentIDs=None, cactusSequencesID = None, checkpoint = False, cactusWorkflowArguments=None, preemptable=True, memPoly=None):
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.flowerNames = flowerNames
        self.flowerSizes = flowerSizes
        self.cactusWorkflowArguments = cactusWorkflowArguments

        #need to do this because the alignment IDs are jobstore promises, and can't 
        #be stored in the config XML until they are respolved into actual IDs, which doesn't
        #happen until the follow-on job after CactusBarWrapperLarge
        self.precomputedAlignmentIDs = precomputedAlignmentIDs

        CactusJob.__init__(self, phaseNode=phaseNode, constantsNode=constantsNode, overlarge=overlarge, 
                           cactusSequencesID=cactusSequencesID, checkpoint=checkpoint, preemptable=preemptable)
        
    def makeFollowOnRecursiveJob(self, job, phaseNode=None):
        """Sets the followon to the given recursive job
        """
        if phaseNode == None:
            phaseNode = self.phaseNode
        return self.addFollowOn(RunAsFollowOn(job, phaseNode=phaseNode, constantsNode=self.constantsNode,
                                              cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                              flowerNames=self.flowerNames, flowerSizes=self.flowerSizes,
                                              overlarge=self.overlarge,
                                              precomputedAlignmentIDs=self.precomputedAlignmentIDs,
                                              cactusSequencesID=self.cactusSequencesID,
                                              cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
        
    def makeChildJobs(self, flowersAndSizes, job, overlargeJob=None, 
                      phaseNode=None, runFlowerStats=False):
        """Make a set of child jobs for a given set of flowers and chosen child job
        """
        cactusSequencesID = None
        if overlargeJob == None:
            overlargeJob = job
        if phaseNode == None:
            phaseNode = self.phaseNode
        
        logger.info("Make wrapper jobs: There are %i flowers" % len(flowersAndSizes))
        for overlarge, flowerNames, flowerSizes in flowersAndSizes:
            if overlarge: #Make sure large flowers are on their own, in their own job
                if runFlowerStats:
                    flowerStatsString = runCactusFlowerStats(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerName=decodeFirstFlowerName(flowerNames))
                    logger.info("Adding an oversize flower for job class %s and stats %s" \
                                             % (overlargeJob, flowerStatsString))
                else:
                    logger.info("Adding an oversize flower %s for job class %s" \
                                             % (decodeFirstFlowerName(flowerNames), overlargeJob))
                cactusSequencesID = self.addChild(overlargeJob(cactusDiskDatabaseString=
                                                               self.cactusDiskDatabaseString,
                                                               phaseNode=phaseNode, 
                                                               constantsNode=self.constantsNode,
                                                               flowerNames=flowerNames,
                                                               flowerSizes=flowerSizes,
                                                               overlarge=True,
                                                               cactusSequencesID=self.cactusSequencesID,
                                                               cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
            else:
                logger.info("Adding recursive flower job")
                cactusSequencesID = self.addChild(job(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                                      phaseNode=phaseNode, constantsNode=self.constantsNode,
                                                      flowerNames=flowerNames,
                                                      flowerSizes=flowerSizes,
                                                      overlarge=False,
                                                      cactusSequencesID=self.cactusSequencesID,
                                                      cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
        if cactusSequencesID:
            return cactusSequencesID
        else:
            return self.cactusSequencesID

    def makeRecursiveJobs(self, fileStore=None, job=None, phaseNode=None, runFlowerStats=False):
        """Make a set of child jobs for a given set of parent flowers.
        """
        if job == None:
            job = self.__class__
        jobNode = getJobNode(self.phaseNode, job)
        logger.info("Flower names: %s" % self.flowerNames)
        flowersAndSizes=runCactusGetFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                            cactusSequencesPath=self.cactusSequencesPath, 
                                            inputSize=sum(self.flowerSizes),
                                            jobName=job.__name__,
                                            fileStore=fileStore,
                                            flowerNames=self.flowerNames, 
                                            minSequenceSizeOfFlower=getOptionalAttrib(jobNode, "minFlowerSize", int, 0), 
                                            maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(jobNode, "maxFlowerGroupSize", int, 
                                            default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault),
                                            maxSequenceSizeOfSecondaryFlowerGrouping=getOptionalAttrib(jobNode, "maxFlowerWrapperGroupSize", int, 
                                            default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault))
        return self.makeChildJobs(flowersAndSizes=flowersAndSizes, 
                              job=job, phaseNode=phaseNode,
                              runFlowerStats=runFlowerStats)
    
    def makeExtendingJobs(self, job, fileStore=None, overlargeJob=None, phaseNode=None,
                          runFlowerStats=False):
        """Make set of child jobs that extend the current cactus tree.
        """
        jobNode = getJobNode(self.phaseNode, job)
        flowersAndSizes=runCactusExtendFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                              cactusSequencesPath=self.cactusSequencesPath, 
                                              inputSize=sum(self.flowerSizes),
                                              jobName=job.__name__,
                                              fileStore=fileStore,
                                              flowerNames=self.flowerNames, 
                                              minSequenceSizeOfFlower=getOptionalAttrib(jobNode, "minFlowerSize", int, 1), 
                                              maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(jobNode, "maxFlowerGroupSize", int, 
                                              default=CactusRecursionJob.maxSequenceSizeOfFlowerGroupingDefault))
        return self.makeChildJobs(flowersAndSizes=flowersAndSizes, 
                              job=job, overlargeJob=overlargeJob, 
                              phaseNode=phaseNode,
                              runFlowerStats=runFlowerStats)

    
    def makeWrapperJobs(self, job, overlargeJob=None, phaseNode=None, runFlowerStats=False):
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
                                  phaseNode=phaseNode, runFlowerStats=runFlowerStats)

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
                tokens = line[1:].split()
                tokens[0] = "id=%d|%s" % (uniqueID, tokens[0])
                out.write(">%s\n" % "".join(tokens))
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

        # Get ingroup and outgroup sequences
        sequenceIDs = self.cactusWorkflowArguments.experimentWrapper.seqIDMap.values()
        sequences = [fileStore.readGlobalFile(seqID) for seqID in sequenceIDs]
        self.cactusWorkflowArguments.totalSequenceSize = sum(os.stat(x).st_size for x in sequences)

        # Prepend unique ID to fasta headers to prevent name collision
        renamedInputSeqDir = os.path.join(fileStore.getLocalTempDir(), "renamedInputs")
        os.mkdir(renamedInputSeqDir)
        uniqueFas = prependUniqueIDs(sequences, renamedInputSeqDir)
        uniqueFaIDs = [fileStore.writeGlobalFile(seq) for seq in uniqueFas]
        
        self.cactusWorkflowArguments.experimentWrapper.seqIDMap = dict(zip(self.cactusWorkflowArguments.experimentWrapper.seqIDMap.keys(), uniqueFaIDs))
        ingroupIDs = map(lambda x: x[1], filter(lambda x: x[0] not in self.cactusWorkflowArguments.experimentWrapper.getOutgroupEvents(), self.cactusWorkflowArguments.experimentWrapper.seqIDMap.items()))
        outgroupIDs = [self.cactusWorkflowArguments.experimentWrapper.seqIDMap[i] for i in self.cactusWorkflowArguments.experimentWrapper.getOutgroupEvents()]
        fileStore.logToMaster("Ingroup sequences: %s" % (ingroupIDs))
        fileStore.logToMaster("Outgroup sequences: %s" % (outgroupIDs))

        # Change the blast arguments depending on the divergence
        setupDivergenceArgs(self.cactusWorkflowArguments)
        setupFilteringByIdentity(self.cactusWorkflowArguments)

        # FIXME: this is really ugly and steals the options from the caf tag
        blastJob = self.addChild(BlastIngroupsAndOutgroups(
                                          BlastOptions(chunkSize=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "chunkSize", int),
                                                        overlapSize=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "overlapSize", int),
                                                        lastzArguments=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "lastzArguments"),
                                                        compressFiles=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "compressFiles", bool),
                                                        realign=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "realign", bool), 
                                                        realignArguments=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "realignArguments"),
                                                        memory=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "lastzMemory", int, sys.maxint),
                                                        smallDisk=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "lastzSmallDisk", int, sys.maxint),
                                                        largeDisk=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "lastzLargeDisk", int, sys.maxint),
                                                        minimumSequenceLength=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "minimumSequenceLengthForBlast", int, 1),
                                                       trimFlanking=self.getOptionalPhaseAttrib("trimFlanking", int, 10),
                                                       trimMinSize=self.getOptionalPhaseAttrib("trimMinSize", int, 0),
                                                       trimThreshold=self.getOptionalPhaseAttrib("trimThreshold", float, 0.8),
                                                       trimWindowSize=self.getOptionalPhaseAttrib("trimWindowSize", int, 10),
                                                       trimOutgroupFlanking=self.getOptionalPhaseAttrib("trimOutgroupFlanking", int, 100),
                                                       trimOutgroupDepth=self.getOptionalPhaseAttrib("trimOutgroupDepth", int, 1),
                                                       keepParalogs=self.getOptionalPhaseAttrib("keepParalogs", bool, False)), ingroupIDs, outgroupIDs))

        self.cactusWorkflowArguments.alignmentsID = blastJob.rv(0)
        self.cactusWorkflowArguments.outgroupFragmentIDs = blastJob.rv(1)
        self.cactusWorkflowArguments.ingroupCoverageIDs = blastJob.rv(2)
        
        return self.makeFollowOnPhaseJob(CactusSetupPhase, "setup")
        
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

class CactusSetupPhase(CactusPhasesJob):
    """Initialises the cactus database and adapts the config file for the run.
    """
    def run(self, fileStore):
        # Point the outgroup sequences to their trimmed versions for
        # phases after this one.
        for i, outgroup in enumerate(self.cactusWorkflowArguments.experimentWrapper.getOutgroupEvents()):
            self.cactusWorkflowArguments.experimentWrapper.seqIDMap[outgroup] = self.cactusWorkflowArguments.outgroupFragmentIDs[i]

        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupDivergenceArgs(self.cactusWorkflowArguments)

        setupJob = CactusSetupPhase2(cactusWorkflowArguments=self.cactusWorkflowArguments,
                                     phaseName='setup', topFlowerName=self.topFlowerName,
                                     index=0)
        loadDBJob = LoadDBPhase(setupJob, ktServerDump=self.ktServerDump,
                                cactusWorkflowArguments=self.cactusWorkflowArguments,
                                phaseName='setup', topFlowerName=self.topFlowerName,
                                index=0, checkpoint=True)
        promise = self.addChild(loadDBJob)
        self.cactusSequencesID, ktServerDump = promise.rv(0), promise.rv(1)
        return self.makeFollowOnPhaseJob(CactusBarPhase, "bar", ktServerDump=ktServerDump)

class LoadDBPhase(CactusPhasesJob):
    def __init__(self, nextJob, *args, **kwargs):
        self.nextJob = nextJob
        super(LoadDBPhase, self).__init__(*args, **kwargs)

    def run(self, fileStore):
        cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)

        if self.cactusWorkflowArguments.experimentWrapper.getDbType() == "kyoto_tycoon":
            memory = self.evaluateResourcePoly([4.10201882, 2.01324291e+08])
            cores = cw.getKtserverCpu(default=getOptionalAttrib(
                    self.constantsNode, "defaultCpu", int, default=sys.maxint))
            dbElem = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
            dbString = self.addService(KtServerService(dbElem=dbElem, isSecondary=False,
                                                       memory=memory, cores=cores))
            self.nextJob.cactusWorkflowArguments.cactusDiskDatabaseString = dbString
            return self.addChild(LoadDBPhase2(self.nextJob, ktServerDump=self.ktServerDump,
                                              cactusWorkflowArguments=self.nextJob.cactusWorkflowArguments,
                                              phaseName='setup', topFlowerName=self.topFlowerName,
                                              index=0)).rv()
        else:
            return self.addFollowOn(self.nextJob).rv()

class LoadDBPhase2(CactusPhasesJob):
    def __init__(self, nextJob, *args, **kwargs):
        self.nextJob = nextJob
        super(LoadDBPhase2, self).__init__(*args, **kwargs)

    def run(self, fileStore):
        if self.ktServerDump is not None:
            dbElem = DbElemWrapper(ET.fromstring(self.cactusWorkflowArguments.cactusDiskDatabaseString))
            fileStore.logToMaster("Restoring cactus DB from %s" % self.ktServerDump)
            restoreKtServer(dbElem, fileStore.readGlobalFile(self.ktServerDump))
        else:
            fileStore.logToMaster("Starting clean cactus DB")
        return self.addFollowOn(self.nextJob).rv()

class CactusSetupPhase2(CactusPhasesJob):   
    def run(self, fileStore):        
        sequenceIDs = []
        tree = self.cactusWorkflowArguments.experimentWrapper.getTree()
        sequenceNames = []
        firstLines = []
        for node in tree.postOrderTraversal():
            if tree.isLeaf(node):
                seqID = self.cactusWorkflowArguments.experimentWrapper.seqIDMap[tree.getName(node)]
                sequenceIDs.append(seqID)
                seq = fileStore.readGlobalFile(seqID)
                sequenceNames.append(tree.getName(node))
                with open(seq, 'r') as fh:
                    firstLines.append(fh.readline())

        sequences = [fileStore.readGlobalFile(fileID) for fileID in sequenceIDs]
        logger.info("Sequences in cactus setup: %s" % sequenceNames)
        logger.info("Sequences in cactus setup filenames: %s" % firstLines)
        cactusSequencesPath = fileStore.getLocalTempFile()
        messages = runCactusSetup(cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString, 
                       cactusSequencesPath = cactusSequencesPath,
                       sequences=sequences,
                       newickTreeString=self.cactusWorkflowArguments.speciesTree, 
                       outgroupEvents=self.cactusWorkflowArguments.outgroupEventNames,
                       makeEventHeadersAlphaNumeric=self.getOptionalPhaseAttrib("makeEventHeadersAlphaNumeric", bool, False))
        self.cactusSequencesID = fileStore.writeGlobalFile(cactusSequencesPath)
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

    def run(self, fileStore):
        if self.cactusWorkflowArguments.ingroupCoverageIDs is not None:
            # Convert the bed files to use 64-bit cactus Names instead
            # of the headers. Ideally this should belong in the bar
            # phase but we run stripUniqueIDs before then.
            bedFiles = [fileStore.readGlobalFile(path) for path in self.cactusWorkflowArguments.ingroupCoverageIDs]
            tempFile = fileStore.getLocalTempFile()
            system("cat %s > %s" % (" ".join(bedFiles), tempFile))
            cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
            ingroupCoverageFile = fileStore.getLocalTempFile()
            runConvertAlignmentsToInternalNames(self.cactusWorkflowArguments.cactusDiskDatabaseString, cactusSequencesPath, tempFile, ingroupCoverageFile, self.topFlowerName, isBedFile=True)
            self.cactusWorkflowArguments.ingroupCoverageID = fileStore.writeGlobalFile(ingroupCoverageFile)

        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupFilteringByIdentity(self.cactusWorkflowArguments)
        #Setup any constraints
        if self.getPhaseIndex() == 0 and self.cactusWorkflowArguments.constraintsID != None: #Setup the constraints arg

            constraintsFile = fileStore.readGlobalFile(self.cactusWorkflowArguments.constraintsID)
            newConstraintsFile = fileStore.getLocalTempFile()
            runCactusConvertAlignmentToCactus(self.cactusWorkflowArguments.cactusDiskDatabaseString,
                                              constraintsFile, newConstraintsFile)
            self.phaseNode.attrib["constraintsID"] = fileStore.writeGlobalFile(newConstraintsFile, cleanup=False)

        assert self.getPhaseNumber() == 1
        alignmentsFile = fileStore.readGlobalFile(self.cactusWorkflowArguments.alignmentsID)
        convertedAlignmentsFile = fileStore.getLocalTempFile()
        # Convert the cigar file to use 64-bit cactus Names instead of the headers.
        cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        runConvertAlignmentsToInternalNames(cactusDiskString=self.cactusWorkflowArguments.cactusDiskDatabaseString, cactusSequencesPath=cactusSequencesPath, alignmentsFile=alignmentsFile, outputFile=convertedAlignmentsFile, flowerName=self.topFlowerName)
        fileStore.logToMaster("Converted headers of cigar file %s to internal names, new file %s" % (self.cactusWorkflowArguments.alignmentsID, convertedAlignmentsFile))
        self.cactusWorkflowArguments.alignmentsID = fileStore.writeGlobalFile(convertedAlignmentsFile, cleanup=False)
        # While we're at it, remove the unique IDs prepended to
        # the headers inside the cactus DB.
        runStripUniqueIDs(self.cactusWorkflowArguments.cactusDiskDatabaseString, cactusSequencesPath)
        self.phaseNode.attrib["alignmentsID"] = self.cactusWorkflowArguments.alignmentsID
        return self.runPhase(CactusCafWrapper, CactusSaveDBPhase, "caf")

class CactusCafWrapper(CactusRecursionJob):
    """Runs cactus_caf on one flower and one alignment file.
    """
    memoryPoly = [1.80395944e+01, 7.96042247e+07]

    def __init__(self, **kwargs):
        # We want to make sure this job is *not* preempted, because it's
        # the longest-running job in the entire process, taking up to
        # 12 hours.
        kwargs['preemptable'] = False
        super(CactusCafWrapper, self).__init__(**kwargs)

    def runCactusCafInWorkflow(self, alignmentFile, constraints=None):
        debugFilePath = self.getOptionalPhaseAttrib("phylogenyDebugPrefix")
        if debugFilePath != None:
            debugFilePath += getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "reference"), "reference")
        messages = runCactusCaf(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                          cactusSequencesPath = self.cactusSequencesPath,
                          alignments=alignmentFile, 
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
                          referenceEventHeader=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "reference"), "reference"),
                          phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce=self.getOptionalPhaseAttrib("phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce"),
                          numTreeBuildingThreads=self.getOptionalPhaseAttrib("numTreeBuildingThreads"),
                          doPhylogeny=self.getOptionalPhaseAttrib("doPhylogeny", bool, True),
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
        alignments = None
        assert "alignmentsID" in self.phaseNode.attrib
        if "alignmentsID" in self.phaseNode.attrib:
            alignments = fileStore.readGlobalFile(self.phaseNode.attrib["alignmentsID"])
        constraints = None
        if "constraintsID" in self.phaseNode.attrib:
            constraints = fileStore.readGlobalFile(self.phaseNode.attrib["constraintsID"])
        logger.info("Alignments file: %s" % alignments)
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        self.runCactusCafInWorkflow(alignmentFile=alignments, constraints=constraints)


class CactusSaveDBPhase(CactusPhasesJob):
    """Saves the DB to a file and clears the DB."""
    def run(self, fileStore):
        tempPath = fileStore.getLocalTempFile()
        dbElem = DbElemWrapper(ET.fromstring(self.cactusWorkflowArguments.cactusDiskDatabaseString))
        dumpKtServer(dbElem, tempPath)
        clearKtServer(dbElem)
        return self.cactusSequencesID, fileStore.writeGlobalFile(tempPath)

############################################################
############################################################
############################################################
#The BAR phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

class CactusBarPhase(CactusPhasesJob):
    def run(self, fileStore):
        barJob = CactusBarPhase2(cactusWorkflowArguments=self.cactusWorkflowArguments,
                                 phaseName='bar', topFlowerName=self.topFlowerName,
                                 index=0, cactusSequencesID=self.cactusSequencesID)
        loadDBJob = LoadDBPhase(barJob, ktServerDump=self.ktServerDump,
                                cactusWorkflowArguments=self.cactusWorkflowArguments,
                                phaseName='bar', topFlowerName=self.topFlowerName,
                                index=0, checkpoint=True)
        return self.addChild(loadDBJob).rv()

class CactusBarPhase2(CactusPhasesJob):
    """Runs bar algorithm."""
    def run(self, fileStore):
        assert self.cactusSequencesID
        logger.info("DatabaseID in BarPhase: %s" % self.cactusSequencesID)
        return self.runPhase(CactusBarRecursion, CactusNormalPhase, "normal", doRecursion=self.getOptionalPhaseAttrib("runBar", bool, False))

class CactusBarRecursion(CactusRecursionJob):
    """This job does the get flowers down pass for the BAR alignment phase."""
    inputSizeFn = lambda self: sum(self.flowerSizes)
    memoryPoly = [1.47507363e-01, 1.34319212e+09]

    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        self.makeRecursiveJobs(fileStore=fileStore)
        self.makeExtendingJobs(fileStore=fileStore,
                               job=CactusBarWrapper, overlargeJob=CactusBarWrapperLarge, runFlowerStats=True)

def runBarForJob(self, fileStore=None, inputSize=None, calculateWhichEndsToComputeSeparately=None, endAlignmentsToPrecomputeOutputFile=None, precomputedAlignments=None):
    return runCactusBar(jobName=self.__class__.__name__,
                 fileStore=fileStore,
                 inputSize=inputSize,
                 cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                 cactusSequencesPath=self.cactusSequencesPath,
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
    inputSizeFn = lambda self: sum(self.flowerSizes)
    memoryPoly = [2.81473430e-01, 2.96245523e+09]

    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        messages = runBarForJob(self, inputSize=sum(self.flowerSizes), fileStore=fileStore)
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
        totalSize = 0
        precomputedAlignmentIDs = []
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        for line in runBarForJob(self, inputSize=sum(self.flowerSizes),
                                 fileStore=fileStore, calculateWhichEndsToComputeSeparately=True):
            endToAlign, sequencesInEndAlignment, basesInEndAlignment = line.split()
            sequencesInEndAlignment = int(sequencesInEndAlignment)
            basesInEndAlignment = int(basesInEndAlignment)
            #If we have a really big end align separately
            if basesInEndAlignment >= veryLargeEndSize:
                alignmentID = self.addChild(CactusBarEndAlignerWrapper(self.phaseNode, self.constantsNode,
                                                        self.cactusDiskDatabaseString, self.flowerNames,
                                                        True, [ endToAlign ], self.cactusSequencesID,
                                                        cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
                precomputedAlignmentIDs.append(alignmentID)
                logger.info("Precomputing very large end alignment for %s with %i caps and %i bases" % \
                             (endToAlign, sequencesInEndAlignment, basesInEndAlignment))
            else:
                endsToAlign.append(endToAlign)
                totalSize += basesInEndAlignment
                if totalSize >= maxFlowerGroupSize:
                    alignmentID = self.addChild(CactusBarEndAlignerWrapper(self.phaseNode,
                                                       self.constantsNode,
                                                       self.cactusDiskDatabaseString,
                                                       self.flowerNames, False,
                                                       endsToAlign,
                                                       self.cactusSequencesID,
                                                       cactusWorkflowArguments=self.cactusWorkflowArguments)).rv()
                    precomputedAlignmentIDs.append(alignmentID)
                    endsToAlign = []
                    totalSize = 0
        if len(endsToAlign) > 0:
            precomputedAlignmentIDs.append(self.addChild(CactusBarEndAlignerWrapper(self.phaseNode, self.constantsNode, self.cactusDiskDatabaseString, self.flowerNames, False, endsToAlign, self.cactusSequencesID,
                                                cactusWorkflowArguments=self.cactusWorkflowArguments)).rv())
        self.precomputedAlignmentIDs = precomputedAlignmentIDs
        self.makeFollowOnRecursiveJob(CactusBarWrapperWithPrecomputedEndAlignments)
        logger.info("Breaking bar job into %i separate jobs" % \
                             (len(precomputedAlignmentIDs)))

class CactusBarEndAlignerWrapper(CactusRecursionJob):
    """Computes an end alignment."""
    memoryPoly = [3.82834548e-01, 1.28869656e+09]
    def __init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, overlarge,
                 endsToAlign, cactusSequencesID, cactusWorkflowArguments):
        self.cactusWorkflowArguments = cactusWorkflowArguments
        CactusRecursionJob.__init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, overlarge, cactusSequencesID=cactusSequencesID, cactusWorkflowArguments=self.cactusWorkflowArguments, preemptable=True)
        self.endsToAlign = endsToAlign

    def run(self, fileStore):
        self.endsToAlign = [ int(i) for i in self.endsToAlign ]
        self.endsToAlign.sort()
        self.flowerNames = encodeFlowerNames((decodeFirstFlowerName(self.flowerNames),) + tuple(self.endsToAlign)) #The ends to align become like extra flower names
        alignmentFile = fileStore.getLocalTempFile()
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        messages = runBarForJob(self,
                                endAlignmentsToPrecomputeOutputFile=alignmentFile)
        for message in messages:
            fileStore.logToMaster(message)
        return fileStore.writeGlobalFile(alignmentFile, cleanup=False)

class CactusBarWrapperWithPrecomputedEndAlignments(CactusRecursionJob):
    """Runs the BAR algorithm implementation with some precomputed end alignments."""
    inputSizeFn = lambda self: sum([fileID.size for fileID in self.precomputedAlignmentIDs])
    memoryPoly = [1.99700749e+00, 3.29659639e+08]

    def run(self, fileStore):
        if self.precomputedAlignmentIDs:
            precomputedAlignments = [readGlobalFileWithoutCache(fileStore, fileID) for fileID in self.precomputedAlignmentIDs]
            self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
            messages = runBarForJob(self, inputSize=self.inputSizeFn(),
                                    fileStore=fileStore,
                                    precomputedAlignments=" ".join(precomputedAlignments))
        else:
            self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
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
        assert self.cactusSequencesID
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
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
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
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        runCactusMakeNormal(self.cactusDiskDatabaseString, self.cactusSequencesPath, flowerNames=self.flowerNames, 
                            maxNumberOfChains=self.getOptionalPhaseAttrib("maxNumberOfChains", int, default=30))
        return self.cactusSequencesID

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
        return self.runPhase(CactusAVGRecursion, CactusReferencePhase, "reference", doRecursion=self.getOptionalPhaseAttrib("buildAvgs", bool, False))

class CactusAVGRecursion(CactusRecursionJob):
    """This job does the recursive pass for the AVG phase.
    """
    def run(self, fileStore):
        self.cactusSequencesID = self.makeWrapperJobs(CactusAVGWrapper)
        return self.makeFollowOnRecursiveJob(CactusAVGRecursion2)
        
class CactusAVGRecursion2(CactusRecursionJob):
    """This job does the recursive pass for the AVG phase.
    """
    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        self.makeRecursiveJobs(job=CactusAVGRecursion, fileStore=fileStore)

class CactusAVGWrapper(CactusRecursionJob):
    """This job runs tree building
    """
    def run(self, fileStore):
        assert self.cactusSequencesID
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        runCactusPhylogeny(self.cactusDiskDatabaseString, self.cactusSequencesPath, flowerNames=self.flowerNames)
        return self.cactusSequencesID

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################
    
class CactusReferencePhase(CactusPhasesJob):
    def run(self, fileStore):
        """Runs the reference problem algorithm
        """
        assert self.cactusSequencesID
        self.setupSecondaryDatabase()
        self.phaseNode.attrib["experimentPath"] = self.cactusWorkflowArguments.experimentFile
        self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
        return self.runPhase(CactusReferenceRecursion, CactusSetReferenceCoordinatesDownPhase, "reference", 
                      doRecursion=self.getOptionalPhaseAttrib("buildReference", bool, False),
                      launchSecondaryKtForRecursiveJob = True, updateDatabase=True)

class CactusReferenceRecursion(CactusRecursionJob):
    """This job creates the wrappers to run the reference problem algorithm, the follow on job then recurses down.
    """
    inputSizeFn = lambda self: sum(self.flowerSizes)
    memoryPoly = [6.15004877e+09]

    def run(self, fileStore):
        assert self.cactusSequencesID
        logger.info("DatabaseID in RefRecursion: %s" % self.cactusSequencesID)
        self.cactusSequencesID = self.makeWrapperJobs(CactusReferenceWrapper, runFlowerStats=True)
        return self.makeFollowOnRecursiveJob(CactusReferenceRecursion2)
        
class CactusReferenceWrapper(CactusRecursionJob):
    """Actually run the reference code.
    """
    inputSizeFn = lambda self: sum(self.flowerSizes)
    memoryPoly = [5.67663286e-01, 6.31948477e+09]

    def run(self, fileStore):
        assert self.cactusSequencesID
        logger.info("Reading database in RefWrapper: %s" % self.cactusSequencesID)
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID, mutable=True)
        runCactusReference(fileStore=fileStore,
                       jobName=self.__class__.__name__,
                       inputSize=sum(self.flowerSizes),
                       cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                       cactusSequencesPath = self.cactusSequencesPath,
                       flowerNames=self.flowerNames, 
                       matchingAlgorithm=self.getOptionalPhaseAttrib("matchingAlgorithm"), 
                       permutations=self.getOptionalPhaseAttrib("permutations", int),
                       referenceEventString=self.getOptionalPhaseAttrib("reference"),
                       useSimulatedAnnealing=self.getOptionalPhaseAttrib("useSimulatedAnnealing", bool),
                       theta=self.getOptionalPhaseAttrib("theta", float),
                       phi=self.getOptionalPhaseAttrib("phi", float),
                       maxWalkForCalculatingZ=self.getOptionalPhaseAttrib("maxWalkForCalculatingZ", int),
                       ignoreUnalignedGaps=self.getOptionalPhaseAttrib("ignoreUnalignedGaps", bool),
                       wiggle=self.getOptionalPhaseAttrib("wiggle", float),
                       numberOfNs=self.getOptionalPhaseAttrib("numberOfNs", int),
                       minNumberOfSequencesToSupportAdjacency=self.getOptionalPhaseAttrib("minNumberOfSequencesToSupportAdjacency", int),
                       makeScaffolds=self.getOptionalPhaseAttrib("makeScaffolds", bool))
        return fileStore.writeGlobalFile(self.cactusSequencesPath)

class CactusReferenceRecursion2(CactusRecursionJob):
    memoryPoly = [9.69305877e-01, 1.67009148e+09]
    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        self.cactusSequencesID = self.makeRecursiveJobs(fileStore=fileStore,
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
    memoryPoly = [1.36245902e+00, 1.36316562e+09]

    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID, mutable=True)
        assert os.path.exists(self.cactusSequencesPath)
        runCactusAddReferenceCoordinates(fileStore=fileStore, jobName=self.__class__.__name__,
                                         inputSize=sum(self.flowerSizes),
                                         cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                         cactusSequencesPath = self.cactusSequencesPath,
                                         secondaryDatabaseString=self.getOptionalPhaseAttrib("secondaryDatabaseString"),
                                         flowerNames=self.flowerNames,
                                         referenceEventString=self.getOptionalPhaseAttrib("reference"),
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"),
                                         bottomUpPhase=True)
        return fileStore.writeGlobalFile(self.cactusSequencesPath)
        
class CactusSetReferenceCoordinatesDownPhase(CactusPhasesJob):
    """This is the second part of the reference coordinate setting, the down pass.
    """
    def run(self, fileStore):
        self.cleanupSecondaryDatabase()
        return self.runPhase(CactusSetReferenceCoordinatesDownRecursion, CactusExtractReferencePhase, "check", doRecursion=self.getOptionalPhaseAttrib("buildReference", bool, False), updateDatabase=True)
        
class CactusSetReferenceCoordinatesDownRecursion(CactusRecursionJob):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """
    memoryPoly = [8.53350165e-01, 1.08914552e+09]

    def run(self, fileStore):
        assert self.cactusSequencesID
        self.cactusSequencesID = self.makeWrapperJobs(CactusSetReferenceCoordinatesDownWrapper)
        return self.makeFollowOnRecursiveJob(CactusSetReferenceCoordinatesDownRecursion2)

class CactusSetReferenceCoordinatesDownRecursion2(CactusRecursionJob):
    memoryPoly = [7.31066727e-01, 1.70083157e+09]

    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        return self.makeRecursiveJobs(fileStore=fileStore,
                                      job=CactusSetReferenceCoordinatesDownRecursion)
        
class CactusSetReferenceCoordinatesDownWrapper(CactusRecursionJob):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """
    memoryPoly = [9.65818693e-01, 1.45719681e+09]

    def run(self, fileStore):
        assert self.cactusSequencesID
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID, mutable=True)
        runCactusAddReferenceCoordinates(fileStore=fileStore, inputSize=sum(self.flowerSizes),
                                         jobName=self.__class__.__name__,
                                         cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                         cactusSequencesPath = self.cactusSequencesPath,
                                         flowerNames=self.flowerNames,
                                         referenceEventString=self.getOptionalPhaseAttrib("reference"),
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"), 
                                         bottomUpPhase=False)
        return fileStore.writeGlobalFile(self.cactusSequencesPath)

class CactusExtractReferencePhase(CactusPhasesJob):
    memoryPoly = [2.24519561e+00, 4.70479486e+08]

    def run(self, fileStore):
        assert self.cactusSequencesID
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        experiment = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
        if hasattr(self.cactusWorkflowArguments, 'buildReference') and\
               self.cactusWorkflowArguments.buildReference:
            fileStore.logToMaster("Starting Reference Extract Phase")
            if experiment.getReferencePath() is not None:
                eventName = os.path.basename(experiment.getReferencePath())
                if eventName.find('.') >= 0:
                    eventName = eventName[:eventName.rfind('.')]
                    referencePath = fileStore.getLocalTempFile()
                    cactus_call(parameters=["cactus_getReferenceSeq"],
                                option_string="--cactusDisk '%s' --cactusSequencesPath '%s' --flowerName 0 --referenceEventString %s --outputFile %s --logLevel %s" % 
                              (self.cactusWorkflowArguments.cactusDiskDatabaseString, os.path.basename(self.cactusSequencesPath), 
                                      eventName, os.path.basename(referencePath), getLogLevelString()))
                    experiment.setReferenceID(fileStore.writeGlobalFile(referencePath))
        self.cactusWorkflowArguments.experimentWrapper = experiment
        return self.makeFollowOnPhaseJob(CactusCheckPhase, "check")

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
        return self.runPhase(CactusCheckRecursion, CactusHalGeneratorPhase, "hal", doRecursion=self.getOptionalPhaseAttrib("runCheck", bool, False))
        
class CactusCheckRecursion(CactusRecursionJob):
    """This job does the recursive pass for the check phase.
    """
    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        self.makeRecursiveJobs(fileStore=fileStore)
        self.makeWrapperJobs(CactusCheckWrapper)
        
class CactusCheckWrapper(CactusRecursionJob):
    """Runs the actual check wrapper
    """
    def run(self, fileStore):
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        runCactusCheck(self.cactusDiskDatabaseString, self.cactusSequencesPath, self.flowerNames, checkNormalised=self.getOptionalPhaseAttrib("checkNormalised", bool, False))
        return self.cactusSequencesID

############################################################
############################################################
############################################################
#Hal generation
############################################################
############################################################
############################################################

class CactusHalGeneratorPhase(CactusPhasesJob):
    def run(self, fileStore):
        referenceNode = findRequiredNode(self.cactusWorkflowArguments.configNode, "reference")
        if referenceNode.attrib.has_key("reference"):
            self.phaseNode.attrib["reference"] = referenceNode.attrib["reference"]
        if self.getOptionalPhaseAttrib("buildFasta", bool, default=False):
            self.phaseNode.attrib["fastaPath"] = self.cactusWorkflowArguments.experimentNode.find("hal").attrib["fastaPath"]
            self.fastaID = self.makeRecursiveChildJob(CactusFastaGenerator)
        return self.makeFollowOnPhaseJob(CactusHalGeneratorPhase2, "hal")

class CactusHalGeneratorPhase2(CactusPhasesJob):
    def run(self, fileStore):
        if self.getOptionalPhaseAttrib("buildHal", bool, default=False):
            self.setupSecondaryDatabase()
            self.phaseNode.attrib["experimentPath"] = self.cactusWorkflowArguments.experimentFile
            self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
            self.phaseNode.attrib["outputFile"]=self.cactusWorkflowArguments.experimentNode.find("hal").attrib["halPath"]

            self.halID = self.makeRecursiveChildJob(CactusHalGeneratorRecursion, launchSecondaryKtForRecursiveJob=True)

        return self.makeFollowOnPhaseJob(CactusHalGeneratorPhase3, "hal")

class CactusHalGeneratorPhase3(CactusPhasesJob):
    def run(self, fileStore):
        self.cactusWorkflowArguments.experimentWrapper.setHalID(self.halID)
        self.cactusWorkflowArguments.experimentWrapper.setHalFastaID(self.fastaID)
        return self.cactusWorkflowArguments.experimentWrapper

class CactusFastaGenerator(CactusRecursionJob):
    memoryPoly = [2.99160856e+00, 4.48507512e+08]

    def run(self, fileStore):
        assert self.cactusSequencesID
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        tmpFasta = fileStore.getLocalTempFile()
        runCactusFastaGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                    cactusSequencesPath=self.cactusSequencesPath,
                                    flowerName=decodeFirstFlowerName(self.flowerNames),
                                    outputFile=tmpFasta,
                                    referenceEventString=self.getOptionalPhaseAttrib("reference"))
        return fileStore.writeGlobalFile(tmpFasta)
            
class CactusHalGeneratorRecursion(CactusRecursionJob):
    """Generate the hal file by merging indexed hal files from the children.
    """
    memoryPoly = [3.19890954e-01, 1.59265705e+09]
    def run(self, fileStore):
        i = extractNode(self.phaseNode)
        if "outputFile" in i.attrib:
            i.attrib.pop("outputFile")

        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        self.makeRecursiveJobs(fileStore=fileStore, phaseNode=i)
        return self.makeFollowOnRecursiveJob(CactusHalGeneratorUpWrapper)

class CactusHalGeneratorUpWrapper(CactusRecursionJob):
    """Generate the .c2h strings for this flower, storing them in the secondary database."""
    memoryPoly = [3.32474798e+00, 1.14948623e+09]

    def run(self, fileStore):
        if self.getOptionalPhaseAttrib("outputFile"):
            tmpHal = fileStore.getLocalTempFile()
        else:
            tmpHal = None
        logger.info("DatabaseID in HalGeneratorUpWrapper: %s" % self.cactusSequencesID)
        self.cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        runCactusHalGenerator(jobName=self.__class__.__name__,
                              inputSize=sum(self.flowerSizes),
                              fileStore=fileStore,
                              cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                              cactusSequencesPath=self.cactusSequencesPath,
                              secondaryDatabaseString=self.getOptionalPhaseAttrib("secondaryDatabaseString"),
                              flowerNames=self.flowerNames,
                              referenceEventString=self.getOptionalPhaseAttrib("reference"),
                              outputFile=tmpHal,
                              showOnlySubstitutionsWithRespectToReference=\
                              self.getOptionalPhaseAttrib("showOnlySubstitutionsWithRespectToReference", bool))
        if tmpHal:
            return fileStore.writeGlobalFile(tmpHal)
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
        self.experimentWrapper.seqIDMap = seqIDMap
        #Get the database string
        self.cactusDiskDatabaseString = ET.tostring(self.experimentNode.find("cactus_disk").find("st_kv_database_conf")).translate(None, '\n')
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
        self.ktServerDump = None

        #Secondary, scratch DB
        secondaryConf = copy.deepcopy(self.experimentNode.find("cactus_disk").find("st_kv_database_conf"))
        secondaryElem = DbElemWrapper(secondaryConf)
        dbPath = secondaryElem.getDbDir()
        assert dbPath is not None
        secondaryDbPath = os.path.join(os.path.dirname(dbPath), "%s_tempSecondaryDatabaseDir_%s" % (
            os.path.basename(dbPath), random.random()))
        secondaryElem.setDbDir(secondaryDbPath)
        if secondaryElem.getDbType() == "kyoto_tycoon":
            secondaryElem.setDbPort(secondaryElem.getDbPort() + 100)
        self.secondaryDatabaseString = secondaryElem.getConfString()
            
        #The config node
        self.configNode = configNode
        self.configWrapper = ConfigWrapper(self.configNode)
        #Now deal with the constants that ned to be added here
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()
        self.configWrapper.setBuildHal(options.buildHal)
        self.configWrapper.setBuildFasta(options.buildFasta)
        
        #Now build the remaining options from the arguments
        if options.buildAvgs:
            findRequiredNode(self.configNode, "avg").attrib["buildAvgs"] = "1"
        if options.buildReference:
            findRequiredNode(self.configNode, "reference").attrib["buildReference"] = "1"
            

def addCactusWorkflowOptions(parser):
    parser.add_argument("--experiment", dest="experimentFile", 
                      help="The file containing a link to the experiment parameters")
    
    parser.add_argument("--buildAvgs", dest="buildAvgs", action="store_true",
                      help="Build trees", default=False)
    
    parser.add_argument("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the flowers", default=False)
    
    parser.add_argument("--buildHal", dest="buildHal", action="store_true",
                      help="Build a hal file", default=False)
    
    parser.add_argument("--buildFasta", dest="buildFasta", action="store_true",
                      help="Build a fasta file of the input sequences (and reference sequence, used with hal output)", 
                      default=False)

    parser.add_argument("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")

class RunCactusPreprocessorThenCactusSetup(RoundedJob):
    def __init__(self, options, cactusWorkflowArguments):
        RoundedJob.__init__(self)
        self.options = options
        self.cactusWorkflowArguments = cactusWorkflowArguments
        
    def run(self, fileStore):
        eW = self.cactusWorkflowArguments.experimentWrapper
        seqIDs = self.addChild(CactusPreprocessor(eW.seqIDMap.values(), self.cactusWorkflowArguments.configNode))
        #Now make the setup, replacing the input sequences with the preprocessed sequences
        eW.seqIDMap = dict(zip(eW.seqIDMap.keys(), [seqIDs.rv(i) for i in range(len(eW.seqIDMap))]))
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
    if options.test:
        _test()
    setLoggingFromOptions(options)
    
    #if len(args) != 0:
     #   raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    experimentWrapper = ExperimentWrapper(ET.parse(options.experimentFile).getroot())
    with Toil(options) as toil:
        seqIDMap = dict()
        seqMap = experimentWrapper.buildSequenceMap()
        for name in seqMap:
            fullSeq = getTempFile()
            if os.path.isdir(seqMap[name]):
                catFiles([os.path.join(seqMap[name], seqFile) for seqFile in os.listdir(seqMap[name])], fullSeq)
            else:
                fullSeq = seqMap[name]
            seqIDMap[name] = toil.importFile(makeURL(fullSeq))


        configNode = ET.parse(experimentWrapper.getConfigPath()).getroot()
        cactusWorkflowArguments = CactusWorkflowArguments(options, experimentFile = options.experimentFile, configNode=configNode, seqIDMap = seqIDMap)

        toil.start(RunCactusPreprocessorThenCactusSetup(options, cactusWorkflowArguments))

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    runCactusWorkflow(sys.argv)
