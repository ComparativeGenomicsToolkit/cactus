#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus package.
"""

import os
import cPickle
import pickle
import sys
import shutil
import subprocess32
import logging
import uuid
import json
import time
import signal

from toil.lib.bioio import logger
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString

from toil.job import Job

from sonLib.bioio import popenCatch

from cactus.shared.version import cactus_commit

_log = logging.getLogger(__name__)

subprocess32._has_poll = False

def makeURL(path):
    if not (path.startswith("file:") or path.startswith("s3:") or path.startswith("http:")):
        return "file://" + os.path.abspath(path)
    else:
        return path

def catFiles(filesToCat, catFile):
    """Cats a bunch of files into one file. Ensures a no more than maxCat files
    are concatenated at each step.
    """
    if len(filesToCat) == 0: #We must handle this case or the cat call will hang waiting for input
        open(catFile, 'w').close()
        return
    maxCat = 25
    system("cat %s > %s" % (" ".join(filesToCat[:maxCat]), catFile))
    filesToCat = filesToCat[maxCat:]
    while len(filesToCat) > 0:
        system("cat %s >> %s" % (" ".join(filesToCat[:maxCat]), catFile))
        filesToCat = filesToCat[maxCat:]

def cactusRootPath():
    """
    function for finding external location
    """
    import cactus
    i = os.path.abspath(cactus.__file__)
    return os.path.split(i)[0]

def getLogLevelString2(logLevelString):
    """Gets the log level string for the binary
    """
    if logLevelString == None:
        return getLogLevelString()
    return logLevelString

def getOptionalAttrib(node, attribName, typeFn=None, default=None):
    """Get an optional attrib, or default if not set or node is None
    """
    if node != None and node.attrib.has_key(attribName):
        if typeFn != None:
            if typeFn == bool:
                return bool(int(node.attrib[attribName]))
            return typeFn(node.attrib[attribName])
        return node.attrib[attribName]
    return default

def findRequiredNode(configNode, nodeName):
    """Retrieve an xml node, complain if it's not there."""
    nodes = configNode.findall(nodeName)
    if nodes == None:
        raise RuntimeError("Could not find any nodes with name %s in %s node" % (nodeName, configNode))
    assert len(nodes) == 1, "More than 1 node for %s in config XML" % nodeName
    return nodes[0]

#############################################
#############################################
#Following used to gather the names of flowers
#in problems
#############################################
#############################################  

def readFlowerNames(flowerStrings):
    ret = []
    for line in flowerStrings.split("\n"):
        if line == '':
            continue
        flowersAndSizes = line[1:].split()
        numFlowers = flowersAndSizes[0]
        flowers = []
        sizes = []
        currentlyAFlower = True
        for token in flowersAndSizes[1:]:
            if token == 'a' or token == 'b':
                flowers += [token]
            elif currentlyAFlower:
                flowers += [token]
                currentlyAFlower = False
            else:
                sizes += [int(token)]
                currentlyAFlower = True
        assert len(sizes) == int(numFlowers)
        ret += [(bool(int(line[0])), " ".join([numFlowers] + flowers), sizes)]
    return ret

def runCactusGetFlowers(cactusDiskDatabaseString, flowerNames,
                        jobName=None, features=None, fileStore=None,
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Gets a list of flowers attached to the given flower. 
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = cactus_call(check_output=True, stdin_string=flowerNames,
                                parameters=["cactus_workflow_getFlowers", logLevel,
                                            cactusDiskDatabaseString,
                                            str(minSequenceSizeOfFlower),
                                            str(maxSequenceSizeOfFlowerGrouping),
                                            str(maxSequenceSizeOfSecondaryFlowerGrouping)],
                                job_name=jobName,
                                features=features,
                                fileStore=fileStore)

    l = readFlowerNames(flowerStrings)
    return l

def runCactusExtendFlowers(cactusDiskDatabaseString, flowerNames, 
                        jobName=None, features=None, fileStore=None,
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = cactus_call(check_output=True, stdin_string=flowerNames,
                                parameters=["cactus_workflow_extendFlowers", logLevel,
                                            cactusDiskDatabaseString,
                                            str(minSequenceSizeOfFlower),
                                            str(maxSequenceSizeOfFlowerGrouping),
                                            str(maxSequenceSizeOfSecondaryFlowerGrouping)],
                                job_name=jobName,
                                features=features,
                                fileStore=fileStore)

    l = readFlowerNames(flowerStrings)
    return l

def encodeFlowerNames(flowerNames):
    if len(flowerNames) == 0:
        return "0"
    return "%i %s" % (len(flowerNames), " ".join([ str(flowerNames[0]) ] + [ str(flowerNames[i] - flowerNames[i-1]) for i in xrange(1, len(flowerNames)) ]))
    
def decodeFirstFlowerName(encodedFlowerNames):
    tokens = encodedFlowerNames.split()
    if int(tokens[0]) == 0:
        return None
    if tokens[1] == 'b':
        return int(tokens[2])
    return int(tokens[1])

def runCactusSplitFlowersBySecondaryGrouping(flowerNames):
    """Splits a list of flowers into smaller lists.
    """
    flowerNames = flowerNames.split()
    flowerGroups = []
    stack = []
    overlarge = False
    name = 0
    for i in flowerNames[1:]:
        if i != '':
            if i in ('a', 'b'):
                if len(stack) > 0:
                    flowerGroups.append((overlarge, encodeFlowerNames(stack))) #b indicates the stack is overlarge
                    stack = []
                overlarge = i == 'b'
            else:
                name = int(i) + name
                stack.append(name)
    if len(stack) > 0:
        flowerGroups.append((overlarge, encodeFlowerNames(stack)))
    return flowerGroups

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################  

def runCactusSetup(cactusDiskDatabaseString, sequences, 
                   newickTreeString, logLevel=None, outgroupEvents=None,
                   makeEventHeadersAlphaNumeric=False):
    logLevel = getLogLevelString2(logLevel)
    args = ["--speciesTree", newickTreeString, "--cactusDisk", cactusDiskDatabaseString,
            "--logLevel", logLevel]
    if makeEventHeadersAlphaNumeric:
        args += ["--makeEventHeadersAlphaNumeric"]
    if outgroupEvents is not None:
        args += ["--outgroupEvents", outgroupEvents]
    masterMessages = cactus_call(check_output=True,
                                 parameters=["cactus_setup"] + args + sequences)

    logger.info("Ran cactus setup okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]

def runConvertAlignmentsToInternalNames(cactusDiskString, alignmentsFile, outputFile, flowerName, isBedFile=False):
    args = [alignmentsFile, outputFile,
            "--cactusDisk", cactusDiskString]
    if isBedFile:
        args += ["--bed"]
    cactus_call(stdin_string=encodeFlowerNames((flowerName,)),
                parameters=["cactus_convertAlignmentsToInternalNames"] + args)

def runStripUniqueIDs(cactusDiskString):
    cactus_call(parameters=["cactus_stripUniqueIDs", "--cactusDisk", cactusDiskString])

def runCactusCaf(cactusDiskDatabaseString, alignments,
                 flowerNames=encodeFlowerNames((0,)),
                 logLevel=None, 
                 writeDebugFiles=False,
                 annealingRounds=None,
                 deannealingRounds=None,
                 trim=None,
                 minimumTreeCoverage=None,
                 blockTrim=None,
                 minimumBlockDegree=None,
                 minimumIngroupDegree=None,
                 minimumOutgroupDegree=None,
                 alignmentFilter=None,
                 lastzArguments=None,
                 minimumSequenceLengthForBlast=None,
                 maxAdjacencyComponentSizeRatio=None,
                 constraints=None,
                 minLengthForChromosome=None,
                 proportionOfUnalignedBasesForNewChromosome=None, 
                 maximumMedianSequenceLengthBetweenLinkedEnds=None,
                 realign=False,
                 realignArguments=None,
                 phylogenyNumTrees=None,
                 phylogenyScoringMethod=None,
                 phylogenyRootingMethod=None,
                 phylogenyBreakpointScalingFactor=None,
                 phylogenySkipSingleCopyBlocks=False,
                 phylogenyMaxBaseDistance=None,
                 phylogenyMaxBlockDistance=None,
                 phylogenyDebugFile=None,
                 phylogenyKeepSingleDegreeBlocks=False,
                 phylogenyTreeBuildingMethod=None,
                 phylogenyCostPerDupPerBase=None,
                 phylogenyCostPerLossPerBase=None,
                 referenceEventHeader=None,
                 phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce=None,
                 numTreeBuildingThreads=None,
                 doPhylogeny=False,
                 removeLargestBlock=None,
                 phylogenyNucleotideScalingFactor=None,
                 minimumBlockDegreeToCheckSupport=None,
                 minimumBlockHomologySupport=None,
                 removeRecoverableChains=None,
                 minimumNumberOfSpecies=None,
                 maxRecoverableChainsIterations=None,
                 maxRecoverableChainLength=None,
                 phylogenyHomologyUnitType=None,
                 phylogenyDistanceCorrectionMethod=None,
                 features=None,
                 jobName=None,
                 fileStore=None):
    logLevel = getLogLevelString2(logLevel)
    args = ["--logLevel", logLevel, "--alignments", alignments, "--cactusDisk", cactusDiskDatabaseString]
    if annealingRounds is not None:
        args += ["--annealingRounds", annealingRounds]
    if deannealingRounds is not None:
        args += ["--deannealingRounds", deannealingRounds]
    if trim is not None:
        args += ["--trim", trim]
    if lastzArguments is not None:
        args += ["--lastzArguments", lastzArguments]
    if minimumTreeCoverage is not None:
        args += ["--minimumTreeCoverage", str(minimumTreeCoverage)]
    if blockTrim is not None:
        args += ["--blockTrim", str(blockTrim)]
    if minimumBlockDegree is not None:
        args += ["--minimumDegree", str(minimumBlockDegree)]
    if minimumSequenceLengthForBlast is not None:
        args += ["--minimumSequenceLengthForBlast", str(minimumSequenceLengthForBlast)]
    if minimumIngroupDegree is not None:
        args += ["--minimumIngroupDegree", str(minimumIngroupDegree)]
    if minimumOutgroupDegree is not None:
        args += ["--minimumOutgroupDegree", str(minimumOutgroupDegree)]
    if alignmentFilter is not None:
        args += ["--alignmentFilter", alignmentFilter]
    if maxAdjacencyComponentSizeRatio is not None:
        args += ["--maxAdjacencyComponentSizeRatio", str(maxAdjacencyComponentSizeRatio)]
    if constraints is not None:
        args += ["--constraints", constraints]
    if realign:
        args += ["--realign"]
    if realignArguments is not None:
        args += ["--realignArguments", realignArguments]
    if phylogenyNumTrees is not None:
        args += ["--phylogenyNumTrees", str(phylogenyNumTrees)]
    if phylogenyRootingMethod is not None:
        args += ["--phylogenyRootingMethod", phylogenyRootingMethod]
    if phylogenyScoringMethod is not None:
        args += ["--phylogenyScoringMethod", phylogenyScoringMethod]
    if phylogenyBreakpointScalingFactor is not None:
        args += ["--phylogenyBreakpointScalingFactor", str(phylogenyBreakpointScalingFactor)]
    if phylogenySkipSingleCopyBlocks:
        args += ["--phylogenySkipSingleCopyBlocks"]
    if phylogenyMaxBaseDistance is not None:
        args += ["--phylogenyMaxBaseDistance", str(phylogenyMaxBaseDistance)]
    if phylogenyMaxBlockDistance is not None:
        args += ["--phylogenyMaxBlockDistance", str(phylogenyMaxBlockDistance)]
    if phylogenyDebugFile is not None:
        args += ["--phylogenyDebugFile", phylogenyDebugFile]
    if phylogenyKeepSingleDegreeBlocks:
        args += ["--phylogenyKeepSingleDegreeBlocks"]
    if phylogenyTreeBuildingMethod is not None:
        args += ["--phylogenyTreeBuildingMethod", phylogenyTreeBuildingMethod]
    if phylogenyCostPerDupPerBase is not None:
        args += ["--phylogenyCostPerDupPerBase", str(phylogenyCostPerDupPerBase)]
    if phylogenyCostPerLossPerBase is not None:
        args += ["--phylogenyCostPerLossPerBase", str(phylogenyCostPerLossPerBase)]
    if referenceEventHeader is not None:
        args += ["--referenceEventHeader", referenceEventHeader]
    if phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce is not None:
        args += ["--phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce", str(phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce)]
    if numTreeBuildingThreads is not None:
        args += ["--numTreeBuildingThreads", str(numTreeBuildingThreads)]
    if doPhylogeny:
        args += ["--phylogeny"]
    if minimumBlockDegreeToCheckSupport is not None:
        args += ["--minimumBlockDegreeToCheckSupport", str(minimumBlockDegreeToCheckSupport)]
    if minimumBlockHomologySupport is not None:
        args += ["--minimumBlockHomologySupport", str(minimumBlockHomologySupport)]
    if phylogenyNucleotideScalingFactor is not None:
        args += ["--phylogenyNucleotideScalingFactor", str(phylogenyNucleotideScalingFactor)]
    if removeRecoverableChains is not None:
        args += ["--removeRecoverableChains", removeRecoverableChains]
    if minimumNumberOfSpecies is not None:
        args += ["--minimumNumberOfSpecies", str(minimumNumberOfSpecies)]
    if maxRecoverableChainsIterations is not None:
        args += ["--maxRecoverableChainsIterations", str(maxRecoverableChainsIterations)]
    if maxRecoverableChainLength is not None:
        args += ["--maxRecoverableChainLength", str(maxRecoverableChainLength)]
    if phylogenyHomologyUnitType is not None:
        args += ["--phylogenyHomologyUnitType", phylogenyHomologyUnitType]
    if phylogenyDistanceCorrectionMethod is not None:
        args += ["--phylogenyDistanceCorrectionMethod", phylogenyDistanceCorrectionMethod]
    if minLengthForChromosome is not None:
        args += ["--minLengthForChromosome", str(minLengthForChromosome)]
    if proportionOfUnalignedBasesForNewChromosome is not None:
        args += ["--proportionOfUnalignedBasesForNewChromosome", str(proportionOfUnalignedBasesForNewChromosome)]
    if maximumMedianSequenceLengthBetweenLinkedEnds is not None:
        args += ["--maximumMedianSequenceLengthBetweenLinkedEnds", str(maximumMedianSequenceLengthBetweenLinkedEnds)]

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 parameters=["cactus_caf"] + args,
                                 features=features, job_name=jobName, fileStore=fileStore)
    logger.info("Ran cactus_caf okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]

def runCactusPhylogeny(cactusDiskDatabaseString,
                       flowerNames=encodeFlowerNames((0,)),
                       logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_phylogeny",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--logLevel", logLevel])
    logger.info("Ran cactus_phylogeny okay")

def runCactusAdjacencies(cactusDiskDatabaseString, flowerNames=encodeFlowerNames((0,)), logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_fillAdjacencies",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--logLevel", logLevel])
    logger.info("Ran cactus_fillAdjacencies OK")

def runCactusConvertAlignmentToCactus(cactusDiskDatabaseString, constraintsFile, newConstraintsFile, logLevel=None):
    """Takes a cigar file and makes an equivalent cigar file using the internal coordinate system format of cactus.
    """
    logLevel = getLogLevelString2(logLevel)
    cactus_call(parameters=["cactus_workflow_convertAlignmentCoordinates",
                            logLevel, cactusDiskDatabaseString,
                            constraintsFile, newConstraintsFile])

def runCactusFlowerStats(cactusDiskDatabaseString, flowerName, logLevel=None):
    """Prints stats for the given flower
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStatsString = cactus_call(check_output=True,
                                    parameters=["cactus_workflow_flowerStats",
                                                logLevel, cactusDiskDatabaseString, str(flowerName)])
    return flowerStatsString

def runCactusMakeNormal(cactusDiskDatabaseString, flowerNames, maxNumberOfChains=0, logLevel=None):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_normalisation",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--maxNumberOfChains", str(maxNumberOfChains),
                            "--logLevel", logLevel])

def runCactusBar(cactusDiskDatabaseString, flowerNames, logLevel=None,
                 spanningTrees=None, maximumLength=None, 
                 gapGamma=None,
                 matchGamma=None,
                 splitMatrixBiggerThanThis=None,
                 anchorMatrixBiggerThanThis=None,
                 repeatMaskMatrixBiggerThanThis=None,
                 diagonalExpansion=None,
                 constraintDiagonalTrim=None,
                 minimumBlockDegree=None,
                 minimumIngroupDegree=None,
                 minimumOutgroupDegree=None,
                 alignAmbiguityCharacters=False,
                 pruneOutStubAlignments=False,
                 useProgressiveMerging=False,
                 calculateWhichEndsToComputeSeparately=False,
                 largeEndSize=None,
                 endAlignmentsToPrecomputeOutputFile=None,
                 precomputedAlignments=None,
                 ingroupCoverageFile=None,
                 minimumSizeToRescue=None,
                 minimumCoverageToRescue=None,
                 minimumNumberOfSpecies=None,
                 jobName=None,
                 fileStore=None,
                 features=None):
    """Runs cactus base aligner."""
    logLevel = getLogLevelString2(logLevel)
    args = ["--logLevel", logLevel, "--cactusDisk", cactusDiskDatabaseString]
    if maximumLength is not None:
        args += ["--maximumLength", str(maximumLength)]
    if spanningTrees is not None:
        args += ["--spanningTrees", str(spanningTrees)]
    if gapGamma is not None:
        args += ["--gapGamma", str(gapGamma)]
    if matchGamma is not None:
        args += ["--matchGamma", str(matchGamma)]
    if splitMatrixBiggerThanThis is not None:
        args += ["--splitMatrixBiggerThanThis", str(splitMatrixBiggerThanThis)]
    if anchorMatrixBiggerThanThis is not None:
        args += ["--anchorMatrixBiggerThanThis", str(anchorMatrixBiggerThanThis)]
    if repeatMaskMatrixBiggerThanThis is not None:
        args += ["--repeatMaskMatrixBiggerThanThis", str(repeatMaskMatrixBiggerThanThis)]
    if diagonalExpansion is not None:
        args += ["--diagonalExpansion", str(diagonalExpansion)]
    if constraintDiagonalTrim is not None:
        args += ["--constraintDiagonalTrim", str(constraintDiagonalTrim)]
    if minimumBlockDegree is not None:
        args += ["--minimumDegree", str(minimumBlockDegree)]
    if minimumIngroupDegree is not None:
        args += ["--minimumIngroupDegree", str(minimumIngroupDegree)]
    if minimumOutgroupDegree is not None:
        args += ["--minimumOutgroupDegree", str(minimumOutgroupDegree)]
    if pruneOutStubAlignments:
        args += ["--pruneOutStubAlignments"]
    if alignAmbiguityCharacters:
        args += ["--alignAmbiguityCharacters"]
    if useProgressiveMerging:
        args += ["--useProgressiveMerging"]
    if calculateWhichEndsToComputeSeparately:
        args += ["--calculateWhichEndsToComputeSeparately"]
    if largeEndSize is not None:
        args += ["--largeEndSize", str(largeEndSize)]
    if endAlignmentsToPrecomputeOutputFile is not None:
        endAlignmentsToPrecomputeOutputFile = os.path.basename(endAlignmentsToPrecomputeOutputFile)
        args += ["--endAlignmentsToPrecomputeOutputFile", endAlignmentsToPrecomputeOutputFile]
    if precomputedAlignments is not None:
        precomputedAlignments = map(os.path.basename, precomputedAlignments)
        precomputedAlignments = " ".join(precomputedAlignments)
        args += ["--precomputedAlignments", precomputedAlignments]
    if ingroupCoverageFile is not None:
        args += ["--ingroupCoverageFile", ingroupCoverageFile]
    if minimumSizeToRescue is not None:
        args += ["--minimumSizeToRescue", str(minimumSizeToRescue)]
    if minimumCoverageToRescue is not None:
        args += ["--minimumCoverageToRescue", str(minimumCoverageToRescue)]
    if minimumNumberOfSpecies is not None:
        args += ["--minimumNumberOfSpecies", str(minimumNumberOfSpecies)]

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 parameters=["cactus_bar"] + args,
                                 job_name=jobName, fileStore=fileStore, features=features)

    logger.info("Ran cactus_bar okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]

def runCactusSecondaryDatabase(secondaryDatabaseString, create=True):
    cactus_call(parameters=["cactus_secondaryDatabase",
                secondaryDatabaseString, create])
            
def runCactusReference(cactusDiskDatabaseString, flowerNames, logLevel=None,
                       jobName=None, features=None, fileStore=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=False,
                       theta=None,
                       phi=None, 
                       maxWalkForCalculatingZ=None,
                       ignoreUnalignedGaps=False,
                       wiggle=None, 
                       numberOfNs=None,
                       minNumberOfSequencesToSupportAdjacency=None,
                       makeScaffolds=False):
    """Runs cactus reference."""
    logLevel = getLogLevelString2(logLevel)
    args = ["--logLevel", logLevel, "--cactusDisk", cactusDiskDatabaseString]
    if matchingAlgorithm is not None:
        args += ["--matchingAlgorithm", matchingAlgorithm]
    if referenceEventString is not None:
        args += ["--referenceEventString", referenceEventString]
    if permutations is not None:
        args += ["--permutations", str(permutations)]
    if useSimulatedAnnealing:
        args += ["--useSimulatedAnnealing"]
    if theta is not None:
        args += ["--theta", str(theta)]
    if phi is not None:
        args += ["--phi", str(phi)]
    if maxWalkForCalculatingZ is not None:
        args += ["--maxWalkForCalculatingZ", str(maxWalkForCalculatingZ)]
    if ignoreUnalignedGaps:
        args += ["--ignoreUnalignedGaps"]
    if wiggle is not None:
        args += ["--wiggle", str(wiggle)]
    if numberOfNs is not None:
        args += ["--numberOfNs", str(numberOfNs)]
    if minNumberOfSequencesToSupportAdjacency is not None:
        args += ["--minNumberOfSequencesToSupportAdjacency", str(minNumberOfSequencesToSupportAdjacency)]
    if makeScaffolds:
        args += ["--makeScaffolds"]

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 parameters=["cactus_reference"] + args,
                                 job_name=jobName,
                                 features=features,
                                 fileStore=fileStore)
    logger.info("Ran cactus_reference okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusAddReferenceCoordinates(cactusDiskDatabaseString, flowerNames,
                                     jobName=None, fileStore=None, features=None,
                                     logLevel=None, referenceEventString=None,
                                     outgroupEventString=None, secondaryDatabaseString=None,
                                     bottomUpPhase=False):
    logLevel = getLogLevelString2(logLevel)
    args = ["--logLevel", logLevel, "--cactusDisk", cactusDiskDatabaseString]
    if bottomUpPhase:
        args += ["--bottomUpPhase"]
    if referenceEventString is not None:
        args += ["--referenceEventString", referenceEventString]
    if outgroupEventString is not None:
        args += ["--outgroupEventString", outgroupEventString]
    if secondaryDatabaseString is not None:
        args += ["--secondaryDisk", secondaryDatabaseString]
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_addReferenceCoordinates"] + args,
                job_name=jobName,
                features=features,
                fileStore=fileStore)

def runCactusCheck(cactusDiskDatabaseString, 
                   flowerNames=encodeFlowerNames((0,)), 
                   logLevel=None, 
                   recursive=False,
                   checkNormalised=False):
    logLevel = getLogLevelString2(logLevel)
    args = ["--cactusDisk", cactusDiskDatabaseString, "--logLevel", logLevel]
    if recursive:
        args += ["--recursive"]
    if checkNormalised:
        args += ["--checkNormalised"]
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_check"] + args)
    logger.info("Ran cactus check")

def _fn(toilDir,
      logLevel=None, retryCount=0, 
      batchSystem="single_machine", 
      rescueJobFrequency=None,
      buildAvgs=False, buildReference=False,
      buildHal=False,
      buildFasta=False,
      toilStats=False,
      maxThreads=None,
      maxCpus=None,
      defaultMemory=None,
      logFile=None):
    logLevel = getLogLevelString2(logLevel)
    args = [toilDir, "--logLevel", logLevel]
    if buildAvgs:
        args += ["--buildAvgs"]
    if buildReference:
        args += ["--buildReference"]
    if buildHal:
        args += ["--buildHal"]
    if buildFasta:
        args += ["--buildFasta"]
    #Jobtree args
    if batchSystem is not None:
        args += ["--batchSystem", batchSystem]
    if retryCount is not None:
        args += ["--retryCount", str(retryCount)]
    if rescueJobFrequency is not None:
        args += ["--rescueJobFrequency", str(rescueJobFrequency)]
    if toilStats:
        args += ["--stats"]
    if maxThreads is not None:
        args += ["--maxThreads", str(maxThreads)]
    if maxCpus is not None:
        args += ["--maxCpus", str(maxCpus)]
    if defaultMemory is not None:
        args += ["--defaultMemory", str(defaultMemory)]
    if logFile is not None:
        args += ["--logFile", logFile]
    return args

def runCactusWorkflow(experimentFile,
                      toilDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildAvgs=False, buildReference=False,
                      buildHal=False,
                      buildFasta=False,
                      toilStats=False,
                      maxThreads=None,
                      maxCpus=None,
                      defaultMemory=None,
                      logFile=None,
                      intermediateResultsUrl=None,
                      extraToilArgumentsString=""):
    args = ["--experiment", experimentFile] + _fn(toilDir,
                      logLevel, retryCount, batchSystem, rescueJobFrequency,
                      buildAvgs, buildReference, buildHal, buildFasta, toilStats, maxThreads, maxCpus, defaultMemory, logFile)
    if intermediateResultsUrl is not None:
        args += ["--intermediateResultsUrl", intermediateResultsUrl]

    import cactus.pipeline.cactus_workflow as cactus_workflow
    cactus_workflow.runCactusWorkflow(args)
    logger.info("Ran the cactus workflow okay")
    
def runCactusProgressive(inputDir,
                         toilDir, 
                         logLevel=None, retryCount=0, 
                         batchSystem="single_machine", 
                         rescueJobFrequency=None,
                         skipAlignments=False,
                         buildHal=None,
                         buildFasta=None,
                         buildAvgs=False, 
                         toilStats=False,
                         maxThreads=None,
                         maxCpus=None,
                         logFile=None,
                         defaultMemory=None):
    command = ["cactus_progressive.py", "--project", inputDir] + _fn(toilDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency,
                      buildAvgs, None,
                      buildHal,
                      buildFasta,
                      toilStats, maxThreads, maxCpus, defaultMemory, logFile)
    system(command)                   
    logger.info("Ran the cactus progressive okay")
    
def runCactusHalGenerator(cactusDiskDatabaseString,
                          secondaryDatabaseString, 
                          flowerNames,
                          referenceEventString, 
                          outputFile=None,
                          showOnlySubstitutionsWithRespectToReference=False,
                          logLevel=None,
                          jobName=None,
                          features=None,
                          fileStore=None):
    logLevel = getLogLevelString2(logLevel)
    if outputFile is not None:
        outputFile = os.path.basename(outputFile)
    args = ["--logLevel", logLevel, "--cactusDisk", cactusDiskDatabaseString,
            "--secondaryDisk", secondaryDatabaseString]
    if referenceEventString is not None:
        args += ["--referenceEventString", referenceEventString]
    if outputFile is not None:
        args += ["--outputFile", outputFile]
    if showOnlySubstitutionsWithRespectToReference:
        args += ["--showOnlySubstitutionsWithRespectToReference"]
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_halGenerator"] + args,
                job_name=jobName, features=features, fileStore=fileStore)

def runCactusFastaGenerator(cactusDiskDatabaseString,
                            flowerName,
                            outputFile,
                            referenceEventString,
                            logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(parameters=["cactus_fastaGenerator",
                            "--flowerName", str(flowerName),
                            "--outputFile", outputFile,
                            "--logLevel", logLevel,
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--referenceEventString", referenceEventString])

def runCactusAnalyseAssembly(sequenceFile):
    return cactus_call(check_output=True,
                parameters=["cactus_analyseAssembly",
                            sequenceFile])[:-1]
    
def runToilStats(toil, outputFile):
    system("toil stats %s --outputFile %s" % (toil, outputFile))
    logger.info("Ran the job-tree stats command apparently okay")

def runToilStatusAndFailIfNotComplete(toilDir):
    command = "toil status %s --failIfNotComplete --verbose" % toilDir
    system(command)

def runLastz(seq1, seq2, alignmentsFile, lastzArguments, work_dir=None):
    #Have to specify the work_dir manually for this, since
    #we're adding arguments to the filename
    assert os.path.dirname(seq1) == os.path.dirname(seq2)
    work_dir = os.path.dirname(seq1)
    cactus_call(work_dir=work_dir, outfile=alignmentsFile,
                parameters=["cPecanLastz",
                            "--format=cigar",
                            "--notrivial"] + lastzArguments.split() +
                           ["%s[multiple][nameparse=darkspace]" % os.path.basename(seq1),
                            "%s[nameparse=darkspace]" % os.path.basename(seq2)],
                soft_timeout=5400)

def runSelfLastz(seq, alignmentsFile, lastzArguments, work_dir=None):
    work_dir = os.path.dirname(seq)
    cactus_call(work_dir=work_dir, outfile=alignmentsFile,
                parameters=["cPecanLastz",
                            "--format=cigar",
                            "--notrivial"] + lastzArguments.split() +
                           ["%s[multiple][nameparse=darkspace]" % os.path.basename(seq),
                            "%s[nameparse=darkspace]" % os.path.basename(seq)],
                soft_timeout=5400)

def runCactusRealign(seq1, seq2, inputAlignmentsFile, outputAlignmentsFile, realignArguments, work_dir=None):
    cactus_call(infile=inputAlignmentsFile, outfile=outputAlignmentsFile, work_dir=work_dir,
                parameters=["cPecanRealign"] + realignArguments.split() + [seq1, seq2])

def runCactusSelfRealign(seq, inputAlignmentsFile, outputAlignmentsFile, realignArguments, work_dir=None):
    cactus_call(infile=inputAlignmentsFile, outfile=outputAlignmentsFile, work_dir=work_dir,
                parameters=["cPecanRealign"] + realignArguments.split() + [seq])

def runCactusCoverage(sequenceFile, alignmentsFile, work_dir=None):
    return cactus_call(check_output=True, work_dir=work_dir,
                parameters=["cactus_coverage", sequenceFile, alignmentsFile])

def runGetChunks(sequenceFiles, chunksDir, chunkSize, overlapSize, work_dir=None):
    chunks = cactus_call(work_dir=work_dir,
                         check_output=True,
                         parameters=["cactus_blast_chunkSequences",
                                     getLogLevelString(),
                                     str(chunkSize),
                                     str(overlapSize),
                         chunksDir] + sequenceFiles)
    return [chunk for chunk in chunks.split("\n") if chunk != ""]

def pullCactusImage():
    """Ensure that the cactus Docker image is pulled."""
    if os.environ.get('CACTUS_DOCKER_MODE') == "0":
        return
    image = getDockerImage()
    call = ["docker", "pull", image]
    process = subprocess32.Popen(call, stdout=subprocess32.PIPE,
                                 stderr=sys.stderr, bufsize=-1)
    output, _ = process.communicate()
    if process.returncode != 0:
        raise RuntimeError("Command %s failed with output: %s" % (call, output))

def getDockerOrg():
    """Get where we should find the cactus containers."""
    if "CACTUS_DOCKER_ORG" in os.environ:
        return os.environ["CACTUS_DOCKER_ORG"]
    else:
        return "quay.io/comparative-genomics-toolkit"

def getDockerTag():
    """Get what docker tag we should use for the cactus image
    (either forced to be latest or the current cactus commit)."""
    if 'CACTUS_USE_LATEST' in os.environ:
        return "latest"
    else:
        return cactus_commit

def getDockerImage():
    """Get fully specified Docker image name."""
    return "%s/cactus:%s" % (getDockerOrg(), getDockerTag())

def maxMemUsageOfContainer(containerInfo):
    """Return the max RSS usage (in bytes) of a container, or None if something failed."""
    if containerInfo['id'] is None:
        # Try to get the internal container ID from the docker name
        try:
            id = popenCatch("docker inspect -f '{{.Id}}' %s" % containerInfo['name']).strip()
            containerInfo['id'] = id
        except:
            # Not yet running
            return None
    # Try to check for the maximum memory usage ever used by that
    # container, in a few different possible locations depending on
    # the distribution
    possibleLocations = ["/sys/fs/cgroup/memory/docker/%s/memory.max_usage_in_bytes",
                         "/sys/fs/cgroup/memory/system.slice.docker-%s.scope/memory.max_usage_in_bytes"]
    possibleLocations = [s % containerInfo['id'] for s in possibleLocations]
    for location in possibleLocations:
        try:
            with open(location) as f:
                return int(f.read())
        except IOError:
            # Not at this location, or sysfs isn't mounted
            continue
    return None

def singularityCommand(tool=None,
                       work_dir=None,
                       parameters=None,
                       port=None):
    base_singularity_call = ["singularity", "--silent", "run", os.environ["CACTUS_SINGULARITY_IMG"]]
    base_singularity_call.extend(parameters)
    return base_singularity_call

def dockerCommand(tool=None,
                  work_dir=None,
                  parameters=None,
                  rm=True,
                  port=None,
                  dockstore=None):
    # This is really dumb, but we have to work around an intersection
    # between two bugs: one in CoreOS where /etc/resolv.conf is
    # sometimes missing temporarily, and one in Docker where it
    # refuses to start without /etc/resolv.conf.
    while not os.path.exists('/etc/resolv.conf'):
        pass

    base_docker_call = ['docker', 'run',
                        '--interactive',
                        '--net=host',
                        '--log-driver=none',
                        '-u', '%s:%s' % (os.getuid(), os.getgid()),
                        '-v', '{}:/data'.format(os.path.abspath(work_dir))]

    if port:
        base_docker_call += ["-p", "%d:%d" % (port, port)]

    containerInfo = { 'name': str(uuid.uuid4()), 'id': None }
    base_docker_call.extend(['--name', containerInfo['name']])
    if rm:
        base_docker_call.append('--rm')

    docker_tag = getDockerTag()
    tool = "%s/%s:%s" % (dockstore, tool, docker_tag)
    call = base_docker_call + [tool] + parameters
    return call, containerInfo

def prepareWorkDir(work_dir, parameters):
    def moveToWorkDir(work_dir, arg):
        if isinstance(arg, str) and os.path.isfile(arg):
            if not os.path.dirname(arg) == work_dir:
                _log.info('Copying file %s to work dir' % arg)
                shutil.copy(arg, work_dir)

    if work_dir:
        for arg in parameters:
            moveToWorkDir(work_dir, arg)

    if not work_dir:
    #Make sure all the paths we're accessing are in the same directory
        files = [par for par in parameters if os.path.isfile(par)]
        folders = [par for par in parameters if os.path.isdir(par)]
        work_dirs = set([os.path.dirname(fileName) for fileName in files] + [os.path.dirname(folder) for folder in folders])
        _log.info("Work dirs: %s" % work_dirs)
        if len(work_dirs) > 1:
            work_dir = os.path.commonprefix(work_dirs)
        elif len(work_dirs) == 1:
            work_dir = work_dirs.pop()

    #If there are no input files, or if their MRCA is '' (when working
    #with relative paths), just set the current directory as the work
    #dir
    if work_dir is None or work_dir == '':
        work_dir = "."
    _log.info("Docker work dir: %s" % work_dir)

    #We'll mount the work_dir containing the paths as /data in the container,
    #so set all the paths to their basenames. The container will access them at
    #/data/<path>
    def adjustPath(path, wd):
        # Hack to relativize paths that are not provided as a
        # single argument (i.e. multiple paths that are
        # space-separated and quoted)
        if wd != '.':
            if not wd.endswith('/'):
                wd = wd + '/'
            return path.replace(wd, '')
        else:
            return path

    if work_dir and os.environ.get('CACTUS_DOCKER_MODE') != "0":
        parameters = [adjustPath(par, work_dir) for par in parameters]
    return work_dir, parameters

def cactus_call(tool=None,
                work_dir=None,
                parameters=None,
                rm=True,
                check_output=False,
                infile=None,
                outfile=None,
                stdin_string=None,
                server=False,
                shell=False,
                port=None,
                check_result=False,
                dockstore=None,
                soft_timeout=None,
                job_name=None,
                features=None,
                fileStore=None,
                swallowStdErr=False):
    mode = os.environ.get("CACTUS_BINARIES_MODE", "docker")

    if dockstore is None:
        dockstore = getDockerOrg()
    if parameters is None:
        parameters = []
    if tool is None:
        tool = "cactus"

    if mode in ("docker", "singularity"):
        work_dir, parameters = prepareWorkDir(work_dir, parameters)

    if mode == "docker":
        call, containerInfo = dockerCommand(tool=tool,
                                            work_dir=work_dir,
                                            parameters=parameters,
                                            rm=rm,
                                            port=port,
                                            dockstore=dockstore)
    elif mode == "singularity":
        call = singularityCommand(tool=tool, work_dir=work_dir,
                                  parameters=parameters, port=port)
    else:
        assert mode == "local"
        call = parameters

    stdinFileHandle = None
    stdoutFileHandle = None
    if stdin_string:
        stdinFileHandle = subprocess32.PIPE
    elif infile:
        stdinFileHandle = open(infile, 'r')
    if outfile:
        stdoutFileHandle = open(outfile, 'w')
    if check_output:
        stdoutFileHandle = subprocess32.PIPE

    _log.info("Running the command %s" % call)
    process = subprocess32.Popen(call, shell=shell,
                                 stdin=stdinFileHandle, stdout=stdoutFileHandle,
                                 stderr=subprocess32.PIPE if swallowStdErr else sys.stderr,
                                 bufsize=-1)

    if server:
        return process

    memUsage = 0
    first_run = True
    start_time = time.time()
    while True:
        try:
            # Wait a bit to see if the process is done
            output, nothing = process.communicate(stdin_string if first_run else None, timeout=10)
        except subprocess32.TimeoutExpired:
            if mode == "docker":
                # Every so often, check the memory usage of the container
                updatedMemUsage = maxMemUsageOfContainer(containerInfo)
                if updatedMemUsage is not None:
                    assert memUsage <= updatedMemUsage, "memory.max_usage_in_bytes should never decrease"
                    memUsage = updatedMemUsage
            first_run = False
            if soft_timeout is not None and time.time() - start_time > soft_timeout:
                # Soft timeout has been triggered. Just return early.
                process.send_signal(signal.SIGINT)
                return None
        else:
            break
    if mode == "docker" and job_name is not None and features is not None and fileStore is not None:
        # Log a datapoint for the memory usage for these features.
        fileStore.logToMaster("Max memory used for job %s (tool %s) "
                              "on JSON features %s: %s" % (job_name, parameters[0],
                                                           json.dumps(features), memUsage))
    if check_result:
        return process.returncode

    if process.returncode != 0:
        raise RuntimeError("Command %s failed with output: %s" % (call, output))

    if check_output:
        return output

class RunAsFollowOn(Job):
    def __init__(self, job, *args, **kwargs):
        Job.__init__(self, memory=100000000, preemptable=True)
        self._args = args
        self._kwargs = kwargs
        self.job = job
    def run(self, fileStore):
        return self.addFollowOn(self.job(*self._args, **self._kwargs)).rv()

class RoundedJob(Job):
    """Thin wrapper around Toil.Job to round up resource requirements.

    Rounding is useful to make Toil's Mesos scheduler more
    efficient--it runs a process that is O(n log n) in the number of
    different resource requirements for every offer received, so
    thousands of slightly different requirements will slow down the
    leader and the workflow.
    """
    # Default rounding amount: 100 MiB
    roundingAmount = 100*1024*1024
    def __init__(self, memory=None, cores=None, disk=None, preemptable=None,
                 unitName=None, checkpoint=False):
        if memory is not None:
            memory = self.roundUp(memory)
        if disk is not None:
            disk = self.roundUp(disk)
        super(RoundedJob, self).__init__(memory=memory, cores=cores, disk=disk,
                                         preemptable=preemptable, unitName=unitName,
                                         checkpoint=checkpoint)

    def roundUp(self, bytesRequirement):
        """
        Round the amount up to the next self.roundingAmount.

        >>> j = RoundedJob()
        >>> j.roundingAmount = 100000000
        >>> j.roundUp(1000)
        10000000
        >>> j.roundUp(200000000)
        200000000
        >>> j.roundUp(200000001)
        300000000
        """
        if bytesRequirement % self.roundingAmount == 0:
            return bytesRequirement
        return (bytesRequirement // self.roundingAmount + 1) * self.roundingAmount

    def _runner(self, jobGraph, jobStore, fileStore):
        if jobStore.config.workDir is not None:
            os.environ['TMPDIR'] = fileStore.getLocalTempDir()
        super(RoundedJob, self)._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore)

def readGlobalFileWithoutCache(fileStore, jobStoreID):
    """Reads a jobStoreID into a file and returns it, without touching
    the cache.

    Works around toil issue #1532.
    """
    f = fileStore.getLocalTempFile()
    fileStore.jobStore.readFile(jobStoreID, f)
    return f

class ChildTreeJob(RoundedJob):
    """Spreads the child-job initialization work among multiple jobs.

    Jobs with many children can often be a bottleneck (because they
    are written serially into the jobStore in a consistent-write
    fashion). Subclasses of this job will automatically spread out
    that work amongst a tree of jobs, increasing the total work done
    slightly, but reducing the wall-clock time taken dramatically.
    """
    def __init__(self, memory=None, cores=None, disk=None, preemptable=None,
                 unitName=None, checkpoint=False, maxChildrenPerJob=20):
        self.queuedChildJobs = []
        self.maxChildrenPerJob = maxChildrenPerJob
        super(ChildTreeJob, self).__init__(memory=memory, cores=cores, disk=disk,
                                           preemptable=preemptable, unitName=unitName,
                                           checkpoint=checkpoint)

    def addChild(self, job):
        self.queuedChildJobs.append(job)
        return job

    def _run(self, jobGraph, fileStore):
        ret = super(ChildTreeJob, self)._run(jobGraph, fileStore)
        if len(self.queuedChildJobs) <= self.maxChildrenPerJob:
            # The number of children is small enough that we can just
            # add them directly.
            for childJob in self.queuedChildJobs:
                super(ChildTreeJob, self).addChild(childJob)
        else:
            # Too many children, so we have to build a tree to avoid
            # bottlenecking on consistently serializing all the jobs.
            for job in self.queuedChildJobs:
                job.prepareForPromiseRegistration(fileStore.jobStore)

            curLevel = self.queuedChildJobs
            while len(curLevel) > self.maxChildrenPerJob:
                curLevel = [curLevel[i:i + self.maxChildrenPerJob] for i in xrange(0, len(curLevel), self.maxChildrenPerJob)]
            # curLevel is now a nested list (of lists, of lists...)
            # representing a tree of out-degree no higher than
            # maxChildrenPerJob. We can pass that to SpawnChildren
            # instances, which will run down the tree and eventually
            # spawn the jobs we're actually interested in running.
            for sublist in curLevel:
                logger.debug(sublist)
                assert isinstance(sublist, list)
                super(ChildTreeJob, self).addChild(SpawnChildren(sublist))
        return ret

class SpawnChildren(RoundedJob):
    """Helper class used only by ChildTreeJob."""
    def __init__(self, childList, *args, **kwargs):
        self.childList = childList

        # cPickle sometimes has major issues pickling childList,
        # crashing and complaining about attempting to pickle a bound
        # instance method. We are certainly *not* doing that, at least
        # not intentionally, and the (Python) pickle library has no
        # issue pickling the same object. So we hack around this and
        # force the use of pickle for this worker process (which
        # should only last as long as the ChildTreeJob that creates
        # this class).
        cPickle.dump = pickle.dump
        cPickle.dumps = pickle.dumps

        super(SpawnChildren, self).__init__(*args, preemptable=True, **kwargs)

    def run(self, fileStore):
        for item in self.childList:
            if isinstance(item, Job):
                # Hit the leaf nodes of the tree, which are the
                # jobs we actually want to run.
                self.addChild(item)
            else:
                # More nested lists of jobs: we need to spawn more
                # SpawnChildren instances to distribute the load.
                self.addChild(SpawnChildren(item))
