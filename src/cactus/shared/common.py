#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus package.
"""

import os
import random
import sys
import shutil
import subprocess
import logging

from toil.lib.bioio import logger
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString

from toil.lib.docker import _fixPermissions

from toil.job import Job

from sonLib.bioio import getTempDirectory
from sonLib.bioio import nameValue
from sonLib.bioio import popenCatch, popenPush

from cactus.shared.version import cactus_commit

_log = logging.getLogger(__name__)


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

def nameValue(name, value, valueType=str, quotes=False):
    """Little function to make it easier to make name value strings for commands.
    """
    if valueType == bool:
        if value:
            return "--%s" % name
        return ""
    if value is None:
        return ""
    if quotes:
        return "--%s '%s'" % (name, valueType(value))
    return "--%s %s" % (name, valueType(value))

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

def findRequiredNode(configNode, nodeName, index=0):
    """Retrieve an xml node, complain if its not there.
    """
    nodes = configNode.findall(nodeName)
    if nodes == None:
        raise RuntimeError("Could not find any nodes with name %s in %s node" % (nodeName, configNode))
    if index >= len(nodes):
        raise RuntimeError("Could not find a node with name %s and index %i in %s node" % (nodeName, index, configNode))
    return nodes[index]

#############################################
#############################################
#Following used to gather the names of flowers
#in problems
#############################################
#############################################  

def readFlowerNames(flowerStrings): 
    return [ (bool(int(line[0])), line[1:]) for line in flowerStrings.split("\n") if line != '' ]
    
def runCactusGetFlowers(cactusDiskDatabaseString, cactusSequencesPath, flowerNames, 
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Gets a list of flowers attached to the given flower. 
    """
    logLevel = getLogLevelString2(logLevel)
    cactusSequencesPath = os.path.basename(cactusSequencesPath)
    flowerStrings = cactus_call(check_output=True, stdin_string=flowerNames,
                                parameters=["cactus_workflow_getFlowers"],
                                option_string="%s '%s' %s %i %i %i" % 
                                            (logLevel, cactusDiskDatabaseString,
                                            cactusSequencesPath, minSequenceSizeOfFlower,
                                            maxSequenceSizeOfFlowerGrouping,
                                            maxSequenceSizeOfSecondaryFlowerGrouping))
                                        
    l = readFlowerNames(flowerStrings)
    return l

def runCactusExtendFlowers(cactusDiskDatabaseString, cactusSequencesPath, flowerNames, 
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    logLevel = getLogLevelString2(logLevel)
    cactusSequencesPath = os.path.basename(cactusSequencesPath)
    flowerStrings = cactus_call(check_output=True, stdin_string=flowerNames,
                                parameters=["cactus_workflow_extendFlowers"],
                                option_string="%s '%s' %s %i %i %i" %
                                            (logLevel,
                                            cactusDiskDatabaseString,
                                            cactusSequencesPath,
                                            minSequenceSizeOfFlower,
                                            maxSequenceSizeOfFlowerGrouping,
                                            maxSequenceSizeOfSecondaryFlowerGrouping))
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

def runCactusSetup(cactusDiskDatabaseString, cactusSequencesPath, sequences, 
                   newickTreeString, logLevel=None, outgroupEvents=None,
                   makeEventHeadersAlphaNumeric=None):
    logLevel = getLogLevelString2(logLevel)
    outgroupEvents = nameValue("outgroupEvents", outgroupEvents, str, quotes=True)
    makeEventHeadersAlphaNumeric=nameValue("makeEventHeadersAlphaNumeric", makeEventHeadersAlphaNumeric, bool)
    masterMessages = cactus_call(check_output=True,
                                 parameters=["cactus_setup"] + sequences,
                                 option_string="--speciesTree '%s' --cactusDisk '%s' --cactusSequencesPath '%s' --logLevel %s %s %s" % (newickTreeString, cactusDiskDatabaseString, os.path.basename(cactusSequencesPath), logLevel, outgroupEvents, makeEventHeadersAlphaNumeric))
    
    logger.info("Ran cactus setup okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    


def runConvertAlignmentsToInternalNames(cactusDiskString, cactusSequencesPath, alignmentsFile, outputFile, flowerName, isBedFile = False):
    isBedFile = nameValue("bed", isBedFile, bool)
    cactus_call(stdin_string=encodeFlowerNames((flowerName,)),
                option_string="--cactusDisk '%s'" % cactusDiskString,
                parameters=["cactus_convertAlignmentsToInternalNames",
                            "--cactusSequencesPath", cactusSequencesPath,
                            alignmentsFile, outputFile,
                            isBedFile])
    
def runStripUniqueIDs(cactusDiskString, cactusSequencesPath):
    cactus_call(option_string="--cactusDisk '%s'" % cactusDiskString,
                parameters=["cactus_stripUniqueIDs",
                            "--cactusSequencesPath", cactusSequencesPath])
    

def runCactusCaf(cactusDiskDatabaseString, cactusSequencesPath, alignments, 
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
                 realign=None,
                 realignArguments=None,
                 phylogenyNumTrees=None,
                 phylogenyScoringMethod=None,
                 phylogenyRootingMethod=None,
                 phylogenyBreakpointScalingFactor=None,
                 phylogenySkipSingleCopyBlocks=None,
                 phylogenyMaxBaseDistance=None,
                 phylogenyMaxBlockDistance=None,
                 phylogenyDebugFile=None,
                 phylogenyKeepSingleDegreeBlocks=None,
                 phylogenyTreeBuildingMethod=None,
                 phylogenyCostPerDupPerBase=None,
                 phylogenyCostPerLossPerBase=None,
                 referenceEventHeader=None,
                 phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce=None,
                 numTreeBuildingThreads=None,
                 doPhylogeny=None,
                 removeLargestBlock=None,
                 phylogenyNucleotideScalingFactor=None,
                 minimumBlockDegreeToCheckSupport=None,
                 minimumBlockHomologySupport=None,
                 removeRecoverableChains=None,
                 minimumNumberOfSpecies=None,
                 maxRecoverableChainsIterations=None,
                 maxRecoverableChainLength=None,
                 phylogenyHomologyUnitType=None,
                 phylogenyDistanceCorrectionMethod=None):
    # remove annoying carriage returns in caf command line.
    cactusDiskDatabaseString = cactusDiskDatabaseString.replace('\n', '')

    if alignments:
        alignments = os.path.basename(alignments)

    logLevel = getLogLevelString2(logLevel)
    annealingRounds = nameValue("annealingRounds", annealingRounds, quotes=True)
    deannealingRounds = nameValue("deannealingRounds", deannealingRounds, quotes=True)
    trim = nameValue("trim", trim, quotes=True)
    alignments = nameValue("alignments", alignments)
    lastzArguments = nameValue("lastzArguments", lastzArguments, quotes=True)
    minimumTreeCoverage = nameValue("minimumTreeCoverage", minimumTreeCoverage, float)
    blockTrim = nameValue("blockTrim", blockTrim, int)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    minimumSequenceLengthForBlast = nameValue("minimumSequenceLengthForBlast", minimumSequenceLengthForBlast, int)
    minimumIngroupDegree = nameValue("minimumIngroupDegree", minimumIngroupDegree, int)
    minimumOutgroupDegree = nameValue("minimumOutgroupDegree", minimumOutgroupDegree, int)
    alignmentFilter = nameValue("alignmentFilter", alignmentFilter)
    maxAdjacencyComponentSizeRatio = nameValue("maxAdjacencyComponentSizeRatio", maxAdjacencyComponentSizeRatio, float)
    constraints = nameValue("constraints", constraints)
    realign = nameValue("realign", realign, bool)
    realignArguments = nameValue("realignArguments", realignArguments, quotes=True)
    phylogenyNumTrees = nameValue("phylogenyNumTrees", phylogenyNumTrees, int)
    phylogenyRootingMethod = nameValue("phylogenyRootingMethod", phylogenyRootingMethod, quotes=True)
    phylogenyScoringMethod = nameValue("phylogenyScoringMethod", phylogenyScoringMethod, quotes=True)
    phylogenyBreakpointScalingFactor = nameValue("phylogenyBreakpointScalingFactor", phylogenyBreakpointScalingFactor)
    phylogenySkipSingleCopyBlocks = nameValue("phylogenySkipSingleCopyBlocks", phylogenySkipSingleCopyBlocks, bool)
    phylogenyMaxBaseDistance = nameValue("phylogenyMaxBaseDistance", phylogenyMaxBaseDistance)
    phylogenyMaxBlockDistance = nameValue("phylogenyMaxBlockDistance", phylogenyMaxBlockDistance)
    phylogenyDebugFile = nameValue("phylogenyDebugFile", phylogenyDebugFile)
    phylogenyKeepSingleDegreeBlocks = nameValue("phylogenyKeepSingleDegreeBlocks", phylogenyKeepSingleDegreeBlocks, bool)
    phylogenyTreeBuildingMethod = nameValue("phylogenyTreeBuildingMethod", phylogenyTreeBuildingMethod)
    phylogenyCostPerDupPerBase = nameValue("phylogenyCostPerDupPerBase", phylogenyCostPerDupPerBase)
    phylogenyCostPerLossPerBase = nameValue("phylogenyCostPerLossPerBase", phylogenyCostPerLossPerBase)
    referenceEventHeader = nameValue("referenceEventHeader", referenceEventHeader, quotes=True)
    phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce = nameValue("phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce", phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce)
    numTreeBuildingThreads = nameValue("numTreeBuildingThreads", numTreeBuildingThreads)
    doPhylogeny = nameValue("phylogeny", doPhylogeny, bool)
    minimumBlockDegreeToCheckSupport = nameValue("minimumBlockDegreeToCheckSupport", minimumBlockDegreeToCheckSupport)
    minimumBlockHomologySupport = nameValue("minimumBlockHomologySupport", minimumBlockHomologySupport)
    phylogenyNucleotideScalingFactor = nameValue("phylogenyNucleotideScalingFactor", phylogenyNucleotideScalingFactor)
    removeRecoverableChains = nameValue("removeRecoverableChains", removeRecoverableChains)
    minimumNumberOfSpecies = nameValue("minimumNumberOfSpecies", minimumNumberOfSpecies, int)
    maxRecoverableChainsIterations = nameValue("maxRecoverableChainsIterations", maxRecoverableChainsIterations, int)
    maxRecoverableChainLength = nameValue("maxRecoverableChainLength", maxRecoverableChainLength, int)
    phylogenyHomologyUnitType = nameValue("phylogenyHomologyUnitType", phylogenyHomologyUnitType, quotes=True)
    phylogenyDistanceCorrectionMethod = nameValue("phylogenyDistanceCorrectionMethod", phylogenyDistanceCorrectionMethod, quotes=True)

    minLengthForChromosome = nameValue("minLengthForChromosome", minLengthForChromosome, int)
    proportionOfUnalignedBasesForNewChromosome = nameValue("proportionOfUnalignedBasesForNewChromosome", proportionOfUnalignedBasesForNewChromosome, float)
    maximumMedianSequenceLengthBetweenLinkedEnds = nameValue("maximumMedianSequenceLengthBetweenLinkedEnds", maximumMedianSequenceLengthBetweenLinkedEnds, int)

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                                 parameters=["cactus_caf",
                                             "--cactusSequencesPath", cactusSequencesPath,
                                             "--logLevel", logLevel, alignments, annealingRounds,
                                             deannealingRounds, 
                                             trim, minimumTreeCoverage, blockTrim, 
                                             minimumBlockDegree, minimumIngroupDegree, minimumOutgroupDegree,  
                                             alignmentFilter, lastzArguments, minimumSequenceLengthForBlast,
                                             maxAdjacencyComponentSizeRatio, constraints,
                                             minLengthForChromosome,
                                             proportionOfUnalignedBasesForNewChromosome,
                                             maximumMedianSequenceLengthBetweenLinkedEnds, realign,
                                             realignArguments, phylogenyNumTrees, phylogenyRootingMethod,
                                             phylogenyScoringMethod, phylogenyBreakpointScalingFactor,
                                             phylogenySkipSingleCopyBlocks, phylogenyMaxBaseDistance,
                                             phylogenyMaxBlockDistance, phylogenyDebugFile,
                                             phylogenyKeepSingleDegreeBlocks, phylogenyTreeBuildingMethod,
                                             phylogenyCostPerDupPerBase, phylogenyCostPerLossPerBase,
                                             referenceEventHeader,
                                             phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce,
                                             numTreeBuildingThreads, doPhylogeny,
                                             minimumBlockDegreeToCheckSupport, minimumBlockHomologySupport,
                                             phylogenyNucleotideScalingFactor, removeRecoverableChains,
                                             minimumNumberOfSpecies, phylogenyHomologyUnitType,
                                             phylogenyDistanceCorrectionMethod,
                                             maxRecoverableChainsIterations, maxRecoverableChainLength])
    logger.info("Ran cactus_core okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusPhylogeny(cactusDiskDatabaseString,
                  cactusSequencesPath,
                  flowerNames=encodeFlowerNames((0,)),
                  logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_phylogeny",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--cactusSequencesPath", cactusSequencesPath,
                            "--logLevel", logLevel])
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(cactusDiskDatabaseString, flowerNames=encodeFlowerNames((0,)), logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_fillAdjacencies",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--logLevel", logLevel])
    logger.info("Ran cactus_fillAdjacencies OK")

def runCactusConvertAlignmentToCactus(cactusDiskDatabaseString, cactusSequencesPath, constraintsFile, newConstraintsFile, logLevel=None):
    """Takes a cigar file and makes an equivalent cigar file using the internal coordinate system format of cactus.
    """
    logLevel = getLogLevelString2(logLevel)
    cactus_call(parameters=["cactus_workflow_convertAlignmentCoordinates",
                            logLevel, "'%s'" % cactusDiskDatabaseString, cactusSequencesPath,
                            constraintsFile, newConstraintsFile])

def runCactusFlowerStats(cactusDiskDatabaseString, flowerName, logLevel=None):
    """Prints stats for the given flower
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStatsString = cactus_call(check_output=True,
                                    parameters=["cactus_workflow_flowerStats",
                                                logLevel, "'%s'" % cactusDiskDatabaseString, flowerName])
    return flowerStatsString.split("\n")[0]

def runCactusMakeNormal(cactusDiskDatabaseString, cactusSequencesPath, flowerNames, maxNumberOfChains=0, logLevel=None):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    logLevel = getLogLevelString2(logLevel)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_normalisation",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--cactusSequencesPath", cactusSequencesPath,
                            "--maxNumberOfChains", maxNumberOfChains,
                            "--logLevel", logLevel])

def runCactusBar(cactusDiskDatabaseString, cactusSequencesPath, flowerNames, logLevel=None,
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
                         alignAmbiguityCharacters=None,
                         pruneOutStubAlignments=None,
                         useProgressiveMerging=None,
                         calculateWhichEndsToComputeSeparately=None,
                         largeEndSize=None,
                         endAlignmentsToPrecomputeOutputFile=None,
                         precomputedAlignments=None,
                         ingroupCoverageFile=None,
                         minimumSizeToRescue=None,
                         minimumCoverageToRescue=None,
                         minimumNumberOfSpecies=None):
    """Runs cactus base aligner.
    """
    logLevel = getLogLevelString2(logLevel)
    maximumLength = nameValue("maximumLength", maximumLength, int)
    spanningTrees = nameValue("spanningTrees", spanningTrees, int)
    gapGamma = nameValue("gapGamma", gapGamma, float)
    matchGamma = nameValue("matchGamma", matchGamma, float)
    splitMatrixBiggerThanThis=nameValue("splitMatrixBiggerThanThis", splitMatrixBiggerThanThis, int)
    anchorMatrixBiggerThanThis=nameValue("anchorMatrixBiggerThanThis", anchorMatrixBiggerThanThis, int)
    repeatMaskMatrixBiggerThanThis=nameValue("repeatMaskMatrixBiggerThanThis", repeatMaskMatrixBiggerThanThis, int)                   
    diagonalExpansion=nameValue("diagonalExpansion", diagonalExpansion, int)
    constraintDiagonalTrim = nameValue("constraintDiagonalTrim", constraintDiagonalTrim, int)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    minimumIngroupDegree = nameValue("minimumIngroupDegree", minimumIngroupDegree, int)
    minimumOutgroupDegree = nameValue("minimumOutgroupDegree", minimumOutgroupDegree, int)
    pruneOutStubAlignments = nameValue("pruneOutStubAlignments", pruneOutStubAlignments, bool)
    alignAmbiguityCharacters = nameValue("alignAmbiguityCharacters", alignAmbiguityCharacters, bool)
    useProgressiveMerging=nameValue("useProgressiveMerging", useProgressiveMerging, bool)
    calculateWhichEndsToComputeSeparately=nameValue("calculateWhichEndsToComputeSeparately", calculateWhichEndsToComputeSeparately, bool)
    largeEndSize=nameValue("largeEndSize", largeEndSize, int)
    endAlignmentsToPrecomputeOutputFile=nameValue("endAlignmentsToPrecomputeOutputFile", endAlignmentsToPrecomputeOutputFile, str)
    precomputedAlignments=nameValue("precomputedAlignments", precomputedAlignments, str, quotes=True)
    ingroupCoverageFile = nameValue("ingroupCoverageFile", ingroupCoverageFile, str, quotes=True)
    minimumSizeToRescue = nameValue("minimumSizeToRescue", minimumSizeToRescue, int)
    minimumCoverageToRescue = nameValue("minimumCoverageToRescue", minimumCoverageToRescue, float)
    minimumNumberOfSpecies = nameValue("minimumNumberOfSpecies", minimumNumberOfSpecies, int)

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                                 option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                                 parameters=["cactus_bar",
                                             "--cactusSequencesPath", cactusSequencesPath,
                                             "--logLevel", logLevel, spanningTrees, maximumLength,
                                             gapGamma, matchGamma,
                                             splitMatrixBiggerThanThis, anchorMatrixBiggerThanThis,
                                             repeatMaskMatrixBiggerThanThis,
                                             constraintDiagonalTrim, minimumBlockDegree, minimumIngroupDegree,
                                             minimumOutgroupDegree,  
                                             alignAmbiguityCharacters, pruneOutStubAlignments,
                                             diagonalExpansion,
                                             useProgressiveMerging, calculateWhichEndsToComputeSeparately,
                                             largeEndSize, endAlignmentsToPrecomputeOutputFile,
                                             precomputedAlignments, ingroupCoverageFile,
                                             minimumSizeToRescue, minimumCoverageToRescue,
                                             minimumNumberOfSpecies])
        
    logger.info("Ran cactus_bar okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]

def runCactusSecondaryDatabase(secondaryDatabaseString, create=True):
    cactus_call(parameters=["cactus_secondaryDatabase",
                secondaryDatabaseString, create])
            
def runCactusReference(cactusDiskDatabaseString, cactusSequencesPath, flowerNames, logLevel=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=None,
                       theta=None,
                       phi=None, 
                       maxWalkForCalculatingZ=None,
                       ignoreUnalignedGaps=None,
                       wiggle=None, 
                       numberOfNs=None,
                       minNumberOfSequencesToSupportAdjacency=None,
                       makeScaffolds=None):
    """Runs cactus reference.
    """
    logLevel = getLogLevelString2(logLevel)
    matchingAlgorithm = nameValue("matchingAlgorithm", matchingAlgorithm)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    permutations = nameValue("permutations", permutations, int)
    useSimulatedAnnealing = nameValue("useSimulatedAnnealing", useSimulatedAnnealing, bool)
    theta = nameValue("theta", theta, float)
    phi = nameValue("phi", phi, float)
    maxWalkForCalculatingZ = nameValue("maxWalkForCalculatingZ", maxWalkForCalculatingZ, int)
    ignoreUnalignedGaps = nameValue("ignoreUnalignedGaps", ignoreUnalignedGaps, bool)
    wiggle = nameValue("wiggle", wiggle, float)
    numberOfNs = nameValue("numberOfNs", numberOfNs, int)
    minNumberOfSequencesToSupportAdjacency = nameValue("minNumberOfSequencesToSupportAdjacency", minNumberOfSequencesToSupportAdjacency, int)
    makeScaffolds = nameValue("makeScaffolds", makeScaffolds, bool)

    masterMessages = cactus_call(stdin_string=flowerNames, check_output=True,
                option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                parameters=["cactus_reference",
                            "--cactusSequencesPath", cactusSequencesPath,
                            "--logLevel", logLevel,
                            matchingAlgorithm, referenceEventString, permutations, 
                            useSimulatedAnnealing, theta, phi, maxWalkForCalculatingZ, ignoreUnalignedGaps,
                            wiggle, numberOfNs, minNumberOfSequencesToSupportAdjacency, makeScaffolds])
    logger.info("Ran cactus_reference okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusAddReferenceCoordinates(cactusDiskDatabaseString, cactusSequencesPath, flowerNames, logLevel=None, referenceEventString=None, outgroupEventString=None, secondaryDatabaseString=None, bottomUpPhase=None):   
    logLevel = getLogLevelString2(logLevel)
    bottomUpPhase = nameValue("bottomUpPhase", bottomUpPhase, bool)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    outgroupEventString = nameValue("outgroupEventString", outgroupEventString)
    secondaryDatabaseString = nameValue("secondaryDisk", secondaryDatabaseString, quotes=True)
    cactus_call(stdin_string=flowerNames,
                option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                parameters=["cactus_addReferenceCoordinates",
                            "--cactusSequencesPath", cactusSequencesPath,
                            secondaryDatabaseString,
                            "--logLevel", logLevel,
                            referenceEventString,
                            outgroupEventString,
                            bottomUpPhase])

def runCactusCheck(cactusDiskDatabaseString, 
                    cactusSequencesPath,
                    flowerNames=encodeFlowerNames((0,)), 
                    logLevel=None, 
                    recursive=None,
                    checkNormalised=None):
    logLevel = getLogLevelString2(logLevel)
    recursive = nameValue("recursive", recursive, bool)
    checkNormalised = nameValue("checkNormalised", checkNormalised, bool)
    cactus_call(stdin_string=flowerNames,
                parameters=["cactus_check",
                            "--cactusDisk", cactusDiskDatabaseString,
                            "--cactusSequencesPath", cactusSequencesPath,
                            "--logLevel", logLevel,
                            recursive, checkNormalised])
    logger.info("Ran cactus check")
    
def _fn(toilDir, 
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
      extraToilArgumentsString=""):
    logLevel = getLogLevelString2(logLevel)
    skipAlignments = nameValue("skipAlignments", skipAlignments, bool)
    buildAvgs = nameValue("buildAvgs", buildAvgs, bool)
    buildReference = nameValue("buildReference", buildReference, bool)
    buildHal = nameValue("buildHal", buildHal, bool)
    buildFasta = nameValue("buildFasta", buildFasta, bool)
    #Jobtree args
    batchSystem = nameValue("batchSystem", batchSystem, str, quotes=True)
    retryCount = nameValue("retryCount", retryCount, int)
    rescueJobFrequency = nameValue("rescueJobsFrequency", rescueJobFrequency, int)
    toilStats = nameValue("stats", toilStats, bool)
    maxThreads = nameValue("maxThreads", maxThreads, int)
    maxCpus = nameValue("maxCpus", maxCpus, int)
    defaultMemory= nameValue("defaultMemory", defaultMemory, int)
    logFile = nameValue("logFile", logFile, str)
    return "%s %s %s %s --logLevel %s %s %s %s %s %s %s %s %s %s %s %s" % (toilDir, skipAlignments, buildAvgs, 
             buildReference, logLevel, buildHal, buildFasta, batchSystem, retryCount, rescueJobFrequency, toilStats, maxThreads, maxCpus, logFile, defaultMemory, extraToilArgumentsString)
     
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
                      extraToilArgumentsString=""):
    arguments = ("--experiment %s" % experimentFile) + " " + _fn(toilDir,
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, buildReference, buildHal, buildFasta, toilStats, maxThreads, maxCpus, defaultMemory, logFile, extraToilArgumentsString=extraToilArgumentsString)

    import cactus.pipeline.cactus_workflow as cactus_workflow
    cactus_workflow.runCactusWorkflow(arguments.split())
    logger.info("Ran the cactus workflow okay")
    
def runCactusCreateMultiCactusProject(experimentFile, outputDir, 
                                      logLevel=None, fixNames=True,
                                      root=None):
    logLevel = getLogLevelString2(logLevel)
    root = nameValue("root", root, str, quotes=True)
    command = "cactus_createMultiCactusProject.py %s %s --fixNames=%s %s" % (experimentFile, outputDir, str(fixNames), root)
    system(command)
    logger.info("Ran the cactus create multi project")
    
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
                      defaultMemory=None,
                      recursive=None,
                      logFile=None,
                      event=None,
                      extraToilArgumentsString="",
                      profileFile=None):
    command = ("cactus_progressive.py --project %s" % inputDir) + " " + _fn(toilDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, None,
                      buildHal,
                      buildFasta,
                      toilStats, maxThreads, maxCpus, defaultMemory, logFile, extraToilArgumentsString=extraToilArgumentsString) + \
                      (" %s %s" % (nameValue("recursive", recursive, bool),
                                      nameValue("event", event)))
    if profileFile != None:
        command = "python -m cProfile -o %s %s/bin/%s" % (profileFile, cactusRootPath(), command)
    system(command)                   
    logger.info("Ran the cactus progressive okay")
    
def runCactusHalGenerator(cactusDiskDatabaseString,
                          cactusSequencesPath,
                          secondaryDatabaseString, 
                          flowerNames,
                          referenceEventString, 
                          outputFile=None,
                          showOnlySubstitutionsWithRespectToReference=None,
                          logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    if outputFile:
        outputFile = os.path.basename(outputFile)
    cactus_call(stdin_string=flowerNames,
                option_string="--cactusDisk '%s' --secondaryDisk '%s'" % (cactusDiskDatabaseString, secondaryDatabaseString),
                parameters=["cactus_halGenerator",
                            "--cactusSequencesPath", cactusSequencesPath,
                            "--logLevel", logLevel,
                            nameValue("referenceEventString", referenceEventString),
                            nameValue("outputFile", outputFile),
                            nameValue("showOnlySubstitutionsWithRespectToReference",
                                      showOnlySubstitutionsWithRespectToReference, bool)])
                            
def runCactusFastaGenerator(cactusDiskDatabaseString,
                          cactusSequencesPath,
                          flowerName,
                          outputFile,
                          referenceEventString=None, 
                          logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    cactus_call(option_string="--cactusDisk '%s'" % cactusDiskDatabaseString,
                parameters=["cactus_fastaGenerator",
                            "--cactusSequencesPath", cactusSequencesPath,
                            "--flowerName", flowerName,
                            "--outputFile", outputFile,
                            "--logLevel", logLevel,
                            nameValue("referenceEventString", referenceEventString)])
    
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
                            "--notrivial",
                            lastzArguments,
                            "%s[multiple][nameparse=darkspace]" % os.path.basename(seq1),
                            "%s[nameparse=darkspace]" % os.path.basename(seq2)])

def runSelfLastz(seq, alignmentsFile, lastzArguments, work_dir=None):
    work_dir = os.path.dirname(seq)
    cactus_call(work_dir=work_dir, outfile=alignmentsFile,
                parameters=["cPecanLastz",
                            "--format=cigar",
                            "--notrivial",
                            lastzArguments,
                            "%s[multiple][nameparse=darkspace]" % os.path.basename(seq),
                            "%s[nameparse=darkspace]" % os.path.basename(seq)])
    
def runCactusRealign(seq1, seq2, inputAlignmentsFile, outputAlignmentsFile, realignArguments, work_dir=None):
    cactus_call(infile=inputAlignmentsFile, outfile=outputAlignmentsFile, work_dir=work_dir,
                parameters=["cPecanRealign", realignArguments, seq1, seq2])

def runCactusSelfRealign(seq, inputAlignmentsFile, outputAlignmentsFile, realignArguments, work_dir=None):
    cactus_call(infile=inputAlignmentsFile, outfile=outputAlignmentsFile, work_dir=work_dir,
                parameters=["cPecanRealign", realignArguments, seq])

def runCactusCoverage(sequenceFile, alignmentsFile, work_dir=None):
    return cactus_call(check_output=True, work_dir=work_dir,
                parameters=["cactus_coverage", sequenceFile, alignmentsFile])

def runGetChunks(sequenceFiles, chunksDir, chunkSize, overlapSize, work_dir=None):
    return [chunk for chunk in cactus_call(work_dir=work_dir,
                                           check_output=True,
                                           parameters=["cactus_blast_chunkSequences",
                                           getLogLevelString(),
                                           chunkSize,
                                           overlapSize,
                                           chunksDir] + sequenceFiles).split("\n") if chunk != ""]

def pullCactusImage():
    """Ensure that the cactus Docker image is pulled."""
    dockerOrg = getDockerOrg()
    dockerTag = getDockerTag()
    image = "%s/cactus:%s" % (dockerOrg, dockerTag)
    call = ["docker", "pull", image]
    process = subprocess.Popen(call, stdout=subprocess.PIPE,
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

#TODO: This function is a mess
def cactus_call(tool=None,
                work_dir=None,
                parameters=None,
                rm=False,
                detached=True,
                check_output=False,
                container_name=None,
                mounts=None,
                infile=None,
                outfile=None,
                stdin_string=None,
                option_string="",
                server=False,
                shell=True,
                port=None,
                check_result=False,
                dockstore=None):
    if dockstore is None:
        dockstore = getDockerOrg()
    if parameters is None:
        parameters = []

    def moveToWorkDir(work_dir, arg):
        if isinstance(arg, str) and os.path.isfile(arg):
            if not os.path.dirname(arg) == work_dir:
                _log.info('Copying file %s to work dir' % arg)
                shutil.copy(arg, work_dir)

    if work_dir:
        for arg in parameters:
            moveToWorkDir(work_dir, arg)

    parameters = [str(par) for par in parameters]
    if not work_dir:
    #Make sure all the paths we're accessing are in the same directory
        files = [par for par in parameters if os.path.isfile(par)]
        folders = [par for par in parameters if os.path.isdir(par)]
        work_dirs = set([os.path.dirname(fileName) for fileName in files] + [os.path.dirname(os.path.dirname(folder)) for folder in folders])
        _log.info("Work dirs: %s" % work_dirs)
        assert len(work_dirs) <= 1
        if len(work_dirs) == 1:
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
        if os.path.isfile(path):
            _log.info('Found file %s' % path)
            return os.path.basename(path)
        if os.path.isdir(path):
            _log.info('Found dir %s' % path)
            return os.path.basename(os.path.dirname(path))
        else:
            _log.info('%s is neither a dir nor a file' % path)
            # Hack to relativize paths that are not provided as a
            # single argument (i.e. multiple paths that are
            # space-separated and quoted)
            if wd != '.':
                if not wd.endswith('/'):
                    wd = wd + '/'
                return path.replace(wd, '')
            else:
                return path

    parameters = [adjustPath(par, work_dir) for par in parameters]

    
    base_docker_call = ['docker', 'run',
                        '--interactive',
                        '--net=host',
                        '--log-driver=none',
                        '-e', 'ST_ABORT=1',
                        '-e', 'ST_ABORT_UNCAUGHT=1',
                        '-v', '{}:/data'.format(os.path.abspath(work_dir))]

    if port:
        base_docker_call += ["-p %d:%d" % (port, port)]

    #base_docker_call.extend(['--name', container_name])
    if rm:
        base_docker_call.append('--rm')


    parameters = [par for par in parameters if par != '']

    parameters = " ".join(parameters)

    if not tool:
        tool = "cactus"

    docker_tag = getDockerTag()

    if os.environ.get('CACTUS_DEVELOPER_MODE'):
        _log.info("Calling tool from local cactus installation.")
        call = parameters
    else:
        tool = "%s/%s:%s" % (dockstore, tool, docker_tag)
        call = " ".join(base_docker_call) + " " + tool + " " + parameters
    if option_string:
        call += " " + option_string
    

    if stdin_string:
        _log.info("Input string: %s" % stdin_string)


    stdinFileHandle = None
    stdoutFileHandle = None
    if stdin_string:
        stdinFileHandle = subprocess.PIPE
    elif infile:
        stdinFileHandle = open(infile, 'r')
    if outfile:
        stdoutFileHandle = open(outfile, 'w')
    if check_output:
        stdoutFileHandle = subprocess.PIPE


    _log.info("Running the command %s" % call)
    if not shell:
        call = call.split()
    process = subprocess.Popen(call, shell=shell,
                               stdin=stdinFileHandle, stdout=stdoutFileHandle, stderr=sys.stderr, bufsize=-1)

    if server:
        return process

    output, nothing = process.communicate(stdin_string)
    if check_result:
        return process.returncode
    
    if process.returncode != 0:
        raise RuntimeError("Command %s failed with output: %s" % (call, output))

    # Fix root ownership of output files
    if not os.getenv('CACTUS_DEVELOPER_MODE'):
        _fixPermissions(tool, work_dir)

    if check_output:
        return output


class RunAsFollowOn(Job):
    def __init__(self, job, *args, **kwargs):
        Job.__init__(self)
        self._args = args
        self._kwargs = kwargs
        self.job = job
    def run(self, fileStore):
        return self.addFollowOn(self.job(*self._args, **self._kwargs)).rv()
        
