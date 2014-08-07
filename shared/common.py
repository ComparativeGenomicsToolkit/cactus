#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus package.
"""

import os
import random
import sys

from sonLib.bioio import logger
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system, popenCatch, popenPush
from sonLib.bioio import nameValue
from sonLib.bioio import getLogLevelString

from jobTree.src.common import runJobTreeStatusAndFailIfNotComplete

def cactusRootPath():
    """
    function for finding external location
    """
    import cactus.shared.test
    i = os.path.abspath(cactus.shared.test.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

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
    
def runCactusGetFlowers(cactusDiskDatabaseString, flowerNames, 
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Gets a list of flowers attached to the given flower. 
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = popenCatch("cactus_workflow_getFlowers %s '%s' %i %i %i" % \
                               (logLevel, cactusDiskDatabaseString, int(minSequenceSizeOfFlower), 
                                int(maxSequenceSizeOfFlowerGrouping), 
                                int(maxSequenceSizeOfSecondaryFlowerGrouping)), 
                                stdinString=flowerNames)
    l = readFlowerNames(flowerStrings)
    return l

def runCactusExtendFlowers(cactusDiskDatabaseString, flowerNames, 
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, 
                        maxSequenceSizeOfSecondaryFlowerGrouping=-1, 
                        logLevel=None):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = popenCatch("cactus_workflow_extendFlowers %s '%s' %i %i %i" % \
                               (logLevel, cactusDiskDatabaseString, int(minSequenceSizeOfFlower), \
                                int(maxSequenceSizeOfFlowerGrouping), 
                                int(maxSequenceSizeOfSecondaryFlowerGrouping)), 
                                stdinString=flowerNames)
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
                   makeEventHeadersAlphaNumeric=None):
    logLevel = getLogLevelString2(logLevel)
    outgroupEvents = nameValue("outgroupEvents", outgroupEvents, str, quotes=True)
    makeEventHeadersAlphaNumeric=nameValue("makeEventHeadersAlphaNumeric", makeEventHeadersAlphaNumeric, bool)
    masterMessages = popenCatch("cactus_setup %s --speciesTree '%s' --cactusDisk '%s' \
--logLevel %s %s %s" \
           % (" ".join(sequences), newickTreeString,
              cactusDiskDatabaseString, logLevel, outgroupEvents, makeEventHeadersAlphaNumeric))
    logger.info("Ran cactus setup okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusBlast(sequenceFiles, outputFile, jobTreeDir,
                   chunkSize=None, overlapSize=None, 
                   logLevel=None, 
                   blastString=None, 
                   selfBlastString=None,
                   compressFiles=None,
                   lastzMemory=None,
                   targetSequenceFiles=None):
    logLevel = getLogLevelString2(logLevel)
    chunkSize = nameValue("chunkSize", chunkSize, int)
    overlapSize = nameValue("overlapSize", overlapSize, int)
    blastString = nameValue("blastString", blastString, str)
    selfBlastString = nameValue("selfBlastString", selfBlastString, str)
    compressFiles = nameValue("compressFiles", compressFiles, bool)
    lastzMemory = nameValue("lastzMemory", lastzMemory, int)
    if targetSequenceFiles != None: 
        targetSequenceFiles = " ".join(targetSequenceFiles)
    targetSequenceFiles = nameValue("targetSequenceFiles", targetSequenceFiles, quotes=True)
    command = "cactus_blast.py %s  --cigars %s %s %s %s %s %s %s %s --jobTree %s --logLevel %s" % \
            (" ".join(sequenceFiles), outputFile,
             chunkSize, overlapSize, blastString, selfBlastString, compressFiles, 
             lastzMemory, targetSequenceFiles, jobTreeDir, logLevel)
    logger.info("Running command : %s" % command)
    system(command)
    logger.info("Ran the cactus_blast command okay")

def runConvertAlignmentsToInternalNames(cactusDiskString, alignmentsFile, outputFile, flowerName):
    popenCatch("cactus_convertAlignmentsToInternalNames --cactusDisk '%s' %s %s" % (cactusDiskString, alignmentsFile, outputFile), stdinString=encodeFlowerNames((flowerName,)))

def runStripUniqueIDs(cactusDiskString):
    system("cactus_stripUniqueIDs --cactusDisk '%s'" % cactusDiskString)

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
                  singleCopyIngroup=None,
                  singleCopyOutgroup=None,
                  lastzArguments=None,
                  minimumSequenceLengthForBlast=None,
                  maxAdjacencyComponentSizeRatio=None,
                  constraints=None,
                  minLengthForChromosome=None,
                  proportionOfUnalignedBasesForNewChromosome=None, 
                  maximumMedianSequenceLengthBetweenLinkedEnds=None,
                  realign=None,
                  realignArguments=None):
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
    singleCopyIngroup = nameValue("singleCopyIngroup", singleCopyIngroup, bool)
    singleCopyOutgroup = nameValue("singleCopyOutgroup", singleCopyOutgroup, bool)
    maxAdjacencyComponentSizeRatio = nameValue("maxAdjacencyComponentSizeRatio", maxAdjacencyComponentSizeRatio, float)
    constraints = nameValue("constraints", constraints)
    realign = nameValue("realign", realign, bool)
    realignArguments = nameValue("realignArguments", realignArguments, quotes=True)
    
    minLengthForChromosome = nameValue("minLengthForChromosome", minLengthForChromosome, int)
    proportionOfUnalignedBasesForNewChromosome = nameValue("proportionOfUnalignedBasesForNewChromosome", proportionOfUnalignedBasesForNewChromosome, float)
    maximumMedianSequenceLengthBetweenLinkedEnds = nameValue("maximumMedianSequenceLengthBetweenLinkedEnds", maximumMedianSequenceLengthBetweenLinkedEnds, int)
    
    command = "cactus_caf --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (cactusDiskDatabaseString, logLevel, alignments, annealingRounds, deannealingRounds, 
     trim, minimumTreeCoverage, blockTrim, 
     minimumBlockDegree, minimumIngroupDegree, minimumOutgroupDegree,  
     singleCopyIngroup, singleCopyOutgroup, lastzArguments, minimumSequenceLengthForBlast, maxAdjacencyComponentSizeRatio, constraints,
     minLengthForChromosome, proportionOfUnalignedBasesForNewChromosome, maximumMedianSequenceLengthBetweenLinkedEnds, realign, realignArguments)
    masterMessages = popenCatch(command, stdinString=flowerNames)
    logger.info("Ran cactus_core okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]
    
def runCactusPhylogeny(cactusDiskDatabaseString,
                  flowerNames=encodeFlowerNames((0,)),
                  logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_phylogeny --cactusDisk '%s' --logLevel %s" % \
    (cactusDiskDatabaseString, logLevel)
    popenPush(command, stdinString=flowerNames)
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(cactusDiskDatabaseString, flowerNames=encodeFlowerNames((0,)), logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_fillAdjacencies --cactusDisk '%s' --logLevel %s" %\
    (cactusDiskDatabaseString, logLevel)
    popenPush(command, stdinString=flowerNames)
    logger.info("Ran cactus_fillAdjacencies OK")

def runCactusConvertAlignmentToCactus(cactusDiskDatabaseString, constraintsFile, newConstraintsFile, logLevel=None):
    """Takes a cigar file and makes an equivalent cigar file using the internal coordinate system format of cactus.
    """
    logLevel = getLogLevelString2(logLevel)
    system("cactus_workflow_convertAlignmentCoordinates %s '%s' %s %s" % \
                               (logLevel, cactusDiskDatabaseString, constraintsFile, newConstraintsFile))

def runCactusFlowerStats(cactusDiskDatabaseString, flowerName, logLevel=None):
    """Prints stats for the given flower
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStatsString = popenCatch("cactus_workflow_flowerStats %s '%s' %s" % 
                              (logLevel, cactusDiskDatabaseString, flowerName))
    return flowerStatsString.split("\n")[0]

def runCactusMakeNormal(cactusDiskDatabaseString, flowerNames, maxNumberOfChains=0, logLevel=None):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    logLevel = getLogLevelString2(logLevel)
    popenPush("cactus_normalisation --cactusDisk '%s' --maxNumberOfChains %i --logLevel %s" % (cactusDiskDatabaseString, maxNumberOfChains, logLevel), stdinString=flowerNames)

def runCactusBar(cactusDiskDatabaseString, flowerNames, logLevel=None,
                         spanningTrees=None, maximumLength=None, 
                         gapGamma=None,
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
                         maximumNumberOfSequencesBeforeSwitchingToFast=None,
                         calculateWhichEndsToComputeSeparately=None,
                         largeEndSize=None,
                         endAlignmentsToPrecomputeOutputFile=None,
                         precomputedAlignments=None):
    """Runs cactus base aligner.
    """
    logLevel = getLogLevelString2(logLevel)
    maximumLength = nameValue("maximumLength", maximumLength, int)
    spanningTrees = nameValue("spanningTrees", spanningTrees, int)
    gapGamma = nameValue("gapGamma", gapGamma, float)
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
    maximumNumberOfSequencesBeforeSwitchingToFast=nameValue("maximumNumberOfSequencesBeforeSwitchingToFast", maximumNumberOfSequencesBeforeSwitchingToFast, int)
    calculateWhichEndsToComputeSeparately=nameValue("calculateWhichEndsToComputeSeparately", calculateWhichEndsToComputeSeparately, bool)
    largeEndSize=nameValue("largeEndSize", largeEndSize, int)
    endAlignmentsToPrecomputeOutputFile=nameValue("endAlignmentsToPrecomputeOutputFile", endAlignmentsToPrecomputeOutputFile, str)
    precomputedAlignments=nameValue("precomputedAlignments", precomputedAlignments, str, quotes=True)
    
    masterMessages = popenCatch("cactus_bar --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % 
           (cactusDiskDatabaseString, logLevel, spanningTrees, maximumLength, gapGamma, 
            splitMatrixBiggerThanThis, anchorMatrixBiggerThanThis, repeatMaskMatrixBiggerThanThis,
            constraintDiagonalTrim, minimumBlockDegree, minimumIngroupDegree, minimumOutgroupDegree,  
            alignAmbiguityCharacters, pruneOutStubAlignments, diagonalExpansion,
            maximumNumberOfSequencesBeforeSwitchingToFast, calculateWhichEndsToComputeSeparately,
            largeEndSize, endAlignmentsToPrecomputeOutputFile, precomputedAlignments), stdinString=flowerNames)
    logger.info("Ran cactus_bar okay")
    return [ i for i in masterMessages.split("\n") if i != '' ]

def runCactusExpectationMaximisation(sequenceFiles, alignmentsFile, outputModelFile,
                                     inputModelFile=None, 
                                     jobTreeDir=None,
                                     iterations=None, randomStart=None, logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    jobTreeDir= nameValue("jobTree", jobTreeDir, str)
    inputModelFile= nameValue("inputModel", inputModelFile, str)
    iterations = nameValue("iterations", iterations, int)
    randomStart = nameValue("randomStart", randomStart, bool)
    system("cactus_expectationMaximisation --sequences '%s' --alignments %s --outputModel %s %s %s %s %s" % \
           (" ".join(sequenceFiles), alignmentsFile, outputModelFile, iterations, randomStart, jobTreeDir, inputModelFile))    

def runCactusSecondaryDatabase(secondaryDatabaseString, create=True):
    command = "cactus_secondaryDatabase '%s' %s" % (secondaryDatabaseString, int(create))
    system(command)
    
def runCactusReference(cactusDiskDatabaseString, flowerNames, logLevel=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=None,
                       theta=None,
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
    maxWalkForCalculatingZ = nameValue("maxWalkForCalculatingZ", maxWalkForCalculatingZ, int)
    ignoreUnalignedGaps = nameValue("ignoreUnalignedGaps", ignoreUnalignedGaps, bool)
    wiggle = nameValue("wiggle", wiggle, float)
    numberOfNs = nameValue("numberOfNs", numberOfNs, int)
    minNumberOfSequencesToSupportAdjacency = nameValue("minNumberOfSequencesToSupportAdjacency", minNumberOfSequencesToSupportAdjacency, int)
    makeScaffolds = nameValue("makeScaffolds", makeScaffolds, bool)
    command = "cactus_reference --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (cactusDiskDatabaseString, logLevel, matchingAlgorithm, referenceEventString, permutations, 
     useSimulatedAnnealing, theta, maxWalkForCalculatingZ, ignoreUnalignedGaps, wiggle, numberOfNs, minNumberOfSequencesToSupportAdjacency, makeScaffolds)
    popenPush(command, stdinString=flowerNames)
    
def runCactusAddReferenceCoordinates(cactusDiskDatabaseString, flowerNames, logLevel=None, referenceEventString=None, outgroupEventString=None, secondaryDatabaseString=None, bottomUpPhase=None):   
    logLevel = getLogLevelString2(logLevel)
    bottomUpPhase = nameValue("bottomUpPhase", bottomUpPhase, bool)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    outgroupEventString = nameValue("outgroupEventString", outgroupEventString)
    secondaryDatabaseString = nameValue("secondaryDisk", secondaryDatabaseString, quotes=True)
    command = "cactus_addReferenceCoordinates --cactusDisk '%s' %s --logLevel %s %s %s %s" % (cactusDiskDatabaseString, secondaryDatabaseString, logLevel, referenceEventString, outgroupEventString, bottomUpPhase)
    popenPush(command, stdinString=flowerNames)

def runCactusCheck(cactusDiskDatabaseString, 
                    flowerNames=encodeFlowerNames((0,)), 
                    logLevel=None, 
                    recursive=None,
                    checkNormalised=None):
    logLevel = getLogLevelString2(logLevel)
    recursive = nameValue("recursive", recursive, bool)
    checkNormalised = nameValue("checkNormalised", checkNormalised, bool)
    popenPush("cactus_check --cactusDisk '%s' --logLevel %s %s %s"  % (cactusDiskDatabaseString, logLevel, recursive, checkNormalised), stdinString=flowerNames)
    logger.info("Ran cactus check")
    
def _fn(jobTreeDir, 
      logLevel=None, retryCount=0, 
      batchSystem="single_machine", 
      rescueJobFrequency=None,
      skipAlignments=False,
      buildAvgs=False, buildReference=False,
      buildHal=False,
      buildFasta=False,
      jobTreeStats=False,
      maxThreads=None,
      maxCpus=None,
      defaultMemory=None,
      logFile=None,
      extraJobTreeArgumentsString=""):
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
    jobTreeStats = nameValue("stats", jobTreeStats, bool)
    maxThreads = nameValue("maxThreads", maxThreads, int)
    maxCpus = nameValue("maxCpus", maxCpus, int)
    defaultMemory= nameValue("defaultMemory", defaultMemory, int)
    logFile = nameValue("logFile", logFile, str)
    return "%s %s %s --jobTree %s --logLevel %s %s %s %s %s %s %s %s %s %s %s %s" % (skipAlignments, buildAvgs, 
             buildReference, jobTreeDir, logLevel, buildHal, buildFasta, batchSystem, retryCount, rescueJobFrequency, jobTreeStats, maxThreads, maxCpus, logFile, defaultMemory, extraJobTreeArgumentsString)
     
def runCactusWorkflow(experimentFile,
                      jobTreeDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildAvgs=False, buildReference=False,
                      buildHal=False,
                      buildFasta=False,
                      jobTreeStats=False,
                      maxThreads=None,
                      maxCpus=None,
                      defaultMemory=None,
                      logFile=None,
                      extraJobTreeArgumentsString=""):
    command = ("cactus_workflow.py --experiment %s" % experimentFile) + " " + _fn(jobTreeDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, buildReference, buildHal, buildFasta, jobTreeStats, maxThreads, maxCpus, defaultMemory, logFile, extraJobTreeArgumentsString=extraJobTreeArgumentsString)
    system(command)
    logger.info("Ran the cactus workflow okay")
    
def runCactusCreateMultiCactusProject(experimentFile, outputDir, 
                                      logLevel=None, fixNames=True,
                                      rootOutgroupPath=None, rootOutgroupDist=None):
    logLevel = getLogLevelString2(logLevel)
    rootOutgroupPath = nameValue("rootOutgroupPath", rootOutgroupPath, str)
    rootOutgroupDist = nameValue("rootOutgroupDist", rootOutgroupDist, float)
    command = "cactus_createMultiCactusProject.py %s %s --fixNames=%s %s %s" % (experimentFile, outputDir, str(fixNames), rootOutgroupPath, rootOutgroupDist)
    system(command)
    logger.info("Ran the cactus create multi project")
    
def runCactusProgressive(inputDir,
                      jobTreeDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildHal=None,
                      buildFasta=None,
                      buildAvgs=False, 
                      jobTreeStats=False,
                      maxThreads=None,
                      maxCpus=None,
                      defaultMemory=None,
                      recursive=None,
                      logFile=None,
                      event=None,
                      extraJobTreeArgumentsString="",
                      profileFile=None):
    command = ("cactus_progressive.py %s" % inputDir) + " " + _fn(jobTreeDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, None,
                      buildHal,
                      buildFasta,
                      jobTreeStats, maxThreads, maxCpus, defaultMemory, logFile, extraJobTreeArgumentsString=extraJobTreeArgumentsString) + \
                      (" %s %s" % (nameValue("recursive", recursive, bool),
                                      nameValue("event", event)))
    if profileFile != None:
        command = "python -m cProfile -o %s %s/bin/%s" % (profileFile, cactusRootPath(), command)
    system(command)                   
    logger.info("Ran the cactus progressive okay")
    
def runCactusHalGenerator(cactusDiskDatabaseString,
                          secondaryDatabaseString, 
                          flowerNames,
                          referenceEventString, 
                          outputFile=None,
                          showOnlySubstitutionsWithRespectToReference=None,
                          logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    popenPush("cactus_halGenerator --cactusDisk '%s' --secondaryDisk '%s' --logLevel %s %s %s %s" % 
           (cactusDiskDatabaseString, secondaryDatabaseString, logLevel, 
            nameValue("referenceEventString", referenceEventString),
            nameValue("outputFile", outputFile),
            nameValue("showOnlySubstitutionsWithRespectToReference", 
                      showOnlySubstitutionsWithRespectToReference, bool)), 
              stdinString=flowerNames)
    
def runCactusFastaGenerator(cactusDiskDatabaseString,
                          flowerName,
                          outputFile,
                          referenceEventString=None, 
                          logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    system("cactus_fastaGenerator --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s %s" % 
           (cactusDiskDatabaseString, flowerName, outputFile, logLevel, 
            nameValue("referenceEventString", referenceEventString)))
    
def runCactusAnalyseAssembly(sequenceFile):
    return popenCatch("cactus_analyseAssembly %s" % sequenceFile)[:-1]
