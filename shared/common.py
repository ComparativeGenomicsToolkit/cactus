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

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete


def cactusRootPath():
    """
    function for finding external location
    """
    import cactus.shared.config
    i = os.path.abspath(cactus.shared.config.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getLogLevelString2(logLevelString):
    """Gets the log level string for the binary
    """
    if logLevelString == None:
        return getLogLevelString()
    return logLevelString

def formatFlowerNames(flowerNames):
    return "%i %s" % (len(flowerNames), " ".join([ str(i) for i in flowerNames ]))

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################  

def runCactusSetup(cactusDiskDatabaseString, sequences, 
                   newickTreeString, logLevel=None, outgroupEvents=None):
    logLevel = getLogLevelString2(logLevel)
    outgroupEvents = nameValue("outgroupEvents", outgroupEvents, str)
    system("cactus_setup %s --speciesTree '%s' --cactusDisk '%s' \
--logLevel %s %s" \
           % (" ".join(sequences), newickTreeString,
              cactusDiskDatabaseString, logLevel, outgroupEvents))
    logger.info("Ran cactus setup okay")
    
def runCactusAligner(cactusDiskDatabaseString, alignmentFile, tempDir, useDummy=True, flowerName=0, logLevel=None):        
    """Runs job tree and fails if not complete.
    """
    logLevel = getLogLevelString2(logLevel)
    tempDir = getTempDirectory(tempDir)
    jobTreeDir = os.path.join(tempDir, "jobTree")
    useDummy = nameValue("useDummy", useDummy, bool)
    command = "cactus_aligner.py --cactusDisk '%s' --flowerName %s \
--resultsFile %s %s --jobTree %s --logLevel %s" % (cactusDiskDatabaseString, flowerName, alignmentFile, useDummy, jobTreeDir, logLevel)
    system(command)
    runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
    system("rm -rf %s" % tempDir)
    logger.info("Ran the cactus aligner okay")
    
def runCactusBatch(sequenceFiles, outputFile, jobTreeDir,
                   chunkSize=None, overlapSize=None, chunksPerJob=None,
                   logLevel=None, 
                   blastString=None, 
                   selfBlastString=None,
                   compressFiles=None):
    
    logLevel = getLogLevelString2(logLevel)
    chunkSize = nameValue("chunkSize", chunkSize, int)
    overlapSize = nameValue("overlapSize", overlapSize, int)
    chunksPerJob = nameValue("chunksPerJob", chunksPerJob, int)
    blastString = nameValue("blastString", blastString, str)
    selfBlastString = nameValue("selfBlastString", selfBlastString, str)
    compressFiles = nameValue("compressFiles", compressFiles, bool)
    command = "cactus_batch.py %s  --cigars %s %s %s %s %s %s %s --jobTree %s --logLevel %s" % \
            (" ".join(sequenceFiles), outputFile,
             chunkSize, overlapSize, chunksPerJob, blastString, selfBlastString, compressFiles, jobTreeDir, logLevel)
    logger.info("Running command : %s" % command)
    system(command)
    logger.info("Ran the cactus_batch command okay")
            
def runCactusCore(cactusDiskDatabaseString, alignments, 
                  flowerNames=(0, 1),
                  logLevel=None, 
                  writeDebugFiles=False,
                  annealingRounds=None,
                  deannealingRounds=None,
                  alignRepeatsAtRound=False,
                  trim=None,
                  minimumTreeCoverage=None,
                  blockTrim=None,
                  minimumBlockDegree=None,
                  requiredIngroupFraction=None,
                  requiredOutgroupFraction=None,
                  requiredAllFraction=None,
                  singleCopyIngroup = None,
                  singleCopyOutgroup = None,
                  lastzArguments = None,
                  minimumSequenceLengthForBlast= None):
    logLevel = getLogLevelString2(logLevel)
    writeDebugFiles = nameValue("writeDebugFiles", writeDebugFiles, bool)
    alignRepeatsAtRound = nameValue("alignRepeatsAtRound", alignRepeatsAtRound, int)
    annealingRounds = nameValue("annealingRounds", annealingRounds, quotes=True)
    deannealingRounds = nameValue("deannealingRounds", deannealingRounds, quotes=True)
    trim = nameValue("trim", trim, quotes=True)
    alignments = nameValue("alignments", alignments)
    lastzArguments = nameValue("lastzArguments", lastzArguments, quotes=True)
    minimumTreeCoverage = nameValue("minimumTreeCoverage", minimumTreeCoverage, float)
    blockTrim = nameValue("blockTrim", blockTrim, int)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    minimumSequenceLengthForBlast = nameValue("minimumSequenceLengthForBlast", minimumSequenceLengthForBlast, int)
    
    requiredIngroupFraction = nameValue("requiredIngroupFraction", requiredIngroupFraction, float)
    requiredOutgroupFraction = nameValue("requiredOutgroupFraction", requiredOutgroupFraction, float)
    requiredAllFraction = nameValue("requiredAllFraction", requiredAllFraction, float)
    singleCopyIngroup = nameValue("singleCopyIngroup", singleCopyIngroup, bool)
    singleCopyOutgroup = nameValue("singleCopyOutgroup", singleCopyOutgroup, bool)
    
    command = "cactus_core --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (cactusDiskDatabaseString, logLevel, alignments, writeDebugFiles, annealingRounds, deannealingRounds, alignRepeatsAtRound,
     trim, minimumTreeCoverage, blockTrim, 
     minimumBlockDegree, requiredIngroupFraction, requiredOutgroupFraction, requiredAllFraction, 
     singleCopyIngroup, singleCopyOutgroup, lastzArguments, minimumSequenceLengthForBlast)
    popenPush(command, stdinString=formatFlowerNames(flowerNames))
    logger.info("Ran cactus_core okay")
    
def runCactusPhylogeny(cactusDiskDatabaseString,
                  flowerNames=(0, 1),
                  logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_phylogeny --cactusDisk '%s' --logLevel %s" % \
    (cactusDiskDatabaseString, logLevel)
    popenPush(command, stdinString=formatFlowerNames(flowerNames))
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(cactusDiskDatabaseString, flowerNames=(0, 1), logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_fillAdjacencies --cactusDisk '%s' --logLevel %s" %\
    (cactusDiskDatabaseString, logLevel)
    popenPush(command, stdinString=formatFlowerNames(flowerNames))
    logger.info("Ran cactus_fillAdjacencies OK")
    
def readFlowerNamesFile(flowerStrings):
    l = []    
    for line in flowerStrings.split("\n"):
        if line != '':
            line = line.split()
            firstFlowerName = int(line[0])
            totalFlowers = int(line[1])
            totalFlowerSizes = int(line[2])
            l.append((totalFlowerSizes, firstFlowerName, totalFlowers))
    #l.sort()
    #l.reverse()
    return l
    
def runCactusGetFlowers(cactusDiskDatabaseString, flowerNames, 
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, logLevel=None):
    """Gets a list of flowers attached to the given flower. 
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = popenCatch("cactus_workflow_getFlowers %s '%s' %i %i %i" % \
                               (logLevel, cactusDiskDatabaseString, int(minSequenceSizeOfFlower), -1, 
                                int(maxSequenceSizeOfFlowerGrouping)), stdinString=formatFlowerNames(flowerNames))
    l = readFlowerNamesFile(flowerStrings)
    return l

def runCactusExtendFlowers(cactusDiskDatabaseString, flowerNames, 
                        minSequenceSizeOfFlower=1,
                        maxSequenceSizeOfFlowerGrouping=-1, logLevel=None):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    logLevel = getLogLevelString2(logLevel)
    flowerStrings = popenCatch("cactus_workflow_extendFlowers %s '%s' %i %i %i" % \
                               (logLevel, cactusDiskDatabaseString, int(minSequenceSizeOfFlower), \
                                -1, int(maxSequenceSizeOfFlowerGrouping)), stdinString=formatFlowerNames(flowerNames))
    l = readFlowerNamesFile(flowerStrings)
    return l

def runCactusMakeNormal(cactusDiskDatabaseString, flowerNames, maxNumberOfChains=0, logLevel=None):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    logLevel = getLogLevelString2(logLevel)
    popenPush("cactus_normalisation --cactusDisk '%s' --maxNumberOfChains %i --logLevel %s" % (cactusDiskDatabaseString, maxNumberOfChains, logLevel), stdinString=formatFlowerNames(flowerNames))

def runCactusBaseAligner(cactusDiskDatabaseString, flowerNames, logLevel=None,
                         spanningTrees=None, maximumLength=None, 
                         gapGamma=None,
                         useBanding=None,
                         maxBandingSize=None,
                         minBandingSize=None,
                         minBandingConstraintDistance=None,
                         minTraceBackDiag=None,
                         minTraceGapDiags=None,
                         constraintDiagonalTrim=None,
                         minimumBlockDegree=None,
                         alignAmbiguityCharacters=None,
                         requiredIngroupFraction=None,
                         requiredOutgroupFraction=None,
                         requiredAllFraction=None,
                         pruneOutStubAlignments=None,
                         numThreads=None):
    """Runs cactus base aligner.
    """
    logLevel = getLogLevelString2(logLevel)
    maximumLength = nameValue("maximumLength", maximumLength, int)
    spanningTrees = nameValue("spanningTrees", spanningTrees, int)
    gapGamma = nameValue("gapGamma", gapGamma, float)
    useBanding = nameValue("useBanding", useBanding, bool)
    maxBandingSize = nameValue("maxBandingSize", maxBandingSize, int)
    minBandingSize = nameValue("minBandingSize", minBandingSize, int)
    minBandingConstraintDistance = nameValue("minBandingConstraintDistance", minBandingConstraintDistance, int)
    minTraceBackDiag = nameValue("minTraceBackDiag", minTraceBackDiag, int)
    minTraceGapDiags = nameValue("minTraceGapDiags", minTraceGapDiags, int)
    constraintDiagonalTrim = nameValue("constraintDiagonalTrim", constraintDiagonalTrim, int)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    pruneOutStubAlignments = nameValue("pruneOutStubAlignments", pruneOutStubAlignments, bool)
    alignAmbiguityCharacters = nameValue("alignAmbiguityCharacters", alignAmbiguityCharacters, bool)
    numThreads = nameValue("numThreads", numThreads, int)
    requiredIngroupFraction = nameValue("requiredIngroupFraction", requiredIngroupFraction, float)
    requiredOutgroupFraction = nameValue("requiredOutgroupFraction", requiredOutgroupFraction, float)
    requiredAllFraction = nameValue("requiredAllFraction", requiredAllFraction, float)
    
    popenPush("cactus_baseAligner --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % 
           (cactusDiskDatabaseString, logLevel, spanningTrees, maximumLength, gapGamma, 
            useBanding, maxBandingSize, minBandingSize, minBandingConstraintDistance, minTraceBackDiag, minTraceGapDiags, 
            constraintDiagonalTrim, minimumBlockDegree, alignAmbiguityCharacters, pruneOutStubAlignments,
            numThreads, requiredIngroupFraction, requiredOutgroupFraction, requiredAllFraction), stdinString=formatFlowerNames(flowerNames))
    
def runCactusReference(cactusDiskDatabaseString, flowerNames, logLevel=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=None,
                       theta=None,
                       maxNumberOfChainsBeforeSwitchingToFast=None):
    """Runs cactus reference.
    """
    logLevel = getLogLevelString2(logLevel)
    matchingAlgorithm = nameValue("matchingAlgorithm", matchingAlgorithm)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    permutations = nameValue("permutations", permutations, int)
    useSimulatedAnnealing = nameValue("useSimulatedAnnealing", useSimulatedAnnealing, bool)
    theta = nameValue("theta", theta, float)
    maxNumberOfChainsBeforeSwitchingToFast = nameValue("maxNumberOfChainsBeforeSwitchingToFast", maxNumberOfChainsBeforeSwitchingToFast, int)
    command = "cactus_reference --cactusDisk '%s' --logLevel %s %s %s %s %s %s" % (cactusDiskDatabaseString, logLevel, matchingAlgorithm, referenceEventString, permutations, useSimulatedAnnealing, theta)
    popenPush(command, stdinString=formatFlowerNames(flowerNames))
    
def runCactusAddReferenceCoordinates(cactusDiskDatabaseString, flowerNames, logLevel=None, referenceEventString=None, outgroupEventString=None, bottomUpPhase=None):   
    logLevel = getLogLevelString2(logLevel)
    bottomUpPhase = nameValue("bottomUpPhase", bottomUpPhase, bool)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    outgroupEventString = nameValue("outgroupEventString", outgroupEventString)
    command = "cactus_addReferenceCoordinates --cactusDisk '%s' --logLevel %s %s %s %s" % (cactusDiskDatabaseString, logLevel, referenceEventString, outgroupEventString, bottomUpPhase)
    popenPush(command, stdinString=formatFlowerNames(flowerNames))

def runCactusCheck(cactusDiskDatabaseString, 
                    flowerNames=(0, 1), 
                    logLevel=None, 
                    recursive=None,
                    checkNormalised=None):
    logLevel = getLogLevelString2(logLevel)
    recursive = nameValue("recursive", recursive, bool)
    checkNormalised = nameValue("checkNormalised", checkNormalised, bool)
    popenPush("cactus_check --cactusDisk '%s' --logLevel %s %s %s"  % (cactusDiskDatabaseString, logLevel, recursive, checkNormalised), stdinString=formatFlowerNames(flowerNames))
    logger.info("Ran cactus check")
    
def _fn(jobTreeDir, 
      logLevel=None, retryCount=0, 
      batchSystem="single_machine", 
      rescueJobFrequency=None,
      setupAndBuildAlignments=True,
      buildTrees=True, buildFaces=True, buildReference=True,
      jobTreeStats=False,
      maxThreads=None,
      maxJobs=None,
      logFile=None):
    logLevel = getLogLevelString2(logLevel)
    buildFaces=False
    setupAndBuildAlignments = nameValue("setupAndBuildAlignments", setupAndBuildAlignments, bool)
    buildTrees = nameValue("buildTrees", buildTrees, bool)
    buildFaces = nameValue("buildFaces", buildFaces, bool)
    buildReference = nameValue("buildReference", buildReference, bool)
    #Jobtree args
    batchSystem = nameValue("batchSystem", batchSystem, str)
    retryCount = nameValue("retryCount", retryCount, int)
    rescueJobFrequency = nameValue("rescueJobsFrequency", rescueJobFrequency, int)
    jobTreeStats = nameValue("stats", jobTreeStats, bool)
    maxThreads = nameValue("maxThreads", maxThreads, int)
    maxJobs = nameValue("maxJobs", maxJobs, int)
    logFile = nameValue("logFile", logFile, str)
    return "%s %s %s %s --jobTree %s --logLevel %s %s %s %s %s %s %s %s" % (setupAndBuildAlignments, buildTrees, buildFaces, 
             buildReference, jobTreeDir, logLevel, batchSystem, retryCount, rescueJobFrequency, jobTreeStats, maxThreads, maxJobs, logFile)
     
def runCactusWorkflow(experimentFile,
                      jobTreeDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      setupAndBuildAlignments=True,
                      buildTrees=True, buildFaces=True, buildReference=True,
                      jobTreeStats=False,
                      maxThreads=None,
                      maxJobs=None,
                      logFile=None):
    command = ("cactus_workflow.py --experiment %s" % experimentFile) + " " + _fn(jobTreeDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, setupAndBuildAlignments,
                      buildTrees, buildFaces, buildReference, jobTreeStats,maxThreads,maxJobs,logFile)
    system(command)
    logger.info("Ran the cactus workflow okay")
    
def runCactusCreateMultiCactusProject(experimentFile, outputDir, 
                                      logLevel=None, fixNames=True):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_createMultiCactusProject.py %s %s --fixNames=%s" % (experimentFile, outputDir, str(fixNames))
    system(command)
    logger.info("Ran the cactus create multi project")
    
def runCactusProgressive(inputDir,
                      jobTreeDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      setupAndBuildAlignments=True,
                      buildMaf=None,
                      joinMaf=None,
                      #buildTrees=True, buildFaces=True, buildReference=True,
                      jobTreeStats=False,
                      maxThreads=None,
                      maxJobs=None,
                      recursive=None,
                      logFile=None,
                      event=None):
    command = ("cactus_progressive.py %s" % inputDir) + " " + _fn(jobTreeDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, setupAndBuildAlignments,
                      None, None, None, #buildTrees, buildFaces, buildReference, 
                      jobTreeStats,maxThreads, maxJobs, logFile) + \
                      (" %s %s %s %s" % (nameValue("recursive", recursive, bool),
                                     nameValue("buildMAF", buildMaf, bool),
                                     nameValue("joinMAF", joinMaf, bool), nameValue("event", event)))
    system(command)
    logger.info("Ran the cactus progressive okay")
    
def runCactusRecursiveMafGenerator(cactusDiskDatabaseString,
                          flowerNames,
                          referenceEventString, 
                          childDir, 
                          parentDir=None, outputFile=None,
                          showOnlySubstitutionsWithRespectToReference=None,
                          logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    print "I am running:", formatFlowerNames(flowerNames)
    popenPush("cactus_recursiveMafGenerator --cactusDisk '%s' --logLevel %s %s %s %s %s %s" % 
           (cactusDiskDatabaseString, logLevel, 
            nameValue("referenceEventString", referenceEventString),
            nameValue("childDir", childDir),
            nameValue("parentDir", parentDir),
            nameValue("outputFile", outputFile),
            nameValue("showOnlySubstitutionsWithRespectToReference", showOnlySubstitutionsWithRespectToReference, bool)), stdinString=formatFlowerNames(flowerNames))
    