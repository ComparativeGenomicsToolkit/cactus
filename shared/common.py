#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus package.
"""

import os
import random

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system
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

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################  

def runCactusSetup(cactusDiskDatabaseString, sequences, 
                   newickTreeString, logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    system("cactus_setup %s --speciesTree '%s' --cactusDisk '%s' \
--logLevel %s" \
           % (" ".join(sequences), newickTreeString,
              cactusDiskDatabaseString, logLevel))
    logger.info("Ran cactus setup okay")
    
def runCactusAligner(cactusDiskDatabaseString, alignmentFile, tempDir, useDummy=True, flowerName="0", logLevel=None):        
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
            
def runCactusCore(cactusDiskDatabaseString, alignmentFile, 
                  flowerName=0,
                  logLevel=None, 
                  writeDebugFiles=False,
                  annealingRounds=None,
                  deannealingRounds=None,
                  minimumChainLength=None,
                  maximumGroupSize=None,
                  alignRepeatsAtRound=False,
                  trim=None,
                  minimumTreeCoverage=None,
                  blockTrim=None,
                  ignoreAllChainsLessThanMinimumTreeCoverage=False,
                  minimumBlockDegree=None,
                  requiredSpecies=None,
                  singleCopySpecies=None):
    logLevel = getLogLevelString2(logLevel)
    writeDebugFiles = nameValue("writeDebugFiles", writeDebugFiles, bool)
    if annealingRounds != None:
        annealingRounds = "--annealingRounds '%s'" % " ".join([ str(i) for i in annealingRounds ])
    if deannealingRounds != None:
        deannealingRounds = "--deannealingRounds '%s'" % " ".join([ str(i) for i in deannealingRounds ])
    alignRepeatsAtRound = nameValue("alignRepeatsAtRound", alignRepeatsAtRound, int)
    if trim != None:
        trim = "--trim '%s'" % " ".join([ str(i) for i in trim ])
    else:
        trim = ""
    minimumTreeCoverage = nameValue("minimumTreeCoverage", minimumTreeCoverage, float)
    blockTrim = nameValue("blockTrim", blockTrim, int)
    ignoreAllChainsLessThanMinimumTreeCoverage = nameValue("ignoreAllChainsLessThanMinimumTreeCoverage", ignoreAllChainsLessThanMinimumTreeCoverage, bool)
    minimumBlockDegree = nameValue("minimumDegree", minimumBlockDegree, int)
    minimumChainLength = nameValue("minimumChainLength", minimumChainLength, int)
    maximumGroupSize = nameValue("maximumGroupSize", maximumGroupSize, int)
    if requiredSpecies != None:
        requiredSpecies = "--requiredSpecies '%s'" % requiredSpecies
    else:
        requiredSpecies = ""
    if singleCopySpecies != None:
        singleCopySpecies = "--singleCopySpecies '%s'" % singleCopySpecies
    else:
        singleCopySpecies = ""
    command = "cactus_core --cactusDisk '%s' --flowerName %s --alignments %s --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (cactusDiskDatabaseString, flowerName, alignmentFile, logLevel, writeDebugFiles, annealingRounds, deannealingRounds, alignRepeatsAtRound,
     trim, minimumTreeCoverage, blockTrim, ignoreAllChainsLessThanMinimumTreeCoverage,
     minimumBlockDegree, requiredSpecies, singleCopySpecies, minimumChainLength, maximumGroupSize)
    #print "command to run", command
    #assert 0
    system(command)
    logger.info("Ran cactus_core okay")
    
def runCactusPhylogeny(cactusDiskDatabaseString,
                  flowerNames=("0",),
                  logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_phylogeny --cactusDisk '%s' --logLevel %s %s" % \
    (cactusDiskDatabaseString, logLevel, " ".join(flowerNames))
    system(command)
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(cactusDiskDatabaseString, flowerNames=("0",), logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    command = "cactus_fillAdjacencies --cactusDisk '%s' --logLevel %s %s" %\
    (cactusDiskDatabaseString, logLevel, " ".join(flowerNames))
    system(command)
    logger.info("Ran cactus_fillAdjacencies OK")
    
def readFlowerNamesFile(flowerNamesFile):
    fileHandle = open(flowerNamesFile, 'r')
    line = fileHandle.readline()
    l = []
    while line != '':
        childFlowerName = line.split()[0]
        childFlowerSize = float(line.split()[1])
        l.append((childFlowerName, childFlowerSize))
        line = fileHandle.readline()
    fileHandle.close()
    return l
    
def runCactusGetFlowers(cactusDiskDatabaseString, flowerNames, tempDir, includeTerminalFlowers=True, logLevel=None):
    """Gets a list of flowers attached to the given flower. 
    """
    logLevel = getLogLevelString2(logLevel)
    flowerNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_getFlowers '%s' %s %s %s %s" % (cactusDiskDatabaseString,  flowerNamesFile, int(includeTerminalFlowers), logLevel, " ".join(flowerNames)))
    l = readFlowerNamesFile(flowerNamesFile)
    os.remove(flowerNamesFile)
    random.shuffle(l) #We shuffle the flowers so we don't end up with an ordering that places all the large problems together.
    return l

def runCactusExtendFlowers(cactusDiskDatabaseString, flowerNames, tempDir,
                        minSizeToExtend=1, logLevel=None):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    logLevel = getLogLevelString2(logLevel)
    flowerNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_extendFlowers %s '%s' %s %i %s" % (logLevel, cactusDiskDatabaseString, flowerNamesFile, int(minSizeToExtend), " ".join(flowerNames)))
    l = readFlowerNamesFile(flowerNamesFile)
    os.remove(flowerNamesFile)
    random.shuffle(l) #We shuffle the flowers so we don't end up with an ordering that places all the large problems together.
    return l

def runCactusMakeNormal(cactusDiskDatabaseString, flowerNames, maxNumberOfChains=0, logLevel=None):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    logLevel = getLogLevelString2(logLevel)
    system("cactus_normalisation --cactusDisk '%s' --maxNumberOfChains %i --logLevel %s %s" % (cactusDiskDatabaseString, maxNumberOfChains, logLevel, " ".join(flowerNames)))

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
                         requiredSpecies=None):
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
    if requiredSpecies != None:
        requiredSpecies = "'%s'" % requiredSpecies
    requiredSpecies = nameValue("requiredSpecies", requiredSpecies, str)
    alignAmbiguityCharacters = nameValue("alignAmbiguityCharacters", alignAmbiguityCharacters, bool)
    
    system("cactus_baseAligner --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % 
           (cactusDiskDatabaseString, logLevel, " ".join(flowerNames), spanningTrees, maximumLength, gapGamma, 
            useBanding, maxBandingSize, minBandingSize, minBandingConstraintDistance, minTraceBackDiag, minTraceGapDiags, 
            constraintDiagonalTrim, minimumBlockDegree, requiredSpecies, alignAmbiguityCharacters))
    
def runCactusReference(cactusDiskDatabaseString, flowerNames, logLevel=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=None,
                       theta=None):
    """Runs cactus reference.
    """
    logLevel = getLogLevelString2(logLevel)
    matchingAlgorithm = nameValue("matchingAlgorithm", matchingAlgorithm)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    permutations = nameValue("permutations", permutations, int)
    useSimulatedAnnealing = nameValue("useSimulatedAnnealing", useSimulatedAnnealing, bool)
    theta = nameValue("theta", theta, float)
    command = "cactus_reference --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s" % (cactusDiskDatabaseString, logLevel, matchingAlgorithm, referenceEventString, permutations, useSimulatedAnnealing, theta, " ".join(flowerNames))
    system(command)
    
def runCactusAddReferenceCoordinates(cactusDiskDatabaseString, logLevel=None, referenceEventString=None):   
    logLevel = getLogLevelString2(logLevel)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    command = "cactus_addReferenceCoordinates --cactusDisk '%s' --logLevel %s %s" % (cactusDiskDatabaseString, logLevel, referenceEventString)
    system(command)

def runCactusCheck(cactusDiskDatabaseString, 
                    flowerNames=("0",), 
                    logLevel=None, 
                    recursive=False):
    logLevel = getLogLevelString2(logLevel)
    recursive = nameValue("recursive", recursive, bool)
    system("cactus_check --cactusDisk '%s' %s --logLevel %s %s"  % (cactusDiskDatabaseString, " ".join(flowerNames), logLevel, recursive))
    logger.info("Ran cactus check")
    
def runCactusWorkflow(experimentFile,
                      jobTreeDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      setupAndBuildAlignments=True,
                      buildTrees=True, buildFaces=True, buildReference=True,
                      jobTreeStats=False,
                      maxThreads=None):
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
    
    command = "cactus_workflow.py --experiment %s %s %s %s %s --jobTree %s --logLevel %s %s %s %s %s %s" % \
            (experimentFile, setupAndBuildAlignments, buildTrees, buildFaces, 
             buildReference, jobTreeDir, logLevel, batchSystem, retryCount, rescueJobFrequency, jobTreeStats, maxThreads)
    #print "going to run the command:", command
    #assert False
    system(command)
    logger.info("Ran the cactus workflow okay")
