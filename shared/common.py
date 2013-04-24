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
                   newickTreeString, logLevel=None, outgroupEvents=None):
    logLevel = getLogLevelString2(logLevel)
    outgroupEvents = nameValue("outgroupEvents", outgroupEvents, str)
    system("cactus_setup %s --speciesTree '%s' --cactusDisk '%s' \
--logLevel %s %s" \
           % (" ".join(sequences), newickTreeString,
              cactusDiskDatabaseString, logLevel, outgroupEvents))
    logger.info("Ran cactus setup okay")
    
def runCactusBlast(sequenceFiles, outputFile, jobTreeDir,
                   chunkSize=None, overlapSize=None, 
                   logLevel=None, 
                   blastString=None, 
                   selfBlastString=None,
                   compressFiles=None,
                   lastzMemory=None):
    logLevel = getLogLevelString2(logLevel)
    chunkSize = nameValue("chunkSize", chunkSize, int)
    overlapSize = nameValue("overlapSize", overlapSize, int)
    blastString = nameValue("blastString", blastString, str)
    selfBlastString = nameValue("selfBlastString", selfBlastString, str)
    compressFiles = nameValue("compressFiles", compressFiles, bool)
    lastzMemory = nameValue("lastzMemory", lastzMemory, int)
    command = "cactus_blast.py %s  --cigars %s %s %s %s %s %s %s --jobTree %s --logLevel %s" % \
            (" ".join(sequenceFiles), outputFile,
             chunkSize, overlapSize, blastString, selfBlastString, compressFiles, lastzMemory, jobTreeDir, logLevel)
    logger.info("Running command : %s" % command)
    system(command)
    logger.info("Ran the cactus_batch command okay")
            
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
                  requiredIngroupFraction=None,
                  requiredOutgroupFraction=None,
                  requiredAllFraction=None,
                  singleCopyIngroup=None,
                  singleCopyOutgroup=None,
                  lastzArguments=None,
                  minimumSequenceLengthForBlast=None,
                  maxAdjacencyComponentSizeRatio=None,
                  constraints=None):
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
    
    requiredIngroupFraction = nameValue("requiredIngroupFraction", requiredIngroupFraction, float)
    requiredOutgroupFraction = nameValue("requiredOutgroupFraction", requiredOutgroupFraction, float)
    requiredAllFraction = nameValue("requiredAllFraction", requiredAllFraction, float)
    singleCopyIngroup = nameValue("singleCopyIngroup", singleCopyIngroup, bool)
    singleCopyOutgroup = nameValue("singleCopyOutgroup", singleCopyOutgroup, bool)
    maxAdjacencyComponentSizeRatio = nameValue("maxAdjacencyComponentSizeRatio", maxAdjacencyComponentSizeRatio, float)
    constraints = nameValue("constraints", constraints)
    
    command = "cactus_caf --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (cactusDiskDatabaseString, logLevel, alignments, annealingRounds, deannealingRounds, 
     trim, minimumTreeCoverage, blockTrim, 
     minimumBlockDegree, requiredIngroupFraction, requiredOutgroupFraction, requiredAllFraction, 
     singleCopyIngroup, singleCopyOutgroup, lastzArguments, minimumSequenceLengthForBlast, maxAdjacencyComponentSizeRatio, constraints)
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
                         alignAmbiguityCharacters=None,
                         requiredIngroupFraction=None,
                         requiredOutgroupFraction=None,
                         requiredAllFraction=None,
                         pruneOutStubAlignments=None,
                         maximumNumberOfSequencesBeforeSwitchingToFast=None,
                         calculateWhichEndsToComputeSeparately=None,
                         largeEndSize=None,
                         alignmentToPrecompute=None,
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
    pruneOutStubAlignments = nameValue("pruneOutStubAlignments", pruneOutStubAlignments, bool)
    alignAmbiguityCharacters = nameValue("alignAmbiguityCharacters", alignAmbiguityCharacters, bool)
    requiredIngroupFraction = nameValue("requiredIngroupFraction", requiredIngroupFraction, float)
    requiredOutgroupFraction = nameValue("requiredOutgroupFraction", requiredOutgroupFraction, float)
    requiredAllFraction = nameValue("requiredAllFraction", requiredAllFraction, float)
    maximumNumberOfSequencesBeforeSwitchingToFast=nameValue("maximumNumberOfSequencesBeforeSwitchingToFast", maximumNumberOfSequencesBeforeSwitchingToFast, int)
    calculateWhichEndsToComputeSeparately=nameValue("calculateWhichEndsToComputeSeparately", calculateWhichEndsToComputeSeparately, bool)
    largeEndSize=nameValue("largeEndSize", largeEndSize, int)
    alignmentToPrecompute=nameValue("alignmentToPrecompute", alignmentToPrecompute, str, quotes=True)
    precomputedAlignments=nameValue("precomputedAlignments", precomputedAlignments, str, quotes=True)
    
    masterMessages = popenCatch("cactus_bar --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % 
           (cactusDiskDatabaseString, logLevel, spanningTrees, maximumLength, gapGamma, 
            splitMatrixBiggerThanThis, anchorMatrixBiggerThanThis, repeatMaskMatrixBiggerThanThis,
            constraintDiagonalTrim, minimumBlockDegree, alignAmbiguityCharacters, pruneOutStubAlignments,
            requiredIngroupFraction, requiredOutgroupFraction, requiredAllFraction, diagonalExpansion,
            maximumNumberOfSequencesBeforeSwitchingToFast, calculateWhichEndsToComputeSeparately,
            largeEndSize, alignmentToPrecompute, precomputedAlignments), stdinString=flowerNames)
    logger.info("Ran cactus_bar okay")
    fn = lambda x : (x[0], int(x[1]), int(x[2]))
    l = [ fn(i.split()) for i in masterMessages.split("\n") if i != '' ]
    assert len(l) != 1, "We got just a single end to precompute, which doesn't achieve any parallelism."

def runCactusSecondaryDatabase(secondaryDatabaseString, create=True):
    command = "cactus_secondaryDatabase '%s' %s" % (secondaryDatabaseString, int(create))
    system(command)
    
def runCactusReference(cactusDiskDatabaseString, flowerNames, logLevel=None,
                       matchingAlgorithm=None, 
                       referenceEventString=None, 
                       permutations=None,
                       useSimulatedAnnealing=None,
                       theta=None,
                       maxWalkForCalculatingZ=None):
    """Runs cactus reference.
    """
    logLevel = getLogLevelString2(logLevel)
    matchingAlgorithm = nameValue("matchingAlgorithm", matchingAlgorithm)
    referenceEventString = nameValue("referenceEventString", referenceEventString)
    permutations = nameValue("permutations", permutations, int)
    useSimulatedAnnealing = nameValue("useSimulatedAnnealing", useSimulatedAnnealing, bool)
    theta = nameValue("theta", theta, float)
    maxWalkForCalculatingZ = nameValue("maxWalkForCalculatingZ", maxWalkForCalculatingZ, int)
    command = "cactus_reference --cactusDisk '%s' --logLevel %s %s %s %s %s %s %s" % (cactusDiskDatabaseString, logLevel, matchingAlgorithm, referenceEventString, permutations, useSimulatedAnnealing, theta, maxWalkForCalculatingZ)
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
      buildMaf=False,
      buildFasta=False,
      jobTreeStats=False,
      maxThreads=None,
      maxJobs=None,
      defaultMemory=None,
      logFile=None,
      extraJobTreeArgumentsString=""):
    logLevel = getLogLevelString2(logLevel)
    skipAlignments = nameValue("skipAlignments", skipAlignments, bool)
    buildAvgs = nameValue("buildAvgs", buildAvgs, bool)
    buildReference = nameValue("buildReference", buildReference, bool)
    buildHal = nameValue("buildHal", buildHal, bool)
    buildMaf= nameValue("buildMaf", buildMaf, bool)
    buildFasta = nameValue("buildFasta", buildFasta, bool)
    #Jobtree args
    batchSystem = nameValue("batchSystem", batchSystem, str, quotes=True)
    retryCount = nameValue("retryCount", retryCount, int)
    rescueJobFrequency = nameValue("rescueJobsFrequency", rescueJobFrequency, int)
    jobTreeStats = nameValue("stats", jobTreeStats, bool)
    maxThreads = nameValue("maxThreads", maxThreads, int)
    maxJobs = nameValue("maxJobs", maxJobs, int)
    defaultMemory= nameValue("defaultMemory", defaultMemory, int)
    logFile = nameValue("logFile", logFile, str)
    return "%s %s %s --jobTree %s --logLevel %s %s %s %s %s %s %s %s %s %s %s %s %s" % (skipAlignments, buildAvgs, 
             buildReference, jobTreeDir, logLevel, buildHal, buildMaf, buildFasta, batchSystem, retryCount, rescueJobFrequency, jobTreeStats, maxThreads, maxJobs, logFile, defaultMemory, extraJobTreeArgumentsString)
     
def runCactusWorkflow(experimentFile,
                      jobTreeDir, 
                      logLevel=None, retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      skipAlignments=False,
                      buildAvgs=False, buildReference=False,
                      buildHal=False,
                      buildMaf=False,
                      buildFasta=False,
                      jobTreeStats=False,
                      maxThreads=None,
                      maxJobs=None,
                      defaultMemory=None,
                      logFile=None,
                      extraJobTreeArgumentsString=""):
    command = ("cactus_workflow.py --experiment %s" % experimentFile) + " " + _fn(jobTreeDir, 
                      logLevel, retryCount, batchSystem, rescueJobFrequency, skipAlignments,
                      buildAvgs, buildReference, buildHal, buildMaf, buildFasta, jobTreeStats, maxThreads, maxJobs, defaultMemory, logFile, extraJobTreeArgumentsString=extraJobTreeArgumentsString)
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
                      skipAlignments=False,
                      buildHal=None,
                      buildMaf=None,
                      buildFasta=None,
                      joinMaf=None,
                      buildAvgs=False, 
                      jobTreeStats=False,
                      maxThreads=None,
                      maxJobs=None,
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
                      False,
                      buildFasta,
                      jobTreeStats, maxThreads, maxJobs, defaultMemory, logFile, extraJobTreeArgumentsString=extraJobTreeArgumentsString) + \
                      (" %s %s %s" % (nameValue("recursive", recursive, bool),
                                     nameValue("joinMAF", joinMaf, bool), nameValue("event", event)))
    if profileFile != None:
        command = "python -m cProfile -o %s %s/bin/%s" % (profileFile, cactusRootPath(), command)
    system(command)
    if buildMaf is True:
        assert buildHal is True
        basePath = os.path.dirname(inputDir)
        assert os.path.isdir(basePath)
        halExportCmd = "cactus2hal.py %s %s" % (inputDir, os.path.join(basePath, 'out.hal'))
        system(halExportCmd)
        mafExportCmd = "hal2maf %s %s --maxRefGap 1000000" % (os.path.join(basePath, 'out.hal'),
                                                         os.path.join(basePath, 'out.maf'))
        system(mafExportCmd)
        
                                           
    logger.info("Ran the cactus progressive okay")
    
def runCactusHalGenerator(cactusDiskDatabaseString,
                          secondaryDatabaseString, 
                          flowerNames,
                          referenceEventString, 
                          outputFile=None,
                          makeMaf=None,
                          showOnlySubstitutionsWithRespectToReference=None,
                          logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    popenPush("cactus_halGenerator --cactusDisk '%s' --secondaryDisk '%s' --logLevel %s %s %s %s %s" % 
           (cactusDiskDatabaseString, secondaryDatabaseString, logLevel, 
            nameValue("referenceEventString", referenceEventString),
            nameValue("outputFile", outputFile),
            nameValue("maf", makeMaf, bool),
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
