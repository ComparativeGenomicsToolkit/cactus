"""Setup functions for assisting in running the various programs of the cactus package.
"""

import os

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system
from sonLib.bioio import nameValue

from workflow.jobTree.jobTree import runJobTree 
from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################  

def runCactusSetup(reconstructionRootDir, sequences, 
                   newickTreeString, logLevel="DEBUG", debug=False):
    debugString = nameValue("debug", debug, bool)
    system("cactus_setup %s --speciesTree '%s' --cactusDisk %s \
--logLevel %s %s" \
           % (" ".join(sequences), newickTreeString,
              reconstructionRootDir, logLevel, debugString))
    logger.info("Ran cactus setup okay")
    
def runCactusAligner(cactusDisk, alignmentFile, tempDir, useDummy=True, flowerName="0", logLevel="DEBUG"):        
    """Runs job tree and fails if not complete.
    """
    tempDir = getTempDirectory(tempDir)
    jobTreeDir = os.path.join(tempDir, "jobTree")
    useDummy = nameValue("useDummy", useDummy, bool)
    command = "cactus_aligner.py --cactusDisk %s --flowerName %s \
--resultsFile %s %s --job JOB_FILE" % (cactusDisk, flowerName, alignmentFile, useDummy)
    runJobTree(command, jobTreeDir, logLevel=logLevel)
    runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
    system("rm -rf %s" % tempDir)
    logger.info("Ran the cactus aligner okay")
            
def runCactusCore(cactusDisk, alignmentFile, 
                  flowerName=0,
                  logLevel="DEBUG", 
                  writeDebugFiles=False,
                  annealingRounds=False,
                  alignRepeatsAtRound=False,
                  trim=None,
                  trimChange=None,
                  minimumTreeCoverage=None,
                  minimumBlockLength=None,
                  minimumBlockLengthChange=None,
                  minimumChainLength=None,
                  minimumChainLengthChange=None,
                  deannealingRounds=None):
    writeDebugFiles = nameValue("writeDebugFiles", writeDebugFiles, bool)
    annealingRounds = nameValue("annealingRounds", annealingRounds, int)
    alignRepeatsAtRound = nameValue("alignRepeatsAtRound", alignRepeatsAtRound, int)
    trim = nameValue("trim", trim, int)
    trimChange = nameValue("trimChange", trimChange, float)
    minimumTreeCoverage = nameValue("minimumTreeCoverage", minimumTreeCoverage, float)
    minimumBlockLength = nameValue("minimumBlockLength", minimumBlockLength, int)
    minimumBlockLengthChange = nameValue("minimumBlockLengthChange", minimumBlockLengthChange, float)
    minimumChainLength = nameValue("minimumChainLength", minimumChainLength, int)
    minimumChainLengthChange = nameValue("minimumChainLengthChange", minimumChainLengthChange, float)
    deannealingRounds = nameValue("deannealingRounds", deannealingRounds, int)
    
    command = "cactus_core --cactusDisk %s --flowerName %s --alignments %s --logLevel %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (cactusDisk, flowerName, alignmentFile, logLevel, writeDebugFiles, annealingRounds, alignRepeatsAtRound,
     trim, trimChange, minimumTreeCoverage, minimumBlockLength,
     minimumBlockLengthChange, minimumChainLength, minimumChainLengthChange, deannealingRounds)

    system(command)
    logger.info("Ran cactus_core okay")
    
def runCactusPhylogeny(cactusDisk,
                  flowerNames=("0",),
                  logLevel="DEBUG"):
    command = "cactus_phylogeny --cactusDisk %s --logLevel %s %s" % \
    (cactusDisk, logLevel, " ".join(flowerNames))
    system(command)
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(cactusDisk, flowerNames=("0",), logLevel="DEBUG"):
    command = "cactus_fillAdjacencies --cactusDisk %s --logLevel %s %s" %\
    (cactusDisk, logLevel, " ".join(flowerNames))
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
    
def runCactusGetFlowers(cactusDisk, flowerNames, tempDir):
    """Gets a list of flowers attached to the given flower. 
    """
    flowerNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_getFlowers %s %s %s" % (cactusDisk,  flowerNamesFile, " ".join(flowerNames)))
    l = readFlowerNamesFile(flowerNamesFile)
    os.remove(flowerNamesFile)
    return l

def runCactusExtendFlowers(cactusDisk, flowerName, tempDir,
                        minSizeToExtend=1):
    """Extends the terminal groups in the cactus and returns the list
    of their child flowers with which to pass to core.
    The order of the flowers is by ascending depth first discovery time.
    """
    flowerNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_extendFlowers %s %s %s %i" % (cactusDisk, flowerName, flowerNamesFile, int(minSizeToExtend)))
    l = readFlowerNamesFile(flowerNamesFile)
    os.remove(flowerNamesFile)
    return l

def runCactusGetUniqueName(cactusDisk, tempDir):
    """Gets a globally unique name.
    """
    uniqueNameFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_getUniqueName %s %s" % (cactusDisk, uniqueNameFile))
    fileHandle = open(uniqueNameFile, 'r')
    nameString = fileHandle.readline()[:-1]
    os.remove(uniqueNameFile)
    return nameString

def runCactusMakeNormal(cactusDisk, flowerNames, maxNumberOfChains=0, logLevel="DEBUG"):
    """Makes the given flowers normal (see normalisation for the various phases)
    """
    system("cactus_normalisation --cactusDisk %s --maxNumberOfChains %i --logLevel %s %s" % (cactusDisk, maxNumberOfChains, logLevel, " ".join(flowerNames)))

def runCactusBaseAligner(cactusDisk, flowerNames, logLevel="DEBUG"):
    """Runs cactus base aligner.
    """
    system("cactus_baseAligner --cactusDisk %s --logLevel %s %s" % (cactusDisk, logLevel, " ".join(flowerNames)))
    
def runCactusReference(cactusDisk, flowerNames, logLevel="DEBUG", bottomUp=False):
    """Runs cactus reference.
    """
    #print "running", "cactus_reference --cactusDisk %s --logLevel %s %s" % (cactusDisk, logLevel, " ".join(flowerNames))
    #assert False
    bottomUp = nameValue("bottomUp", bottomUp, bool)
    system("cactus_reference --cactusDisk %s --logLevel %s %s %s" % (cactusDisk, logLevel, bottomUp, " ".join(flowerNames)))

def runCactusCheck(cactusDisk, 
                    flowerNames=("0",), 
                    logLevel="DEBUG", 
                    recursive=False):
    recursive = nameValue("recursive", recursive, bool)
    system("cactus_check --cactusDisk %s %s --logLevel %s %s"  % (cactusDisk, " ".join(flowerNames), logLevel, recursive))
    logger.info("Ran cactus check")
    
def runCactusWorkflow(cactusDisk, sequenceFiles, 
                      newickTreeString, 
                      jobTreeDir, 
                      logLevel="DEBUG", retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      setupAndBuildAlignments=True,
                      buildTrees=True, buildFaces=True, buildReference=True,
                      configFile=None):
    buildFaces=False
    setupAndBuildAlignments = nameValue("setupAndBuildAlignments", setupAndBuildAlignments, bool)
    buildTrees = nameValue("buildTrees", buildTrees, bool)
    buildFaces = nameValue("buildFaces", buildFaces, bool)
    buildReference = nameValue("buildReference", buildReference, bool)
    configFile = nameValue("configFile", configFile)
    command = "cactus_workflow.py %s --speciesTree '%s' \
--cactusDisk %s %s %s %s %s %s --job JOB_FILE" % \
            (" ".join(sequenceFiles), newickTreeString,
             cactusDisk, setupAndBuildAlignments, buildTrees, buildFaces, buildReference, configFile)
    #print "going to run the command:", command
    #assert False
    runJobTree(command, jobTreeDir, logLevel, retryCount, batchSystem, rescueJobFrequency)
    logger.info("Ran the cactus workflow for %s okay" % cactusDisk)
    
#############################################
#############################################
#Runs cactus utilities.
#############################################
#############################################    
    
def runCactusTreeStats(outputFile, cactusDisk, flowerName='0'):
    command = "cactus_treeStats --cactusDisk %s --flowerName %s --outputFile %s" % (cactusDisk, flowerName, outputFile)
    system(command)
    logger.info("Ran the cactus tree stats command apprently okay")

def runCactusTreeStatsToLatexTables(inputFiles, regionNames, outputFile):
    assert len(regionNames) == len(inputFiles)
    k = " ".join([ "%s %s" % (i, j) for i, j in zip(inputFiles, regionNames) ])
    command = "cactus_treeStatsToLatexTables.py --outputFile %s %s" % (outputFile, k)
    system(command)
    logger.info("Ran cactus_treeStatsToLatexTables okay")
    
def runCactusTreeViewer(graphFile,
                        cactusDisk, 
                        flowerName="0", 
                        logLevel="DEBUG"):
    system("cactus_treeViewer --cactusDisk %s --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDisk, flowerName, graphFile, logLevel))
    logger.info("Created a cactus tree graph")
    
def runCactusAdjacencyGraphViewer(graphFile,
                             cactusDisk, flowerName="0",
                             logLevel="DEBUG", includeInternalAdjacencies=False):
    includeInternalAdjacencies = nameValue("includeInternalAdjacencies", includeInternalAdjacencies, bool)
    system("cactus_adjacencyGraphViewer --cactusDisk %s --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDisk, flowerName, graphFile, logLevel))
    logger.info("Created a break point graph of the problem")
    
def runCactusReferenceGraphViewer(graphFile,
                                  cactusDisk, flowerName="0",
                                  logLevel="DEBUG"):
    system("cactus_referenceViewer --cactusDisk %s --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDisk, flowerName, graphFile, logLevel))
    logger.info("Created a cactus reference graph")
    

def runCactusMAFGenerator(mAFFile, cactusDisk, flowerName="0",
                          logLevel="DEBUG"):
    system("cactus_MAFGenerator --cactusDisk %s --flowerName %s --outputFile %s --logLevel %s" \
            % (cactusDisk, flowerName, mAFFile, logLevel))
    logger.info("Created a MAF for the given cactusDisk")
