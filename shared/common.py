"""Setup functions for assisting in testing the various modules of the cactus package.
"""

import os

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system

from workflow.jobTree.jobTree import runJobTree 
from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def nameValue(name, value, valueType=str):
    """Little function to make it easier to make name value strings for commands.
    """
    if valueType == bool:
        if value:
            return "--%s" % name
        return ""
    if value == None:
        return ""
    return "--%s %s" % (name, valueType(value))

#############################################
#############################################
#All the following provide command line wrappers
#for core programs in the cactus pipeline.
#############################################
#############################################  

def runCactusSetup(reconstructionRootDir, sequences, 
                   newickTreeString, logLevel="DEBUG", debug=False):
    debugString = nameValue("debug", debug, bool)
    system("cactus_setup %s --speciesTree '%s' --netDisk %s \
--logLevel %s %s" \
           % (" ".join(sequences), newickTreeString,
              reconstructionRootDir, logLevel, debugString))
    logger.info("Ran cactus setup okay")
    
def runCactusAligner(netDisk, alignmentFile, tempDir, useDummy=True, netName="0", logLevel="DEBUG"):        
    """Runs job tree and fails if not complete.
    """
    tempDir = getTempDirectory(tempDir)
    jobTreeDir = os.path.join(tempDir, "jobTree")
    useDummy = nameValue("useDummy", useDummy, bool)
    command = "cactus_aligner.py --netDisk %s --netName %s \
--resultsFile %s %s --job JOB_FILE" % (netDisk, netName, alignmentFile, useDummy)
    runJobTree(command, jobTreeDir, logLevel=logLevel)
    runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
    system("rm -rf %s" % tempDir)
    logger.info("Ran the cactus aligner okay")
            
def runCactusCore(netDisk, alignmentFile, 
                  netName=0,
                  logLevel="DEBUG", 
                  writeDebugFiles=False,
                  alignUndoLoops=False,
                  alignRepeatsAtLoop=False,
                  maximumEdgeDegree=None,
                  extensionSteps=None,
                  extensionStepsReduction=None,
                  trim=None,
                  trimReduction=None,
                  maximumTreeCoverageUndo = None,
                  maximumTreeCoverageUndoReduction = None,
                  maximumChainLengthUndo = None,
                  maximumChainLengthUndoReduction = None,
                  minimumTreeCoverage=None,
                  minimumTreeCoverageReduction=None,
                  minimumBlockLength=None,
                  minimumChainLength=None,
                  minimumChainLengthReduction=None):
    writeDebugFiles = nameValue("writeDebugFiles", writeDebugFiles, bool)
    alignUndoLoops = nameValue("alignUndoLoops", alignUndoLoops, int)
    alignRepeatsAtLoop = nameValue("alignRepeatsAtLoop", alignRepeatsAtLoop, int)
    maximumEdgeDegree = nameValue("maxEdgeDegree", maximumEdgeDegree, int)
    extensionSteps = nameValue("extensionSteps", extensionSteps, int)
    extensionStepsReduction = nameValue("extensionStepsReduction", extensionStepsReduction, int)
    trim = nameValue("trim", trim, int)
    trimReduction = nameValue("trimReduction", trimReduction, int)
    
    maximumTreeCoverageUndo = nameValue("maximumTreeCoverageUndo", maximumTreeCoverageUndo, int)
    maximumTreeCoverageUndoReduction = nameValue("maximumTreeCoverageUndoReduction", maximumTreeCoverageUndoReduction, float)
    maximumChainLengthUndo = nameValue("maximumChainLengthUndo", maximumChainLengthUndo, int)
    maximumChainLengthUndoReduction = nameValue("maximumChainLengthUndoReduction", maximumChainLengthUndoReduction, int)
    
    minimumTreeCoverage = nameValue("minimumTreeCoverage", minimumTreeCoverage, float)
    minimumTreeCoverageReduction = nameValue("minimumTreeCoverageReduction", minimumTreeCoverageReduction, float)
    minimumBlockLength = nameValue("minimumBlockLength", minimumBlockLength, int)
    minimumChainLength = nameValue("minimumChainLength", minimumChainLength, int)
    minimumChainLengthReduction = nameValue("minimumChainLengthReduction", minimumChainLengthReduction, int)
    
    command = "cactus_core --netDisk %s --netName %s --alignments %s --logLevel %s \
%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % \
    (netDisk, netName, alignmentFile, 
     logLevel, 
     writeDebugFiles,
     alignUndoLoops,
     alignRepeatsAtLoop,
     maximumEdgeDegree,
     extensionSteps,
     extensionStepsReduction,
     trim,
     trimReduction,
     maximumTreeCoverageUndo,
     maximumTreeCoverageUndoReduction,
     maximumChainLengthUndo,
     maximumChainLengthUndoReduction,
     minimumTreeCoverage,
     minimumTreeCoverageReduction,
     minimumBlockLength,
     minimumChainLength,
     minimumChainLengthReduction)

    system(command)
    logger.info("Ran cactus_core okay")
    
def runCactusPhylogeny(netDisk,
                  netNames=("0",),
                  logLevel="DEBUG"):
    command = "cactus_phylogeny --netDisk %s --logLevel %s %s" % \
    (netDisk, logLevel, " ".join(netNames))
    system(command)
    logger.info("Ran cactus_phylogeny okay")
    
def runCactusAdjacencies(netDisk, netNames=("0",), logLevel="DEBUG"):
    command = "cactus_fillAdjacencies --netDisk %s --logLevel %s %s" %\
    (netDisk, logLevel, " ".join(netNames))
    system(command)
    logger.info("Ran cactus_fillAdjacencies OK")
    
def readNetNamesFile(netNamesFile):
    fileHandle = open(netNamesFile, 'r')
    line = fileHandle.readline()
    l = []
    while line != '':
        childNetName = line.split()[0]
        childNetSize = float(line.split()[1])
        l.append((childNetName, childNetSize))
        line = fileHandle.readline()
    fileHandle.close()
    return l
    
def runCactusGetNets(netDisk, netNames, tempDir):
    """Gets a list of nets attached to the given net. 
    """
    netNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_getNets %s %s %s" % (netDisk,  netNamesFile, " ".join(netNames)))
    l = readNetNamesFile(netNamesFile)
    os.remove(netNamesFile)
    return l

def runCactusExtendNets(netDisk, netName, tempDir,
                        minSizeToExtend=1):
    """Extends the terminal groups in the cactus and returns the list
    of their child nets with which to pass to core.
    The order of the nets is by ascending depth first discovery time.
    """
    netNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_extendNets %s %s %s %i" % (netDisk, netName, netNamesFile, int(minSizeToExtend)))
    l = readNetNamesFile(netNamesFile)
    os.remove(netNamesFile)
    return l

def runCactusGetUniqueName(netDisk, tempDir):
    """Gets a globally unique name.
    """
    uniqueNameFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_getUniqueName %s %s" % (netDisk, uniqueNameFile))
    fileHandle = open(uniqueNameFile, 'r')
    nameString = fileHandle.readline()[:-1]
    os.remove(uniqueNameFile)
    return nameString

def runCactusMakeTerminalNormal(netDisk, netNames):
    """Makes the given nets terminal-normal (not a mix of terminal and non-terminal children).
    """
    system("cactus_workflow_makeTerminalNormal %s %s" % (netDisk, " ".join(netNames)))

def runCactusBaseAligner(netDisk, netNames, logLevel="DEBUG"):
    """Runs cactus base aligner.
    """
    system("cactus_baseAligner --netDisk %s --logLevel %s %s" % (netDisk, logLevel, " ".join(netNames)))
    
def runCactusReference(netDisk, netNames, referenceName, logLevel="DEBUG"):
    """Runs cactus reference.
    """
    system("cactus_reference --netDisk %s --logLevel %s --referenceName %s %s" % (netDisk, logLevel, referenceName, " ".join(netNames)))

def runCactusCheck(netDisk, 
                    netNames=("0",), 
                    logLevel="DEBUG", 
                    recursive=False):
    recursive = nameValue("recursive", recursive, bool)
    system("cactus_check --netDisk %s %s --logLevel %s %s"  % (netDisk, " ".join(netNames), logLevel, recursive))
    logger.info("Ran cactus check")
    
def runCactusWorkflow(netDisk, sequenceFiles, 
                      newickTreeString, 
                      jobTreeDir, 
                      logLevel="DEBUG", retryCount=0, 
                      batchSystem="single_machine", 
                      rescueJobFrequency=None,
                      setupAndBuildAlignments=True,
                      buildTrees=True, buildFaces=True, buildReference=True):
    buildFaces = False
    setupAndBuildAlignments = nameValue("setupAndBuildAlignments", setupAndBuildAlignments, bool)
    buildTrees = nameValue("buildTrees", buildTrees, bool)
    buildFaces = nameValue("buildFaces", buildFaces, bool)
    buildReference = nameValue("buildReference", buildReference, bool)
    command = "cactus_workflow.py %s --speciesTree '%s' \
--netDisk %s %s %s %s %s --job JOB_FILE" % \
            (" ".join(sequenceFiles), newickTreeString,
             netDisk, setupAndBuildAlignments, buildTrees, buildFaces, buildReference)
    #print "going to run the command:", command
    #assert False
    runJobTree(command, jobTreeDir, logLevel, retryCount, batchSystem, rescueJobFrequency)
    logger.info("Ran the cactus workflow for %s okay" % netDisk)
    
#############################################
#############################################
#Runs cactus utilities.
#############################################
#############################################    
    
def runCactusTreeStats(outputFile, netDisk, netName='0'):
    command = "cactus_treeStats --netDisk %s --netName %s --outputFile %s" % (netDisk, netName, outputFile)
    system(command)
    logger.info("Ran the cactus tree stats command apprently okay")

def runCactusTreeStatsToLatexTables(inputFiles, regionNames, outputFile):
    assert len(regionNames) == len(inputFiles)
    k = " ".join([ "%s %s" % (i, j) for i, j in zip(inputFiles, regionNames) ])
    command = "cactus_treeStatsToLatexTables.py --outputFile %s %s" % (outputFile, k)
    system(command)
    logger.info("Ran cactus_treeStatsToLatexTables okay")
    
def runCactusTreeViewer(graphFile,
                        netDisk, 
                        netName="0", 
                        logLevel="DEBUG"):
    system("cactus_treeViewer --netDisk %s --netName %s --outputFile %s --logLevel %s" \
                    % (netDisk, netName, graphFile, logLevel))
    logger.info("Created a cactus tree graph")
    
def runCactusAdjacencyGraphViewer(graphFile,
                             netDisk, netName="0",
                             logLevel="DEBUG", includeInternalAdjacencies=False):
    includeInternalAdjacencies = nameValue("includeInternalAdjacencies", includeInternalAdjacencies, bool)
    system("cactus_adjacencyGraphViewer --netDisk %s --netName %s --outputFile %s --logLevel %s" \
                    % (netDisk, netName, graphFile, logLevel))
    logger.info("Created a break point graph of the problem")
    
def runCactusReferenceGraphViewer(graphFile,
                                  netDisk, netName="0",
                                  logLevel="DEBUG"):
    system("cactus_referenceViewer --netDisk %s --netName %s --outputFile %s --logLevel %s" \
                    % (netDisk, netName, graphFile, logLevel))
    logger.info("Created a cactus reference graph")
    

def runCactusMAFGenerator(mAFFile, netDisk, netName="0",
                          logLevel="DEBUG"):
    system("cactus_MAFGenerator --netDisk %s --netName %s --outputFile %s --logLevel %s" \
            % (netDisk, netName, mAFFile, logLevel))
    logger.info("Created a MAF for the given netDisk")
