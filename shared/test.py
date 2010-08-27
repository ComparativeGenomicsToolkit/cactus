"""Common test functions used for generating inputs to run cactus workflow and running
cactus workflow and the various utilities.
"""

import random
import os

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getRandomSequence
from sonLib.bioio import mutateSequence
from sonLib.bioio import reverseComplement
from sonLib.bioio import fastaWrite
from sonLib.bioio import printBinaryTree
from sonLib.bioio import system
from sonLib.bioio import getRandomAlphaNumericString
from sonLib.bioio import runGraphViz

from cactus.shared.common import runCactusWorkflow
from cactus.shared.common import runCactusTreeViewer
from cactus.shared.common import runCactusAdjacencyGraphViewer
from cactus.shared.common import runCactusReferenceGraphViewer 
from cactus.shared.common import runCactusCheck
from cactus.shared.common import runCactusTreeStats
from cactus.shared.common import runCactusMAFGenerator

from sonLib.bioio import TestStatus

from sonLib.tree import makeRandomBinaryTree

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def getCactusInputs_random(regionNumber=0, tempDir=None,
                          sequenceNumber=random.choice(xrange(100)), 
                          avgSequenceLength=random.choice(xrange(1, 5000)), 
                          treeLeafNumber=random.choice(xrange(1, 10))):
    """Gets a random set of sequences, each of length given, and a species
    tree relating them. Each sequence is a assigned an event in this tree.
    """
    #Make tree
    binaryTree = makeRandomBinaryTree(treeLeafNumber)
    newickTreeString = printBinaryTree(binaryTree, includeDistances=True)
    newickTreeLeafNames = []
    def fn(tree):
        if tree.internal:
            fn(tree.left)
            fn(tree.right)
        else:
            newickTreeLeafNames.append(tree.iD)
    fn(binaryTree)
    logger.info("Made random binary tree: %s" % newickTreeString)
    
    sequenceDirs = []
    for i in xrange(len(newickTreeLeafNames)):
        seqDir = getTempDirectory(rootDir=tempDir)
        sequenceDirs.append(seqDir)
        
    logger.info("Made a set of random directories: %s" % " ".join(sequenceDirs))
    
    #Random sequences and species labelling
    sequenceFile = None
    fileHandle = None
    parentSequence = getRandomSequence(length=random.choice(xrange(1, 2*avgSequenceLength)))[1]
    for i in xrange(sequenceNumber):
        if sequenceFile == None:
            if random.random() > 0.5: #Randomly choose the files to be attached or not
                suffix = ".fa.complete"
            else:
                suffix = ".fa"
            sequenceFile = getTempFile(rootDir=random.choice(sequenceDirs), suffix=suffix)
            fileHandle = open(sequenceFile, 'w')
        if random.random() > 0.8: #Get a new root sequence
            parentSequence = getRandomSequence(length=random.choice(xrange(1, 2*avgSequenceLength)))[1]
        sequence = mutateSequence(parentSequence, distance=random.random()*0.5)
        name = getRandomAlphaNumericString(15)
        if random.random() > 0.5:
            sequence = reverseComplement(sequence)
        fastaWrite(fileHandle, name, sequence)
        if random.random() > 0.5:
            fileHandle.close()
            fileHandle = None
            sequenceFile = None
    if fileHandle != None:
        fileHandle.close()
        
    logger.info("Made %s sequences in %s directories" % (sequenceNumber, len(sequenceDirs)))
    
    return sequenceDirs, newickTreeString

def parseNewickTreeFile(newickTreeFile):
    fileHandle = open(newickTreeFile, 'r')
    newickTreeString = fileHandle.readline()
    fileHandle.close()
    return newickTreeString

def getCactusInputs_blanchette(regionNumber=0, tempDir=None):
    """Gets the inputs for running cactus_workflow using a blanchette simulated region
    (0 <= regionNumber < 50).
    
    Requires setting SON_TRACE_DATASETS variable and having access to datasets.
    """
    assert regionNumber >= 0
    assert regionNumber < 50
    blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
    sequences = [ os.path.join(blanchettePath, ("%.2i.job" % regionNumber), species) \
                 for species in ("HUMAN", "CHIMP", "BABOON", "MOUSE", "RAT", "DOG", "CAT", "PIG", "COW") ] #Same order as tree
    newickTreeString = parseNewickTreeFile(os.path.join(blanchettePath, "tree.newick"))
    return sequences, newickTreeString
    
    
def getCactusInputs_encode(regionNumber=0, tempDir=None):
    """Gets the inputs for running cactus_workflow using an Encode pilot project region.
     (0 <= regionNumber < 15).
    
    Requires setting SON_TRACE_DATASETS variable and having access to datasets.
    """
    assert regionNumber >= 0
    assert regionNumber < 14
    encodeRegionString = "ENm%03i" % (regionNumber+1)
    encodeDatasetPath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
    sequences = [ os.path.join(encodeDatasetPath, encodeRegionString, ("%s.%s.fa" % (species, encodeRegionString))) for\
                species in ("human", "chimp", "baboon", "mouse", "rat", "dog", "cow") ]
    newickTreeString = parseNewickTreeFile(os.path.join(encodeDatasetPath, "reducedTree.newick"))
    return sequences, newickTreeString

def getCactusInputs_chromosomeX(regionNumber=0, tempDir=None):
    """Gets the inputs for running cactus_workflow using an some mammlian chromosome
    X's.
    
    Requires setting SON_TRACE_DATASETS variable and having access to datasets.
    """
    chrXPath = os.path.join(TestStatus.getPathToDataSets(), "chr_x")
    sequences = [ os.path.join(chrXPath, seqFile) for seqFile in ("hg18.fa", "panTro2.fa", "mouse_chrX.fa", "dog_chrX.fa") ]
    newickTreeString = parseNewickTreeFile(os.path.join(chrXPath, "newickTree.txt"))
    return sequences, newickTreeString

def runWorkflow_TestScript(sequences, newickTreeString, 
                           tempDir=None,
                           outputDir=None, 
                           batchSystem="single_machine",
                           buildTrees=True, buildFaces=True, buildReference=True,
                           buildCactusPDF=False,
                           buildAdjacencyPDF=False,
                           buildReferencePDF=False,
                           makeCactusTreeStats=False,
                           makeMAFs=False, 
                           cleanup=True,
                           configFile=None):
    """Runs the workflow and various downstream utilities.
    """
    logger.info("Running cactus workflow test script")
    logger.info("Got the following sequence dirs/files: %s" % " ".join(sequences))
    logger.info("Got the following tree %s" % newickTreeString)
    
    #Setup the temp dir
    if tempDir == None:
        tempDir = getTempDirectory(".")
    logger.info("Using the temp dir: %s" % tempDir)
        
    #Setup the output dir
    if outputDir == None: 
        outputDir = getTempDirectory(tempDir)
    logger.info("Using the output dir: %s" % outputDir)
    
    #Setup the flower disk.
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    cactusDisk = os.path.join(outputDir, "cactusDisk")
    logger.info("Using the cactusDisk: %s" % cactusDisk)
    system("rm -rf %s" % cactusDisk)
    logger.info("Cleaned up any previous flower disk: %s" % cactusDisk)
    
    #Setup the job tree dir.
    jobTreeDir = os.path.join(getTempDirectory(tempDir), "jobTree")
    logger.info("Got a job tree dir for the test: %s" % jobTreeDir)
    
    #Run the actual workflow
    runCactusWorkflow(cactusDisk, sequences, newickTreeString, jobTreeDir, 
                      batchSystem=batchSystem, buildTrees=buildTrees, 
                      buildFaces=buildFaces, buildReference=buildReference,
                      configFile=configFile)
    logger.info("Ran the the workflow")
    
    #Check if the jobtree completed sucessively.
    runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
    logger.info("Checked the job tree dir")
    
    #Check if the cactusDisk is okay..
    runCactusCheck(cactusDisk, recursive=True) #This should also occur during the workflow, so this
    #is redundant, but defensive
    logger.info("Checked the cactus tree")
    
    #Now run various utilities..
    
    #Run the cactus tree graph-viz plot
    if buildCactusPDF:
        cactusTreeDotFile = os.path.join(outputDir, "cactusTree.dot")
        cactusTreePDFFile = os.path.join(outputDir, "cactusTree.pdf")
        runCactusTreeViewer(cactusTreeDotFile, cactusDisk)
        runGraphViz(cactusTreeDotFile, cactusTreePDFFile)
        logger.info("Ran the cactus tree plot script")
    else:
        logger.info("Not building a cactus tree plot")
    
    #Run the cactus tree graph-viz plot
    if buildAdjacencyPDF:
        adjacencyGraphDotFile = os.path.join(outputDir, "adjacencyGraph.dot")
        adjacencyGraphPDFFile = os.path.join(outputDir, "adjacencyGraph.pdf")
        runCactusAdjacencyGraphViewer(adjacencyGraphDotFile, cactusDisk)
        runGraphViz(adjacencyGraphDotFile, adjacencyGraphPDFFile)
        logger.info("Ran the adjacency graph plot script")
    else:
        logger.info("Not building a adjacency graph plot")
    
    #Run the cactus tree graph-viz plot
    if buildReferencePDF:
        referenceGraphDotFile = os.path.join(outputDir, "referenceGraph.dot")
        referenceGraphPDFFile = os.path.join(outputDir, "referenceGraph.pdf")
        runCactusReferenceGraphViewer(referenceGraphDotFile, cactusDisk)
        runGraphViz(referenceGraphDotFile, referenceGraphPDFFile)
        logger.info("Ran the reference graph plot script")
    else:
        logger.info("Not building a reference graph plot")
    
    if makeCactusTreeStats:
        cactusTreeFile = os.path.join(outputDir, "cactusStats.xml")
        runCactusTreeStats(cactusTreeFile, cactusDisk)
        #Then run the latex script?
        logger.info("Ran the tree stats script")
    else:
        logger.info("Not running cactus tree stats")
    
    if makeMAFs:
        mAFFile = os.path.join(outputDir, "cactus.maf")
        runCactusMAFGenerator(mAFFile, cactusDisk)
        logger.info("Ran the MAF building script")
    else:
        logger.info("Not building the MAFs")
        
    if cleanup:
        #Now remove everything
        system("rm -rf %s" % tempDir)
        logger.info("Cleaned everything up")
    else:
        logger.info("Not cleaning up")
        
testRestrictions_NotShort = ()
        
def runWorkflow_multipleExamples(inputGenFunction,
                                 testNumber=1, 
                                 testRestrictions=(TestStatus.TEST_SHORT, TestStatus.TEST_MEDIUM, \
                                                   TestStatus.TEST_LONG, TestStatus.TEST_VERY_LONG,),
                               inverseTestRestrictions=False,
                               outputDir=None,
                               batchSystem="single_machine",
                               buildTrees=True, buildFaces=True, buildReference=True,
                               buildCactusPDF=False, buildAdjacencyPDF=False,
                               buildReferencePDF=False,
                               makeCactusTreeStats=False, makeMAFs=False,
                               configFile=None):
    """A wrapper to run a number of examples.
    """
    if (inverseTestRestrictions and TestStatus.getTestStatus() not in testRestrictions) or \
        (not inverseTestRestrictions and TestStatus.getTestStatus() in testRestrictions):
        for test in xrange(testNumber): 
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = inputGenFunction(regionNumber=test, tempDir=tempDir)
            if outputDir != None:
                out = os.path.join(outputDir, str(test))
            else:
                out = None
            runWorkflow_TestScript(sequences, newickTreeString, tempDir=tempDir,
                                   outputDir=out, batchSystem=batchSystem,
                                   buildTrees=buildTrees, buildFaces=buildFaces, buildReference=buildReference, 
                                   buildCactusPDF=buildCactusPDF, buildAdjacencyPDF=buildAdjacencyPDF,
                                   buildReferencePDF=buildReferencePDF,
                                   makeCactusTreeStats=makeCactusTreeStats, makeMAFs=makeMAFs, configFile=configFile)
            system("rm -rf %s" % tempDir)
            logger.info("Finished random test %i" % test)
    
