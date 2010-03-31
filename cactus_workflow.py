#!/usr/bin/env python

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework so structure all the related wrappers.

There are four high level wrappers, a SetupPhase, DownPassPhase, UpPassPhase, VerificationPahse. 

In the setup phase the system sets up the files needed for the reconstruction problem.

In the down pass phase alignments and trees are built.

In the up pass phase the adjacencies are added.

In the verification phase the reconstruction tree is checked against the expected spec.

"""
import xml.etree.ElementTree as ET
import os
import math

from sonLib.bioio import logger
from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import getTempFile
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from workflow.jobTree.scriptTree.target import Target

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import runCactusCore
from cactus.cactus_common import runCactusAdjacencyBuilder
from cactus.cactus_common import runCactusCheckReconstructionTree

from cactus.cactus_aligner import MakeSequences
from cactus.cactus_aligner import MakeSequencesOptions

from pecan2.pecan2_batch import pecan2BatchWrapperMiddleLevel
from pecan2.pecan2_batch import pecan2BatchWrapperTopLevel

############################################################
############################################################
############################################################
##Common functions
############################################################
############################################################
############################################################

alphaNumericChars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

def getUniquePrefix(jobNumber):
    """Gets a unique string prefix for a given job number.
    """
    jobNumber = int(jobNumber) #defensive, to protect against failure to convert from string
    return "".join([ alphaNumericChars[int(jobNumber / math.pow(len(alphaNumericChars), 4-i))\
                                         % len(alphaNumericChars)] for i in xrange(5) ])

def getUniqueJobNumber(prefix):
    """The inverse function to getUniquePrefix.
    
    >>> for i in xrange(1000):
    ...     import random
    ...     j = random.choice(xrange(10000000))
    ...     k = getUniquePrefix(j)
    ...     l = getUniqueJobNumber(k)
    ...     assert j == int(l)
    ...
    """
    return sum([ alphaNumericChars.index(prefix[i]) * math.pow(len(alphaNumericChars), 4-i)\
                 for i in xrange(5) ])
    
############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################

class SetupPhase(Target):
    def __init__(self, job, options, sequences):
        self.options = options 
        self.sequences = sequences 
        Target.__init__(self, job, None)

    def run(self, job):
        logger.info("Starting setup phase target")
        #Make the child setup job.
        self.addChildTarget(CactusSetupWrapper(job, self.options, self.sequences))
        #initialise the down pass as the follow on..
        DownPassPhase(job, self, self.options)
        logger.info("Created child target cactus_setup job, and follow on down pass job")

class CactusSetupWrapper(Target):
    def __init__(self, job, options, sequences):
        self.options = options
        self.sequences = sequences
        Target.__init__(self, job, None)
        
    def run(self, job):
        logger.info("Starting cactus setup target")
        runCactusSetup(self.options.reconstructionTree, self.sequences, 
                       self.options.speciesTree,
                       uniqueNamePrefix=getUniquePrefix(int(job.attrib["job_number"])),
                       tempDir=job.attrib["local_temp_dir"], logLevel=getLogLevelString())
        logger.info("Finished the setup phase target")

############################################################
############################################################
############################################################
#The down-pass phase.
#
#Creates the reconstruction structure with atoms + trees
#in a downward sweep.
############################################################
############################################################
############################################################
    
class DownPassPhase(Target):
    def __init__(self, job, previousTarget, options):
        self.options = options
        Target.__init__(self, job, previousTarget)

    def run(self, job):
        logger.info("Starting the down pass target")
        #Setup call to cactus aligner wrapper as child
        childTarget = CactusAlignerWrapper(job, self.options, "reconstructionProblem.xml", 0)
        self.addChildTarget(childTarget)
        UpPassPhase(job, self, self.options)
        logger.info("Created child target aligner/core job, and follow on cleanup job")
    
class CactusAlignerWrapper(Target):
    def __init__(self, job, options, reconstructionProblem, iteration):
        self.options = options
        self.reconstructionProblem = reconstructionProblem
        self.iteration = int(iteration)
        Target.__init__(self, job, None)
    
    def run(self, job):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", job.attrib["global_temp_dir"])
        logger.info("Got an alignments file")
        
        #Now make the child aligner target
        if self.iteration == 0:
            makeSequencesOptions = MakeSequencesOptions(pecan2BatchWrapperTopLevel)
        else:
            makeSequencesOptions = MakeSequencesOptions(pecan2BatchWrapperMiddleLevel)
        self.addChildTarget(MakeSequences(job, self.options.reconstructionTree, 
                                          self.reconstructionProblem, 
                                          alignmentFile, makeSequencesOptions))
        logger.info("Created the cactus_aligner child target")
    
        #First copy the reconstruction problem, so that we can rerun cactus core target if it fails.
        inputReconstructionProblem = getTempFile(suffix=".xml", rootDir=job.attrib["global_temp_dir"])
        system("cp %s %s" % (os.path.join(self.options.reconstructionTree, self.reconstructionProblem), inputReconstructionProblem))
        
        #Now setup a call to cactus core wrapper as a follow on
        CactusCoreWrapper(job, self, self.options, inputReconstructionProblem, self.reconstructionProblem, alignmentFile, self.iteration)
        logger.info("Setup the follow on cactus_core target")
        
def cactusCoreParameters1():
    return { "maximumEdgeDegree":50, "proportionToKeep":1.0, "discardRatio":0.0, "minimumTreeCoverage":0.6, "minimumChainLength":10 }
    
def cactusCoreParameters2():
    return { "maximumEdgeDegree":50, "proportionToKeep":1.0, "discardRatio":0.5, "minimumTreeCoverage":0.5, "minimumChainLength":1 }

class CactusCoreWrapper(Target):
    def __init__(self, job, previousTarget, options, inputReconstructionProblem,
                 outputReconstructionProblem, alignmentFile, iteration):
        self.options = options
        self.inputReconstructionProblem = inputReconstructionProblem
        self.outputReconstructionProblem = outputReconstructionProblem
        self.alignmentFile = alignmentFile
        self.iteration = iteration
        
        #Debug check that directory of reconstruction problem contains initially no directories (which will contain the child problems)
        absReconstructionProblem = os.path.join(self.options.reconstructionTree, self.outputReconstructionProblem)
        parentDir = os.path.split(absReconstructionProblem)[0]
        for fileName in list(os.listdir(parentDir)):
            absFileName = os.path.join(parentDir, fileName)
            assert not os.path.isdir(absFileName)
            
        Target.__init__(self, job, previousTarget)
    
    def run(self, job):
        logger.info("Starting the core wrapper target")
        #Remove the child reconstruction directories from the directory.
        absReconstructionProblem = os.path.join(self.options.reconstructionTree, self.outputReconstructionProblem)
        parentDir = os.path.split(absReconstructionProblem)[0]
        for fileName in list(os.listdir(parentDir)):
            absFileName = os.path.join(parentDir, fileName)
            if os.path.isdir(absFileName):
                system("rm -rf %s" % absFileName)
        #Copy the input reconstruction problem to the output copy (This restores the input file if it was previously broken).
        system("cp %s %s" % (self.inputReconstructionProblem, absReconstructionProblem))
        if self.iteration == 0: 
            coreParameters = cactusCoreParameters1()
        else:
            coreParameters = cactusCoreParameters2()
        runCactusCore(reconstructionRootDir=self.options.reconstructionTree, 
                      alignmentFile=self.alignmentFile, 
                      tempDir=job.attrib["local_temp_dir"], 
                      reconstructionProblem=self.outputReconstructionProblem,
                      treeProgram=self.options.treeBuilder, 
                      uniqueNamePrefix=getUniquePrefix(int(job.attrib["job_number"])),
                      logLevel=getLogLevelString(), 
                      maximumEdgeDegree=coreParameters["maximumEdgeDegree"],
                      proportionOfAtomsToKeep=coreParameters["proportionToKeep"],
                      discardRatio=coreParameters["discardRatio"],
                      minimumTreeCoverage=coreParameters["minimumTreeCoverage"],
                      minimumChainLength=coreParameters["minimumChainLength"])
        logger.info("Ran the cactus core program okay")
        #Setup call to core and aligner recursive as follow on.
        CactusDownPass(job, self, self.options, self.inputReconstructionProblem,
                       self.outputReconstructionProblem, self.alignmentFile, self.iteration)
        logger.info("Issued a the recursive/cleanup wrapper as a follow on job")

class CactusDownPass(Target):
    def __init__(self, job, previousTarget, options, inputReconstructionProblem,
                 outputReconstructionProblem, alignmentFile, iteration):
        self.options = options
        self.inputReconstructionProblem = inputReconstructionProblem
        self.outputReconstructionProblem = outputReconstructionProblem
        self.alignmentFile = alignmentFile
        self.iteration = iteration    
        Target.__init__(self, job, previousTarget)
    
    def cleanup(self, job):
        #Cleans up from a round
        system("rm -rf %s" % self.inputReconstructionProblem) #Clean up the alignments file
        system("rm -rf %s" % self.alignmentFile) #Clean up the alignments file
    
    def run(self, job):
        logger.info("Starting the cactus down pass (recursive) target")
        #Traverses leaf jobs and create aligner wrapper targets as children.
        if self.iteration+1 < int(self.options.alignmentIterations):
            def fn(reconstructionProblem):
                tree = ET.parse(os.path.join(self.options.reconstructionTree, 
                                            reconstructionProblem)).getroot()
                adjacencyComponent = tree.find("adjacency_components")    
                adjacencyComponents = adjacencyComponent.findall("adjacency_component")     
                if len(adjacencyComponents) > 0:      
                    for adjacencyComponent in adjacencyComponents:
                        fn(adjacencyComponent.attrib["child_file"])
                else:
                    self.addChildTarget(CactusAlignerWrapper(job, self.options, reconstructionProblem, self.iteration+1))
            fn(self.outputReconstructionProblem)
            logger.info("Created child targets for all the recursive reconstruction jobs")
    
############################################################
############################################################
############################################################
# The up-pass phase.
#
# This phase adds the adjacencies.
############################################################
############################################################
############################################################
    
class UpPassPhase(Target):
    def __init__(self, job, previousTarget, options):
        self.options = options
        Target.__init__(self, job, previousTarget)

    def run(self, job):
        logger.info("Starting the cactus up pass phase target")
        #Create child
        self.addChildTarget(CactusUpPass(job, self.options, "reconstructionProblem.xml"))
        
        ValidationPhase(job, self, self.options)
        logger.info("Created the adjacency building child targets and the follow on verification target")
   
def createChildTargets(job, options, reconstructionProblem, addChildTarget, childTarget):     
    tree = ET.parse(os.path.join(options.reconstructionTree,
                                 reconstructionProblem)).getroot()
    #Do the child jobs first.
    adjacencyComponent = tree.find("adjacency_components")
    adjacencyComponents = adjacencyComponent.findall("adjacency_component")
    for adjacencyComponent in adjacencyComponents:
        addChildTarget(childTarget(job, options, adjacencyComponent.attrib["child_file"]))
    
class CactusUpPass(Target):
    def __init__(self, job, options, reconstructionProblem):
        self.options = options
        self.reconstructionProblem = reconstructionProblem
        Target.__init__(self, job, None)
    
    def run(self, job):
        logger.info("Starting the cactus up pass target")
        createChildTargets(job, self.options, self.reconstructionProblem, self.addChildTarget, CactusUpPass)
        #Now do the follow on jobs.
        CactusAdjacencyBuilder(job, self, self.options, self.reconstructionProblem)
        logger.info("Issued the cactus up pass follow on jobs")
        
class CactusAdjacencyBuilder(Target):
    def __init__(self, job, previousTarget, options, reconstructionProblem):
        self.options = options
        self.reconstructionProblem = reconstructionProblem
        self.inputReconstructionProblem = getTempFile(suffix=".xml", rootDir=job.attrib["global_temp_dir"])
        system("cp %s %s" % (os.path.join(self.options.reconstructionTree, self.reconstructionProblem), self.inputReconstructionProblem))
        Target.__init__(self, job, previousTarget)
    
    def run(self, job):
        logger.info("Starting the cactus up pass follow on target, to build the actual adjacencies")
        #Copy the input reconstruction problem to the output copy (This restores the input file if it was previously broken).
        system("cp %s %s" % (self.inputReconstructionProblem, os.path.join(self.options.reconstructionTree, self.reconstructionProblem)))    
        #Run the adjacency builder.. 
        runCactusAdjacencyBuilder(self.options.reconstructionTree, self.reconstructionProblem, 
                                  tempDir=job.attrib["local_temp_dir"], 
                                  uniqueNamePrefix=getUniquePrefix(int(job.attrib["job_number"])),
                                  adjacencyProgram=self.options.adjacencyBuilder, 
                                  logLevel=job.attrib["log_level"])
        #Make a follow on job to clean up the adjacency builder.
        CactusAdjacencyBuilderFollowOn(job, self, self.inputReconstructionProblem)

class CactusAdjacencyBuilderFollowOn(Target):
    def __init__(self, job, previousTarget, inputReconstructionProblem):
        self.inputReconstructionProblem = inputReconstructionProblem
        Target.__init__(self, job, previousTarget)
    
    def cleanup(self, job):
        system("rm -rf %s" % self.inputReconstructionProblem)

############################################################
############################################################
############################################################
# The validation phase.
#
# This phase checks the reconstruction tree.
############################################################
############################################################
############################################################  
   
class ValidationPhase(UpPassPhase):
    def run(self, job):
        logger.info("Starting the validation phase target")
        #Create child
        self.addChildTarget(CactusCheckReconstructionTree(job, self.options, "reconstructionProblem.xml"))
        #Currently has no follow on.. will have verification target
    
class CactusCheckReconstructionTree(CactusUpPass):
    def run(self, job):
        logger.info("Starting the cactus check reconstruction target")
        
        createChildTargets(job, self.options, self.reconstructionProblem, self.addChildTarget, CactusCheckReconstructionTree)
        #Now do the follow on jobs.
        CactusCheckReconstructionTreeFollowOn(job, self, self.options, self.reconstructionProblem)
        logger.info("Issued the cactus up pass follow on jobs")
        
class CactusCheckReconstructionTreeFollowOn(Target):
    def __init__(self, job, previousTarget, options, reconstructionProblem):
        self.options = options
        self.reconstructionProblem = reconstructionProblem
        Target.__init__(self, job, previousTarget)
    
    def run(self, job):
        logger.info("Starting the cactus check reconstruction tree follow on target, to verify the reconstruction tree.")
        #Run the adjacency builder.. 
        runCactusCheckReconstructionTree(self.options.reconstructionTree, self.reconstructionProblem, 
                                         logLevel=job.attrib["log_level"], recursive=False)
      
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = getBasicOptionParser("usage: %prog [options] [sequence files]", "%prog 0.1")
    
    parser.add_option("--job", dest="jobFile", 
                      help="Job file containing command to run")

    parser.add_option("--speciesTree", dest="speciesTree", help="The species tree relating the input sequences")
    
    parser.add_option("--reconstructionTree", dest="reconstructionTree", help="Top level directory that will be created to write the reconstruction tree structure in") 
    
    parser.add_option("--aligner", dest="aligner", help="The program to build alignments from (is used as prefix string, and so may contain additional arguments to the program (such as a configuration file)")  
    
    parser.add_option("--alignmentIterations", dest="alignmentIterations", help="The number of recursive alignment iterations to emply", default="2")  
    
    parser.add_option("--treeBuilder", dest="treeBuilder", help="The program to build trees from", default="cactus_coreTestTreeBuilder.py")  
    
    parser.add_option("--adjacencyBuilder", dest="adjacencyBuilder", help="The program to build adjacencies from", default="cactus_adjacencyTestAdjacencyBuilder.py")
    
    options, args = parseBasicOptions(parser)

    logger.info("Parsed arguments")
    
    job = ET.parse(options.jobFile).getroot()
    setupPhase = SetupPhase(job, options, args)
    setupPhase.execute(options.jobFile)
    
    logger.info("Done with first target")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
