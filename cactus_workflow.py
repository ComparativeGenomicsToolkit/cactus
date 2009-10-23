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

from sonLib.bioio import logger
from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import getTempFile
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from workflow.jobTree.scriptTree.target import Target

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import runCactusCore

from cactus.cactus_aligner import MakeSequences
from cactus.cactus_aligner import MakeSequencesOptions

from pecan2.pecan2_batch import pecan2BatchWrapperMiddleLevel
from pecan2.pecan2_batch import pecan2BatchWrapperTopLevel
    
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
        DownPassPhase(job, self, '0', self.options)
        logger.info("Created child target cactus_setup job, and follow on down pass job")

class CactusSetupWrapper(Target):
    def __init__(self, job, options, sequences):
        self.options = options
        self.sequences = sequences
        Target.__init__(self, job, None)
        
    def run(self, job):
        logger.info("Starting cactus setup target")
        runCactusSetup(self.options.netDisk, self.sequences, 
                       self.options.speciesTree,
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
    def __init__(self, job, previousTarget, netName, options):
        self.netName = netName
        self.options = options
        Target.__init__(self, job, previousTarget)

    def run(self, job):
        logger.info("Starting the down pass target")
        #Setup call to cactus aligner wrapper as child
        childTarget = CactusAlignerWrapper(job, self.options, self.netName, 0)
        self.addChildTarget(childTarget)
        #UpPassPhase(job, self, self.options)
        logger.info("Created child target aligner/core job, and follow on cleanup job")
    
class CactusAlignerWrapper(Target):
    def __init__(self, job, options, netName, iteration):
        self.options = options
        self.netName = netName
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
        self.addChildTarget(MakeSequences(job, self.options.netDisk, 
                                          self.netName, alignmentFile, makeSequencesOptions))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        CactusCoreWrapper(job, self, self.options, self.netName, alignmentFile, self.iteration)
        logger.info("Setup the follow on cactus_core target")
        
def cactusCoreParameters1():
    return { "maximumEdgeDegree":50, "proportionToKeep":1.0, "discardRatio":0.0, "minimumTreeCoverage":0.6, "minimumChainLength":10 }
    
def cactusCoreParameters2():
    return { "maximumEdgeDegree":50, "proportionToKeep":1.0, "discardRatio":0.5, "minimumTreeCoverage":0.5, "minimumChainLength":1 }

class CactusCoreWrapper(Target):
    def __init__(self, job, previousTarget, options, netName, alignmentFile, iteration):
        self.options = options
        self.netName = netName
        self.alignmentFile = alignmentFile
        self.iteration = iteration     
        Target.__init__(self, job, previousTarget)
    
    def run(self, job):
        logger.info("Starting the core wrapper target")

        if self.iteration == 0: 
            coreParameters = cactusCoreParameters1()
        else:
            coreParameters = cactusCoreParameters2()
            
        runCactusCore(netDisk=self.options.netDisk, 
                      alignmentFile=self.alignmentFile, 
                      tempDir=job.attrib["local_temp_dir"], 
                      netName=self.netName,
                      logLevel=getLogLevelString(), 
                      maximumEdgeDegree=coreParameters["maximumEdgeDegree"],
                      proportionOfAtomsToKeep=coreParameters["proportionToKeep"],
                      discardRatio=coreParameters["discardRatio"],
                      minimumTreeCoverage=coreParameters["minimumTreeCoverage"],
                      minimumChainLength=coreParameters["minimumChainLength"])
        logger.info("Ran the cactus core program okay")
        #Setup call to core and aligner recursive as follow on.
        CactusDownPass(job, self, self.options, self.netName, self.alignmentFile, self.iteration)
        logger.info("Issued a the recursive/cleanup wrapper as a follow on job")

class CactusDownPass(Target):
    def __init__(self, job, previousTarget, options, netName, alignmentFile, iteration):
        self.options = options
        self.iteration = iteration
        self.netName = netName
        self.alignmentFile = alignmentFile
        Target.__init__(self, job, previousTarget)
    
    def cleanup(self, job):
        #Cleans up from a round
        system("rm -rf %s" % self.alignmentFile) #Clean up the alignments file
    
    def run(self, job):
        logger.info("Starting the cactus down pass (recursive) target")
        #Traverses leaf jobs and create aligner wrapper targets as children.
        if self.iteration+1 < int(self.options.alignmentIterations):
            netNamesFile = getTempFile(".txt", job.attrib["global_temp_dir"])
            system("cactus_workflow_getNets %s %s %s" % (self.options.netDisk, self.netName, netNamesFile))
            fileHandle = open(netNamesFile, 'r')
            line = fileHandle.readline()
            while line != '':
                childNetName = line.split()[0]
                self.addChildTarget(CactusAlignerWrapper(job, self.options, childNetName, self.iteration+1))
                line = fileHandle.readline()
            fileHandle.close()
            logger.info("Created child targets for all the recursive reconstruction jobs")
      
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = getBasicOptionParser("usage: %prog [options] [sequence files]", "%prog 0.1")
    
    parser.add_option("--job", dest="jobFile", 
                      help="Job file containing command to run")

    parser.add_option("--speciesTree", dest="speciesTree", help="The species tree relating the input sequences")
    
    parser.add_option("--netDisk", dest="netDisk", help="The location of the net disk.") 
    
    parser.add_option("--alignmentIterations", dest="alignmentIterations", help="The number of recursive alignment iterations to emply", default="2")  
    
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
