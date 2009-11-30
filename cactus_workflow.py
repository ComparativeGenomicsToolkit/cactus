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

from pecan2.pecan2_batch import makeTopLevelBlastOptions
from pecan2.pecan2_batch import makeUpperMiddleLevelBlastOptions
from pecan2.pecan2_batch import makeMiddleLevelBlastOptions
from pecan2.pecan2_batch import makeLowLevelBlastOptions
from pecan2.pecan2_batch import makeBlastFromOptions

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

def getChildNets(netDisk, netName, tempDir):
    """Gets a list of leaf nets attached to the given net. If the net has no children,
    as is therefore a leaf, it will also be returned. 
    
    The net names are returned in a list of tuples with the size (in terms of total bases pairs 
    of sequence contained within the net).
    """
    netNamesFile = getTempFile(".txt", tempDir)
    system("cactus_workflow_getNets %s %s %s" % (netDisk, netName, netNamesFile))
    fileHandle = open(netNamesFile, 'r')
    line = fileHandle.readline()
    l = []
    while line != '':
        childNetName = line.split()[0]
        childNetSize = float(line.split()[1])
        l.append((childNetName, childNetSize))
        line = fileHandle.readline()
    fileHandle.close()
    os.remove(netNamesFile)
    return l
    
class DownPassPhase(Target):
    def __init__(self, job, previousTarget, netName, options):
        self.netName = netName
        self.options = options
        Target.__init__(self, job, previousTarget)

    def run(self, job):
        logger.info("Starting the down pass target")
        #Setup call to cactus aligner wrapper as child
        #Calculate the size of the child.
        l = getChildNets(self.options.netDisk, self.netName, job.attrib["local_temp_dir"])
        assert(len(l) == 1)
        assert(l[0][0] == self.netName)
        assert(l[0][1] >= 0)
        childTarget = CactusAlignerWrapper(job, self.options, self.netName, l[0][1], 0)
        self.addChildTarget(childTarget)
        #UpPassPhase(job, self, self.options)
        logger.info("Created child target aligner/core job, and follow on cleanup job")

#Each level around 30x larger than the last
#BASE_LEVEL_SIZE = 10000
GENE_LEVEL_SIZE = 500000
LOCI_LEVEL_SIZE = 100000000
CHR_LEVEL_SIZE =  3000000000

class CactusAlignerWrapper(Target):
    def __init__(self, job, options, netName, netSize, iteration):
        self.options = options
        self.netName = netName
        self.netSize = netSize
        self.iteration = int(iteration)
        Target.__init__(self, job, None)
    
    def run(self, job):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", job.attrib["global_temp_dir"])
        logger.info("Got an alignments file")
        
        #Now make the child aligner target
        assert self.iteration <= 3
        #if self.iteration == 4 or self.netSize < BASE_LEVEL_SIZE:
        #    self.iteration = 4
        #    blastOptions = makeBlastFromOptions(makeLowLevelBlastOptions())
        if self.iteration == 3 or self.netSize < GENE_LEVEL_SIZE:
            self.iteration = 3
            blastOptions = makeBlastFromOptions(makeLowLevelBlastOptions())
        elif self.iteration == 2 or self.netSize < LOCI_LEVEL_SIZE:
            self.iteration = 2
            blastOptions = makeBlastFromOptions(makeMiddleLevelBlastOptions())
        elif self.iteration == 1 or self.netSize < CHR_LEVEL_SIZE:
            self.iteration = 1
            blastOptions = makeBlastFromOptions(makeUpperMiddleLevelBlastOptions())
        else:
            self.iteration = 0
            blastOptions = makeBlastFromOptions(makeTopLevelBlastOptions())
            
        self.addChildTarget(MakeSequences(job, self.options.netDisk, 
                                          self.netName, alignmentFile, blastOptions))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        CactusCoreWrapper(job, self, self.options, self.netName, self.netSize, alignmentFile, self.iteration)
        logger.info("Setup the follow on cactus_core target")
        
cactusCoreParameters = { 
    0:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.6, "minimumAtomLength":4, "minimumChainLength":12, "trim":3, "alignRepeats":False },
    1:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.5, "minimumAtomLength":4, "minimumChainLength":8, "trim":3, "alignRepeats":False },
    2:{ "maximumEdgeDegree":30, "minimumTreeCoverage":0.2, "minimumAtomLength":4, "minimumChainLength":8, "trim":3, "alignRepeats":False },
    3:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.0, "minimumAtomLength":0, "minimumChainLength":0, "trim":0, "alignRepeats":False },
    #4:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.0, "minimumAtomLength":0, "minimumChainLength":0, "trim":0, "alignRepeats":False }
}

class CactusCoreWrapper(Target):
    def __init__(self, job, previousTarget, options, netName, netSize, alignmentFile, iteration):
        self.options = options
        self.netName = netName
        self.netSize = netSize
        self.alignmentFile = alignmentFile
        self.iteration = iteration     
        Target.__init__(self, job, previousTarget)
    
    def run(self, job):
        logger.info("Starting the core wrapper target")
    
        coreParameters = cactusCoreParameters[self.iteration]
            
        runCactusCore(netDisk=self.options.netDisk, 
                      alignmentFile=self.alignmentFile, 
                      tempDir=job.attrib["local_temp_dir"], 
                      netName=self.netName,
                      logLevel=getLogLevelString(), 
                      maximumEdgeDegree=coreParameters["maximumEdgeDegree"],
                      minimumTreeCoverage=coreParameters["minimumTreeCoverage"],
                      minimumAtomLength=coreParameters["minimumAtomLength"],
                      minimumChainLength=coreParameters["minimumChainLength"],
                      trim=coreParameters["trim"],
                      alignRepeats=coreParameters["alignRepeats"])
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
        if self.iteration+1 <= 3:
            for childNetName, childNetSize in getChildNets(self.options.netDisk, self.netName, job.attrib["local_temp_dir"]):
                if childNetSize > 0: #Does not do any refinement if the net is completely specified.
                    self.addChildTarget(CactusAlignerWrapper(job, self.options, childNetName, childNetSize, self.iteration+1))
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
