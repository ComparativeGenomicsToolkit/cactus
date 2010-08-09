#!/usr/bin/env python

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework so structure all the related wrappers.

There are four high level wrappers, a SetupPhase, DownPassPhase, UpPassPhase, VerificationPahse. 

In the setup phase the system sets up the files needed for the reconstruction problem.

In the down pass phase alignments and trees are built.

In the up pass phase the adjacencies are added.

In the verification phase the reconstruction tree is checked against the expected spec.

"""

import os
import xml.etree.ElementTree as ET

from sonLib.bioio import logger
from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import getTempFile
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from sonLib.misc import sonTraceRootPath

from workflow.jobTree.scriptTree.target import Target

from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusCore
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import runCactusPhylogeny
from cactus.shared.common import runCactusAdjacencies
from cactus.shared.common import runCactusBaseAligner
from cactus.shared.common import runCactusMakeNormal 
from cactus.shared.common import runCactusReference
from cactus.shared.common import runCactusCheck

from cactus.blastAlignment.cactus_aligner import MakeSequences
from cactus.blastAlignment.cactus_batch import MakeBlastOptions
from cactus.blastAlignment.cactus_batch import makeBlastFromOptions

############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################

class CactusSetupPhase(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0)
        self.options = options 
        self.sequences = sequences 

    def run(self, localTempDir, globalTempDir):
        logger.info("Starting setup phase target")
        #Make the child setup job.
        self.addChildTarget(CactusSetupWrapper(self.options, self.sequences))
        #initialise the down pass as the follow on.. using special '0'
        self.setFollowOnTarget(CactusAlignmentPhase('0', self.options))
        logger.info("Created child target cactus_setup job, and follow on down pass job")

class CactusSetupWrapper(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0)
        self.options = options
        self.sequences = sequences
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting cactus setup target")
        runCactusSetup(self.options.cactusDisk, self.sequences, 
                       self.options.speciesTree, logLevel=getLogLevelString())
        logger.info("Finished the setup phase target")

############################################################
############################################################
############################################################
#The alignment phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################
    
class CactusAlignmentPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the down pass target")
        #Setup call to cactus aligner wrapper as child
        #Calculate the size of the child.
        self.addChildTarget(CactusAlignmentWrapper(self.options, self.flowerName, None, 0))
        self.setFollowOnTarget(CactusNormalPhase(self.flowerName, self.options))

def getAlignmentIteration(iterations, iterationNumber, flowerSize):
    assert len(iterations) > 0
    i = len(iterations)-1
    iterations = iterations[:]
    iterations.reverse() #Ensure we navigate backwards
    for iteration in iterations:
        assert int(iteration.attrib["number"]) == i
        if flowerSize < float(iteration.attrib["max_sequence_size"]) or iterationNumber >= i:
            assert i >= iterationNumber
            return i
        i -= 1

class CactusAlignmentWrapper(Target):
    def __init__(self, options, flowerName, alignmentFile, iteration):
        Target.__init__(self, time=0)
        self.options = options
        self.iteration = iteration
        self.flowerName = flowerName
        self.alignmentFile = alignmentFile
    
    def cleanup(self, localTempDir, globalTempDir):
        #Cleans up from a round
        if self.alignmentFile != None:
            system("rm -rf %s" % self.alignmentFile) #Clean up the alignments file
     
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the cactus down pass (recursive) target")
        #Traverses leaf jobs and create aligner wrapper targets as children.
        iterations = self.options.config.find("alignment").find("iterations").findall("iteration")
        if self.iteration < len(iterations):
            idealJobRuntime = float(self.options.config.attrib["ideal_job_runtime"])
            #base level flowers.
            baseLevelTime = 0.0
            baseLevelFlowers = []
            #This loop is properly atomic, because if it is run twice it will return the same
            #set of flowernames
            for childFlowerName, childFlowerSize in runCactusExtendFlowers(self.options.cactusDisk, self.flowerName, 
                                                                  localTempDir):
                assert childFlowerSize > 0
                nextIteration = getAlignmentIteration(iterations, self.iteration, childFlowerSize)
                iterationTime = float(iterations[nextIteration].attrib["time"])
                if iterations[nextIteration].attrib["type"] == "blast":
                    self.addChildTarget(CactusBlastWrapper(self.options, childFlowerName, nextIteration, iterationTime))
                else:
                    assert(iterations[nextIteration].attrib["type"] == "base")
                    baseLevelFlowers.append(childFlowerName)
                    baseLevelTime += iterationTime
                    if baseLevelTime >= idealJobRuntime:
                        self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelFlowers, baseLevelTime))
                        baseLevelFlowers = []
                        baseLevelTime = 0.0
            if len(baseLevelFlowers) > 0:
                self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelFlowers, baseLevelTime))        
            logger.info("Created child targets for all the recursive reconstruction jobs")
        else:
            logger.info("We've no more iterations to consider")

class CactusBlastWrapper(Target):
    def __init__(self, options, flowerName, iteration, time):
        Target.__init__(self, time=time)
        self.options = options
        self.flowerName = flowerName
        self.iteration = int(iteration)
    
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", globalTempDir)
        logger.info("Got an alignments file")
        
        #Now make the child aligner target
        alignmentNode = self.options.config.find("alignment")
        blastNode = alignmentNode.find("iterations").findall("iteration")[self.iteration].find("blast")
        blastMiscNode = alignmentNode.find("blast_misc")
        blastOptions =  \
        makeBlastFromOptions(MakeBlastOptions(int(blastNode.attrib["chunkSize"]),
                                              int(blastMiscNode.attrib["overlapSize"]), 
                                              blastNode.attrib["blastString"], 
                                              blastNode.attrib["selfBlastString"], 
                                              int(blastMiscNode.attrib["chunksPerJob"]), 
                                              bool(blastMiscNode.attrib["compressFiles"])))
            
        self.addChildTarget(MakeSequences(self.options.cactusDisk, 
                                          self.flowerName, alignmentFile, blastOptions))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        self.setFollowOnTarget(CactusCoreWrapper(self.options, self.flowerName, alignmentFile, self.iteration))
        logger.info("Setup the follow on cactus_core target")

class CactusCoreWrapper(Target):
    def __init__(self, options, flowerName, alignmentFile, iteration):
        Target.__init__(self, time=0)
        self.options = options
        self.flowerName = flowerName
        self.alignmentFile = alignmentFile
        self.iteration = iteration
    
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the core wrapper target")
    
        coreParameters = self.options.config.find("alignment").find("iterations").findall("iteration")[self.iteration].find("core")
        
        runCactusCore(cactusDisk=self.options.cactusDisk,
                      alignmentFile=self.alignmentFile, 
                      flowerName=self.flowerName,
                      logLevel=getLogLevelString(), 
                      annealingRounds=float(coreParameters.attrib["annealingRounds"]),
                      alignRepeatsAtRound=float(coreParameters.attrib["alignRepeatsAtRound"]),
                      trim=float(coreParameters.attrib["trim"]),
                      trimChange=float(coreParameters.attrib["trimChange"]),
                      minimumTreeCoverage=float(coreParameters.attrib["minimumTreeCoverage"]),
                      minimumBlockLength=float(coreParameters.attrib["minimumBlockLength"]),
                      minimumBlockLengthChange=float(coreParameters.attrib["minimumBlockLengthChange"]),
                      minimumChainLength=float(coreParameters.attrib["minimumChainLength"]),
                      minimumChainLengthChange=float(coreParameters.attrib["minimumChainLengthChange"]),
                      deannealingRounds=float(coreParameters.attrib["deannealingRounds"]))
        logger.info("Ran the cactus core program okay")
        
        #Setup call to core and aligner recursive as follow on.
        self.setFollowOnTarget(CactusAlignmentWrapper(self.options, self.flowerName, self.alignmentFile, self.iteration+1))
        logger.info("Issued a the recursive/cleanup wrapper as a follow on job")
     
class CactusBaseLevelAlignerWrapper(Target):
    #We split, to deal with cleaning up the alignment file
    def __init__(self, options, flowerNames, time):
        Target.__init__(self, time=time)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self, localTempDir, globalTempDir):
        #return
        runCactusBaseAligner(self.options.cactusDisk, self.flowerNames, getLogLevelString())
        logger.info("Run the cactus base aligner")
        
############################################################
############################################################
############################################################
#Normalisation pass
############################################################
############################################################
############################################################
    
class CactusNormalPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the normalisation phase")
        time = float(self.options.config.find("normal").attrib["time"])
        self.addChildTarget(CactusExtensionWrapper(self.options, [ self.flowerName ], MAKE_NORMAL, time))
        self.setFollowOnTarget(CactusPhylogenyPhase(self.flowerName, self.options))
        
class CactusNormalRunnable(Target):
    """This targets run the normalisation script.
    """
    def __init__(self, flowerNames, options):
        Target.__init__(self, time=0)
        self.flowerNames = flowerNames
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        runCactusMakeNormal(self.options.cactusDisk, flowerNames=self.flowerNames)

############################################################
############################################################
############################################################
#Phylogeny pass
############################################################
############################################################
############################################################
    
class CactusPhylogenyPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the phylogeny phase")
        if self.options.buildTrees:
            time = float(self.options.config.find("phylogeny").attrib["time"])
            self.addChildTarget(CactusExtensionWrapper(self.options, [ self.flowerName ], BUILD_TREES, time))
        self.setFollowOnTarget(CactusFacesPhase(self.flowerName, self.options))

############################################################
############################################################
############################################################
#Faces pass
############################################################
############################################################
############################################################
    
class CactusFacesPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the faces phase")
        if self.options.buildFaces:
            time = float(self.options.config.find("faces").attrib["time"])
            self.addChildTarget(CactusExtensionWrapper(self.options, [ self.flowerName ], BUILD_FACES, time))
        self.setFollowOnTarget(CactusReferencePhase(self.flowerName, self.options))

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################
    
class CactusReferencePhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the faces phase")
        if self.options.buildReference:
            time = float(self.options.config.find("reference").attrib["time"])
            self.addChildTarget(CactusExtensionWrapper(self.options, [ self.flowerName ], BUILD_REFERENCE, time))
        self.setFollowOnTarget(CactusCheckPhase(self.flowerName, self.options))
            
############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################
    
class CactusCheckPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the verification phase")
        time = float(self.options.config.find("check").attrib["time"])
        self.addChildTarget(CactusExtensionWrapper(self.options, [ self.flowerName ], CHECK, time))
        
############################################################
############################################################
############################################################
#The extension phase, used to modify the cactus tree.
############################################################
############################################################
############################################################
            
MAKE_NORMAL = 0
BUILD_TREES = 1
BUILD_FACES = 2
BUILD_REFERENCE = 3
CHECK = 4

class CactusExtensionWrapper(Target):
    def __init__(self, options, flowerNames, switch, unitTime):
        Target.__init__(self, time=unitTime*len(flowerNames))
        self.options = options
        self.flowerNames = flowerNames
        self.switch = switch
        self.unitTime = unitTime
    
    def run(self, localTempDir, globalTempDir):
        #The following are atomic, in that we check if they have already been run successfully.
        #This ensures things end up terminal normal.. which we need for face building.
        if self.switch == MAKE_NORMAL: #We set this as a follow on, as it is run in bottom up order (currently the only one, so it's on its own as a target)
            self.setFollowOnTarget(CactusNormalRunnable(options=self.options, flowerNames=self.flowerNames))
            #runCactusMakeNormal(self.options.cactusDisk, flowerNames=self.flowerNames)
        elif self.switch == BUILD_TREES:
            runCactusPhylogeny(self.options.cactusDisk, flowerNames=self.flowerNames)
            #Not atomic!
        elif self.switch == BUILD_FACES:
            runCactusAdjacencies(self.options.cactusDisk, flowerNames=self.flowerNames)
        elif self.switch == BUILD_REFERENCE:
            runCactusReference(self.options.cactusDisk, flowerNames=self.flowerNames)
        elif self.switch == CHECK:
            runCactusCheck(self.options.cactusDisk, self.flowerNames)
        #Make child jobs
        childFlowerNames = []
        idealJobRuntime = float(self.options.config.attrib["ideal_job_runtime"])
        for childFlowerName, childFlowerSize in runCactusGetFlowers(self.options.cactusDisk, self.flowerNames, localTempDir):
            assert(childFlowerSize) >= 0
            childFlowerNames.append(childFlowerName)
            if self.unitTime*len(childFlowerNames) >= idealJobRuntime:
                self.addChildTarget(CactusExtensionWrapper(self.options, childFlowerNames, self.switch, self.unitTime))
                childFlowerNames = []
        if len(childFlowerNames) > 0:
            self.addChildTarget(CactusExtensionWrapper(self.options, childFlowerNames, self.switch, self.unitTime))
      
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = getBasicOptionParser("usage: %prog [options] [sequence files]", "%prog 0.1")
    
    parser.add_option("--job", dest="jobFile", 
                      help="Job file containing command to run")

    parser.add_option("--speciesTree", dest="speciesTree", help="The species tree relating the input sequences")
    
    parser.add_option("--cactusDisk", dest="cactusDisk", help="The location of the flower disk.", default="cactusDisk") 
    
    parser.add_option("--configFile", dest="config", help="The file XML file containing the parameters for the pipeline", 
                      default=os.path.join(sonTraceRootPath(), "src", "cactus", "pipeline", "cactus_workflow_config.xml"))
    
    parser.add_option("--setupAndBuildAlignments", dest="setupAndBuildAlignments", action="store_true",
                      help="Setup and build alignments", default=False)
    
    parser.add_option("--buildTrees", dest="buildTrees", action="store_true",
                      help="Build trees", default=False) 
    
    parser.add_option("--buildFaces", dest="buildFaces", action="store_true",
                      help="Build adjacencies", default=False)
    
    parser.add_option("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the flowers", default=False)
    
    options, args = parseBasicOptions(parser)

    logger.info("Parsed arguments")
    
    options.config = ET.parse(options.config).getroot()
    logger.info("Parsed the XML options file")
    
    if options.setupAndBuildAlignments:
        baseTarget = CactusSetupPhase(options, args)
        logger.info("Going to create alignments and define the cactus tree")
    elif options.buildTrees or options.buildFaces or options.buildReference:
        baseTarget = CactusNormalPhase('0', options)
        logger.info("Starting from extension phase")
    else:
        logger.info("Nothing to do!")
        return
    baseTarget.execute(options.jobFile) 
    
    
    logger.info("Done with first target")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
