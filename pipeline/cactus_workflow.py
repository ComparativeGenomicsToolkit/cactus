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
from cactus.shared.common import runCactusGetNets
from cactus.shared.common import runCactusExtendNets
from cactus.shared.common import runCactusPhylogeny
from cactus.shared.common import runCactusAdjacencies
from cactus.shared.common import runCactusBaseAligner
from cactus.shared.common import runCactusMakeTerminalNormal
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
        runCactusSetup(self.options.netDisk, self.sequences, 
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
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the down pass target")
        #Setup call to cactus aligner wrapper as child
        #Calculate the size of the child.
        self.addChildTarget(CactusAlignmentWrapper(self.options, self.netName, None, 0))
        self.setFollowOnTarget(CactusTerminalNormalPhase(self.netName, self.options))

def getAlignmentIteration(iterations, iterationNumber, netSize):
    assert len(iterations) > 0
    i = len(iterations)-1
    iterations = iterations[:]
    iterations.reverse() #Ensure we navigate backwards
    for iteration in iterations:
        assert int(iteration.attrib["number"]) == i
        if netSize < float(iteration.attrib["max_sequence_size"]) or iterationNumber >= i:
            assert i >= iterationNumber
            return i
        i -= 1

class CactusAlignmentWrapper(Target):
    def __init__(self, options, netName, alignmentFile, iteration):
        Target.__init__(self, time=0)
        self.options = options
        self.iteration = iteration
        self.netName = netName
        self.alignmentFile = alignmentFile
    
    def cleanup(self, localTempDir, globalTempDir):
        #Cleans up from a round
        if self.alignmentFile != None:
            system("rm -rf %s" % self.alignmentFile) #Clean up the alignments file
     
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the cactus down pass (recursive) target")
        #Traverses leaf jobs and create aligner wrapper targets as children.
        #This loop is properly atomic, because if it is run twice it will return the same
        #set of netnames!
        iterations = self.options.config.find("alignment").find("iterations").findall("iteration")
        idealJobRuntime = float(self.options.config.attrib["ideal_job_runtime"])
        #base level nets.
        baseLevelTime = 0.0
        baseLevelNets = []
        for childNetName, childNetSize in runCactusExtendNets(self.options.netDisk, self.netName, 
                                                              localTempDir):
            assert childNetSize > 0
            nextIteration = getAlignmentIteration(iterations, self.iteration, childNetSize)
            iterationTime = float(iterations[nextIteration].attrib["time"])
            if iterations[nextIteration].attrib["type"] == "blast":
                self.addChildTarget(CactusBlastWrapper(self.options, childNetName, nextIteration, iterationTime))
            else:
                assert(iterations[nextIteration].attrib["type"] == "base")
                baseLevelNets.append(childNetName)
                baseLevelTime += iterationTime
                if baseLevelTime >= idealJobRuntime:
                    self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelNets, baseLevelTime))
                    baseLevelNets = []
                    baseLevelTime = 0.0
        if len(baseLevelNets) > 0:
            self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelNets, baseLevelTime))        
        logger.info("Created child targets for all the recursive reconstruction jobs")

class CactusBlastWrapper(Target):
    def __init__(self, options, netName, iteration, time):
        Target.__init__(self, time=time)
        self.options = options
        self.netName = netName
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
            
        self.addChildTarget(MakeSequences(self.options.netDisk, 
                                          self.netName, alignmentFile, blastOptions))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        self.setFollowOnTarget(CactusCoreWrapper(self.options, self.netName, alignmentFile, self.iteration))
        logger.info("Setup the follow on cactus_core target")
        

class CactusCoreWrapper(Target):
    def __init__(self, options, netName, alignmentFile, iteration):
        Target.__init__(self, time=0)
        self.options = options
        self.netName = netName
        self.alignmentFile = alignmentFile
        self.iteration = iteration
    
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the core wrapper target")
    
        coreParameters = self.options.config.find("alignment").find("iterations").findall("iteration")[self.iteration].find("core")
        
        runCactusCore(netDisk=self.options.netDisk, 
                      alignmentFile=self.alignmentFile, 
                      netName=self.netName,
                      logLevel=getLogLevelString(), 
                      alignUndoLoops=float(coreParameters.attrib["alignUndoLoops"]),
                      alignRepeatsAtLoop=float(coreParameters.attrib["alignRepeatsAtLoop"]),
                      maximumEdgeDegree=float(coreParameters.attrib["maximumEdgeDegree"]),
                      extensionSteps=float(coreParameters.attrib["extensionSteps"]),
                      extensionStepsReduction=float(coreParameters.attrib["extensionStepsReduction"]),
                      trim=float(coreParameters.attrib["trim"]),
                      trimReduction=float(coreParameters.attrib["trimReduction"]),
                      maximumTreeCoverageUndo=float(coreParameters.attrib["maximumTreeCoverageUndo"]),
                      maximumTreeCoverageUndoReduction=float(coreParameters.attrib["maximumTreeCoverageUndoReduction"]),
                      maximumChainLengthUndo=float(coreParameters.attrib["maximumChainLengthUndo"]),
                      maximumChainLengthUndoReduction=float(coreParameters.attrib["maximumChainLengthUndoReduction"]),
                      minimumTreeCoverage=float(coreParameters.attrib["minimumTreeCoverage"]),
                      minimumTreeCoverageReduction=float(coreParameters.attrib["minimumTreeCoverageReduction"]),
                      minimumBlockLength=float(coreParameters.attrib["minimumBlockLength"]),
                      minimumChainLength=float(coreParameters.attrib["minimumChainLength"]),
                      minimumChainLengthReduction=float(coreParameters.attrib["minimumChainLengthReduction"]))
        logger.info("Ran the cactus core program okay")
        
        #Setup call to core and aligner recursive as follow on.
        self.setFollowOnTarget(CactusAlignmentWrapper(self.options, self.netName, self.alignmentFile, self.iteration+1))
        logger.info("Issued a the recursive/cleanup wrapper as a follow on job")
     
class CactusBaseLevelAlignerWrapper(Target):
    #We split, to deal with cleaning up the alignment file
    def __init__(self, options, netNames, time):
        Target.__init__(self, time=time)
        self.options = options
        self.netNames = netNames
    
    def run(self, localTempDir, globalTempDir):
        #return
        runCactusBaseAligner(self.options.netDisk, self.netNames, getLogLevelString())
        logger.info("Run the cactus base aligner")
        
############################################################
############################################################
############################################################
#Terminal normal pass
############################################################
############################################################
############################################################
    
class CactusTerminalNormalPhase(Target):
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the terminal normal phase")
        time = float(self.options.config.find("terminal_normal").attrib["time"])
        self.addChildTarget(CactusExtensionWrapper(self.options, [ self.netName ], MAKE_TERMINAL_NORMAL, time))
        self.setFollowOnTarget(CactusPhylogenyPhase(self.netName, self.options))

############################################################
############################################################
############################################################
#Phylogeny pass
############################################################
############################################################
############################################################
    
class CactusPhylogenyPhase(Target):
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the phylogeny phase")
        if self.options.buildTrees:
            time = float(self.options.config.find("phylogeny").attrib["time"])
            self.addChildTarget(CactusExtensionWrapper(self.options, [ self.netName ], BUILD_TREES, time))
        self.setFollowOnTarget(CactusFacesPhase(self.netName, self.options))

############################################################
############################################################
############################################################
#Faces pass
############################################################
############################################################
############################################################
    
class CactusFacesPhase(Target):
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the faces phase")
        if self.options.buildFaces:
            time = float(self.options.config.find("faces").attrib["time"])
            self.addChildTarget(CactusExtensionWrapper(self.options, [ self.netName ], BUILD_FACES, time))
        self.setFollowOnTarget(CactusReferencePhase(self.netName, self.options))

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################
    
class CactusReferencePhase(Target):
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the faces phase")
        if self.options.buildReference:
            time = float(self.options.config.find("reference").attrib["time"])
            self.addChildTarget(CactusExtensionWrapper(self.options, [ self.netName ], BUILD_REFERENCE, time))
        self.setFollowOnTarget(CactusCheckPhase(self.netName, self.options))
            
############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################
    
class CactusCheckPhase(Target):
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the verification phase")
        time = float(self.options.config.find("check").attrib["time"])
        self.addChildTarget(CactusExtensionWrapper(self.options, [ self.netName ], CHECK, time))
        
############################################################
############################################################
############################################################
#The extension phase, used to modify the cactus tree.
############################################################
############################################################
############################################################
            
MAKE_TERMINAL_NORMAL = 0
BUILD_TREES = 1
BUILD_FACES = 2
BUILD_REFERENCE = 3
CHECK = 4

class CactusExtensionWrapper(Target):
    def __init__(self, options, netNames, switch, unitTime):
        Target.__init__(self, time=unitTime*len(netNames))
        self.options = options
        self.netNames = netNames
        self.switch = switch
        self.unitTime = unitTime
    
    def run(self, localTempDir, globalTempDir):
        #The following are atomic, in that we check if they have already been run successfully.
        #This ensures things end up terminal normal.. which we need for face building.
        if self.switch == MAKE_TERMINAL_NORMAL:
            runCactusMakeTerminalNormal(self.options.netDisk, netNames=self.netNames)
        elif self.switch == BUILD_TREES:
            runCactusPhylogeny(self.options.netDisk, netNames=self.netNames)
            #Not atomic!
        elif self.switch == BUILD_FACES:
            runCactusAdjacencies(self.options.netDisk, netNames=self.netNames)
        elif self.switch == BUILD_REFERENCE:
            runCactusReference(self.options.netDisk, netNames=self.netNames)
        elif self.switch == CHECK:
            runCactusCheck(self.options.netDisk, self.netNames)
        #Make child jobs
        childNetNames = []
        idealJobRuntime = float(self.options.config.attrib["ideal_job_runtime"])
        for childNetName, childNetSize in runCactusGetNets(self.options.netDisk, self.netNames, localTempDir):
            assert(childNetSize) >= 0
            childNetNames.append(childNetName)
            if self.unitTime*len(childNetNames) >= idealJobRuntime:
                self.addChildTarget(CactusExtensionWrapper(self.options, childNetNames, self.switch, self.unitTime))
                childNetNames = []
        if len(childNetNames) > 0:
            self.addChildTarget(CactusExtensionWrapper(self.options, childNetNames, self.switch, self.unitTime))
      
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = getBasicOptionParser("usage: %prog [options] [sequence files]", "%prog 0.1")
    
    parser.add_option("--job", dest="jobFile", 
                      help="Job file containing command to run")

    parser.add_option("--speciesTree", dest="speciesTree", help="The species tree relating the input sequences")
    
    parser.add_option("--netDisk", dest="netDisk", help="The location of the net disk.", default="netDisk") 
    
    parser.add_option("--configFile", dest="config", help="The file XML file containing the parameters for the pipeline", 
                      default=os.path.join(sonTraceRootPath(), "src", "cactus", "pipeline", "cactus_workflow_config.xml"))
    
    parser.add_option("--setupAndBuildAlignments", dest="setupAndBuildAlignments", action="store_true",
                      help="Setup and build alignments", default=False)
    
    parser.add_option("--buildTrees", dest="buildTrees", action="store_true",
                      help="Build trees", default=False) 
    
    parser.add_option("--buildFaces", dest="buildFaces", action="store_true",
                      help="Build adjacencies", default=False)
    
    parser.add_option("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the nets", default=False)
    
    options, args = parseBasicOptions(parser)

    logger.info("Parsed arguments")
    
    options.config = ET.parse(options.config).getroot()
    logger.info("Parsed the XML options file")
    
    if options.setupAndBuildAlignments:
        baseTarget = CactusSetupPhase(options, args)
        logger.info("Going to create alignments and define the cactus tree")
    elif options.buildTrees or options.buildFaces or options.buildReference:
        baseTarget = CactusTerminalNormalPhase('0', options)
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
