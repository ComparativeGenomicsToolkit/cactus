#!/usr/bin/env python

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework so structure all the related wrappers.

There are four high level wrappers, a SetupPhase, DownPassPhase, UpPassPhase, VerificationPahse. 

In the setup phase the system sets up the files needed for the reconstruction problem.

In the down pass phase alignments and trees are built.

In the up pass phase the adjacencies are added.

In the verification phase the reconstruction tree is checked against the expected spec.

"""

from sonLib.bioio import logger
from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import getTempFile
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from workflow.jobTree.scriptTree.target import Target

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import runCactusCore
from cactus.cactus_common import runCactusGetNets
from cactus.cactus_common import runCactusPhylogeny
from cactus.cactus_common import runCactusAdjacencies
from cactus.cactus_common import runCactusBaseAligner

from cactus.cactus_aligner import MakeSequences

from cactus.cactus_batch import makeTopLevelBlastOptions
from cactus.cactus_batch import makeUpperMiddleLevelBlastOptions
from cactus.cactus_batch import makeMiddleLevelBlastOptions
from cactus.cactus_batch import makeLowLevelBlastOptions
from cactus.cactus_batch import makeBlastFromOptions

IDEAL_JOB_RUNTIME = 200

############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################

class SetupPhase(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0)
        self.options = options 
        self.sequences = sequences 

    def run(self, localTempDir, globalTempDir):
        logger.info("Starting setup phase target")
        #Make the child setup job.
        self.addChildTarget(CactusSetupWrapper(self.options, self.sequences))
        #initialise the down pass as the follow on..
        self.setFollowOnTarget(AlignmentPhase('0', self.options))
        logger.info("Created child target cactus_setup job, and follow on down pass job")

class CactusSetupWrapper(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0)
        self.options = options
        self.sequences = sequences
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting cactus setup target")
        runCactusSetup(self.options.netDisk, self.sequences, 
                       self.options.speciesTree,
                       tempDir=localTempDir, logLevel=getLogLevelString())
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
    
class AlignmentPhase(Target):
    def __init__(self, netName, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.options = options
        
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the down pass target")
        #Setup call to cactus aligner wrapper as child
        #Calculate the size of the child.
        self.addChildTarget(CactusCoreWrapper2(self.options, '0', None, -1))
        self.setFollowOnTarget(HistoryPhase('0', None, self.options))

#Each level around 30x larger than the last
BASE_LEVEL_SIZE = 2000 #30000
GENE_LEVEL_SIZE = 1000000
LOCI_LEVEL_SIZE = 30000000
CHR_LEVEL_SIZE =  1000000000

iterationMaxSize = { 4:BASE_LEVEL_SIZE, 3:GENE_LEVEL_SIZE, 2:LOCI_LEVEL_SIZE, 1:CHR_LEVEL_SIZE }

def getIteration(iteration, netSize):
    if iteration == 4 or netSize < BASE_LEVEL_SIZE:
        return 4
    if iteration == 3 or netSize < GENE_LEVEL_SIZE:
        return 3
    elif iteration == 2 or netSize < LOCI_LEVEL_SIZE:
        return 2
    elif iteration == 1 or netSize < CHR_LEVEL_SIZE:
        return 1
    return 0

timeParameters = { 0:10000000, 1:10000000, 2:100, 3:20, 4:1 }

blastParameters = { 3:makeLowLevelBlastOptions, 2:makeMiddleLevelBlastOptions, 1:makeUpperMiddleLevelBlastOptions, 0:makeTopLevelBlastOptions }

cactusCoreParameters = { 
    0:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.7, "minimumTreeCoverageForBlocks":0.9, 
       "minimumBlockLength":4, "minimumChainLength":8, "alignRepeats":True, "alignUndoLoops":10, 
       "trim":20, "trimReduction":2, "extensionSteps":31500, "extensionStepsReduction":3500, 
       "minimumTreeCoverageForAlignUndoBlock":0.95, "minimumTreeCoverageForAlignUndoBlockReduction":0.1 },
    
    1:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.7, "minimumTreeCoverageForBlocks":0.9, 
       "minimumBlockLength":4, "minimumChainLength":8, "alignRepeats":True, "alignUndoLoops":10, 
       "trim":20, "trimReduction":2, "extensionSteps":31500, "extensionStepsReduction":3500, 
       "minimumTreeCoverageForAlignUndoBlock":0.95, "minimumTreeCoverageForAlignUndoBlockReduction":0.1 },
    
    2:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.7, "minimumTreeCoverageForBlocks":0.7, 
       "minimumBlockLength":0, "minimumChainLength":8, "alignRepeats":True, "alignUndoLoops":10, 
       "trim":20, "trimReduction":2, "extensionSteps":360, "extensionStepsReduction":40, 
       "minimumTreeCoverageForAlignUndoBlock":0.95, "minimumTreeCoverageForAlignUndoBlockReduction":0.1 },
    
    3:{ "maximumEdgeDegree":50, "minimumTreeCoverage":0.01, "minimumTreeCoverageForBlocks":0.01, 
       "minimumBlockLength":0, "minimumChainLength":2, "alignRepeats":True, "alignUndoLoops":10, 
       "trim":9, "trimReduction":1, "extensionSteps":20, "extensionStepsReduction":3, 
       "minimumTreeCoverageForAlignUndoBlock":0.95, "minimumTreeCoverageForAlignUndoBlockReduction":0.1 },
}

class CactusAlignerWrapper(Target):
    def __init__(self, options, netName, iteration):
        Target.__init__(self, time=timeParameters[iteration])
        self.options = options
        self.netName = netName
        self.iteration = int(iteration)
    
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", globalTempDir)
        logger.info("Got an alignments file")
        
        #Now make the child aligner target
        blastOptions = makeBlastFromOptions(blastParameters[self.iteration]())
            
        self.addChildTarget(MakeSequences(self.options.netDisk, 
                                          self.netName, alignmentFile, blastOptions))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        self.setFollowOnTarget(CactusCoreWrapper(self.options, self.netName, alignmentFile, self.iteration))
        logger.info("Setup the follow on cactus_core target")

class CactusCoreWrapper(Target):
    def __init__(self, options, netName, alignmentFile, iteration):
        Target.__init__(self, time=5)
        self.options = options
        self.netName = netName
        self.alignmentFile = alignmentFile
        self.iteration = iteration     
    
    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the core wrapper target")
    
        coreParameters = cactusCoreParameters[self.iteration]
        
        runCactusCore(netDisk=self.options.netDisk, 
                      alignmentFile=self.alignmentFile, 
                      netName=self.netName,
                      logLevel=getLogLevelString(), 
                      maximumEdgeDegree=coreParameters["maximumEdgeDegree"],
                      minimumTreeCoverage=coreParameters["minimumTreeCoverage"],
                      minimumTreeCoverageForBlocks=coreParameters["minimumTreeCoverageForBlocks"],
                      minimumBlockLength=coreParameters["minimumBlockLength"],
                      minimumChainLength=coreParameters["minimumChainLength"],
                      alignRepeats=coreParameters["alignRepeats"],
                      alignUndoLoops=coreParameters["alignUndoLoops"],
                      trim=coreParameters["trim"],
                      trimReduction=coreParameters["trimReduction"],
                      extensionSteps=coreParameters["extensionSteps"],
                      extensionStepsReduction=coreParameters["extensionStepsReduction"],
                      minimumTreeCoverageForAlignUndoBlock=coreParameters["minimumTreeCoverageForAlignUndoBlock"],
                      minimumTreeCoverageForAlignUndoBlockReduction=coreParameters["minimumTreeCoverageForAlignUndoBlockReduction"])
        logger.info("Ran the cactus core program okay")
        
        #Setup call to core and aligner recursive as follow on.
        self.setFollowOnTarget(CactusCoreWrapper2(self.options, self.netName, self.alignmentFile, self.iteration))
        logger.info("Issued a the recursive/cleanup wrapper as a follow on job")

class CactusCoreWrapper2(Target):
    #We split, to deal with cleaning up the alignment file
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
        assert self.iteration+1 <= 4
        baseLevelNets = []
        assert int(self.options.maxIteration) >= 0
        if int(self.options.maxIteration) >= 4:
            minSizeToExtend = 1
        else:
            minSizeToExtend = iterationMaxSize[int(self.options.maxIteration)+1]
        for childNetName, childNetSize in runCactusGetNets(self.options.netDisk, self.netName, localTempDir,
                                                           includeInternalNodes=False, 
                                                           recursive=True,
                                                           extendNonZeroTrivialGroups=True, 
                                                           minSizeToExtend=minSizeToExtend):
            assert childNetSize > 0
            nextIteration = getIteration(self.iteration+1, childNetSize)
            if nextIteration == 4:
                baseLevelNets.append(childNetName)
                if timeParameters[4]*len(baseLevelNets) > IDEAL_JOB_RUNTIME:
                    self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelNets))
                    baseLevelNets = []
            else: #Does not do any refinement if the net is completely specified.
                self.addChildTarget(CactusAlignerWrapper(self.options, childNetName, nextIteration))
        if len(baseLevelNets) > 0:
            self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelNets))        
        
        logger.info("Created child targets for all the recursive reconstruction jobs")
     
class CactusBaseLevelAlignerWrapper(Target):
    #We split, to deal with cleaning up the alignment file
    def __init__(self, options, netNames):
        Target.__init__(self, timeParameters[4]*len(netNames))
        self.options = options
        self.netNames = netNames
    
    def run(self, localTempDir, globalTempDir):
        runCactusBaseAligner(self.options.netDisk, self.netNames, getLogLevelString())
        logger.info("Run the cactus base aligner")
        
############################################################
############################################################
############################################################
#The history building phase
#
#Adds trees and adjacencies to the AVG component of the datastructure.
############################################################
############################################################
############################################################

class HistoryPhase(Target):
    def __init__(self, netName, netSize, options):
        Target.__init__(self, time=0)
        self.netName = netName
        self.netSize = netSize
        self.options = options

    def run(self, localTempDir, globalTempDir):
        logger.info("Starting the down pass target")
        if self.options.buildTrees:
            childTarget = CactusHistoryWrapper(self.options, [ self.netName ], self.netSize)
            self.addChildTarget(childTarget)

class CactusHistoryWrapper(Target):
    def __init__(self, options, netNames, cummulativeNetSize):
        Target.__init__(self, time=timeParameters[getIteration(0, cummulativeNetSize)])
        self.options = options
        self.netNames = netNames
    
    def run(self, localTempDir, globalTempDir):
        #Run cactus phylogeny
        if self.options.buildTrees:
            runCactusPhylogeny(self.options.netDisk, tempDir=localTempDir, netNames=self.netNames)
        if self.options.buildAdjacencies:
            runCactusAdjacencies(self.options.netDisk, tempDir=localTempDir, netNames=self.netNames)
        #Make child jobs
        childNetNames = []
        cummulativeNetSize = 0
        for netName in self.netNames:
            for childNetName, childNetSize in runCactusGetNets(self.options.netDisk, netName, localTempDir, includeInternalNodes=True, 
                                                               recursive=False, extendNonZeroTrivialGroups=False):
                if childNetName != netName: #Avoids running again for leaf without children
                    childNetNames.append(childNetName)
                    cummulativeNetSize += 1 #childNetSize
                    if timeParameters[getIteration(0, cummulativeNetSize)] > IDEAL_JOB_RUNTIME:
                        self.addChildTarget(CactusHistoryWrapper(self.options, childNetNames, cummulativeNetSize))
                        childNetNames = []
                        cummulativeNetSize = 0
        if len(childNetNames) > 0:
            self.addChildTarget(CactusHistoryWrapper(self.options, childNetNames, cummulativeNetSize))
      
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = getBasicOptionParser("usage: %prog [options] [sequence files]", "%prog 0.1")
    
    parser.add_option("--job", dest="jobFile", 
                      help="Job file containing command to run")

    parser.add_option("--speciesTree", dest="speciesTree", help="The species tree relating the input sequences")
    
    parser.add_option("--netDisk", dest="netDisk", help="The location of the net disk.", default="netDisk") 
    
    parser.add_option("--maxIteraton", dest="maxIteration", help="The maximum iteration to align to (0-4)..", 
                      default="4")
    
    parser.add_option("--setupAndBuildAlignments", dest="setupAndBuildAlignments", action="store_true",
                      help="Setup and build alignments", default=False)
    
    parser.add_option("--buildTrees", dest="buildTrees", action="store_true",
                      help="Build trees", default=False) 
    
    parser.add_option("--buildAdjacencies", dest="buildAdjacencies", action="store_true",
                      help="Build adjacencies", default=False) 
    
    options, args = parseBasicOptions(parser)

    logger.info("Parsed arguments")
    
    if options.setupAndBuildAlignments:
        baseTarget = SetupPhase(options, args)
    elif options.buildTrees or options.buildAdjacencies:
        baseTarget = HistoryPhase('0', None, options)
        
    baseTarget.execute(options.jobFile) 
    
    logger.info("Done with first target")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
