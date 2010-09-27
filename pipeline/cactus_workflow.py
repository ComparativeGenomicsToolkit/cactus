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
from optparse import OptionParser

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from sonLib.misc import sonTraceRootPath

from workflow.jobTree.scriptTree.target import Target
from workflow.jobTree.scriptTree.stack import Stack

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
        Target.__init__(self, time=0.0002)
        self.options = options 
        self.sequences = sequences 

    def run(self):
        logger.info("Starting setup phase target")
        #Make the child setup job.
        self.addChildTarget(CactusSetupWrapper(self.options, self.sequences))
        #initialise the down pass as the follow on.. using special '0'
        self.setFollowOnTarget(CactusAlignmentPhase('0', self.options))
        logger.info("Created child target cactus_setup job, and follow on down pass job")

class CactusSetupWrapper(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=1.0)
        self.options = options
        self.sequences = sequences
        
    def run(self):
        logger.info("Starting cactus setup target")
        runCactusSetup(self.options.cactusDiskDatabaseString, self.sequences, 
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
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
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
        Target.__init__(self, time=1.5)
        self.options = options
        self.iteration = iteration
        self.flowerName = flowerName
        self.alignmentFile = alignmentFile
    
    def run(self):
        #Cleans up from a round
        if self.alignmentFile != None:
            system("rm -rf %s" % self.alignmentFile) #Clean up the alignments file
        
        logger.info("Starting the cactus down pass (recursive) target")
        #Traverses leaf jobs and create aligner wrapper targets as children.
        iterations = self.options.config.find("alignment").find("iterations").findall("iteration")
        if self.iteration < len(iterations):
            #base level flowers.
            baseLevelFlowers = []
            #This loop is properly atomic, because if it is run twice it will return the same
            #set of flowernames
            for childFlowerName, childFlowerSize in runCactusExtendFlowers(self.options.cactusDiskDatabaseString, self.flowerName, 
                                                                  self.getLocalTempDir()):
                assert childFlowerSize >= 0
                nextIteration = getAlignmentIteration(iterations, self.iteration, childFlowerSize)
                if iterations[nextIteration].attrib["type"] == "blast":
                    self.addChildTarget(CactusBlastWrapper(self.options, childFlowerName, nextIteration))
                else:
                    assert(iterations[nextIteration].attrib["type"] == "base")
                    baseLevelFlowers.append(childFlowerName)
                    if len(baseLevelFlowers) > 50:
                        self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelFlowers))
                        baseLevelFlowers = []
            if len(baseLevelFlowers) > 0:
                self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, baseLevelFlowers))        
            logger.info("Created child targets for all the recursive reconstruction jobs")
        else:
            logger.info("We've no more iterations to consider")

class CactusBlastWrapper(Target):
    def __init__(self, options, flowerName, iteration):
        Target.__init__(self, time=0.01)
        self.options = options
        self.flowerName = flowerName
        self.iteration = int(iteration)
    
    def run(self):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", self.getGlobalTempDir())
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
            
        self.addChildTarget(MakeSequences(self.options.cactusDiskDatabaseString, 
                                          self.flowerName, alignmentFile, blastOptions))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        self.setFollowOnTarget(CactusCoreWrapper(self.options, self.flowerName, alignmentFile, self.iteration))
        logger.info("Setup the follow on cactus_core target")

class CactusCoreWrapper(Target):
    def __init__(self, options, flowerName, alignmentFile, iteration):
        Target.__init__(self, time=150, memory=4294967295) #Request 2^32 (4 gigs of ram)
        self.options = options
        self.flowerName = flowerName
        self.alignmentFile = alignmentFile
        self.iteration = iteration
    
    def run(self):
        logger.info("Starting the core wrapper target")
    
        coreParameters = self.options.config.find("alignment").find("iterations").findall("iteration")[self.iteration].find("core")
        
        runCactusCore(cactusDiskDatabaseString=self.options.cactusDiskDatabaseString,
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
                      deannealingRounds=float(coreParameters.attrib["deannealingRounds"]),
                      adjacencyComponentOverlap=int(coreParameters.attrib["adjacencyComponentOverlap"]))
        logger.info("Ran the cactus core program okay")
        
        #Setup call to core and aligner recursive as follow on.
        self.setFollowOnTarget(CactusAlignmentWrapper(self.options, self.flowerName, self.alignmentFile, self.iteration+1))
        logger.info("Issued a the recursive/cleanup wrapper as a follow on job")
     
class CactusBaseLevelAlignerWrapper(Target):
    #We split, to deal with cleaning up the alignment file
    def __init__(self, options, flowerNames):
        Target.__init__(self, time=30) #time)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        #return
        runCactusBaseAligner(self.options.cactusDiskDatabaseString, self.flowerNames, getLogLevelString())
        logger.info("Run the cactus base aligner")
        
############################################################
############################################################
############################################################
#Normalisation pass
############################################################
############################################################
############################################################

def makeChildTargets(options, flowerNames, target, childTarget, jobNumber=100):
    #Make child jobs
    childFlowerNames = []
    for childFlowerName, childFlowerSize in runCactusGetFlowers(options.cactusDiskDatabaseString, flowerNames, target.getLocalTempDir()):
        assert(childFlowerSize) >= 0
        childFlowerNames.append(childFlowerName)
        if len(childFlowerNames) >= jobNumber:
            target.addChildTarget(childTarget(options, childFlowerNames))
            childFlowerNames = []
    if len(childFlowerNames) > 0:
        target.addChildTarget(childTarget(options, childFlowerNames))
    
class CactusNormalPhase(Target):
    def __init__(self, flowerName, options, normalisationRounds=-1):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        if(normalisationRounds < 0):
            normalisationRounds =  int(self.options.config.find("normal").attrib["rounds"])
        assert(normalisationRounds > 0)
        self.normalisationRounds=normalisationRounds
        
    def run(self):
        logger.info("Starting the normalisation phase")
        self.addChildTarget(CactusNormalDown(self.options, [ self.flowerName ]))
        if self.normalisationRounds-1 > 0:
            self.setFollowOnTarget(CactusNormalPhase(self.flowerName, self.options, self.normalisationRounds-1))
        else:
            self.setFollowOnTarget(CactusPhylogenyPhase(self.flowerName, self.options))
     
class CactusNormalDown(Target):
    """This target does the down pass for the normal phase.
    """
    def __init__(self, options, flowerNames):
        Target.__init__(self, time=0.1)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        self.setFollowOnTarget(CactusNormalRunnable(options=self.options, flowerNames=self.flowerNames))
        makeChildTargets(self.options, self.flowerNames, self, CactusNormalDown)
        
class CactusNormalRunnable(Target):
    """This targets run the normalisation script.
    """
    def __init__(self, flowerNames, options):
        Target.__init__(self, time=1.0)
        self.flowerNames = flowerNames
        self.options = options
        
    def run(self):
        maxNumberOfChains = int(self.options.config.find("normal").attrib["max_number_of_chains"])
        runCactusMakeNormal(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames, maxNumberOfChains=maxNumberOfChains)

############################################################
############################################################
############################################################
#Phylogeny pass
############################################################
############################################################
############################################################
    
class CactusPhylogenyPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
        logger.info("Starting the phylogeny phase")
        if self.options.buildTrees:
            self.addChildTarget(CactusPhylogeny(self.options, [ self.flowerName ]))
        self.setFollowOnTarget(CactusReferencePhase(self.flowerName, self.options))

class CactusPhylogeny(Target):
    """This target does the down pass for the phylogeny phase.
    """
    def __init__(self, options, flowerNames):
        Target.__init__(self, time=1.0)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        runCactusPhylogeny(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames)
        makeChildTargets(self.options, self.flowerNames, self, CactusPhylogeny)

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################
    
class CactusReferencePhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
        logger.info("Starting the reference phase")
        if self.options.buildReference:
            self.addChildTarget(CactusReferenceDown(self.options, [ self.flowerName ]))
        self.setFollowOnTarget(CactusFacesPhase(self.flowerName, self.options))
        
class CactusReferenceDown(Target):
    """This target does the down pass for the reference phase.
    """
    def __init__(self, options, flowerNames):
        Target.__init__(self, time=1.0)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        matchingAlgorithm = self.options.config.find("reference").attrib["matching_algorithm"]
        runCactusReference(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames, matchingAlgorithm=matchingAlgorithm) #We first run the top down phase
        self.setFollowOnTarget(CactusReferenceRunnable(options=self.options, flowerNames=self.flowerNames)) #We second run a bottom up phase
        makeChildTargets(self.options, self.flowerNames, self, CactusReferenceDown)

class CactusReferenceRunnable(Target):
    """This target runs the reference script bottom up (second phase).
    """
    def __init__(self, flowerNames, options):
        Target.__init__(self, time=1.0)
        self.flowerNames = flowerNames
        self.options = options
        
    def run(self):
        runCactusReference(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames, bottomUp=True)
            
############################################################
############################################################
############################################################
#Faces pass
############################################################
############################################################
############################################################
    
class CactusFacesPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
        logger.info("Starting the faces phase")
        if self.options.buildFaces:
            self.addChildTarget(CactusFaces(self.options, [ self.flowerName ]))
        self.setFollowOnTarget(CactusCheckPhase(self.flowerName, self.options))
        
class CactusFaces(Target):
    """This target does the down pass for the faces phase.
    """
    def __init__(self, options, flowerNames):
        Target.__init__(self, time=0.0)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        runCactusAdjacencies(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames)
        makeChildTargets(self.options, self.flowerNames, self, CactusFaces)

############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################
    
class CactusCheckPhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
        logger.info("Starting the verification phase")
        self.addChildTarget(CactusCheck(self.options, [ self.flowerName ]))
        
class CactusCheck(Target):
    """This target does the down pass for the check phase.
    """
    def __init__(self, options, flowerNames):
        Target.__init__(self, time=1.0)
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        runCactusCheck(self.options.cactusDiskDatabaseString, self.flowerNames)
        makeChildTargets(self.options, self.flowerNames, self, CactusCheck)   
      
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    
    parser.add_option("--experiment", dest="experimentFile", help="The file containing a link to the experiment parameters")
    
    parser.add_option("--setupAndBuildAlignments", dest="setupAndBuildAlignments", action="store_true",
                      help="Setup and build alignments then normalise the resulting structure", default=False)
    
    parser.add_option("--buildTrees", dest="buildTrees", action="store_true",
                      help="Build trees", default=False) 
    
    parser.add_option("--buildFaces", dest="buildFaces", action="store_true",
                      help="Build adjacencies", default=False)
    
    parser.add_option("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the flowers", default=False)
    
    options, args = parser.parse_args()
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    options.experimentFile = ET.parse(options.experimentFile).getroot()
    #Get the database string
    options.cactusDiskDatabaseString = ET.tostring(options.experimentFile.find("cactus_disk").find("st_kv_database_conf"))
    #Get the species tree
    options.speciesTree = options.experimentFile.attrib["species_tree"]
    #Parse the config file which contains all the program options
    if options.experimentFile.attrib["config"] == "default":
        options.experimentFile.attrib["config"]=os.path.join(sonTraceRootPath(), "src", "cactus", "pipeline", "cactus_workflow_config.xml")
    #Get the config file for the experiment
    options.config = ET.parse(options.experimentFile.attrib["config"]).getroot()
    #Get the sequences
    sequences = options.experimentFile.attrib["sequences"].split()
    
    logger.info("Parsed the XML options file")
    
    if options.setupAndBuildAlignments:
        baseTarget = CactusSetupPhase(options, sequences)
        logger.info("Going to create alignments and define the cactus tree")
    elif options.buildTrees:
        baseTarget = CactusPhylogenyPhase('0', options)
        logger.info("Starting from phylogeny phase")
    elif options.buildReference:
        baseTarget = CactusReferencePhase('0', options)
        logger.info("Starting from reference phase")
    elif options.buildFaces:
        baseTarget = CactusFacesPhase('0', options)
        logger.info("Starting from faces phase")
    else:
        logger.info("Nothing to do!")
        return
    
    Stack(baseTarget).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.pipeline.cactus_workflow import *
    _test()
    main()
