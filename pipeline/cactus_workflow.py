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

from sonLib.bioio import getTempFile
from sonLib.bioio import system

from workflow.jobTree.lib.bioio import getLogLevelString
from workflow.jobTree.lib.bioio import logger
from workflow.jobTree.lib.bioio import setLoggingFromOptions

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
#The alignment phases, split into the Caf and Bar phases.
############################################################
############################################################
############################################################
    
class CactusAlignmentPhase(Target):
    def __init__(self, flowerName, options, iteration=0):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        self.iteration = iteration
        
    def run(self):
        logger.info("Starting the alignment phase for iteration %i", self.iteration)
        iterations = self.options.config.find("alignment").find("iterations").findall("iteration")
        if self.iteration < len(iterations):
            iterationNode = iterations[self.iteration]
            if iterationNode.attrib["type"] == "blast":
                self.addChildTarget(CactusCafDown(self.options, iterationNode, [ self.flowerName ]))
            else:
                assert iterationNode.attrib["type"] == "base"
                self.addChildTarget(CactusBarDown(self.options, iterationNode, [ self.flowerName ]))
            self.setFollowOnTarget(CactusAlignmentPhase(self.flowerName, self.options, self.iteration+1))
        else:
            self.setFollowOnTarget(CactusNormalPhase(self.flowerName, self.options))
        logger.info("Finished the alignment phase for this iteration")

############################################################
############################################################
############################################################
#The CAF phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

def makeChildTargets(options, extraArgs, flowerNames, target, childTarget, maxSequenceSize=100000, jobNumber=200):
    #Make child jobs
    childFlowerNames = []
    totalSequenceSize = 0
    
    minChildSize = max(1, maxSequenceSize/jobNumber)
    childFlowers = runCactusGetFlowers(options.cactusDiskDatabaseString, flowerNames, target.getLocalTempDir(), logLevel=getLogLevelString())
    totalChildSize = sum([ max(childFlowerSize, minChildSize) for childFlowerName, childFlowerSize in childFlowers])
    avgChildSize = totalChildSize / (totalChildSize / maxSequenceSize + 1)
    
    for childFlowerName, childFlowerSize, in childFlowers:
        assert(childFlowerSize) >= 0
        totalSequenceSize += max(childFlowerSize, minChildSize)
        childFlowerNames.append(childFlowerName)
        if totalSequenceSize >= avgChildSize:
            target.addChildTarget(childTarget(options, extraArgs, childFlowerNames))
            childFlowerNames = []
            totalSequenceSize = 0
    if len(childFlowerNames) > 0:
        target.addChildTarget(childTarget(options, extraArgs, childFlowerNames))

class CactusCafDown(Target):
    """This target does the down pass for the CAF alignment phase.
    """
    def __init__(self, options, iteration, flowerNames):
        Target.__init__(self, time=0.2)
        self.options = options
        self.iteration = iteration
        assert self.iteration.attrib["type"] == "blast"
        self.flowerNames = flowerNames
    
    def run(self):
        makeChildTargets(self.options, self.iteration, self.flowerNames, self, CactusCafDown)
        minSize = int(self.iteration.attrib["min_sequence_size"])
        for childFlowerName, childFlowerSize in runCactusExtendFlowers(self.options.cactusDiskDatabaseString, self.flowerNames, 
                                                                       self.getLocalTempDir(),  logLevel=getLogLevelString()):
            if childFlowerSize > minSize:
                self.addChildTarget(CactusBlastWrapper(self.options, self.iteration, childFlowerName)) 
        
class CactusBlastWrapper(Target):
    """Runs blast on the given flower and passes the resulting alignment to cactus core.
    """
    def __init__(self, options, iteration, flowerName):
        Target.__init__(self, time=0.01)
        self.options = options
        self.iteration = iteration
        self.flowerName = flowerName
    
    def run(self):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", self.getGlobalTempDir())
        logger.info("Got an alignments file")
        
        #Now make the child aligner target
        alignmentNode = self.options.config.find("alignment")
        blastNode = self.iteration.find("blast")
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
        self.setFollowOnTarget(CactusCoreWrapper(self.options, self.iteration, self.flowerName, alignmentFile))
        logger.info("Setup the follow on cactus_core target")

class CactusCoreWrapper(Target):
    """Runs cactus_core upon a given flower and alignment file.
    """
    def __init__(self, options, iteration, flowerName, alignmentFile,):
        Target.__init__(self, time=100, memory=4294967295) #Request 2^32 (4 gigs of ram)
        self.options = options
        self.iteration = iteration
        self.flowerName = flowerName
        self.alignmentFile = alignmentFile
    
    def run(self):
        logger.info("Starting the core wrapper target")
        coreParameters = self.iteration.find("core")
        runCactusCore(cactusDiskDatabaseString=self.options.cactusDiskDatabaseString,
                      alignmentFile=self.alignmentFile, 
                      flowerName=self.flowerName,
                      logLevel=getLogLevelString(), 
                      annealingRounds=[ int(i) for i in coreParameters.attrib["annealingRounds"].split() ],
                      deannealingRounds=[ int(i) for i in coreParameters.attrib["deannealingRounds"].split() ],
                      alignRepeatsAtRound=float(coreParameters.attrib["alignRepeatsAtRound"]),
                      trim=[ int(i) for i in coreParameters.attrib["trim"].split() ],
                      minimumTreeCoverage=float(coreParameters.attrib["minimumTreeCoverage"]),
                      minimumBlockLength=float(coreParameters.attrib["minimumBlockLength"]),
                      adjacencyComponentOverlap=int(coreParameters.attrib["adjacencyComponentOverlap"]))
        logger.info("Ran the cactus core program okay")
        
############################################################
############################################################
############################################################
#The BAR phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

class CactusBarDown(Target):
    """This target does the down pass for the BAR alignment phase.
    """
    def __init__(self, options, iteration, flowerNames):
        Target.__init__(self, time=0.2)
        self.options = options
        self.flowerNames = flowerNames
        self.iteration = iteration
    
    def run(self):
        children = []
        makeChildTargets(self.options, self.iteration, self.flowerNames, self, CactusBarDown)
        for childFlowerName, childFlowerSize in runCactusExtendFlowers(self.options.cactusDiskDatabaseString, self.flowerNames, 
                                                              self.getLocalTempDir(),  logLevel=getLogLevelString()):
            assert childFlowerSize >= 0
            children.append(childFlowerName)
            if len(children) > 50:
                self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, self.iteration, children)) 
                children = []
        if len(children) > 0:
            self.addChildTarget(CactusBaseLevelAlignerWrapper(self.options, self.iteration, children)) 
     
class CactusBaseLevelAlignerWrapper(Target):
    """Runs cactus_baseAligner (the BAR algorithm implementation.
    """
    #We split, to deal with cleaning up the alignment file
    def __init__(self, options, iteration, flowerNames):
        Target.__init__(self, time=10)
        self.options = options
        self.iteration = iteration
        self.flowerNames = flowerNames
    
    def run(self):
        assert self.iteration.attrib["type"] =="base"
        runCactusBaseAligner(self.options.cactusDiskDatabaseString, self.flowerNames, getLogLevelString(),
                             maximumLength=float(self.iteration.attrib["banding_limit"]),
                             spanningTrees=float(self.iteration.attrib["spanning_trees"]),
                             gapGamma=float(self.iteration.attrib["gap_gamma"]),
                             useBanding=bool(int(self.iteration.attrib["use_banding"])),
                             bandingSize=int(self.iteration.attrib["banding_size"]))
        
############################################################
############################################################
############################################################
#Normalisation pass
############################################################
############################################################
############################################################

    
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
        self.addChildTarget(CactusNormalDown(self.options, None, [ self.flowerName ]))
        if self.normalisationRounds-1 > 0:
            self.setFollowOnTarget(CactusNormalPhase(self.flowerName, self.options, self.normalisationRounds-1))
        else:
            self.setFollowOnTarget(CactusPhylogenyPhase(self.flowerName, self.options))
     
class CactusNormalDown(Target):
    """This target does the down pass for the normal phase.
    """
    def __init__(self, options, extras, flowerNames):
        Target.__init__(self, time=0.2)
        assert extras == None #We currently don't use this argument
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        self.setFollowOnTarget(CactusNormalRunnable(options=self.options, flowerNames=self.flowerNames))
        makeChildTargets(self.options, None, self.flowerNames, self, CactusNormalDown)
        
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
            self.addChildTarget(CactusPhylogeny(self.options, None, [ self.flowerName ]))
        self.setFollowOnTarget(CactusReferencePhase(self.flowerName, self.options))

class CactusPhylogeny(Target):
    """This target does the down pass for the phylogeny phase.
    """
    def __init__(self, options, extras, flowerNames):
        Target.__init__(self, time=1.0)
        assert extras == None #We currently don't use this argument
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        runCactusPhylogeny(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames)
        makeChildTargets(self.options, None, self.flowerNames, self, CactusPhylogeny)

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
            self.addChildTarget(CactusReferenceDown(self.options, None, [ self.flowerName ]))
        self.setFollowOnTarget(CactusFacesPhase(self.flowerName, self.options))
        
class CactusReferenceDown(Target):
    """This target does the down pass for the reference phase.
    """
    def __init__(self, options, extras, flowerNames):
        Target.__init__(self, time=1.0)
        assert extras == None #We currently don't use this argument
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        matchingAlgorithm = self.options.config.find("reference").attrib["matching_algorithm"]
        runCactusReference(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames, matchingAlgorithm=matchingAlgorithm) #We first run the top down phase
        self.setFollowOnTarget(CactusReferenceRunnable(options=self.options, flowerNames=self.flowerNames)) #We second run a bottom up phase
        makeChildTargets(self.options, None, self.flowerNames, self, CactusReferenceDown)

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
            self.addChildTarget(CactusFaces(self.options, None, [ self.flowerName ]))
        self.setFollowOnTarget(CactusCheckPhase(self.flowerName, self.options))
        
class CactusFaces(Target):
    """This target does the down pass for the faces phase.
    """
    def __init__(self, options, extras, flowerNames):
        Target.__init__(self, time=0.0)
        assert extras == None #We currently don't use this argument
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        runCactusAdjacencies(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames)
        makeChildTargets(self.options, None, self.flowerNames, self, CactusFaces)

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
        self.addChildTarget(CactusCheck(self.options, None, [ self.flowerName ]))
        
class CactusCheck(Target):
    """This target does the down pass for the check phase.
    """
    def __init__(self, options, extras, flowerNames):
        Target.__init__(self, time=1.0)
        assert extras == None #We currently don't use this argument
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        runCactusCheck(self.options.cactusDiskDatabaseString, self.flowerNames)
        makeChildTargets(self.options, None, self.flowerNames, self, CactusCheck)   
      
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
    setLoggingFromOptions(options)
    
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
