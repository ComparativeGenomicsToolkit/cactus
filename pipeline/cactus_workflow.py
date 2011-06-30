#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
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
import sys
import xml.etree.ElementTree as ET
import math
from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser

from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions

from cactus.shared.common import cactusRootPath
  
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 

from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusCore
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import runCactusPhylogeny
from cactus.shared.common import runCactusAdjacencies
from cactus.shared.common import runCactusBaseAligner
from cactus.shared.common import runCactusMakeNormal 
from cactus.shared.common import runCactusReference
from cactus.shared.common import runCactusAddReferenceCoordinates
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

def getLongestPath(node, distance=0.0):
    """Identify the longest path from root to leaves of the species tree
    and add the min-distance.
    """
    i, j = distance, distance
    if node.left != None:
        i = getLongestPath(node.left, node.left.distance) + distance
    if node.right != None:  
        j = getLongestPath(node.right, node.right.distance) + distance
    return max(i, j)

def inverseJukesCantor(d):
    """Takes a substitution distance and calculates the number of expected changes per site (inverse jukes cantor)
    
    d = -3/4 * log(1 - 4/3 * p)
    exp(-4/3 * d) = 1 - 4/3 * p
    4/3 * p = 1 - exp(-4/3 * d)
    p = 3/4 * (1 - exp(-4/3 * d))
    
    >>> inverseJukesCantor(0.5)
    0.36493716072555599
    >>> inverseJukesCantor(1.0)
    0.55230214641320496
    >>> inverseJukesCantor(10.0)
    0.74999878530240571
    >>> inverseJukesCantor(100000.0)
    0.75
    """
    assert d >= 0.0
    return 0.75 * (1 - math.exp(-d * 4.0/3.0))

def getOptionalAttrib(node, attribName, default=None):
    """Get an optional attrib, or None, if not set.
    """
    if node.attrib.has_key(attribName):
        return node.attrib[attribName]
    return default

class CactusSetupPhase(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0.0002)
        self.options = options 
        self.sequences = sequences 
            
    def modifyConfig(self):
        #Add the identity clause into the blast strings
        alignmentNode = self.options.config.find("alignment")
        if int(alignmentNode.find("blast_misc").attrib["filterByIdentity"]):
            longestPath = getLongestPath(newickTreeParser(self.options.speciesTree))
            adjustedPath = float(alignmentNode.find("blast_misc").attrib["identityRatio"]) * longestPath + float(alignmentNode.find("blast_misc").attrib["minimumDistance"])
            identity = str(100 - int(100 * inverseJukesCantor(adjustedPath)))
            logger.info("The blast stage will filter by identity, the calculated minimum identity is %s from a longest path of %s and an adjusted path of %s" % (identity, longestPath, adjustedPath))
            for iterationNode in alignmentNode.find("iterations").findall("iteration"):
                if iterationNode.attrib["type"] == "blast":
                    blastNode = iterationNode.find("blast")
                    assert "IDENTITY" in blastNode.attrib["blastString"]
                    blastNode.attrib["blastString"] = blastNode.attrib["blastString"].replace("IDENTITY", identity)
                    assert "IDENTITY" in blastNode.attrib["selfBlastString"]
                    blastNode.attrib["selfBlastString"] = blastNode.attrib["selfBlastString"].replace("IDENTITY", identity)
                else:
                    assert iterationNode.attrib["type"] == "base"

    def run(self):
        logger.info("Starting setup phase target")
        #Modify the config options
        self.modifyConfig()
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
                       self.options.speciesTree)
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
            assert int(iterationNode.attrib["number"]) == self.iteration
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

MAX_SEQUENCE_SIZE=1000000
MAX_JOB_NUMBER=1000

def makeTargets(options, extraArgs, flowersAndSizes, parentTarget, target, 
                maxSequenceSize=MAX_SEQUENCE_SIZE, jobNumber=MAX_JOB_NUMBER,
                ignoreFlowersLessThanThisSize=0, ignoreFlowersGreaterThanThisSize=sys.maxint):
    """Make a set of targets for a given set of flowers.
    """
    #Make child jobs
    flowerNames = []
    totalSequenceSize = 0.0
    
    minChildSize = max(1, float(maxSequenceSize)/jobNumber)
    totalChildSize = sum([ max(flowerSize, minChildSize) for flowerName, flowerSize in flowersAndSizes ])
    
    for flowerName, flowerSize, in flowersAndSizes:
        assert(flowerSize) >= 0
        if flowerSize >= ignoreFlowersLessThanThisSize:
            if flowerSize <= ignoreFlowersGreaterThanThisSize:
                totalSequenceSize += max(flowerSize, minChildSize)
                flowerNames.append(flowerName)
                if totalSequenceSize >= maxSequenceSize: 
                    parentTarget.addChildTarget(target(options, extraArgs, flowerNames))
                    flowerNames = []
                    totalSequenceSize = 0.0
            else:
                parentTarget.logToMaster("Ignoring a flower %s of size %s" % (flowerName, flowerSize))
    if len(flowerNames) > 0:
        parentTarget.addChildTarget(target(options, extraArgs, flowerNames))

def makeChildTargets(options, extraArgs, flowerNames, target, childTarget, maxSequenceSize=MAX_SEQUENCE_SIZE, jobNumber=MAX_JOB_NUMBER,
                     ignoreFlowersLessThanThisSize=0, ignoreFlowersGreaterThanThisSize=sys.maxint):
    """Make a set of child targets for a given set of parent flowers.
    """
    childFlowers = runCactusGetFlowers(options.cactusDiskDatabaseString, flowerNames, target.getLocalTempDir())
    makeTargets(options, extraArgs, childFlowers, target, childTarget, maxSequenceSize, jobNumber,
                ignoreFlowersLessThanThisSize, ignoreFlowersGreaterThanThisSize)

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
        ignoreFlowersLessThanThisSize = int(self.iteration.attrib["min_sequence_size"])
        ignoreFlowersGreaterThanThisSize = int(getOptionalAttrib(self.iteration, "max_sequence_size", sys.maxint))
        makeChildTargets(self.options, self.iteration, self.flowerNames, self, CactusCafDown, 
                         ignoreFlowersLessThanThisSize=ignoreFlowersLessThanThisSize)
        for childFlowerName, childFlowerSize in runCactusExtendFlowers(self.options.cactusDiskDatabaseString, self.flowerNames, 
                                                                       self.getLocalTempDir()):
            if childFlowerSize >= ignoreFlowersLessThanThisSize:
                if childFlowerSize <= ignoreFlowersGreaterThanThisSize:
                    self.addChildTarget(CactusBlastWrapper(self.options, self.iteration, childFlowerName))
                else:
                    self.logToMaster("Ignoring a flower from CAF alignment %s of size %s" % (childFlowerName, childFlowerSize))

def getOption(node, attribName, default):
    if node.attrib.has_key(attribName):
        return node.attrib[attribName]
    return default
        
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
                                          self.flowerName, alignmentFile, blastOptions,
                                          minimumSequenceLength=int(getOption(blastMiscNode, "minimumSequenceLength", 0))))
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
                      annealingRounds=[ int(i) for i in coreParameters.attrib["annealingRounds"].split() ],
                      deannealingRounds=[ int(i) for i in coreParameters.attrib["deannealingRounds"].split() ],
                      alignRepeatsAtRound=float(coreParameters.attrib["alignRepeatsAtRound"]),
                      trim=[ int(i) for i in coreParameters.attrib["trim"].split() ],
                      minimumTreeCoverage=float(coreParameters.attrib["minimumTreeCoverage"]),
                      blockTrim=float(coreParameters.attrib["blockTrim"]),
                      ignoreAllChainsLessThanMinimumTreeCoverage=bool(coreParameters.attrib["ignoreAllChainsLessThanMinimumTreeCoverage"]),
                      minimumBlockDegree=int(coreParameters.attrib["minimumBlockDegree"]),
                      requiredSpecies=self.options.requiredSpecies)
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
        Target.__init__(self, time=1.0)
        self.options = options
        self.flowerNames = flowerNames
        self.iteration = iteration
    
    def run(self):
        children = []
        makeChildTargets(self.options, self.iteration, self.flowerNames, self, CactusBarDown)
        ignoreFlowersGreaterThanThisSize = int(getOptionalAttrib(self.iteration, "max_sequence_size", sys.maxint))
        childFlowersAndSizes = runCactusExtendFlowers(self.options.cactusDiskDatabaseString, self.flowerNames, 
                                                              self.getLocalTempDir())
        makeTargets(self.options, self.iteration, childFlowersAndSizes, self, CactusBaseLevelAlignerWrapper, maxSequenceSize=10000, 
                    ignoreFlowersGreaterThanThisSize=ignoreFlowersGreaterThanThisSize)
     
class CactusBaseLevelAlignerWrapper(Target):
    """Runs cactus_baseAligner (the BAR algorithm implementation.
    """
    #We split, to deal with cleaning up the alignment file
    def __init__(self, options, iteration, flowerNames):
        Target.__init__(self, time=30)
        self.options = options
        self.iteration = iteration
        self.flowerNames = flowerNames
    
    def run(self):
        assert self.iteration.attrib["type"] == "base"
        runCactusBaseAligner(self.options.cactusDiskDatabaseString, self.flowerNames, 
                             maximumLength=float(self.iteration.attrib["banding_limit"]),
                             spanningTrees=float(self.iteration.attrib["spanning_trees"]),
                             gapGamma=float(self.iteration.attrib["gap_gamma"]),
                             useBanding=bool(int(self.iteration.attrib["use_banding"])),
                             maxBandingSize=int(self.iteration.attrib["max_banding_size"]),
                             minBandingSize=int(self.iteration.attrib["min_banding_size"]),
                             minBandingConstraintDistance=int(self.iteration.attrib["min_banding_constraint_distance"]),
                             minTraceBackDiag=int(self.iteration.attrib["min_trace_back_diag"]),
                             minTraceGapDiags=int(self.iteration.attrib["min_trace_gap_diags"]),
                             constraintDiagonalTrim=int(self.iteration.attrib["constraint_diagonal_trim"]),
                             minimumBlockDegree=int(self.iteration.attrib["minimumBlockDegree"]),
                             alignAmbiguityCharacters=bool(int(self.iteration.attrib["alignAmbiguityCharacters"])),
                             requiredSpecies=self.options.requiredSpecies)
        
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
        Target.__init__(self, time=2.0)
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
        Target.__init__(self, time=3.0)
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
        Target.__init__(self, time=5.0)
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

def getReferenceEventString(config):
    referenceNode = config.find("reference")
    if referenceNode.attrib.has_key("reference"):
        return referenceNode.attrib["reference"]
    return None
    
class CactusReferencePhase(Target):
    def __init__(self, flowerName, options):
        Target.__init__(self, time=0.0002)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
        logger.info("Starting the reference phase")
        if self.options.buildReference:
            self.addChildTarget(CactusReferenceDown(self.options, None, [ self.flowerName ]))
            self.setFollowOnTarget(CactusSetReferenceCoordinates(self.flowerName, self.options))
        else:
            self.setFollowOnTarget(CactusFacesPhase(self.flowerName, self.options))
        
class CactusReferenceDown(Target):
    """This target does the down pass for the reference phase.
    """
    def __init__(self, options, extras, flowerNames):
        Target.__init__(self, time=2.0)
        assert extras == None #We currently don't use this argument
        self.options = options
        self.flowerNames = flowerNames
    
    def run(self):
        matchingAlgorithm = self.options.config.find("reference").attrib["matching_algorithm"]
        referenceEventString = getReferenceEventString(self.options.config)
        runCactusReference(self.options.cactusDiskDatabaseString, flowerNames=self.flowerNames, matchingAlgorithm=matchingAlgorithm, 
                           maxNumberOfChainsToSolvePerRound=getOptionalAttrib(self.options.config.find("reference"), "maxNumberOfChainsToSolvePerRound"),
                           referenceEventString=referenceEventString)
        makeChildTargets(self.options, None, self.flowerNames, self, CactusReferenceDown)

class CactusSetReferenceCoordinates(Target):
    """Fills in the coordinates, once a reference is added.
    """
    def __init__(self, flowerName, options):
        Target.__init__(self, time=100.0)
        self.flowerName = flowerName
        self.options = options
        
    def run(self):
        referenceEventString = getReferenceEventString(self.options.config)
        runCactusAddReferenceCoordinates(self.options.cactusDiskDatabaseString, referenceEventString=referenceEventString)
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
        Target.__init__(self, time=5.0)
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
        options.experimentFile.attrib["config"] = os.path.join(cactusRootPath(), "pipeline", "cactus_workflow_config.xml")
    else:
        logger.info("Using user specified experiment file")
    #Get the config file for the experiment
    options.config = ET.parse(options.experimentFile.attrib["config"]).getroot()
    #Get the sequences
    sequences = options.experimentFile.attrib["sequences"].split()
    #Get any list of 'required species' for the blocks of the cactus.
    options.requiredSpecies = None
    if options.experimentFile.attrib.has_key("required_species"):
        options.requiredSpecies = options.experimentFile.attrib["required_species"]
    
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
