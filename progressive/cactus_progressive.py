#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Wrapper to run the cactus_workflow progressively, using the input species tree as a guide

tree.  The input and options are the same as cactus_worklfow, with a few addtions:  The --cladeSize

option (which should get moved to config.xml?), specifies the maximum number of sequences to 

align at once;  The --buildMAF and --joinMAF options can be used to create and merge MAFs from the

cacti.   The --buildReference option is removed, as references *are always* built when 

--setupAndBuildAlignments is specified.  If the kyoto tycoon DB is used, it is best to use the helper

script, preKtserverDbs.py to set up the server.
"""

import os
import xml.etree.ElementTree as ET
import math
from optparse import OptionParser
from collections import deque
import random
from itertools import izip
from shutil import move

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree

from jobTree.src.bioio import getLogLevelString
from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions

from cactus.shared.common import cactusRootPath
  
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 

from cactus.pipeline.cactus_workflow import CactusSetupPhase
from cactus.pipeline.cactus_workflow import CactusPhylogenyPhase
from cactus.pipeline.cactus_workflow import CactusReferencePhase
from cactus.pipeline.cactus_workflow import CactusFacesPhase

from cactus.progressive.progressiveSplitUtils import getCladeLeaves
from cactus.progressive.progressiveSplitUtils import nameUnlabeledInternalNodes
from cactus.progressive.progressiveSplitUtils import createSeqeunceMap
from cactus.progressive.progressiveSplitUtils import getCladeSequences
from cactus.progressive.progressiveSplitUtils import createCladeOptions
from cactus.progressive.progressiveSplitUtils import createCladeFileStructure
from cactus.progressive.progressiveSplitUtils import getMAFGeneratorOptions
from cactus.progressive.progressiveSplitUtils import getMAFJoinOptions
from cactus.progressive.progressiveSplitUtils import getReferenceSeqOptions
from cactus.progressive.progressiveSplitUtils import createProgWorkDir
from cactus.progressive.progressiveSplitUtils import addDotsToMAF
from cactus.progressive.progressiveSplitUtils import getCladeMAFJoinTempPath
from cactus.progressive.progressiveSplitUtils import getCladeMAFPath

from cactus.progressive.progressiveKtserver import isKyotoTycoon
from cactus.progressive.progressiveKtserver import spawnLocalKtserver
from cactus.progressive.progressiveKtserver import killLocalKtserver

class ProgressiveSetup(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0.0002)
        self.options = options
        self.sequences = sequences
       
    def run(self):
        logger.info("Starting Progressive setup phase target")
        
        # get the full phylogeny and label internal nodes (default: Anc0,Anc1....)
        root = newickTreeParser(self.options.speciesTree)
        nameUnlabeledInternalNodes(root)
        createProgWorkDir(self.options, root)
              
        # create map between node names and fasta file paths
        self.options.lookup = createSeqeunceMap(root, self.options, self.sequences)
        
        # start recursive alignment at root
        self.setFollowOnTarget(ProgressiveAlignmentDown(root, self.options, self.sequences))
        
class ProgressiveAlignmentDown(Target):
    def __init__(self, root, options, sequences):
        Target.__init__(self)
        self.root = root
        self.options = options
        self.sequences = sequences
    
    def run(self):
        logger.info("Progressive Down: " + self.root.iD)
        
        # get the bottom nodes of the clade
        leaves = getCladeLeaves(self.root, self.options)
        cladeSequences = getCladeSequences(leaves, self.options)
        
        # NOTE TO SELF:  CONSIDER CASE OF INTERNAL NODE WITH DEGREE 2
        if len(leaves) > 1:
            for leaf in leaves:
                self.addChildTarget(ProgressiveAlignmentDown(leaf, self.options, cladeSequences))
            self.setFollowOnTarget(ProgressiveAlignmentUp(self.root, leaves, self.options, cladeSequences))
        
class ProgressiveAlignmentUp(Target):
    def __init__(self, root, leaves, options, sequences):
        Target.__init__(self)
        self.root = root
        self.leaves = leaves
        self.options = options
        self.sequences = sequences
    
    def run(self):
        logger.info("Doing progressive Up on " + self.root.iD)
        
        # create options object for clade alignment
        cladeOptions = createCladeOptions(self.root, self.leaves, self.options)
       
        # spawn ktserver and update cladeOptions 
        if isKyotoTycoon(cladeOptions) and cladeOptions.autoKtserver:
            spawnLocalKtserver(cladeOptions)
        
        if self.options.setupAndBuildAlignments:
            createCladeFileStructure(self.root, self.leaves, self.options, cladeOptions)       
            self.addChildTarget(CactusSetupPhase(cladeOptions, self.sequences))
            logger.info("Going to create alignments and define the cactus tree")
        elif self.options.buildTrees:
            self.addChildTarget(CactusPhylogenyPhase('0', cladeOptions))
            logger.info("Starting from phylogeny phase")
        #elif self.options.buildReference:
        #    self.addChildTarget(CactusReferencePhase('0', cladeOptions))
        #    logger.info("Starting from reference phase")
        elif self.options.buildFaces:
            self.addChildTarget(CactusFacesPhase('0', cladeOptions))
            logger.info("Starting from faces phase")
        else:
            logger.info("Nothing to do!")
        
        self.setFollowOnTarget(ProgressiveExtractReference(self.root, self.leaves, self.options,
                                                           cladeOptions, self.sequences))
        
class ProgressiveExtractReference(Target):
    def __init__(self, root, leaves, options, cladeOptions, sequences):
        Target.__init__(self)
        self.root = root
        self.leaves = leaves
        self.options = options
        self.cladeOptions = cladeOptions
        self.sequences = sequences

    def run(self):        
        assert(self.root.left is not None or self.root.right is not None)
        
        if self.options.buildReference:
            logger.info("Starting Reference Extract Phase")
            refOptions = getReferenceSeqOptions(self.root, self.options, self.cladeOptions)
            refStatus = os.system("cactus_getReferenceSeq " + refOptions)
            assert refStatus == 0
            
        self.setFollowOnTarget(ProgressiveBuildMAF(self.root, self.leaves, self.options,
                                                      self.cladeOptions, self.sequences))
        
                
class ProgressiveBuildMAF(Target):
    def __init__(self, root, leaves, options, cladeOptions, sequences):
        Target.__init__(self)
        self.root = root
        self.leaves = leaves
        self.options = options
        self.cladeOptions = cladeOptions
        self.sequences = sequences
    
    def run(self):
        
        assert(self.root.left is not None or self.root.right is not None)
        
        if self.options.buildMAF:
            logger.info("Starting MAF Build phase")
            mafOptions = getMAFGeneratorOptions(self.root, self.options, self.cladeOptions)
            mafStatus = os.system("cactus_MAFGenerator " + mafOptions)
            assert mafStatus == 0
            # mafJoin requires periods. add a .1 to all sequence names without a period
            # addDotsToMAF(self.root, self.options)
             
        
        self.setFollowOnTarget(ProgressiveJoinMAF(self.root, self.leaves, self.options,
                                                  self.cladeOptions, self.sequences))


class ProgressiveJoinMAF(Target):
    def __init__(self, root, leaves, options, cladeOptions, sequences):
        Target.__init__(self)
        self.root = root
        self.leaves = leaves
        self.options = options
        self.cladeOptions = cladeOptions
        self.sequences = sequences
    
    def run(self):
        
        # don't need the ktserver anymore, so we kill it
        if isKyotoTycoon(cladeOptions) and options.autoKtserver:
            killLocalKtserver(cladeOptions)
        
        if self.options.joinMAF:
            assert(self.root.left is not None or self.root.right is not Note)
            assert(len(self.leaves) == len(self.sequences))
            logger.info("Starting MAF Join phase")
            
            rootIsTreeMAF = False
            for leaf in self.leaves:
                if leaf.left is not None or leaf.right is not None:
                    mafOptions = getMAFJoinOptions(self.root, leaf, 
                                                   self.options, self.cladeOptions,
                                                   rootIsTreeMAF)
                                 
                    mafStatus = os.system("mafJoin " + mafOptions)
                    assert mafStatus == 0
                    move(getCladeMAFJoinTempPath(self.root, leaf, self.options), 
                         getCladeMAFPath(self.root, self.options, True))
                    rootIsTreeMAF = True    
                    
                            
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
                      help="Deprecated", default=False)
    
    parser.add_option("--buildMAF", dest="buildMAF", action="store_true",
                     help="Create a MAF file from the cactus and reference", default=False)
    
    parser.add_option("--joinMAF", dest="joinMAF", action="store_true",
                     help="Progressively join all cactus MAFs", default=False)
    
    parser.add_option("--cladeSize", dest="cladeSize", type="int", 
                      help="Max number of sequences to align at a time", default=2)
    
    parser.add_option("--autoKtserver", dest="autoKtserver", action="store_true",
                      help="Autmatically manage ktservers (if cactus disk is KT) [default=True]",
                      default=True)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    # cannot align without building reference and vice versa in progressive mode
    options.buildReference = options.setupAndBuildAlignments
    
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
    
    baseTarget = ProgressiveSetup(options, sequences)
    Stack(baseTarget).startJobTree(options)
    
   
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.progressive.cactus_progressive import *
    _test()
    main()