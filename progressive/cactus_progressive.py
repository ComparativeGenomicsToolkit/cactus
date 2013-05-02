#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Wrapper to run the cactus_workflow progressively, using the input species tree as a guide

tree.  
"""

import os
import xml.etree.ElementTree as ET
import math
from optparse import OptionParser
from collections import deque
import random
from itertools import izip
from shutil import move
import copy
from time import sleep

from sonLib.bioio import getTempFile
from sonLib.bioio import printBinaryTree
from sonLib.bioio import system

from jobTree.src.bioio import getLogLevelString
from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions

from cactus.shared.common import cactusRootPath
  
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 

from cactus.pipeline.cactus_workflow import CactusPreprocessorPhase
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import getOptionalAttrib
from cactus.pipeline.cactus_workflow import findRequiredNode

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
        
class ProgressiveDown(Target):
    def __init__(self, options, project, event, schedule):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
    
    def run(self):
        logger.info("Progressive Down: " + self.event)
        
        if not self.options.nonRecursive:
            deps = self.schedule.deps(self.event)
            for child in deps:
                if child in self.project.expMap:
                    self.addChildTarget(ProgressiveDown(self.options,
                                                        self.project, child, 
                                                        self.schedule))
        
        self.setFollowOnTarget(ProgressiveUp(self.options, self.project, self.event))
        
class ProgressiveUp(Target):
    def __init__(self, options, project, event):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.event = event
    
    def run(self):
        logger.info("Progressive Up: " + self.event)

        # open up the experiment
        # note that we copy the path into the options here
        self.options.experimentFile = self.project.expMap[self.event]
        expXml = ET.parse(self.options.experimentFile).getroot()
        experiment = ExperimentWrapper(expXml)
        configXml = ET.parse(experiment.getConfigPath()).getroot()
        configWrapper = ConfigWrapper(configXml)

        # need at least 3 processes for every event when using ktserver:
        # 1 proc to run jobs, 1 proc to run server, 1 proc to run 2ndary server
        if experiment.getDbType() == "kyoto_tycoon":            
            maxParallel = min(len(self.project.expMap),
                             configWrapper.getMaxParallelSubtrees()) 
            if self.options.batchSystem == "singleMachine":
                if int(self.options.maxThreads) < maxParallel * 3:
                    raise RuntimeError("At least %d threads are required to handle up to %d events using kyoto tycoon. Either increase the number of threads using the --maxThreads option or decrease the number of parallel jobs (currently %d) by adjusting max_parallel_subtrees in the config file" % (maxParallel * 3, maxParallel, configWrapper.getMaxParallelSubtrees()))
            else:
                if int(self.options.maxJobs) < maxParallel * 3:
                    raise RuntimeError("At least %d concurrent jobs are required to handle up to %d events using kyoto tycoon. Either increase the number of jobs using the --maxJobs option or decrease the number of parallel jobs (currently %d) by adjusting max_parallel_subtrees in the config file" % (maxParallel * 3, maxParallel, configWrapper.getMaxParallelSubtrees()))
                    
        # take union of command line options and config options for hal and reference
        if self.options.buildReference == False:
            refNode = findRequiredNode(configXml, "reference")
            self.options.buildReference = getOptionalAttrib(refNode, "buildReference", bool, False)
        halNode = findRequiredNode(configXml, "hal")
        if self.options.buildHal == False:
            self.options.buildHal = getOptionalAttrib(halNode, "buildHal", bool, False)
        if self.options.buildFasta == False:
            self.options.buildFasta = getOptionalAttrib(halNode, "buildFasta", bool, False)

        # get parameters that cactus_workflow stuff wants
        workFlowArgs = CactusWorkflowArguments(self.options)
        # copy over the options so we don't trail them around
        workFlowArgs.skipAlignments = self.options.skipAlignments
        workFlowArgs.buildReference = self.options.buildReference
        workFlowArgs.buildHal = self.options.buildHal
        workFlowArgs.buildFasta = self.options.buildFasta
        workFlowArgs.overwrite = self.options.overwrite
        workFlowArgs.globalLeafEventSet = self.options.globalLeafEventSet
        
        experiment = ExperimentWrapper(workFlowArgs.experimentNode)

        donePath = os.path.join(os.path.dirname(workFlowArgs.experimentFile), "DONE")
        doneDone = os.path.isfile(donePath)
        refDone = not workFlowArgs.buildReference or os.path.isfile(experiment.getReferencePath())
        halDone = not workFlowArgs.buildHal or (os.path.isfile(experiment.getHALFastaPath()) and
                                                os.path.isfile(experiment.getHALPath()))
                                                               
        if not workFlowArgs.overwrite and doneDone and refDone and halDone:
            self.logToMaster("Skipping %s because it is already done and overwrite is disabled" %
                             self.event)
        else:
            system("rm -f %s" % donePath)
            # delete database files if --setupAndBuildAlignments
            # and overwrite specified (or if reference not present)
            if not self.options.skipAlignments:
                dbPath = os.path.join(experiment.getDbDir(), 
                                      experiment.getDbName())
                seqPath = os.path.join(experiment.getDbDir(), "sequences")
                system("rm -f %s* %s %s" % (dbPath, seqPath, 
                                            experiment.getReferencePath()))

            if not workFlowArgs.skipAlignments:
                self.addChildTarget(CactusPreprocessorPhase(cactusWorkflowArguments=workFlowArgs,
                                                            phaseName="preprocessor"))
                logger.info("Going to create alignments and define the cactus tree")

def main():    
    usage = "usage: %prog [options] <multicactus project>"
    description = "Progressive version of cactus_workflow"
    parser = OptionParser(usage=usage, description=description)
    Stack.addJobTreeOptions(parser)
    addCactusWorkflowOptions(parser)
    
    parser.add_option("--nonRecursive", dest="nonRecursive", action="store_true",
                      help="Only process given event (not children) [default=False]", 
                      default=False)
    
    parser.add_option("--event", dest="event", 
                      help="Target event to process [default=root]", default=None)
    
    parser.add_option("--overwrite", dest="overwrite", action="store_true",
                      help="Recompute and overwrite output files if they exist [default=False]",
                      default=False)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)

    if len(args) != 1:
        parser.print_help()
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    project = MultiCactusProject()
    project.readXML(args[0])
    schedule = Schedule()
    schedule.loadProject(project)
    schedule.compute()
    if options.event == None:
        options.event = project.mcTree.getRootName()
    assert options.event in project.expMap
    leafNames = [project.mcTree.getName(i) for i in project.mcTree.getLeaves()]
    options.globalLeafEventSet = set(leafNames)
    
    baseTarget = ProgressiveDown(options, project, options.event, schedule)
    Stack(baseTarget).startJobTree(options)
    
   
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.progressive.cactus_progressive import *
    _test()
    main()
