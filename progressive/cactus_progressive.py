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
from argparse import ArgumentParser
from collections import deque
import random
from itertools import izip
from shutil import move
import copy
from time import sleep

from sonLib.bioio import getTempFile
from sonLib.bioio import printBinaryTree
from sonLib.bioio import system

from toil.lib.bioio import getLogLevelString
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions

from cactus.shared.common import cactusRootPath
from cactus.shared.common import getOptionalAttrib
  
from toil.job import Job

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import findRequiredNode
from cactus.pipeline.cactus_workflow import CactusSetupPhase
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
        
class ProgressiveDown(Job):
    def __init__(self, options, project, event, schedule):
        Job.__init__(self)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
    
    def run(self, fileStore):
        logger.info("Progressive Down: " + self.event)
        
        if not self.options.nonRecursive:
            deps = self.schedule.deps(self.event)
            for child in deps:
                self.addChild(ProgressiveDown(self.options,
                                                    self.project, child, 
                                                    self.schedule))
        
        self.addFollowOn(ProgressiveNext(self.options, self.project, self.event,
                                               self.schedule))

class ProgressiveNext(Job):
    def __init__(self, options, project, event, schedule):
        Job.__init__(self)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
    
    def run(self, fileStore):
        logger.info("Progressive Next: " + self.event)

        if not self.schedule.isVirtual(self.event):
            self.addChild(ProgressiveUp(self.options, self.project, self.event))
        followOnEvent = self.schedule.followOn(self.event)
        if followOnEvent is not None:
            self.addChild(ProgressiveDown(self.options, self.project, followOnEvent,
                                                self.schedule))
    
class ProgressiveUp(Job):
    def __init__(self, options, project, event):
        Job.__init__(self)
        self.options = options
        self.project = project
        self.event = event
    
    def run(self, fileStore):
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
                    raise RuntimeError("At least %d threads are required (only %d were specified) to handle up to %d events using kyoto tycoon. Either increase the number of threads using the --maxThreads option or decrease the number of parallel jobs (currently %d) by adjusting max_parallel_subtrees in the config file" % (maxParallel * 3, self.options.maxThreads, maxParallel, configWrapper.getMaxParallelSubtrees()))
            else:
                if int(self.options.maxCpus) < maxParallel * 3:
                    raise RuntimeError("At least %d concurrent cpus are required to handle up to %d events using kyoto tycoon. Either increase the number of cpus using the --maxCpus option or decrease the number of parallel jobs (currently %d) by adjusting max_parallel_subtrees in the config file" % (maxParallel * 3, maxParallel, configWrapper.getMaxParallelSubtrees()))
                    
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
            # delete database 
            # and overwrite specified (or if reference not present)
            dbPath = os.path.join(experiment.getDbDir(), 
                                  experiment.getDbName())
            seqPath = os.path.join(experiment.getDbDir(), "sequences")
            system("rm -f %s* %s %s" % (dbPath, seqPath, 
                                        experiment.getReferencePath()))

            if workFlowArgs.configWrapper.getDoTrimStrategy() and workFlowArgs.outgroupEventNames is not None:
                # Use the trimming strategy to blast ingroups vs outgroups.
                self.addChild(CactusTrimmingBlastPhase(cactusWorkflowArguments=workFlowArgs, phaseName="trimBlast"))
            else:
                self.addChild(CactusSetupPhase(cactusWorkflowArguments=workFlowArgs,
                                                     phaseName="setup"))
        logger.info("Going to create alignments and define the cactus tree")

        self.addFollowOn(FinishUp(workFlowArgs, self.project))
                               
class FinishUp(Job):
    def __init__(self, workFlowArgs, project,):
        Job.__init__(self)
        self.workFlowArgs = workFlowArgs
        self.project = project
    
    def run(self, fileStore):
        donePath = os.path.join(os.path.dirname(self.workFlowArgs.experimentFile), "DONE")
        doneFile = open(donePath, "w")
        doneFile.write("")
        doneFile.close()

class RunCactusPreprocessorThenProgressiveDown(Job):
    def __init__(self, options, inputDir):
        Job.__init__(self)
        self.options = options
        self.inputDir = inputDir
        
    def run(self, fileStore):
        #Load the multi-cactus project
        project = MultiCactusProject()
        project.readXML(self.inputDir)
        #Create jobs to create the output sequences
        configNode = ET.parse(project.getConfigPath()).getroot()
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals() #This is necessary..
        #Create the preprocessor
        self.addChild(CactusPreprocessor(project.getInputSequencePaths(), 
                                               CactusPreprocessor.getOutputSequenceFiles(project.getInputSequencePaths(), project.getOutputSequenceDir()),
                                               configNode))
        #Now build the progressive-down job
        schedule = Schedule()
        schedule.loadProject(project)
        schedule.compute()
        if self.options.event == None:
            self.options.event = project.mcTree.getRootName()
        assert self.options.event in project.expMap
        leafNames = [ project.mcTree.getName(i) for i in project.mcTree.getLeaves() ]
        self.options.globalLeafEventSet = set(leafNames)
        self.addFollowOn(ProgressiveDown(self.options, project, self.options.event, schedule))

def main():
    usage = "usage: prog [options] <multicactus project>"
    description = "Progressive version of cactus_workflow"
    parser = ArgumentParser(usage=usage, description=description)
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)
    
    parser.add_argument("--nonRecursive", dest="nonRecursive", action="store_true",
                      help="Only process given event (not children) [default=False]", 
                      default=False)
    
    parser.add_argument("--event", dest="event", 
                      help="Target event to process [default=root]", default=None)
    
    parser.add_argument("--overwrite", dest="overwrite", action="store_true",
                      help="Recompute and overwrite output files if they exist [default=False]",
                      default=False)
    parser.add_argument("--inputDir", dest="inputDir", help="Directory of multicactus project.")
    
    options = parser.parse_args()
    setLoggingFromOptions(options)

    #if len(args) != 1:
    #    parser.print_help()
    #    raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    Job.Runner.startToil(RunCactusPreprocessorThenProgressiveDown(options, options.inputDir), options)

if __name__ == '__main__':
    from cactus.progressive.cactus_progressive import *
    main()
