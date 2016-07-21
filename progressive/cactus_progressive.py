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
from sonLib.bioio import catFiles

from toil.lib.bioio import getLogLevelString
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions

from cactus.shared.common import cactusRootPath
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL
  
from toil.job import Job
from toil.common import Toil

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
from cactus2hal.cactus2hal import exportHal

class ProgressiveDown(Job):
    def __init__(self, options, project, event, schedule, memory=None, cores=None, disk=None):
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
    
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()
        logger.info("Progressive Down: " + self.event)

        depProjects = dict()
        if not self.options.nonRecursive:
            deps = self.schedule.deps(self.event)
            for child in deps:
                depProjects[child] = self.addChild(ProgressiveDown(self.options,
                                                    self.project, child, 
                                                                   self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()
        
        return self.addFollowOn(ProgressiveNext(self.options, self.project, self.event,
                                                              self.schedule, depProjects, memory=self.configWrapper.getDefaultMemory())).rv()
class ProgressiveNext(Job):
    def __init__(self, options, project, event, schedule, depProjects, memory=None, cores=None, disk=None):
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
        self.depProjects = depProjects
    
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        for projName in self.depProjects:
            depProject = self.depProjects[projName]
            for expName in depProject.expIDMap: 
                expID = depProject.expIDMap[expName]
                logger.info("ExpID %s" % expID)
                experiment = ExperimentWrapper(ET.parse(fileStore.readGlobalFile(expID)).getroot())
                logger.info("Reference id: %s, reference path: %s" % (experiment.getReferenceID(), experiment.getReferencePath()))
                logger.info("Getting the reference sequence for %s" % expName)
                if experiment.getReferenceID():
                    self.project.expIDMap[expName] = expID
                    self.project.outputSequenceIDMap[expName] = experiment.getReferenceID()
                        
        eventExpWrapper = None
        logger.info("Progressive Next: " + self.event)
        if not self.schedule.isVirtual(self.event):
            eventExpWrapper = self.addChild(ProgressiveUp(self.options, self.project, self.event, memory=self.configWrapper.getDefaultMemory())).rv()
        return self.addFollowOn(ProgressiveOut(self.options, self.project, self.event, eventExpWrapper, self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

class ProgressiveOut(Job):
    def __init__(self, options, project, event, eventExpWrapper, schedule, memory=None, cores=None, disk=None):
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.options = options
        self.project = project
        self.event = event
        self.eventExpWrapper = eventExpWrapper
        self.schedule = schedule
        
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        tmpExp = fileStore.getLocalTempFile()
        self.eventExpWrapper.writeXML(tmpExp)
        self.project.expIDMap[self.event] = fileStore.writeGlobalFile(tmpExp)
        followOnEvent = self.schedule.followOn(self.event)
        if followOnEvent is not None:
            logger.info("Adding follow-on event %s" % followOnEvent)
            return self.addFollowOn(ProgressiveDown(self.options, self.project, followOnEvent,
                                                    self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

        return self.project
    
class ProgressiveUp(Job):
    def __init__(self, options, project, event, memory=None, cores=None, disk=None):
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.options = options
        self.project = project
        self.event = event
    
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        logger.info("Progressive Up: " + self.event)

        # open up the experiment
        # note that we copy the path into the options here
        experimentFile = fileStore.readGlobalFile(self.project.expIDMap[self.event])
        expXml = ET.parse(experimentFile).getroot()
        experiment = ExperimentWrapper(expXml)
        configPath = fileStore.readGlobalFile(experiment.getConfigID())
        configXml = ET.parse(configPath).getroot()
        configWrapper = ConfigWrapper(configXml)

        seqMap = experiment.buildSequenceMap()
        seqIDMap = dict()
        tree = experiment.getTree()
        seqNames = []
        for node in tree.postOrderTraversal():
            if tree.isLeaf(node):
                name = tree.getName(node)
                seqIDMap[name] = self.project.outputSequenceIDMap[name]
                seqNames.append(name)
        logger.info("Sequences in progressive, %s: %s" % (self.event, seqNames))
            
        experimentFile = fileStore.getLocalTempFile()
        experiment.writeXML(experimentFile)
        self.options.experimentFileID = fileStore.writeGlobalFile(experimentFile)

        # need at least 3 processes for every event when using ktserver:
        # 1 proc to run jobs, 1 proc to run server, 1 proc to run 2ndary server
        if experiment.getDbType() == "kyoto_tycoon":            
            maxParallel = min(len(self.project.expMap),
                             configWrapper.getMaxParallelSubtrees()) 
            if self.options.batchSystem == "singleMachine":
                pass
                #if int(self.options.maxThreads) < maxParallel * 3:
                    #raise RuntimeError("At least %d threads are required (only %d were specified) to handle up to %d events using kyoto tycoon. Either increase the number of threads using the --maxThreads option or decrease the number of parallel jobs (currently %d) by adjusting max_parallel_subtrees in the config file" % (maxParallel * 3, self.options.maxThreads, maxParallel, configWrapper.getMaxParallelSubtrees()))
            else:
                pass
                #if int(self.options.maxCores) < maxParallel * 3:
                    #raise RuntimeError("At least %d concurrent cpus are required to handle up to %d events using kyoto tycoon. Either increase the number of cpus using the --maxCpus option or decrease the number of parallel jobs (currently %d) by adjusting max_parallel_subtrees in the config file" % (maxParallel * 3, maxParallel, configWrapper.getMaxParallelSubtrees()))
                    
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
        configFile = fileStore.readGlobalFile(experiment.getConfigID())
        configNode = ET.parse(configFile).getroot()
        workFlowArgs = CactusWorkflowArguments(self.options, experimentFile=experimentFile, configNode=configNode, seqIDMap = seqIDMap)
        # copy over the options so we don't trail them around
        workFlowArgs.buildReference = self.options.buildReference
        workFlowArgs.buildHal = self.options.buildHal
        workFlowArgs.buildFasta = self.options.buildFasta
        workFlowArgs.overwrite = self.options.overwrite
        workFlowArgs.globalLeafEventSet = self.options.globalLeafEventSet
        

        donePath = os.path.join(os.path.dirname(workFlowArgs.experimentFile), "DONE")
        doneDone = os.path.isfile(donePath)
        refDone = not workFlowArgs.buildReference or os.path.isfile(experiment.getReferencePath())
        halDone = not workFlowArgs.buildHal or (os.path.isfile(experiment.getHALFastaPath()) and
                                                os.path.isfile(experiment.getHALPath()))
        finalExpID = None
                                                               
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
                finalExpWrapper = self.addChild(CactusTrimmingBlastPhase(cactusWorkflowArguments=workFlowArgs, phaseName="trimBlast")).rv()
            else:
                finalExpWrapper = self.addChild(CactusSetupPhase(cactusWorkflowArguments=workFlowArgs,
                                                     phaseName="setup")).rv()
        logger.info("Going to create alignments and define the cactus tree")

        return finalExpWrapper

class RunCactusPreprocessorThenProgressiveDown(Job):
    def __init__(self, options, project, memory=None, cores=None, disk=None):
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.options = options
        self.project = project
        
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        #Create jobs to create the output sequences
        logger.info("Reading config file from: %s" % self.project.getConfigID())
        configFile = fileStore.readGlobalFile(self.project.getConfigID())
        configNode = ET.parse(configFile).getroot()
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals() #This is necessary..
        #Add the preprocessor child job. The output is a job promise value that will be
        #converted into a list of the IDs of the preprocessed sequences in the follow on job.
        preprocessorJob = self.addChild(CactusPreprocessor(self.project.getInputSequenceIDs(), configNode))
        self.project.setOutputSequenceIDs([preprocessorJob.rv(i) for i in range(len(self.project.getInputSequenceIDs()))])

        #Now build the progressive-down job
        schedule = Schedule()
        schedule.loadProject(self.project, fileStore = fileStore)
        schedule.compute()
        if self.options.event == None:
            self.options.event = self.project.mcTree.getRootName()
        assert self.options.event in self.project.expMap
        leafNames = [ self.project.mcTree.getName(i) for i in self.project.mcTree.getLeaves() ]
        self.options.globalLeafEventSet = set(leafNames)

        return self.addFollowOn(RunCactusPreprocessorThenProgressiveDown2(self.options, self.project, self.options.event, schedule, memory=self.configWrapper.getDefaultMemory())).rv()

class RunCactusPreprocessorThenProgressiveDown2(Job):
    def __init__(self, options, project, event, schedule, memory=None, cores=None, disk=None):
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        project = self.addChild(ProgressiveDown(self.options, self.project, self.event, self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

        #Combine the smaller HAL files from each experiment
        return self.addFollowOnJobFn(exportHal, project=project, memory=self.configWrapper.getDefaultMemory(),
                disk=self.configWrapper.getExportHalDisk()).rv()

        
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
    parser.add_argument("--project", dest="project", help="Directory of multicactus project.")
    parser.add_argument('--outputHal', dest="outputHAL", help="Output HAL file from the alignment, given \
            as either a path starting with file:// or an S3 bucket.")
    parser.add_argument('--resume', dest='resume', help="Resume the workflow", action="store_true")
    
    options = parser.parse_args()
    setLoggingFromOptions(options)
    project = MultiCactusProject()

    with Toil(options) as toil:
        project.readXML(options.project)
        #import the sequences
        seqIDs = []
        for seq in project.getInputSequencePaths():
            if os.path.isdir(seq):
                tmpSeq = getTempFile()
                catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                seq = tmpSeq
            seq = makeURL(seq)
            seqIDs.append(toil.importFile(seq))
        project.setInputSequenceIDs(seqIDs)

        
        #import cactus config
        cactusConfigID = toil.importFile(makeURL(project.getConfigPath()))
        logger.info("Setting config id to: %s" % cactusConfigID)
        project.setConfigID(cactusConfigID)

        project.syncToFileStore(toil)
        configNode = ET.parse(project.getConfigPath()).getroot()
        configWrapper = ConfigWrapper(configNode)
        configWrapper.substituteAllPredefinedConstantsWithLiterals()


        project.writeXML(options.project)

        #Run the workflow
        if options.resume:
            halID = toil.restart()

        else:
            halID = toil.start(RunCactusPreprocessorThenProgressiveDown(options, project, memory=configWrapper.getDefaultMemory()))

        if options.outputHAL:
            toil.exportFile(halID, makeURL(options.outputHAL))


if __name__ == '__main__':
    from cactus.progressive.cactus_progressive import *
    main()
