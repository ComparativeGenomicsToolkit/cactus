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
from argparse import ArgumentParser

from toil.lib.bioio import getTempFile

from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions

from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL
from cactus.shared.common import catFiles
from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob

from toil.job import Job
from toil.common import Toil

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper

from sonLib.nxnewick import NXNewick

class ProgressiveDown(RoundedJob):
    def __init__(self, options, project, event, schedule, memory=None, cores=None):
        RoundedJob.__init__(self, memory=memory, cores=cores, preemptable=True)
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
        deps = self.schedule.deps(self.event)
        fileStore.logToMaster("There are %i dependent projects" % len(deps))
        for child in deps:
            fileStore.logToMaster("Adding dependent project %s" % child)
            depProjects[child] = self.addChild(ProgressiveDown(self.options,
                                                               self.project, child,
                                                               self.schedule)).rv()

        return self.addFollowOn(ProgressiveNext(self.options, self.project, self.event,
                                                              self.schedule, depProjects, memory=self.configWrapper.getDefaultMemory())).rv()
class ProgressiveNext(RoundedJob):
    def __init__(self, options, project, event, schedule, depProjects, memory=None, cores=None):
        RoundedJob.__init__(self, memory=memory, cores=cores, preemptable=True)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
        self.depProjects = depProjects
    
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        fileStore.logToMaster("Project has %i dependencies" % len(self.depProjects))
        for projName in self.depProjects:
            depProject = self.depProjects[projName]
            for expName in depProject.expIDMap: 
                expID = depProject.expIDMap[expName]
                experiment = ExperimentWrapper(ET.parse(fileStore.readGlobalFile(expID)).getroot())
                fileStore.logToMaster("Reference ID for experiment %s: %s" % (expName, experiment.getReferenceID()))
                if experiment.getReferenceID():
                    self.project.expIDMap[expName] = expID
                    self.project.outputSequenceIDMap[expName] = experiment.getReferenceID()
                        
        eventExpWrapper = None
        logger.info("Progressive Next: " + self.event)
        if not self.schedule.isVirtual(self.event):
            eventExpWrapper = self.addChild(ProgressiveUp(self.options, self.project, self.event, memory=self.configWrapper.getDefaultMemory())).rv()
        return self.addFollowOn(ProgressiveOut(self.options, self.project, self.event, eventExpWrapper, self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

class ProgressiveOut(RoundedJob):
    def __init__(self, options, project, event, eventExpWrapper, schedule, memory=None, cores=None):
        RoundedJob.__init__(self, memory=memory, cores=cores, preemptable=True)
        self.options = options
        self.project = project
        self.event = event
        self.eventExpWrapper = eventExpWrapper
        self.schedule = schedule
        
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        if not self.schedule.isVirtual(self.event):
            tmpExp = fileStore.getLocalTempFile()
            self.eventExpWrapper.writeXML(tmpExp)
            self.project.expIDMap[self.event] = fileStore.writeGlobalFile(tmpExp)
        followOnEvent = self.schedule.followOn(self.event)
        if followOnEvent is not None:
            logger.info("Adding follow-on event %s" % followOnEvent)
            return self.addFollowOn(ProgressiveDown(self.options, self.project, followOnEvent,
                                                    self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

        return self.project
    
class ProgressiveUp(RoundedJob):
    def __init__(self, options, project, event, memory=None, cores=None):
        RoundedJob.__init__(self, memory=memory, cores=cores, preemptable=True)
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
        workFlowArgs.globalLeafEventSet = self.options.globalLeafEventSet
        if self.options.intermediateResultsUrl is not None:
            # Give the URL prefix a special name for this particular
            # subproblem (by suffixing it with the name of the
            # internal node in the guide tree)
            workFlowArgs.intermediateResultsUrl = self.options.intermediateResultsUrl + '-' + self.event

        # Use the trimming strategy to blast ingroups vs outgroups.
        finalExpWrapper = self.addChild(CactusTrimmingBlastPhase(cactusWorkflowArguments=workFlowArgs, phaseName="trimBlast")).rv()
        logger.info("Going to create alignments and define the cactus tree")

        return finalExpWrapper

def logAssemblyStats(job, message, name, sequenceID):
    sequenceFile = job.fileStore.readGlobalFile(sequenceID)
    analysisString = cactus_call(parameters=["cactus_analyseAssembly", sequenceFile], check_output=True)
    job.fileStore.logToMaster("%s, got assembly stats for genome %s: %s" % (message, name, analysisString))

class RunCactusPreprocessorThenProgressiveDown(RoundedJob):
    def __init__(self, options, project, memory=None, cores=None):
        RoundedJob.__init__(self, memory=memory, cores=cores, preemptable=True)
        self.options = options
        self.project = project
        
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        # Log the stats for the un-preprocessed assemblies
        for name, sequence in self.project.getInputSequenceIDMap().items():
            self.addChildJobFn(logAssemblyStats, "Before preprocessing", name, sequence)

        # Create jobs to create the output sequences
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
        schedule.loadProject(self.project, fileStore=fileStore)
        schedule.compute()
        self.options.event = self.project.mcTree.getRootName()
        leafNames = [ self.project.mcTree.getName(i) for i in self.project.mcTree.getLeaves() ]
        fileStore.logToMaster("Leaf names = %s" % leafNames)
        self.options.globalLeafEventSet = set(leafNames)

        return self.addFollowOn(RunCactusPreprocessorThenProgressiveDown2(options=self.options, project=self.project, event=self.options.event, schedule=schedule, memory=self.configWrapper.getDefaultMemory())).rv()


class RunCactusPreprocessorThenProgressiveDown2(RoundedJob):
    def __init__(self, options, project, event, schedule, memory=None, cores=None):
        RoundedJob.__init__(self, memory=memory, cores=cores, preemptable=True)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule

    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()

        # Save preprocessed sequences
        if self.options.intermediateResultsUrl is not None:
            preprocessedSequences = self.project.getOutputSequenceIDMap()
            for genome, seqID in preprocessedSequences.items():
                fileStore.exportFile(seqID, self.options.intermediateResultsUrl + '-preprocessed-' + genome)

        # Log the stats for the preprocessed assemblies
        for name, sequence in self.project.getOutputSequenceIDMap().items():
            self.addChildJobFn(logAssemblyStats, "After preprocessing", name, sequence)

        project = self.addChild(ProgressiveDown(options=self.options, project=self.project, event=self.event, schedule=self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

        #Combine the smaller HAL files from each experiment
        return self.addFollowOnJobFn(exportHal, project=project, memory=self.configWrapper.getDefaultMemory(),
                                     disk=self.configWrapper.getExportHalDisk(),
                                     preemptable=False).rv()

def exportHal(job, project, event=None, cacheBytes=None, cacheMDC=None, cacheRDC=None, cacheW0=None, chunk=None, deflate=None, inMemory=False):

    HALPath = "tmp_alignment.hal"

    # traverse tree to make sure we are going breadth-first
    tree = project.mcTree

    # find subtree if event specified
    rootNode = None
    if event is not None:
        assert event in tree.nameToId and not tree.isLeaf(tree.nameToId[event])
        rootNode = tree.nameToId[event]

    for node in tree.breadthFirstTraversal(rootNode):
        genomeName = tree.getName(node)
        if genomeName in project.expMap:
            experimentFilePath = job.fileStore.readGlobalFile(project.expIDMap[genomeName])
            experiment = ExperimentWrapper(ET.parse(experimentFilePath).getroot())

            outgroups = experiment.getOutgroupEvents()
            experiment.setConfigPath(job.fileStore.readGlobalFile(experiment.getConfigID()))
            expTreeString = NXNewick().writeString(experiment.getTree(onlyThisSubtree=True))
            assert len(expTreeString) > 1
            assert experiment.getHalID() is not None
            assert experiment.getHalFastaID() is not None
            subHALPath = job.fileStore.readGlobalFile(experiment.getHalID())
            halFastaPath = job.fileStore.readGlobalFile(experiment.getHalFastaID())


            opts = "\'{0}\' \'{1}\' \'{2}\' \'{3}\'".format(os.path.basename(subHALPath), os.path.basename(halFastaPath), expTreeString, os.path.basename(HALPath))
            
            if len(outgroups) > 0:
                opts += " --outgroups {0}".format(",".join(outgroups))
            if cacheBytes is not None:
                opts += " --cacheBytes {0}".format(cacheBytes)
            if cacheMDC is not None:
                opts += " --cacheMDC {0}".format(cacheMDC)
            if cacheRDC is not None:
                opts += " --cacheRDC {0}".format(cacheRDC)
            if cacheW0 is not None:
                opts += " --cacheW0 {0}".format(cacheW0)
            if chunk is not None:
                opts += " --chunk {0}".format(chunk)
            if deflate is not None:
                opts += " --deflate {0}".format(deflate)
            if inMemory is True:
                opts += " --inMemory"

            cactus_call(parameters=["halAppendCactusSubtree"],
                        option_string=opts)

    return job.fileStore.writeGlobalFile(HALPath)

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("cactusDir", help = "Work dir")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--database", dest="database",
                      help="Database type: tokyo_cabinet or kyoto_tycoon"
                      " [default: %(default)s]",
                      default="kyoto_tycoon")
    parser.add_argument("--configFile", dest="configFile",
                      help="Specify cactus configuration file",
                      default=None)
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                      " must appear in NEWICK tree in <seqfile>) to use as a "
                      "root for the alignment.  Any genomes not below this node "
                      "in the tree may be used as outgroups but will never appear"
                      " in the output.  If no root is specifed then the root"
                      " of the tree is used. ", default=None)   
    parser.add_argument("--latest", dest="latest", action="store_true",
            help="Use the latest, locally-built docker container rather than \
                    pulling from quay.io")

    options = parser.parse_args()

    if options.latest:
        os.environ["CACTUS_USE_LATEST"] = "1"

    setLoggingFromOptions(options)

    options.cactusDir = os.path.abspath(options.cactusDir)

    # Mess with some toil options to create useful defaults.

    # Caching generally slows down the cactus workflow, plus some
    # methods like readGlobalFileStream don't support forced
    # reads directly from the job store rather than from cache.
    options.disableCaching = True

    if options.retryCount is None:
        # If the user didn't specify a retryCount value, make it 5
        # instead of Toil's default (1).
        options.retryCount = 5

    #Create the progressive cactus project 
    projWrapper = ProjectWrapper(options)
    projWrapper.writeXml()

    pjPath = os.path.join(options.cactusDir, ProjectWrapper.alignmentDirName,
                          '%s_project.xml' % ProjectWrapper.alignmentDirName)
    assert os.path.exists(pjPath)

    project = MultiCactusProject()

    if not os.path.isdir(options.cactusDir):
        os.makedirs(options.cactusDir)

    with Toil(options) as toil:
        #Run the workflow
        if options.restart:
            halID = toil.restart()
        else:
            project.readXML(pjPath)
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
            if options.configFile:
                cactusConfigID = toil.importFile(makeURL(options.configFile))
            else:
                cactusConfigID = toil.importFile(makeURL(project.getConfigPath()))
            logger.info("Setting config id to: %s" % cactusConfigID)
            project.setConfigID(cactusConfigID)

            project.syncToFileStore(toil)
            configNode = ET.parse(project.getConfigPath()).getroot()
            configWrapper = ConfigWrapper(configNode)
            configWrapper.substituteAllPredefinedConstantsWithLiterals()


            project.writeXML(pjPath)
            halID = toil.start(RunCactusPreprocessorThenProgressiveDown(options, project, memory=configWrapper.getDefaultMemory()))

        toil.exportFile(halID, makeURL(options.outputHal))

if __name__ == '__main__':
    main()
