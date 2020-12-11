#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Wrapper to run the cactus_workflow progressively, using the input species tree as a guide

tree.
"""

import logging
import os
import xml.etree.ElementTree as ET
import timeit
import multiprocessing
from argparse import ArgumentParser
from base64 import b64encode
from subprocess import check_call
from subprocess import CalledProcessError

from toil.lib.bioio import getTempFile

from toil.lib.bioio import setLoggingFromOptions
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count

from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL
from cactus.shared.common import catFiles
from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import getDockerImage
from cactus.shared.version import cactus_commit
from cactus.shared.common import cactusRootPath
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options

from toil.job import Job
from toil.common import Toil

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
from cactus.preprocessor.dnabrnnMasking import loadDnaBrnnModel
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.shared.common import setupBinaries, importSingularityImage

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

logger = logging.getLogger(__name__)

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
            name = tree.getName(node)
            if tree.isLeaf(node) or (name == experiment.getRootGenome() and experiment.isRootReconstructed() == False):
                seqIDMap[name] = self.project.outputSequenceIDMap[name]
                seqNames.append(name)
        logger.info("Sequences in progressive, %s: %s" % (self.event, seqNames))

        experimentFile = fileStore.getLocalTempFile()
        experiment.writeXML(experimentFile)
        self.options.experimentFileID = fileStore.writeGlobalFile(experimentFile)

        # take union of command line options and config options for hal and reference
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

def logAssemblyStats(job, message, name, sequenceID, preemptable=True):
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

        fileStore.logToMaster("Using the following configuration:\n%s" % ET.tostring(self.configNode, encoding='unicode'))

        # Log the stats for the un-preprocessed assemblies
        for name, sequence in list(self.project.inputSequenceIDMap.items()):
            self.addChildJobFn(logAssemblyStats, "Before preprocessing", name, sequence)

        # Create jobs to create the output sequences
        logger.info("Reading config file from: %s" % self.project.getConfigID())
        configFile = fileStore.readGlobalFile(self.project.getConfigID())
        configNode = ET.parse(configFile).getroot()
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals() #This is necessary..
        #Add the preprocessor child job. The output is a job promise value that will be
        #converted into a list of the IDs of the preprocessed sequences in the follow on job.
        preprocessorJob = self.addChild(CactusPreprocessor(list(self.project.inputSequenceIDMap.values()), configNode))
        rvs = [preprocessorJob.rv(i) for i in range(len(self.project.inputSequenceIDMap))]
        fileStore.logToMaster('input sequence IDs: %s' % self.project.inputSequenceIDMap)
        for genome, rv in zip(list(self.project.inputSequenceIDMap.keys()), rvs):
            self.project.outputSequenceIDMap[genome] = rv

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
            preprocessedSequences = self.project.outputSequenceIDMap
            for genome, seqID in list(preprocessedSequences.items()):
                fileStore.exportFile(seqID, self.options.intermediateResultsUrl + '-preprocessed-' + genome)

        # Log the stats for the preprocessed assemblies
        for name, sequence in list(self.project.outputSequenceIDMap.items()):
            self.addChildJobFn(logAssemblyStats, "After preprocessing", name, sequence)

        project = self.addChild(ProgressiveDown(options=self.options, project=self.project, event=self.event, schedule=self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()

        #Combine the smaller HAL files from each experiment
        return self.addFollowOnJobFn(exportHal, project=project, memory=self.configWrapper.getDefaultMemory(),
                                     disk=self.configWrapper.getExportHalDisk(),
                                     preemptable=False).rv()

def exportHal(job, project, event=None, cacheBytes=None, cacheMDC=None, cacheRDC=None, cacheW0=None, chunk=None, deflate=None, inMemory=True):

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

            outgroups = experiment.getOutgroupGenomes()
            experiment.setConfigPath(job.fileStore.readGlobalFile(experiment.getConfigID()))
            expTreeString = NXNewick().writeString(experiment.getTree(onlyThisSubtree=True))
            assert len(expTreeString) > 1
            assert experiment.getHalID() is not None
            assert experiment.getHalFastaID() is not None
            subHALPath = job.fileStore.readGlobalFile(experiment.getHalID())
            halFastaPath = job.fileStore.readGlobalFile(experiment.getHalFastaID())

            args = [os.path.basename(subHALPath), os.path.basename(halFastaPath), expTreeString, os.path.basename(HALPath)]

            if len(outgroups) > 0:
                args += ["--outgroups", ",".join(outgroups)]
            if cacheBytes is not None:
                args += ["--cacheBytes", cacheBytes]
            if cacheMDC is not None:
                args += ["--cacheMDC", cacheMDC]
            if cacheRDC is not None:
                args += ["--cacheRDC", cacheRDC]
            if cacheW0 is not None:
                args += ["--cacheW0", cacheW0]
            if chunk is not None:
                args += ["--chunk", chunk]
            if deflate is not None:
                args += ["--deflate", deflate]
            if inMemory is True:
                args += ["--inMemory"]

            cactus_call(parameters=["halAppendCactusSubtree"] + args)

    cactus_call(parameters=["halSetMetadata", HALPath, "CACTUS_COMMIT", cactus_commit])
    with job.fileStore.readGlobalFileStream(project.configID) as configFile:
        cactus_call(parameters=["halSetMetadata", HALPath, "CACTUS_CONFIG", b64encode(configFile.read()).decode()])

    return job.fileStore.writeGlobalFile(HALPath)
        
def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                      " must appear in NEWICK tree in <seqfile>) to use as a "
                      "root for the alignment.  Any genomes not below this node "
                      "in the tree may be used as outgroups but will never appear"
                      " in the output.  If no root is specifed then the root"
                      " of the tree is used. ", default=None)
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)
    parser.add_argument("--database", choices=["kyoto_tycoon", "redis"],
                        help="The type of database", default="kyoto_tycoon")

    options = parser.parse_args()

    setupBinaries(options)
    setLoggingFromOptions(options)
    enableDumpStack()

    # cactus doesn't run with 1 core
    if options.batchSystem == 'singleMachine':
        if options.maxCores is not None:
            if int(options.maxCores) < 2:
                raise RuntimeError('Cactus requires --maxCores > 1')
        else:
            # is there a way to get this out of Toil?  That would be more consistent
            if cpu_count() < 2:
                raise RuntimeError('Only 1 CPU detected.  Cactus requires at least 2')

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()
    runCactusProgressive(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("Cactus has finished after {} seconds".format(run_time))

def runCactusProgressive(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            halID = toil.restart()
        else:

            options.cactusDir = getTempDirectory()
            #Create the progressive cactus project
            projWrapper = ProjectWrapper(options, options.configFile)
            projWrapper.writeXml()

            pjPath = os.path.join(options.cactusDir, ProjectWrapper.alignmentDirName,
                                  '%s_project.xml' % ProjectWrapper.alignmentDirName)
            assert os.path.exists(pjPath)

            project = MultiCactusProject()

            if not os.path.isdir(options.cactusDir):
                os.makedirs(options.cactusDir)

            project.readXML(pjPath)
            #import the sequences
            for genome, seq in list(project.inputSequenceMap.items()):
                if os.path.isdir(seq):
                    tmpSeq = getTempFile()
                    catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                    seq = tmpSeq
                seq = makeURL(seq)
                project.inputSequenceIDMap[genome] = toil.importFile(seq)
                
            #import cactus config
            cactusConfigID = toil.importFile(makeURL(options.configFile))
            project.setConfigID(cactusConfigID)

            project.syncToFileStore(toil)
            configNode = ET.parse(project.getConfigPath()).getroot()
            configWrapper = ConfigWrapper(configNode)
            configWrapper.substituteAllPredefinedConstantsWithLiterals()

            # Make sure we have the dna-brnn model in the filestore if we need it
            loadDnaBrnnModel(toil, configNode)

            project.writeXML(pjPath)
            halID = toil.start(RunCactusPreprocessorThenProgressiveDown(options, project, memory=configWrapper.getDefaultMemory()))

        toil.exportFile(halID, makeURL(options.outputHal))

if __name__ == '__main__':
    main()
