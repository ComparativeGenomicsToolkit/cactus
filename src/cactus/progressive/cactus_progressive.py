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

from toil.lib.bioio import getTempFile
from toil.lib.bioio import system


from toil.lib.bioio import getLogLevelString
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions

from cactus.shared.common import cactusRootPath
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL
from cactus.shared.common import catFiles
  
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
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.progressive.seqFile import SeqFile

class ProgressiveDown(Job):
    def __init__(self, options, project, event, schedule, memory=None, cores=None):
        Job.__init__(self, memory=memory, cores=cores)
        self.options = options
        self.project = project
        self.event = event
        self.schedule = schedule
    
    def run(self, fileStore):
        self.configNode = ET.parse(fileStore.readGlobalFile(self.project.getConfigID())).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()
        logger.info("Progressive Down: " + self.event)
        assert self.event == "Anc2"

        depProjects = dict()
        if not self.options.nonRecursive:
            deps = self.schedule.deps(self.event)
            fileStore.logToMaster("There are %i dependent projects" % len(deps))
            for child in deps:
                fileStore.logToMaster("Adding dependent project %s" % child)
                depProjects[child] = self.addChild(ProgressiveDown(self.options,
                                                    self.project, child, 
                                                                   self.schedule, memory=self.configWrapper.getDefaultMemory())).rv()
        
        return self.addFollowOn(ProgressiveNext(self.options, self.project, self.event,
                                                              self.schedule, depProjects, memory=self.configWrapper.getDefaultMemory())).rv()
class ProgressiveNext(Job):
    def __init__(self, options, project, event, schedule, depProjects, memory=None, cores=None):
        Job.__init__(self, memory=memory, cores=cores)
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

class ProgressiveOut(Job):
    def __init__(self, options, project, event, eventExpWrapper, schedule, memory=None, cores=None):
        Job.__init__(self, memory=memory, cores=cores)
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
    def __init__(self, options, project, event, memory=None, cores=None):
        Job.__init__(self, memory=memory, cores=cores)
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
    def __init__(self, options, project, memory=None, cores=None):
        Job.__init__(self, memory=memory, cores=cores)
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

def exportHal(job, project, event=None, cacheBytes=None, cacheMDC=None, cacheRDC=None, cacheW0=None, chunk=None, deflate=None, inMemory=False):

    # some quick stats
    totalTime = time.time()
    totalAppendTime = 0

    HALPath = job.fileStore.getLocalTempFile()

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
            expTreeString = NXNewick().writeString(experiment.getTree())
            assert len(expTreeString) > 1
            assert experiment.getHalID() is not None
            assert experiment.getHalFastaID() is not None
            subHALPath = job.fileStore.readGlobalFile(experiment.getHalID())
            halFastaPath = job.fileStore.readGlobalFile(experiment.getHalFastaID())


            opts = "\'{0}\' \'{1}\' \'{2}\' \'{3}\'".format(os.path.basename(subHALPath), os.path.basename(halFastaPath), expTreeStrinqg, HALPath)
            
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

            cactus_call(tool="halAppendCactusSubtree",
                        option_string=opts)

            
            #appendTime = time.time() - appendTime
            #totalAppendTime += appendTime
#            print "time of above command: {0:.2f}".format(appendTime)
 
    #totalTime = time.time() - totalTime
    #print "total time: {0:.2f}  total halAppendCactusSubtree time: {1:.2f}".format(totalTime, totalAppendTime)
    return job.fileStore.writeGlobalFile(HALPath)

class RunCactusPreprocessorThenProgressiveDown2(Job):
    def __init__(self, options, project, event, schedule, memory=None, cores=None):
        Job.__init__(self, memory=memory, cores=cores)
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
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("workDir", help = "Work dir")
    parser.add_argument("outputHal", help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--jobStore", dest="jobStore",
                      help="JobStore to use for Toil. If not given,\
                      the FileJobStore will be used.", default=None)
    parser.add_argument("--optionsFile", dest="optionsFile",
                      help="Text file containing command line options to use as"\
                      " defaults", default=None)
    parser.add_argument("--database", dest="database",
                      help="Database type: tokyo_cabinet or kyoto_tycoon"
                      " [default: %default]",
                      default="kyoto_tycoon")
    parser.add_argument("--outputMaf", dest="outputMaf",
                      help="[DEPRECATED use hal2maf on the ouput file instead] Path of output alignment in .maf format.  This option should be avoided and will soon be removed.  It may cause sequence names to be mangled, and use a tremendous amount of memory. ",
                      default=None)
    parser.add_argument("--configFile", dest="configFile",
                      help="Specify cactus configuration file",
                      default=None)
    parser.add_argument("--legacy", dest="legacy", action="store_true", help=
                      "Run cactus directly on all input sequences "
                      "without any progressive decomposition (ie how it "
                      "was originally published in 2011)",
                      default=False)
    parser.add_argument("--autoAbortOnDeadlock", dest="autoAbortOnDeadlock",
                      action="store_true",
                      help="Abort automatically when jobTree monitor" +
                      " suspects a deadlock by deleting the jobTree folder." +
                      " Will guarantee no trailing ktservers but still " +
                      " dangerous to use until we can more robustly detect " +
                      " deadlocks.",
                      default=False)
    parser.add_argument("--overwrite", dest="overwrite", action="store_true",
                      help="Re-align nodes in the tree that have already" +
                      " been successfully aligned.",
                      default=False)
    parser.add_argument("--rootOutgroupDists", dest="rootOutgroupDists",
                      help="root outgroup distance (--rootOutgroupPaths must " +
                      "be given as well)", default=None)
    parser.add_argument("--rootOutgroupPaths", dest="rootOutgroupPaths", type=str,
                      help="root outgroup path (--rootOutgroup must be given " +
                      "as well)", default=None)
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                      " must appear in NEWICK tree in <seqfile>) to use as a "
                      "root for the alignment.  Any genomes not below this node "
                      "in the tree may be used as outgroups but will never appear"
                      " in the output.  If no root is specifed then the root"
                      " of the tree is used. ", default=None)
    parser.add_argument("--resume", dest="resume", 
                      help="Resume the workflow", action="store_true")

    #Kyoto Tycoon Options
    ktGroup = parser.add_argument_group("kyoto_tycoon Options",
                          "Kyoto tycoon provides a client/server framework "
                          "for large in-memory hash tables and is available "
                          "via the --database option.")
    ktGroup.add_argument("--ktPort", dest="ktPort",
                       help="starting port (lower bound of range) of ktservers"
                       " [default: %default]",
                       default=1978)
    ktGroup.add_argument("--ktHost", dest="ktHost",
                       help="The hostname to use for connections to the "
                       "ktserver (this just specifies where nodes will attempt"
                       " to find the server, *not* where the ktserver will be"
                       " run)",
                       default=None)
    ktGroup.add_argument("--ktType", dest="ktType",
                       help="Kyoto Tycoon server type "
                       "(memory, snapshot, or disk)"
                       " [default: %default]",
                       default='memory')
    # sonlib doesn't allow for spaces in attributes in the db conf
    # which renders this options useless
    #ktGroup.add_argument("--ktOpts", dest="ktOpts",
    #                   help="Command line ktserver options",
    #                   default=None)
    ktGroup.add_argument("--ktCreateTuning", dest="ktCreateTuning",
                       help="ktserver options when creating db "\
                            "(ex #bnum=30m#msiz=50g)",
                       default=None)
    ktGroup.add_argument("--ktOpenTuning", dest="ktOpenTuning",
                       help="ktserver options when opening existing db "\
                            "(ex #opts=ls#ktopts=p)",
                       default=None)
    parser.add_argument_group(ktGroup)
   
    parser.add_argument("--nonRecursive", dest="nonRecursive", action="store_true",
                      help="Only process given event (not children) [default=False]", 
                      default=False)
    
    parser.add_argument("--event", dest="event", 
                      help="Target event to process [default=root]", default=None)
    


    options = parser.parse_args()
    setLoggingFromOptions(options)

    #Create the progressive cactus project 
    projWrapper = ProjectWrapper(options)
    projWrapper.writeXml()

    pjPath = os.path.join(options.workDir, ProjectWrapper.alignmentDirName,
                          '%s_project.xml' % ProjectWrapper.alignmentDirName)

    project = MultiCactusProject()

    if not os.path.isdir(options.workDir):
        os.makedirs(options.workDir)

    with Toil(options) as toil:
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
        cactusConfigID = toil.importFile(makeURL(project.getConfigPath()))
        logger.info("Setting config id to: %s" % cactusConfigID)
        project.setConfigID(cactusConfigID)

        project.syncToFileStore(toil)
        configNode = ET.parse(project.getConfigPath()).getroot()
        configWrapper = ConfigWrapper(configNode)
        configWrapper.substituteAllPredefinedConstantsWithLiterals()


        project.writeXML(pjPath)

        #Run the workflow
        if options.resume:
            halID = toil.restart()

        else:
            halID = toil.start(RunCactusPreprocessorThenProgressiveDown(options, project, memory=configWrapper.getDefaultMemory()))

        if options.outputHAL:
            toil.exportFile(halID, makeURL(options.outputHAL))


if __name__ == '__main__':
    main()
