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

from cactus.pipeline.cactus_workflow import CactusSetupPhase
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import getOptionalAttrib
from cactus.pipeline.cactus_workflow import findRequiredNode

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.ktserverLauncher import KtserverLauncher
from cactus.progressive.mafFilter import mafFilterOutgroup
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

        # take union of command line options and config options for hal and reference
        if self.options.buildReference == False:
            refNode = findRequiredNode(configXml, "reference")
            self.options.buildReference = getOptionalAttrib(refNode, "buildReference", bool, False)
        if self.options.buildHal == False:
            halNode = findRequiredNode(configXml, "hal")
            self.options.buildHal = getOptionalAttrib(refNode, "buildHal", bool, False)
        halNode = findRequiredNode(configXml, "hal")
        if self.options.buildHal == False:
            self.options.buildHal = getOptionalAttrib(halNode, "buildHal", bool, False)
        self.options.makeMaf = getOptionalAttrib(halNode, "makeMaf", bool, False)
        self.options.joinMaf = getOptionalAttrib(halNode, "joinMaf", bool, False)
        
        # delete database files if --setupAndBuildAlignments
        # and overwrite specified (or if reference not present)
        if self.options.skipAlignments is False and\
         (self.options.overwrite or\
          not os.path.exists(experiment.getReferencePath())):
            dbPath = os.path.join(experiment.getDbDir(), 
                                  experiment.getDbName())
            seqPath = os.path.join(experiment.getDbDir(), "sequences")
            system("rm -f %s* %s %s" % (dbPath, seqPath, 
                                        experiment.getReferencePath()))
            
        ktserver = None
        if experiment.getDbType() == "kyoto_tycoon" and \
        (self.options.skipAlignments is False or \
         self.options.buildReference is True or \
         self.options.buildHal is True):
            ktserver = KtserverLauncher()
            ktserver.spawnServer(experiment)
            # port may be updated.  let's save the experiment back to disk
            experiment.writeXML(self.options.experimentFile)
            
        self.addChildTarget(StartWorkflow(self.options,
                                          self.project,
                                          self.event,
                                          ktserver))
        
        if ktserver is not None:
            self.setFollowOnTarget(EndWorkflow(experiment, self.event, 
                                               ktserver))

class StartWorkflow(Target):
    def __init__(self, options, project, event, ktserver):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.event = event
        self.ktserver = ktserver
    
    def run(self):
        logger.info("StartWorkflow: " + self.event)
        
        # get parameters that cactus_workflow stuff wants
        workFlowArgs = CactusWorkflowArguments(self.options)
        # copy over the options so we don't trail them around
        workFlowArgs.skipAlignments = self.options.skipAlignments
        workFlowArgs.buildReference = self.options.buildReference
        workFlowArgs.buildHal = self.options.buildHal
        workFlowArgs.makeMaf = self.options.makeMaf
        workFlowArgs.joinMaf = self.options.joinMaf
        workFlowArgs.overwrite = self.options.overwrite
        experiment = ExperimentWrapper(workFlowArgs.experimentNode)

        if workFlowArgs.skipAlignments is False and \
        not os.path.exists(experiment.getReferencePath()):       
            self.addChildTarget(CactusSetupPhase(cactusWorkflowArguments=workFlowArgs, phaseName="setup"))
            logger.info("Going to create alignments and define the cactus tree")
        
        self.setFollowOnTarget(ExtractReference(workFlowArgs, self.project, self.event))

class EndWorkflow(Target):
    def __init__(self, workFlowArgs, event, ktserver):
        Target.__init__(self)
        self.workFlowArgs = workFlowArgs
        self.event = event
        self.ktserver = ktserver
    
    def run(self):
        logger.info("EndWorkflow: " + self.event)
        
        # don't need the ktserver anymore, so we kill it
        if self.ktserver is not None:
             experiment = ExperimentWrapper(self.workFlowArgs.experimentNode)
             self.ktserver.killServer(experiment)   
                         
class ExtractReference(Target):
    def __init__(self, workFlowArgs, project, event):
        Target.__init__(self)
        self.workFlowArgs = workFlowArgs
        self.project = project
        self.event = event

    def run(self):                
        if self.workFlowArgs.buildReference:
            logger.info("Starting Reference Extract Phase")
            experiment = ExperimentWrapper(self.workFlowArgs.experimentNode)
            cmdLine = "cactus_getReferenceSeq --cactusDisk \'%s\' --flowerName 0 --referenceEventString %s --outputFile %s --logLevel %s" % \
            (experiment.getDiskDatabaseString(), self.event,
             experiment.getReferencePath(), getLogLevelString())            
            
            system(cmdLine)
          
        self.setFollowOnTarget(BuildMAF(self.workFlowArgs, self.project))
                
class BuildMAF(Target):
    def __init__(self, workFlowArgs, project,):
        Target.__init__(self)
        self.workFlowArgs = workFlowArgs
        self.project = project
    
    def run(self):
        experiment =  ExperimentWrapper(self.workFlowArgs.experimentNode)
        if self.workFlowArgs.buildHal and self.workFlowArgs.makeMaf and\
               os.path.exists(experiment.getMAFPath()) and \
               os.path.splitext(experiment.getMAFPath())[1] == ".maf":
            logger.info("Filtering outgroup from MAF")

            mafFilterOutgroup(experiment) 
        
        self.setFollowOnTarget(JoinMAF(self.workFlowArgs, self.project))

class JoinMAF(Target):
    def __init__(self, workFlowArgs, project):
        Target.__init__(self)
        self.workFlowArgs = workFlowArgs
        self.project = project
    
    # hack to see if a maf file was created by cactus
    # (by searching for the word cactus in the first few comments).
    # used to know if we rerun mafjoin when overwrite is false.
    # if it's cactus we do, otherwise we assume it's been 
    # made by a succesfsul mafjoin in the past and we don't.
    def isMafFromCactus(self, path):
        file = open(path, "r")
        isCactus = False
        for line in file.readlines()[:10]:
            if line.lstrip().find('cactus') >= 0:
                isCactus = True
                break
        file.close()
        return isCactus
            
    def run(self):
        experiment = ExperimentWrapper(self.workFlowArgs.experimentNode)
        if self.workFlowArgs.joinMaf and\
        (self.workFlowArgs.overwrite or not os.path.exists(experiment.getMAFPath())
         or self.isMafFromCactus(experiment.getMAFPath())):
            logger.info("Starting MAF Join")
            
            rootRefName = experiment.getReferenceNameFromConfig()
            rootIsTreeMAF = False
            for child, path in experiment.seqMap.items():
                if child in self.project.expMap and\
                child not in experiment.getOutgroupEvents():
                    childExpPath = self.project.expMap[child]
                    childExpXML = ET.parse(childExpPath).getroot()
                    childExp = ExperimentWrapper(childExpXML)
                    childRefName = childExp.getReferenceNameFromConfig()
                    cmdline = "mafJoin -maxBlkWidth=128 -maxInputBlkWidth=1000 "
                    if not rootIsTreeMAF:
                        cmdline += "-treelessRoot1=%s " % rootRefName
                    cmdline += "-treelessRoot2=%s " % childRefName
                    cmdline = "%s %s %s %s %s_tmp.maf" % (cmdline, childRefName, 
                                                  experiment.getMAFPath(),
                                                  childExp.getMAFPath(),
                                                  experiment.getMAFPath())
                    system(cmdline)
                    
                    move("%s_tmp.maf" % experiment.getMAFPath(), 
                         experiment.getMAFPath())
                    rootIsTreeMAF = True

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
<<<<<<< HEAD
=======
    
    # cannot align without building reference and vice versa in progressive mode
    options.buildReference = options.setupAndBuildAlignments
    # forget about this stuff for now
    options.buildTrees = False
>>>>>>> 50d70bd4b64548ff6a6bd3240fd7a0d680cfe136

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
