#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Wrapper to run the cactus_workflow progressively, using the input species tree as a guide

tree.   The --buildMAF and --joinMAF options can be used to create and merge MAFs from the

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
from cactus.pipeline.cactus_workflow import CactusPhylogenyPhase
from cactus.pipeline.cactus_workflow import CactusReferencePhase
from cactus.pipeline.cactus_workflow import CactusFacesPhase
from cactus.pipeline.cactus_workflow import expandWorkflowOptions

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.ktserverLauncher import KtserverLauncher
from cactus.progressive.mafFilter import removeOutgroupFromMaf
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
        
        expXml = ET.parse(self.project.expMap[self.event]).getroot()
        experiment = ExperimentWrapper(expXml)
        
        if not self.options.nonRecursive:
            deps = self.schedule.deps(self.event)
            for child in deps:
                if child in self.project.expMap:
                    self.addChildTarget(ProgressiveDown(self.options,
                                                        self.project, child, 
                                                        self.schedule))
        
        self.setFollowOnTarget(ProgressiveUp(self.options, self.project, 
                               experiment, self.event))
        
class ProgressiveUp(Target):
    def __init__(self, options, project, experiment, event):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.experiment = experiment
        self.event = event
    
    def run(self):
        logger.info("Progressive Up: " + self.event)
        
        # delete database files if --setupAndBuildAlignments
        if self.options.setupAndBuildAlignments and self.options.overwrite:
            dbPath = os.path.join(self.experiment.getDbDir(), 
                                  self.experiment.getDbName())
            seqPath = os.path.join(self.experiment.getDbDir(), "sequences")
            system("rm -f %s* %s %s" % (dbPath, seqPath, 
                                        self.experiment.getReferencePath()))
            
        ktserver = None
        if self.experiment.getDbType() == "kyoto_tycoon" and \
        (self.options.setupAndBuildAlignments is True or \
         self.options.buildReference is True or \
         self.options.buildMAF is True):
            ktserver = KtserverLauncher()
            ktserver.spawnServer(self.experiment)
            
        self.addChildTarget(StartWorkflow(self.options,
                                          self.project,
                                          self.experiment,
                                          self.event,
                                          ktserver))
        
        if ktserver is not None:
            self.setFollowOnTarget(EndWorkflow(self.experiment, self.event, 
                                               ktserver))

class StartWorkflow(Target):
    def __init__(self, options, project, experiment, event, ktserver):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.experiment = experiment
        self.event = event
        self.ktserver = ktserver
    
    def run(self):
        logger.info("StartWorkflow: " + self.event)
        
        # get parameters that cactus_workflow stuff wants 
        wfOpts, sequences = getWorkflowParams(self.options, self.experiment)
         
        if self.options.setupAndBuildAlignments and \
        not os.path.exists(self.experiment.getReferencePath()):       
            self.addChildTarget(CactusSetupPhase(wfOpts, sequences))
            logger.info("Going to create alignments and define the cactus tree")
        
        self.setFollowOnTarget(ExtractReference(self.options,
                                                self.project,
                                                self.experiment,
                                                self.event,
                                                self.ktserver))

class EndWorkflow(Target):
    def __init__(self, experiment, event, ktserver):
        Target.__init__(self)
        self.experiment = experiment
        self.event = event
        self.ktserver = ktserver
    
    def run(self):
        logger.info("EndWorkflow: " + self.event)
        
        # don't need the ktserver anymore, so we kill it
        if self.ktserver is not None:
            self.ktserver.killServer(self.experiment)   
                         
class ExtractReference(Target):
    def __init__(self, options, project, experiment, event, ktserver):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.experiment = experiment
        self.event = event
        self.ktserver = ktserver

    def run(self):                
        if self.options.buildReference:
            logger.info("Starting Reference Extract Phase")
            
            cmdLine = "cactus_getReferenceSeq --cactusDisk \'%s\' --flowerName 0 --referenceEventString %s --outputFile %s --logLevel %s" % \
            (self.experiment.getDiskDatabaseString(), self.event,
             self.experiment.getReferencePath(), getLogLevelString())            
            
            system(cmdLine)
            
        self.setFollowOnTarget(BuildMAF(self.options, self.project,
                                        self.experiment,
                                        self.event, self.ktserver))
                
class BuildMAF(Target):
    def __init__(self, options, project, experiment, event, ktserver):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.experiment = experiment
        self.event = event
        self.ktserver = ktserver
    
    def run(self):
        if self.options.buildMAF and\
        (self.options.overwrite or not os.path.exists(self.experiment.getMAFPath())):             
            logger.info("Starting MAF Build phase")
            
            cmdLine = "cactus_MAFGenerator --cactusDisk \'%s\' --flowerName 0 --outputFile %s --logLevel %s" % \
            (self.experiment.getDiskDatabaseString(),
             self.experiment.getMAFPath(), getLogLevelString())            

            system(cmdLine) 
            removeOutgroupFromMaf(self.experiment.getMAFPath(), 
                                  self.experiment.getOutgroupName()) 
        
        self.setFollowOnTarget(JoinMAF(self.options, self.project,
                                       self.experiment,
                                       self.event, self.ktserver))

class JoinMAF(Target):
    def __init__(self, options, project, experiment, event, ktserver):
        Target.__init__(self)
        self.options = options
        self.project = project
        self.experiment = experiment
        self.event = event
        self.ktserver = ktserver
    
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
        if self.options.joinMAF and\
        (self.options.overwrite or not os.path.exists(self.experiment.getMAFPath())
         or self.isMafFromCactus(self.experiment.getMAFPath())):
            logger.info("Starting MAF Join phase")
            
            rootRefName = self.experiment.getReferenceNameFromConfig()
            rootIsTreeMAF = False
            for child, path in self.experiment.seqMap.items():
                if child in self.project.expMap and\
                child != self.experiment.getOutgroupName():
                    childExpPath = self.project.expMap[child]
                    childExpXML = ET.parse(childExpPath).getroot()
                    childExp = ExperimentWrapper(childExpXML)
                    childRefName = childExp.getReferenceNameFromConfig()
                    cmdline = "mafJoin -maxBlkWidth=128 -maxInputBlkWidth=1000 "
                    if not rootIsTreeMAF:
                        cmdline += "-treelessRoot1=%s " % rootRefName
                    cmdline += "-treelessRoot2=%s " % childRefName
                    cmdline = "%s %s %s %s %s_tmp.maf" % (cmdline, childRefName, 
                                                  self.experiment.getMAFPath(),
                                                  childExp.getMAFPath(),
                                                  self.experiment.getMAFPath())
                    system(cmdline)
                    
                    move("%s_tmp.maf" % self.experiment.getMAFPath(), 
                         self.experiment.getMAFPath())
                    rootIsTreeMAF = True
    
                    
def getWorkflowParams(baseOptions, experiment):
    options = copy.deepcopy(baseOptions)
    expandWorkflowOptions(options, experiment.xmlRoot)
    #Get the sequences
    sequences = options.experimentFile.attrib["sequences"].split()
    return options, sequences
                           
def main():    
    usage = "usage: %prog [options] <multicactus project>"
    description = "Progressive version of cactus_workflow"
    parser = OptionParser(usage=usage, description=description)
    Stack.addJobTreeOptions(parser)
    
    parser.add_option("--setupAndBuildAlignments", dest="setupAndBuildAlignments", action="store_true",
                      help="Setup and build alignments then normalise the resulting structure", default=False)
    
    parser.add_option("--buildMAF", dest="buildMAF", action="store_true",
                     help="Create a MAF file from the cactus [default=false]", default=False)
    
    parser.add_option("--joinMAF", dest="joinMAF", action="store_true",
                     help="Progressively join all cactus MAFs[default=false]", default=False)
    
    parser.add_option("--skipCheck", dest="skipCheck", action="store_true",
                      help="Don't run cactus_check", default=False)
    
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
    
    # cannot align without building reference and vice versa in progressive mode
    options.buildReference = options.setupAndBuildAlignments
    # forget about this stuff for now
    options.buildTrees = False
    options.buildFaces = False

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
    
    baseTarget = ProgressiveDown(options, project, options.event, schedule)
    Stack(baseTarget).startJobTree(options)
    
   
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.progressive.cactus_progressive import *
    _test()
    main()