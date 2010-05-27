#!/usr/bin/env python

"""Wrapper to run cactus on different combinations of cactus_workflow_config.xml (cactus
parameters) and simulations, followed by evaluations steps.
"""
#nknguyen@soe.ucsc.edu
#05/26/2010

import os, re, sys, time
from optparse import OptionParser

from workflow.jobTree.jobTree import runJobTree
from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete
from workflow.jobTree.scriptTree.target import Target

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue

from cactus.shared.common import runCactusWorkflow
from cactus.shared.common import runCactusMAFGenerator

class CactusTuningWrapper(Target):
   """Creates children, each of which (is a childTarget and) corresponds to a cactus run 
   """
   def __init__(self, options):
      Target.__init__(self)
      self.options = options
   def run(self, localTempDir, globalTempDir):
      #---------------------------------------------------------------------------------
      #Find combinations of simulation and configuration & issue children to jobTree:#
      #---------------------------------------------------------------------------------
      for sim in self.options.sim:
         sim = modify_dirname(sim)
         #convert mfa file of current simulation into MAF format:
         simRoot = getRootDir(sim)
         simDir = self.options.dir + simRoot + "/"
         check_dir(simDir)
         runEvalMFAToMAF(sim + "true.mfa", simDir + "true.maf")
         logger.info("Converted true.mfa of simulation %s to %strue.maf\n" % (sim, simDir))
         #Get path to sequence file of each species
         sequenceFiles = ''
         for spc in self.options.species:
            sequenceFiles = sequenceFiles + sim + spc + " "
         for cf in self.options.config:
      	    self.addChildTarget(CactusWorkflowWrapper(sim, cf, self.options.dir, sequenceFiles, self.options.tree))
            logger.info("Added child CactusWorkflowWrapper for sim %s and confi %s\n" % (sim, cf))

class CactusWorkflowWrapper(Target):
   """runCactusWorkFlow and issue child Target to generate MAF for the cactus results
   """
   def __init__(self, simulation, config, dir, sequenceFiles, tree):
      Target.__init__(self)
      self.simulation = simulation
      self.config = config
      self.dir = dir
      self.sequenceFiles = sequenceFiles
      self.tree = tree
   def run(self, localTempDir, globalTempDir):
      #------------------- Run cactus_workflow.py ---------------------#
      simRoot = getRootDir(self.simulation)
      cfRoot = getRoot(self.config)
      outDirName = self.dir + simRoot + "/" + cfRoot + "/"
      check_dir(outDirName)
      netdisk = outDirName + "netDisk"
      jobtreeDir = outDirName + "jobTree"
      batchSystem = "parasol"
      command = "cactus_workflow.py --speciesTree='" + self.tree + "' " + self.sequenceFiles + "--configFile " + self.config + " --buildTrees --setupAndBuildAlignments --netDisk " + netdisk + " --logDebug --job=JOB_FILE"
      starttime = time.clock()
      runJobTree(command, jobtreeDir, "DEBUG", 0, 'parasol', None)
      #runCactusWorkflow(netdisk, self.sequenceFiles, self.tree, jobtreeDir, "DEBUG", 0, batchSystem, None, True, True, False, False, self.config)
      runtime = time.clock() - starttime
      f = open(outDirName + "time", "w")
      f.write("%s\n" % (str(runtime)))
      f.close()
      logger.info("Done cactus_workflow for simulation %s, config %s\n" %(self.simulation, self.config))
      #self.addChildCommand(command)
      #------------------- Adding child ------------------------#
      self.addChildTarget(CactusMAFGeneratorWrapper(outDirName))
      logger.info("Added child CactusMAFGeneratorWrapper at %s\n" % outDirName)

class CactusMAFGeneratorWrapper(Target):
   """run cactus_MAFGenerator and issue child EvalMafComparatorWrapper
   """
   def __init__(self, dir):
      Target.__init__(self)
      self.dir = dir
   def run(self, localTempDir, globalTempDir):
      netdisk = self.dir + "netDisk"
      maffile = self.dir + "cactus.maf"
      runCactusMAFGenerator(mAFFile = maffile, netDisk = netdisk)
      truemaffile = self.dir + "../true.maf"
      self.addChildTarget(EvalMafComparatorWrapper(truemaffile, maffile, self.dir))

class EvalMafComparatorWrapper(Target):
   def __init__(self, maf1, maf2, dir):
      Target.__init__(self)
      self.maf1 = maf1
      self.maf2 = maf2
      self.outdir = dir
   def run(self, localTempDir, globalTempDir):
      sampleNumber = "1000000"
      runEvalMAFComparator(self.maf1, self.maf2, self.outdir + "sensitivity.xml", sampleNumber)
      runEvalMAFComparator(self.maf2, self.maf1, self.outdir + "selectivity.xml", sampleNumber)

def runEvalMAFComparator(mafFile1, mafFile2, outputFile, sampleNumber):
   command = "eval_MAFComparator -b %s -c %s -d %s -e %s" %(mafFile1, mafFile2, outputFile, sampleNumber)
   system(command)
   logger.info("Compared MAF %s with MAF %s\n" %(mafFile1, mafFile2))

def runEvalMFAToMAF(mfa, maf):
   command = "eval_MFAToMAF -b %s -d %s" %(mfa, maf)
   system(command)
   logger.info("Converted MFA %s to MAF %s\n" %(mfa, maf))

def modify_dirname(dir):
   """Add slash / at the end of the directory name if it doesnt have yet"""
   if (not re.search('/$', dir)): #not end with /
      dir = dir + '/'
   return dir

def check_dir(path):
   """Check if directories on the path, and create them if not."""
   if not os.path.exists(path):
      os.makedirs(path)

def getList(file):
   f = open(file, 'r')
   list = []
   for line in f.readlines():
      list.append(line.rstrip()) 
   f.close()
   return list

def getFirstLine(file):
   f = open(file, 'r')
   line = f.readline().rstrip()
   f.close()
   return line

def getRoot(path):
   pathLi = path.split('/')
   if len(pathLi) < 1:
      return ''
   else:
      li = pathLi[len(pathLi) -1].split('.')
      return li[0]

def getRootDir(path):
   pathLi = path.split('/')
   if len(pathLi) < 2:
      return ''
   else:
      li = pathLi[len(pathLi) -2].split('.')
      return li[0]

def main():
   usg = "Usage: %prog [options]\n"
   parser = OptionParser(usage=usg)
   parser.add_option("-d", "--simList", dest="sim", help="List of simulation directories. Default: simulations.lst", default="simulations.lst")
   parser.add_option("-c", "--configList", dest="config", help="List of cactus_workflow_config files. Default: configs.lst", default="configs.lst")
   parser.add_option("-o", "--outputDir", dest="dir", help="Directory for the outputs of the runs. Default: out", default="out/")
   parser.add_option("-t", "--tree", dest="tree", help="Phylogeny tree of the species of interest, in Newick format.Default: tree", default="tree")
   parser.add_option("-s", "--species", dest="species", help="List of species in the order as they appear in the  Newick tree. Default: species.lst", default="species.lst")
   parser.add_option("-j", "--job", dest="jobFile", help="Job file containing command to run.", default=None)
   (options, args) = parser.parse_args()
   #Process options:
   options.dir = modify_dirname(options.dir)
   options.tree = getFirstLine(options.tree)
   #assert options.tree == ''
   options.species = getFirstLine(options.species).split()
   #assert len(options.species) == 0
   options.sim = getList(options.sim)
   #assert len(options.sim) == 0
   options.config = getList(options.config)
   #assert len(options.config) == 0
   logger.info("Processed options\n")
   #Tuning
   cactusTuningWrapper = CactusTuningWrapper(options)
   cactusTuningWrapper.execute(options.jobFile)

if __name__ == "__main__":
   main()
