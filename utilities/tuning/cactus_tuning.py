#!/usr/bin/env python

"""Wrapper to run cactus on different combinations of cactus_workflow_config.xml (cactus
parameters) and simulations, followed by evaluations steps.
"""
#nknguyen@soe.ucsc.edu
#05/26/2010

import os, re, sys, time
from optparse import OptionParser
import xml.etree.ElementTree as ET


from workflow.jobTree.jobTree import runJobTree
from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete
from workflow.jobTree.scriptTree.target import Target

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel
from sonLib.misc import sonTraceRootPath

from cactus.shared.common import runCactusWorkflow
from cactus.shared.common import runCactusMAFGenerator
from cactus.shared.common import runCactusTreeStats

class CactusTuningWrapper(Target):
   """Wrapper to run cactus on different sets of parameters and different simulation data
   """
   def __init__(self, options):
      Target.__init__(self)
      self.options = options

   def run(self, localTempDir, globalTempDir):
      #--------------------------------------------------------------------------------------
      #Get parameter sets. For each set, issue job to run cactus on different simulation data
      #--------------------------------------------------------------------------------------
      setLogLevel("DEBUG")
      system("rm -rf %s/*" % self.options.outputDir)
      for parameterFile, parameterName in getParameters(self.options.config):
         outDir = os.path.join(self.options.outputDir, parameterName)
	 #system("rm -rf %s" % outDir)
	 os.mkdir(outDir)
	 system("cp %s %s/" % (parameterFile, outDir))
         paraFile = os.path.join(outDir, parameterFile)
         statsDir = os.path.join(outDir, "stats")
	 os.mkdir(statsDir)
         self.addChildTarget(CactusTuningSimulationsWrapper(self.options, paraFile, outDir))
      #Summarize results
      self.setFollowOnTarget(CactusTuningSummary(self.options.outputDir))

class CactusTuningSimulationsWrapper(Target):
   """Run cactus for a set of different simulation data and report results
   """
   def __init__(self, options, paraFile, outDir):
      Target.__init__(self)
      self.options = options
      self.paraFile = paraFile
      self.outDir = outDir

   def run(self, localTempDir, globalTempDir):
      #--------------------------------------------
      #Run cactus & evaluations for each simulation
      #--------------------------------------------
      simNum = 0
      for sim in self.options.sim:
         sim = modify_dirname(sim)

         #convert mfa file of current simulation into MAF format:
	 trueMFA = os.path.join(sim, "true.mfa")
	 trueMAF = os.path.join(self.outDir, "true.maf")
         runEvalMFAToMAF(trueMFA, trueMAF)
         logger.info("Converted true.mfa of simulation %s to %s\n" % (sim, trueMAF))

         #Get path to sequence file of each species
	 sequenceFiles = " ".join([ os.path.join(sim, spc) for spc in self.options.species ])
         logger.info("Got sequence files: %s\n" % (sequenceFiles))

	 #add child
      	 self.addChildTarget(CactusWorkflowWrapper(sim, simNum, self.paraFile, self.outDir, sequenceFiles, self.options.tree))
         logger.info("Added child CactusWorkflowWrapper for sim %s and confi %s\n" % (sim, self.paraFile))
         simNum += 1
	
      #----------------------------------------------------------------
      #Done running cactus & evaluations steps for all the simulations. 
      #Now Merge results & clean up.
      #----------------------------------------------------------------
      logger.info("Done running cactus & evaluations for parameter %s. Now merge results and clean up.\n" %(self.paraFile))
      self.setFollowOnTarget(CactusMergeResultsAndCleanup(simNum, self.outDir))

class CactusWorkflowWrapper(Target):
   """runCactusWorkFlow and issue child Target to generate MAF for the cactus results
   """
   def __init__(self, simulation, simNum, paraFile, outDir, sequenceFiles, tree):
      Target.__init__(self)
      self.simulation = simulation
      self.simNum = str(simNum)
      self.paraFile = paraFile
      self.outDir = outDir
      self.sequenceFiles = sequenceFiles
      self.tree = tree

   def run(self, localTempDir, globalTempDir):
      #----------------------------------------
      # Run cactus_workflow.py and report time#
      #----------------------------------------
      tempDir = getTempDirectory(self.outDir)
      netdisk = os.path.join(tempDir, "netDisk")
      jobtreeDir = os.path.join(tempDir, "jobTree")
      batchSystem = "single_machine"
      command = "cactus_workflow.py --speciesTree='%s' %s --configFile %s --buildTrees --setupAndBuildAlignments --netDisk %s --logDebug --job=JOB_FILE" %(self.tree, self.sequenceFiles, self.paraFile, netdisk)
      starttime = time.time()
      runJobTree(command, jobtreeDir, "DEBUG", 0, batchSystem, None)
      #command = "cactus_workflow.py --speciesTree='" + self.tree + "' " + self.sequenceFiles + "--configFile " + self.parameterFile + " --buildTrees --setupAndBuildAlignments --netDisk " + netdisk + " --logDebug --job=JOB_FILE"
      #runCactusWorkflow(netdisk, self.sequenceFiles, self.tree, jobtreeDir, "DEBUG", 0, batchSystem, None, True, True, False, False, self.config)
      runtime = time.time() - starttime
      timeFile = os.path.join(self.outDir, "time")
      if not os.path.isfile(timeFile):
         f = open(timeFile, "w")
      else:
         f = open(timeFile, "a")
      f.write("%s\n" % (str(runtime)))
      f.close()
      logger.info("Done cactus_workflow for simulation %s, config %s\n" %(self.simulation, self.paraFile))

      #-----------------------
      # Run cactus_treeStats #
      #-----------------------
      statsFile = os.path.join(self.outDir, "stats", "%s.xml" % self.simNum)
      runCactusTreeStats(outputFile=statsFile, netDisk=netdisk)
      #self.addChildCommand(command)

      #------------------- Adding child ------------------------#
      self.addChildTarget(CactusMAFGeneratorWrapper(self.outDir, tempDir, self.simNum))
      logger.info("Added child CactusMAFGeneratorWrapper at %s\n" % self.outDir)

      #------------------- Cleaning up -------------------------#
      self.setFollowOnTarget(CactusWorkflowWrapperCleanup(tempDir))


class CactusMAFGeneratorWrapper(Target):
   """run cactus_MAFGenerator and issue child EvalMafComparatorWrapper
   """
   def __init__(self, outDir, resultsDir, simNum):
      Target.__init__(self)
      self.outDir = outDir
      self.resultsDir = resultsDir #Directory contains cactus netDisk and jobTree
      self.simNum = simNum
   def run(self, localTempDir, globalTempDir):
      netdisk = os.path.join(self.resultsDir, "netDisk")
      maffile = os.path.join(self.resultsDir, "cactus.maf")
      runCactusMAFGenerator(mAFFile = maffile, netDisk = netdisk)
      truemaffile = os.path.join(self.outDir, "true.maf")
      mafCompareFile = os.path.join(self.outDir, "mafCompare%s.xml" %self.simNum)
      self.addChildTarget(EvalMafComparatorWrapper(truemaffile, maffile, mafCompareFile))

class EvalMafComparatorWrapper(Target):
   def __init__(self, maf1, maf2, outputFile):
      Target.__init__(self)
      self.maf1 = maf1
      self.maf2 = maf2
      self.outputFile = outputFile
   def run(self, localTempDir, globalTempDir):
      sampleNumber = "1000000"
      runEvalMAFComparator(self.maf1, self.maf2, self.outputFile, sampleNumber)

class CactusWorkflowWrapperCleanup(Target):
   def __init__(self, dir):
      Target.__init__(self)
      self.dir = dir
   def run(self, localTempDir, globalTempDir):
      system("rm -rf %s" % self.dir)
      logger.info("Clean up tempDir for next run\n")

class CactusMergeResultsAndCleanup(Target):
   """
   """
   def __init__(self, count, outDir):
      Target.__init__(self)
      self.count = count #number of files to merge
      self.outDir = outDir
      
   def run(self, localTempDir, globalTempDir):
      mergedFile = os.path.join(self.outDir, "mafCompare.xml")
      currentFile = os.path.join(self.outDir, "mafCompare0.xml")
      system("mv %s %s" % (currentFile, mergedFile))
      logger.info("Moved %s to %s\n" %(currentFile, mergedFile))
      for i in range(1, self.count):
         currentFile = os.path.join(self.outDir, "mafCompare%s.xml" %(str(i)))
	 system("eval_mergeMAFComparatorResults.py --logLevel DEBUG --results1 %s --results2 %s --outputFile %s" % (mergedFile, currentFile, mergedFile))
         logger.info("Merged %s to %s\n" %(currentFile, mergedFile))
         system("rm -f %s" % currentFile)
	 logger.info("Removed %s\n" %(currentFile))

class CactusTuningSummary(Target):
   """
   """
   def __init__(self, outDir):
      Target.__init__(self)
      self.outDir = outDir
      
   def run(self, localTempDir, globalTempDir):
      getCactusTuningSummary(self.outDir)

#============================ Getting parameters =======================================#
def getParameters(startFile):
   #for minimumTreeCoverage in [0.0]:
   #   for alignUndoLoops in [1, 10]:
   #      for (trim, trimChange) in [(0, 0)]:
   #         for alignRepeatsAtLoop in [0]:
   #            for minimumChainLength in [20]:
   #               for minimumChainLengthChange in [0]:
   for minimumTreeCoverage in (0.0, 0.01, 0.1):
      for alignUndoLoops in (1, 10):
         for (trim, trimChange) in set(((0, 0), (alignUndoLoops-1, -1), (alignUndoLoops-1, 0))):
            for alignRepeatsAtLoop in set((0, alignUndoLoops, alignUndoLoops/2)):
               for minimumChainLength in (20, 40, 80):
                  for minimumChainLengthChange in (0, minimumChainLength/5):
                     minimumChainLengthCactusUndoLoopStepSize=minimumChainLength/10
                     #startFile = os.path.join(sonTraceRootPath(), "src", "cactus", "pipeline", "cactus_workflow_config.xml")
                     config = ET.parse(startFile).getroot()
                     node = config.find("alignment").find("iterations").findall("iteration")[-2].find("core")

                     node.attrib["minimumTreeCoverage"] = str(minimumTreeCoverage)
                     node.attrib["alignUndoLoops"] = str(alignUndoLoops)
                     node.attrib["trim"] = str(trim)
                     node.attrib["trimChange"] = str(trimChange)
                     node.attrib["alignRepeatsAtLoop"] = str(alignRepeatsAtLoop)
                     node.attrib["minimumChainLength"] = str(minimumChainLength)
                     node.attrib["minimumChainLengthChange"] = str(minimumChainLengthChange)
                     node.attrib["minimumChainLengthCactusUndoLoopStepSize"] = str(minimumChainLengthCactusUndoLoopStepSize)

                     paramFile = os.path.join(os.getcwd(), "param.xml")
                     fileHandle = open(paramFile, 'w')
                     tree = ET.ElementTree(config)
                     tree.write(fileHandle)
                     fileHandle.close()
                     yield (paramFile, ("results_%s_%s_%s_%s_%s_%s_%s" % (minimumTreeCoverage, alignUndoLoops, trim, trimChange, alignRepeatsAtLoop, minimumChainLength, minimumChainLengthChange)))

#============================ Results summary ==========================================#
def getCactusTuningSummary(dir):
   l = []
   for resultsName in os.listdir(dir):
      results = os.path.join(dir, resultsName)
      statsDir = os.path.join(results, "stats")
      avgBlockMaxDegree = str(getStats(statsDir))
      try:
         config = ET.parse(os.path.join(results, "param.xml")).getroot()
         scores = ET.parse(os.path.join(results, "mafCompare.xml")).getroot()
      except IOError:
         continue

      sensNode = scores.findall("homology_tests")[0]
      specNode = scores.findall("homology_tests")[1]
      node = config.find("alignment").find("iterations").findall("iteration")[-2].find("core")
      l.append((sensNode.attrib["average"], specNode.attrib["average"], node.attrib["minimumTreeCoverage"], node.attrib
["alignUndoLoops"], \
      node.attrib["trim"], node.attrib["trimChange"], \
      node.attrib["alignRepeatsAtLoop"], node.attrib["minimumChainLength"], \
      node.attrib["minimumChainLengthChange"], node.attrib["minimumChainLengthCactusUndoLoopStepSize"], \
      avgBlockMaxDegree))

   l.sort(cmpFn)
   l2 = ("sens\t\tspec\t\tmTCo", "aULo", "trim", "trimR", "ARaL", "mCL", "mCLI", "mCLUSS", "deg")
   outFile = os.path.join(dir, "summary")
   f = open(outFile, 'w')
   f.write("\t".join(l2) + "\n")
   for i in l:
      f.write("\t".join(i) + "\n")
   f.write("\t".join(l2) + "\n")
   f.close()

def getStats(statsDir):
   sumBlockMaxDegree = 0
   statsList = os.listdir(statsDir)
   for stats in statsList:
      statsFile = os.path.join(statsDir, stats)
      try:
         stats = ET.parse(statsFile).getroot()
      except IOError:
         continue

      blocksNode = stats.find("blocks")
      sumBlockMaxDegree += float(blocksNode.find("degrees").attrib["max"])
   
   return sumBlockMaxDegree/len(statsList)   

#============================ Utilities functions ======================================#
def runEvalMAFComparator(mafFile1, mafFile2, outputFile, sampleNumber):
   command = "eval_MAFComparator -b %s -c %s -d %s -e %s" %(mafFile1, mafFile2, outputFile, sampleNumber)
   system(command)
   logger.info("Compared MAF %s with MAF %s\n" %(mafFile1, mafFile2))

def runEvalMFAToMAF(mfa, maf):
   command = "eval_MFAToMAF -b %s -d %s --logLevel DEBUG" %(mfa, maf)
   system(command)
   logger.info("Converted MFA %s to MAF %s\n" %(mfa, maf))

def cmpFn(a, b):
   i = float(a[0])
   j = float(b[0])
   return cmp(i, j)

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
   parser.add_option("-c", "--configStartFile", dest="config", help="cactus_workflow_config.xml", default="cactus_workflow_config.xml")
   parser.add_option("-o", "--outputDir", dest="outputDir", help="Directory for the outputs of the runs. Default: out", default="out/")
   parser.add_option("-t", "--tree", dest="tree", help="Phylogeny tree of the species of interest, in Newick format.Default: tree", default="tree")
   parser.add_option("-s", "--species", dest="species", help="List of species in the order as they appear in the  Newick tree. Default: species.lst", default="species.lst")
   parser.add_option("-j", "--job", dest="jobFile", help="Job file containing command to run.", default=None)
   (options, args) = parser.parse_args()
   #Process options:
   options.outputDir = modify_dirname(options.outputDir)
   check_dir(options.outputDir)
   options.tree = getFirstLine(options.tree)
   #assert options.tree == ''
   options.species = getFirstLine(options.species).split()
   #assert len(options.species) == 0
   options.sim = getList(options.sim)
   #assert len(options.sim) == 0
   #options.config = getList(options.config)
   #assert len(options.config) == 0
   logger.info("Processed options\n")
   #Tuning
   cactusTuningWrapper = CactusTuningWrapper(options)
   cactusTuningWrapper.execute(options.jobFile)

if __name__ == "__main__":
   main()

