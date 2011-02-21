#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Wrapper to run cactus on different combinations of cactus_workflow_config.xml (cactus
parameters) and simulations, followed by evaluations steps.
"""
#nknguyen@soe.ucsc.edu
#05/26/2010

import os, re, sys, time
from optparse import OptionParser
import xml.etree.ElementTree as ET

from jobTree.src.jobTree import runJobTree  
from jobTree.scriptTree.target import Target

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel
from cactus.shared.common import cactusRootPath

from cactus.shared.common import runCactusWorkflow
from cactus.shared.common import runCactusMAFGenerator 
from cactus.shared.common import runCactusTreeStats

class CactusTuningWrapper(Target):
   """Wrapper to run cactus on different sets of parameters and different simulation data
   """
   def __init__(self, options):
      Target.__init__(self)
      self.options = options

   def run(self):
      #--------------------------------------------------------------------------------------
      #Get parameter sets. For each set, issue job to run cactus on different simulation data
      #--------------------------------------------------------------------------------------
      setLogLevel("DEBUG")
      system("rm -rf %s*" % self.options.outputDir)
      logger.info("Remove output directory if exists\n")
      
      #Convert true.mfa of each simulation to maf format
      #simTrueMafDir = os.path.join(self.options.outputDir, "sim")
      simTrueMafDir = self.options.simTrueMafDir
      check_dir(simTrueMafDir)
      for sim in self.options.sim:
         #convert mfa file of current simulation into MAF format:
         sim = modify_dirname(sim)
         simName = getRootDir(sim)
         trueMAF = os.path.join(simTrueMafDir, "%s_true.maf" %(simName))
         if not os.path.exists(trueMAF):
            trueMFA = os.path.join(sim, "true.mfa")
            runEvalMFAToMAF(trueMFA, trueMAF)
            logger.info("Converted true.mfa of simulation %s to %s\n" % (sim, trueMAF))
         else:
            logger.info("TrueMAF already exists: %s\n" %(trueMAF))

      for parameterFile, parameterName in getParameters(self.options.config):
         outDir = os.path.join(self.options.outputDir, parameterName)
	 #system("rm -rf %s" % outDir)
	 os.mkdir(outDir)
	 system("mv %s %s/" % (parameterFile, outDir))
         logger.info("Created output directory %s for parameter set %s and moved config file to that directory\n" % (outDir, parameterName))
         paraFile = os.path.join(outDir, 'param.xml')
         statsDir = os.path.join(outDir, "stats")
	 os.mkdir(statsDir)
         logger.info("Created directory for stats files: %s\n" % (statsDir))
         self.addChildTarget(CactusTuningSimulationsWrapper(self.options, paraFile, outDir))
         logger.info("Added CactusTuningSimulationsWrapper as child for parameter %s\n" %(parameterName))
      #Summarize results
      #self.setFollowOnTarget(CactusTuningSummary(self.options))
      logger.info("Added CactusTuningSummary\n")

class CactusTuningSimulationsWrapper(Target):
   """Run cactus for a set of different simulation data and report results
   """
   def __init__(self, options, paraFile, outDir):
      Target.__init__(self)
      self.options = options
      self.paraFile = paraFile
      self.outDir = outDir

   def run(self):
      #--------------------------------------------
      #Run cactus & evaluations for each simulation
      #--------------------------------------------
      logger.info("CactusTuningSimulationsWrapper: going to issue cactus runs for all simulations for parameter %s\n" %(self.paraFile))
      simNum = 0
      for sim in self.options.sim:
         sim = modify_dirname(sim)
         simName = getRootDir(sim)
         
         #Get path to sequence file of each species
	 sequenceFiles = " ".join([ os.path.join(sim, spc) for spc in self.options.species ])
         logger.info("Got sequence files: %s\n" % (sequenceFiles))

	 #add child
      	 #self.addChildTarget(CactusWorkflowWrapper(sim, simNum, self.paraFile, self.outDir, sequenceFiles, self.options.tree))
      	 self.addChildTarget(CactusWorkflowWrapper(sim, simName, self.options.simTrueMafDir, self.paraFile, self.outDir, sequenceFiles, self.options.tree))
         logger.info("Added child CactusWorkflowWrapper for sim %s and confi %s\n" % (sim, self.paraFile))
         simNum += 1
	
      #----------------------------------------------------------------
      #Done running cactus & evaluations steps for all the simulations. 
      #Now Merge results & clean up.
      #----------------------------------------------------------------
      logger.info("Done running cactus & evaluations for parameter %s. Now merge results and clean up.\n" %(self.paraFile))
      self.setFollowOnTarget(CactusMergeResultsAndCleanup(simNum, self.outDir, self.options))
      logger.info("Added CactusMergeResultsAndCleanup as FollowOnTarget for %s\n" %(self.outDir))

class CactusWorkflowWrapper(Target):
   """runCactusWorkFlow and issue child Target to generate MAF for the cactus results
   """
   #def __init__(self, simulation, simNum, paraFile, outDir, sequenceFiles, tree):
   def __init__(self, simulation, simName, simTrueMafDir, paraFile, outDir, sequenceFiles, tree):
      Target.__init__(self)
      self.simulation = simulation
      #self.simNum = str(simNum)
      self.simName = simName
      self.simTrueMafDir = simTrueMafDir
      self.paraFile = paraFile
      self.outDir = outDir
      self.sequenceFiles = sequenceFiles
      self.tree = tree

   def run(self):
      #----------------------------------------
      # Run cactus_workflow.py and report time#
      #----------------------------------------
      logger.info("CactusWorkflowWrapper: going to issue cactus run for simulation %s, parameter %s\n" %(self.simulation, self.paraFile))
      tempDir = getTempDirectory(self.outDir)
      flowerdisk = os.path.join(tempDir, "cactusDisk")
      jobtreeDir = os.path.join(tempDir, "jobTree")
      #batchSystem = "single_machine"
      batchSystem = "parasol"
      retryCount = 0
      command = "cactus_workflow.py --speciesTree='%s' %s --configFile %s --buildTrees --setupAndBuildAlignments --cactusDisk %s --logDebug --job=JOB_FILE" %(self.tree, self.sequenceFiles, self.paraFile, flowerdisk)
      starttime = time.time()
      runJobTree(command, jobtreeDir, "DEBUG", retryCount, batchSystem, None)
      #runCactusWorkflow(flowerdisk, self.sequenceFiles, self.tree, jobtreeDir, "DEBUG", 0, batchSystem, None, True, True, False, False, self.config)
      runtime = time.time() - starttime
      logger.info("Done cactus_workflow for simulation %s, config %s\n" %(self.simulation, self.paraFile))

      #-----------------------
      # Run cactus_treeStats #
      #-----------------------
      #statsFile = os.path.join(self.outDir, "stats", "%s.xml" % self.simNum)
      statsFile = os.path.join(self.outDir, "stats", "%s.xml" % self.simName)
      runCactusTreeStats(outputFile=statsFile, cactusDisk=flowerdisk)
      #self.addChildCommand(command)

      #------------------- Adding child ------------------------#
      #self.addChildTarget(CactusMAFGeneratorWrapper(self.outDir, tempDir, self.simNum, runtime))
      self.addChildTarget(CactusMAFGeneratorWrapper(self.outDir, tempDir, self.simTrueMafDir, self.simName, runtime))
      logger.info("Added child CactusMAFGeneratorWrapper at %s\n" % self.outDir)

      #------------------- Cleaning up -------------------------#
      self.setFollowOnTarget(CactusWorkflowWrapperCleanup(tempDir))


class CactusMAFGeneratorWrapper(Target):
   """run cactus_MAFGenerator and issue child EvalMafComparatorWrapper
   """
   #def __init__(self, outDir, resultsDir, simNum, cactusRunTime):
   def __init__(self, outDir, resultsDir, simTrueMafDir, simName, cactusRunTime):
      Target.__init__(self)
      self.outDir = outDir
      self.resultsDir = resultsDir #Directory contains cactus cactusDisk and jobTree
      #self.simNum = simNum
      self.simTrueMafDir = simTrueMafDir
      self.simName = simName
      self.cactusRunTime = cactusRunTime

   def run(self):
      flowerdisk = os.path.join(self.resultsDir, "cactusDisk")
      maffile = os.path.join(self.resultsDir, "cactus.maf")
      runCactusMAFGenerator(mAFFile = maffile, cactusDisk = flowerdisk)
      #truemaffile = os.path.join(self.outDir,"..","sim", "%s_true.maf" %(self.simNum))
      #mafCompareFile = os.path.join(self.outDir, "mafCompare%s.xml" %self.simNum)
      truemaffile = os.path.join(self.simTrueMafDir, "%s_true.maf" %(self.simName))
      mafCompareFile = os.path.join(self.outDir, "mafCompare%s.xml" %self.simName)
      self.addChildTarget(EvalMafComparatorWrapper(truemaffile, maffile, mafCompareFile, self.cactusRunTime))

class EvalMafComparatorWrapper(Target):
   def __init__(self, maf1, maf2, outputFile, time):
      Target.__init__(self)
      self.maf1 = maf1
      self.maf2 = maf2
      self.outputFile = outputFile
      self.time = time
   def run(self):
      sampleNumber = "1000000"
      runEvalMAFComparator(self.maf1, self.maf2, self.outputFile, sampleNumber)
      #Add the the run time to the results
      resultsNode = ET.parse(self.outputFile).getroot()
      resultsNode.attrib["time"] = str(self.time)
      fileHandle = open(self.outputFile, 'w')
      ET.ElementTree(resultsNode).write(fileHandle)
      fileHandle.close()


class CactusWorkflowWrapperCleanup(Target):
   def __init__(self, dir):
      Target.__init__(self)
      self.dir = dir
   def run(self):
      system("rm -rf %s" % self.dir)
      logger.info("Clean up tempDir for next run\n")

class CactusMergeResultsAndCleanup(Target):
   """
   """
   def __init__(self, count, outDir, options):
      Target.__init__(self)
      self.count = count #number of files to merge
      self.outDir = outDir
      self.options = options
      
   def run(self):
      mergedFile = os.path.join(self.outDir, "mafCompare.xml")
      count = 0
      for sim in self.options.sim:
         simName = getRootDir(modify_dirname(sim))
         currentFile = os.path.join(self.outDir, "mafCompare%s.xml" %simName)
         if count == 0:
            system("mv %s %s" % (currentFile, mergedFile))
            logger.info("Moved %s to %s\n" %(currentFile, mergedFile))
         else:
	    system("mergeMafComparatorResults.py --logLevel DEBUG --results1 %s --results2 %s --outputFile %s" % (mergedFile, currentFile, mergedFile))
            logger.info("Merged %s to %s\n" %(currentFile, mergedFile))
         count += 1
         #system("rm -f %s" % currentFile)
	 #logger.info("Removed %s\n" %(currentFile))

class CactusTuningSummary(Target):
   """
   """
   def __init__(self, options):
      Target.__init__(self)
      self.options = options
      
   def run(self):
      getCactusTuningSummary(self.options.outputDir, self.options.species, self.options.sim)

#============================ Getting parameters =======================================#
def fn(min, max, loops):
   if loops == 0 or loops == 1:
      return 0
   return float(min - max)/(loops - 1) #The value must be zero

def getParameters(startFile):
   for minimumTreeCoverage in (0.0,):
      for annealingRounds in (5,): #1+trim+minimumChainLength/2):
         for lastzThreshold in (1800,):
            for (minimumTrim, maximumTrim) in ((0, 3),):
               trimChange = fn(minimumTrim, maximumTrim, annealingRounds)
               for minimumChainLength, maximumChainLength in ((5,30),):
                  minimumChainLengthChange = fn(minimumChainLength, maximumChainLength, annealingRounds)
                  for minimumBlockLength , maximumBlockLength in ((0, 0),):
                     minimumBlockLengthChange = fn(minimumBlockLength, maximumBlockLength, annealingRounds)
                     for alignRepeatsAtRound in set((0,)):
                        for deannealingRounds in set((10,)):
                           for baseLevel in (True,):
   #for minimumTreeCoverage in (0.0,0.5,1.0):
   #   for annealingRounds in (5,): #1+trim+minimumChainLength/2):
   #      for lastzThreshold in (1800,2200,2600,3000):
   #         for (minimumTrim, maximumTrim) in ((0, 3),(1, 3),(2, 3),(3, 5)):
   #            trimChange = fn(minimumTrim, maximumTrim, annealingRounds)
   #            for minimumChainLength, maximumChainLength in ((5,30), (5,100), (10,30),(10, 100),(20,30),(20, 100)):
   #               minimumChainLengthChange = fn(minimumChainLength, maximumChainLength, annealingRounds)
   #               for minimumBlockLength , maximumBlockLength in ((0, 0),(2, 2),):
   #                  minimumBlockLengthChange = fn(minimumBlockLength, maximumBlockLength, annealingRounds)
   #                  for alignRepeatsAtRound in set((0,)):
   #                     for deannealingRounds in set((10,)):
   #                        for baseLevel in (True, False):
                               config = ET.parse(startFile).getroot()
                               iterationNode = config.find("alignment").find("iterations").findall("iteration")[-2]
                               #node = config.find("alignment").find("iterations").findall("iteration")[-2].find("core")

                               blastNode = iterationNode.find("blast")
                               blastNode.attrib["blastString"] = "lastz --format=cigar --hspthresh=%s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE" % lastzThreshold
                               blastNode.attrib["selfBlastString"]="lastz --format=cigar --hspthresh=%s SEQ_FILE[nameparse=darkspace] --self > CIGARS_FILE" % lastzThreshold

                               node = iterationNode.find("core")
                               node.attrib["minimumTreeCoverage"] = str(minimumTreeCoverage)
                               node.attrib["annealingRounds"] = str(annealingRounds)
                               node.attrib["trim"] = str(maximumTrim)
                               node.attrib["trimChange"] = str(trimChange)                                        
                               node.attrib["alignRepeatsAtRound"] = str(alignRepeatsAtRound)
                               node.attrib["minimumBlockLength"] = str(maximumBlockLength)
                               node.attrib["minimumBlockLengthChange"] = str(minimumBlockLengthChange)
                               node.attrib["minimumChainLength"] = str(maximumChainLength)
                               node.attrib["minimumChainLengthChange"] = str(minimumChainLengthChange)
                               node.attrib["deannealingRounds"] = str(deannealingRounds)

                               #Remove the base alignment stage:
                               if not baseLevel:
                                  config.find("alignment").find("iterations").remove(config.find("alignment").find("iterations").findall("iteration")[-1])

                               paramFile = os.path.join(os.getcwd(), "param.xml")
                               fileHandle = open(paramFile, 'w')
                               tree = ET.ElementTree(config)
                               tree.write(fileHandle)
                               fileHandle.close()
                               yield (paramFile, ("results_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s" % (minimumTreeCoverage, annealingRounds, maximumTrim, trimChange, alignRepeatsAtRound, maximumChainLength, minimumChainLengthChange, maximumBlockLength, minimumBlockLengthChange, deannealingRounds, lastzThreshold, baseLevel)))
                               #yield (paramFile, ("results_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s" % (minimumTreeCoverage, annealingRounds, maximumTrim, trimChange, alignRepeatsAtRound, maximumChainLength, minimumChainLengthChange, maximumBlockLength, minimumBlockLengthChange, deannealingRounds)))

#============================ Results summary ==========================================#
def getCactusTuningSummary(dir, species, sim):
   l = []
   #Use the stat file of first simulation as an estimate for stats of other simulations:
   firstSimName = getRootDir(modify_dirname(sim[0]))
   for resultsName in os.listdir(dir):
      results = os.path.join(dir, resultsName)
      statsDir = os.path.join(results, "stats")
      #maxBlockMaxDegree = str(getStats(statsDir))
      try:
         stats = ET.parse(os.path.join(statsDir, "%s.xml" %(firstSimName))).getroot()
         config = ET.parse(os.path.join(results, "param.xml")).getroot()
         scores = ET.parse(os.path.join(results, "mafCompare.xml")).getroot()
      except IOError:
         continue

      blocksNode = stats.find("blocks")
      sensNode = scores.findall("homology_tests")[0]
      specNode = scores.findall("homology_tests")[1]
      if len(config.find("alignment").find("iterations").findall("iteration")) == 4:
         baseLevel = True
         iterationNode = config.find("alignment").find("iterations").findall("iteration")[-2]
      else:
         baseLevel = False
         iterationNode = config.find("alignment").find("iterations").findall("iteration")[-1]

      node = iterationNode.find("core")
      #node = config.find("alignment").find("iterations").findall("iteration")[-2].find("core")
      #node = config.find("alignment").find("iterations").findall("iteration")[-1].find("core")#Because obmit the base level

      l.append((sensNode.attrib["average"], specNode.attrib["average"], node.attrib["minimumTreeCoverage"], node.attrib["annealingRounds"], \
      node.attrib["trim"], node.attrib["trimChange"], \
      node.attrib["alignRepeatsAtRound"], node.attrib["minimumChainLength"], \
      node.attrib["minimumChainLengthChange"], node.attrib["deannealingRounds"], \
      node.attrib["minimumBlockLength"], node.attrib["minimumBlockLengthChange"], \
      node.attrib["deannealingRounds"],
      str(baseLevel),
      blocksNode.find("degrees").attrib["max"], scores.attrib["time"], \
      fn("HUMAN", "MOUSE", sensNode), fn("HUMAN", "MOUSE", specNode),\
      fn("HUMAN", "DOG", sensNode), fn("HUMAN", "DOG", specNode),\
      fn("HUMAN", "CHIMP", sensNode), fn("HUMAN", "CHIMP", specNode),\
      fn("HUMAN", "BABOON", sensNode), fn("HUMAN", "BABOON", specNode),\
      fn("MOUSE", "RAT", sensNode), fn("MOUSE", "RAT", specNode),\
      fn("DOG", "COW", sensNode), fn("DOG", "COW", specNode),\
      fn("MOUSE", "DOG", sensNode), fn("MOUSE", "DOG", specNode),
      stats.find("terminal_group_sizes").attrib["max"],
      stats.find("chains").find("base_block_lengths").attrib["median"],
      stats.find("chains").find("base_block_lengths").attrib["avg"],
      stats.find("chains").find("base_block_lengths").attrib["max"],
      fn2("tangles", stats), fn2("links", stats), fn3(iterationNode.find("blast"))))

      
      #currList = [(sensNode.attrib["average"], specNode.attrib["average"], node.attrib["minimumTreeCoverage"], node.attrib["annealingRounds"], \
      #node.attrib["trim"], node.attrib["trimChange"], \
      #node.attrib["alignRepeatsAtRound"], node.attrib["minimumChainLength"], \
      #node.attrib["minimumChainLengthChange"], node.attrib["deannealingRounds"], \
      #node.attrib["minimumBlockLength"], node.attrib["minimumBlockLengthChange"], \
      #node.attrib["deannealingRounds"], \
      #blocksNode.find("degrees").attrib["max"], scores.attrib["time"], \
      #stats.find("terminal_group_sizes").attrib["max"], \
      #stats.find("chains").find("base_block_lengths").attrib["median"], \
      #stats.find("chains").find("base_block_lengths").attrib["avg"], \
      #stats.find("chains").find("base_block_lengths").attrib["max"], \
      #fn2("tangles", stats), fn2("links", stats))]


      #for i in range(len(species) -1):
      #   for j in range(i+1, len(species)):
      #      currList.extend(fn1(species[i], species[j], sensNode))
      #      currList.extend(fn1(species[i], species[j], specNode))
      #l.append(currList)

      stats = None
      config = None
      scores = None

   l.sort(cmpFn)
   l2 = ("sens\t\tspec\t\tmTCo", "AR", "trim", "trimR", "ARaL", "mCL", "mCLI", "mCLUSS", "mBL", "mBLC", "DAR", "BASE", "deg", "time", \
      "HMSE", "HMSP", "HDSE", "HDSP", "HCSE", "HCSP", "HBSE", "HBSP", "MRSE", "MRSP", "DCSE", "DCSP", "MDSE", "MDSP", \
      "TGS", "CML", "CAL", "CMXL", "NTTN", "NTLN")

   #l2 = ["sens", "spec", "mTCo", "AR", "trim", "trimR", "ARaL", "mCL", "mCLI", "mCLUSS", "mBL", "mBLC", "DAR", "deg", "time", \
   #   "TGS", "CML", "CAL", "CMXL", "NTTN", "NTLN"]
   #for i in range(len(species) -1):
   #   for j in range(i+1, len(species)):
   #      sensCol = "%s_%s_SE" %(species[i], species[j])
   #      specCol = "%s_%s_SP" %(species[i], species[j])
   #      l2.extend((sensCol, specCol))
         #"HMSE", "HMSP", "HDSE", "HDSP", "HCSE", "HCSP", "HBSE", "HBSP", "MRSE", "MRSP", "DCSE", "DCSP", "MDSE", "MDSP", \

   outFile = os.path.join(dir, "summary")
   f = open(outFile, 'w')
   f.write("\t".join(l2) + "\n")
   for i in l:
      f.write("\t".join(i) + "\n")
      #for c in i:
      #   f.write(c + "\t")
      #f.write("\n")
   f.write("\t".join(l2) + "\n")
   f.close()

def fn1(speciesA, speciesB, node):
    for hTest in node.findall("homology_test"):
        #print hTest.attrib, speciesA, speciesB
        if (hTest.attrib["sequenceA"] == speciesA and hTest.attrib["sequenceB"] == speciesB) or \
            (hTest.attrib["sequenceA"] == speciesB and hTest.attrib["sequenceB"] == speciesA):
            return hTest.attrib["average"]
    assert False

def fn2(type, node):
    for ends in node.findall("ends"):
        if ends.attrib["include_terminal_groups"] == '1' and ends.attrib["include_non_terminal_groups"] == '0':
            if type == "tangles" and ends.attrib["include_tangle_groups"] == '1' and ends.attrib["include_link_groups"] == '0':
                return ends.find("counts").attrib["total"]
            if type == "links" and ends.attrib["include_tangle_groups"] == '0' and ends.attrib["include_link_groups"] == '1':
                return ends.find("counts").attrib["total"]
    assert False

def fn3(node):
    return node.attrib["blastString"].split()[2]

def getStats(statsDir):
   maxBlockMaxDegree = 0
   statsList = os.listdir(statsDir)
   for s in statsList:
      statsFile = os.path.join(statsDir, s)
      try:
         stats = ET.parse(statsFile).getroot()
      except IOError:
         continue

      blocksNode = stats.find("blocks")
      blockMaxDegree = int(blocksNode.find("degrees").attrib["max"])
      if maxBlockMaxDegree < blockMaxDegree:
         maxBlockMaxDegree = blockMaxDegree
   
   return maxBlockMaxDegree   

#============================ Utilities functions ======================================#
def runEvalMAFComparator(mafFile1, mafFile2, outputFile, sampleNumber):
   command = "mafComparator -b %s -c %s -d %s -e %s" %(mafFile1, mafFile2, outputFile, sampleNumber)
   system(command)
   logger.info("Compared MAF %s with MAF %s\n" %(mafFile1, mafFile2))

def runEvalMFAToMAF(mfa, maf):
   command = "mfaToMaf -b %s -d %s --logLevel DEBUG" %(mfa, maf)
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
   parser.add_option("-m", "--simTrueMafDir", dest="simTrueMafDir", help="Directory for 'true' mafs of the simulations. Default: sim/", default="sim/")
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

