#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#04/17/2011

"""
Wrapper to run cactus on multiple regions (same set of species/input sequences), and aggregate the genemap results.
"""

import os, sys, re, time
from optparse import OptionParser
import xml.etree.ElementTree as ET

#from jobTree.src.jobTree import runJobTree
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getTempFile

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

#from mhc.genemap.lib import *

class Setup(Target):
    """Set up cactus runs, and subsequent analysis for all regions
    """
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
        self.options = options

    def run(self):
        setLogLevel("DEBUG")
        system("mkdir -p %s" %self.options.outdir)
        system("chmod ug+xrw %s" %self.options.outdir)

        regions = getList(self.options.regions)
        logger.info("Regions: %s\n" %(" ".join(regions)))

        for r in regions:
            self.addChildTarget(RunRegion(r, self.options))

	self.setFollowOnTarget(AggregateResults(self.options))

class RunRegion(Target):
    """Get experiment.xml, runCactus, cactus_genemapChain, cactus_genemapHomolog
    """
    def __init__(self, region, options):
        #Target.__init__(self, time=18000, memory=50000000)
        Target.__init__(self, time=18000)
        self.region = region
        self.options = options
        self.seqdir = os.path.join(os.getcwd(), options.seqDir, region)
        
        self.config = options.config
        if re.search(".xml", options.config) and not os.path.isabs(options.config):
            self.config = os.path.join(os.getcwd(), options.config)

        self.species = getFirstLine(options.species)
        self.tree = getFirstLine(options.tree)
        self.dbStr = getDatabaseString(options.db, options.dbTabBase + self.region)

        self.genedir = os.path.join(options.geneDir, region)

    def run(self):
        #localTempDir = getTempFile(rootDir=self.getGlobalTempDir())
        localTempDir = self.getLocalTempDir()
        config = os.path.join(localTempDir, "cactus_workflow_config.xml")
        system("cp %s %s" %(self.config, config)) #Copy the config file to local disk
        
        #Copy sequences to localTempDir:
        localSeqdir = os.path.join(localTempDir, "data")
        system("mkdir -p %s" %localSeqdir)
        for spc in self.species.split():
            currseqdir = os.path.join(self.seqdir, spc)
            system("cp -r %s %s" %(currseqdir, localSeqdir))
            

        #Make dir for this region if not already existed
        #system("rm -fR %s" %self.region)
        system("mkdir -p %s" % os.path.join(os.getcwd(), self.region))

        #Write experiment.xml for this region:
        experimentFile = os.path.join(localTempDir, "experiment.xml")
        writeExpCommand = "cactus_writeExperimentXml.py --species \"%s\" --tree \"%s\" --output %s --sequenceDir %s --config %s --databaseString %s"\
                          %(self.species, self.tree, experimentFile, localSeqdir, config, self.dbStr)
        system("%s" %writeExpCommand )
        system("cp %s %s" %( experimentFile, os.path.join(os.getcwd(), self.region, "experiment.xml") ))
        logger.info("Got experiment.xml file for %s with command: %s\n" %(self.region, writeExpCommand))
        
        #Now ready to runCactus:
        batchSystem = "singleMachine"
        jobTree = os.path.join(localTempDir, "jobTree")
        cactusCommand = "cactus_workflow.py --batchSystem %s --experiment %s --buildReference --setupAndBuildAlignments --logDebug --jobTree %s" \
                        %(batchSystem, experimentFile, jobTree)
        logger.info("Going to run cactus now, the command is %s" %cactusCommand)
        system("%s" %cactusCommand)
        system("cp -r %s %s" %(jobTree, os.path.join(os.getcwd(), self.region, "jobTree")) )
        logger.info("Done cactusRun for %s\n" %self.region)
        
        #Run genemapChain:
        self.addChildTarget(RunGenemapChain(self.region, self.dbStr, self.options.outdir, self.options.refSpecies, self.genedir))
        self.addChildTarget(RunGenemapHomolog(self.region, self.dbStr, self.options.outdir, self.options.refSpecies, self.genedir))

class RunGenemapChain(Target):
    """
    """
    def __init__(self, region, databaseString, outputDir, refSpecies, geneDir):
        Target.__init__(self, time=120)
        self.region = region
        self.dbStr = databaseString
        self.refSpecies = refSpecies
        self.outputFile = os.path.join(outputDir, "%s-%s.xml"%("genemapChain", region) )
        
        self.geneFile = os.path.join(geneDir, "refGeneConverted.bed")

    def run(self):
        geneFile = os.path.join(self.getLocalTempDir(), "refgene.bed")
        system("cp %s %s" %(self.geneFile, geneFile))

        system("cactus_genemapChain -c %s -o \"%s\" -s \"%s\" -g \"%s\"" \
                %(self.dbStr, self.outputFile, self.refSpecies, geneFile))
        logger.info("Done genemapChain for %s\n" %self.region)

class RunGenemapHomolog(Target):
    """
    """
    def __init__(self, region, databaseString, outputDir, refSpecies, geneDir):
        Target.__init__(self, time=120)
        self.region = region
        self.dbStr = databaseString
        self.output1 = os.path.join(outputDir, "%s-%s.txt" %("genemapHomolog", region))
        self.output2 = os.path.join(outputDir, "%s-%s.xml" %("genemapHomolog", region))
        self.refSpecies = refSpecies

        self.geneFile = os.path.join(geneDir, "refGeneSorted.bed")

    def run(self):
        geneFile = os.path.join(self.getLocalTempDir(), "refgene.bed")
        system("cp %s %s" %(self.geneFile, geneFile))
        
        command = "cactus_genemapHomolog -c %s -o \"%s\" -s \"%s\" -g \"%s\" > %s" \
                  %(self.dbStr, self.output1, self.refSpecies, geneFile, self.output2)
        system("%s" %command)
        logger.info("Done genemapHomolog for %s, command: %s\n" %(self.region, command))

class AggregateResults(Target):
    def __init__(self, options):
        Target.__init__(self, time=3600)
        self.options = options
        self.output = options.outdir
        self.species = getFirstLine(options.species)

    def run(self):
	regions = getList(self.options.regions)
        genemapChainXmls = [] #list of all genemapChain output Xmls
        genemapHomologXmls = [] #list of all genemapHomology output Xmls
        for r in regions:
            genemapChainXmls.append( os.path.join(self.output, "%s-%s.xml" %("genemapChain", r)) )
            genemapHomologXmls.append( os.path.join(self.output, "%s-%s.xml" %("genemapHomolog", r)) )
            
        #Directory of more details information if interested
        extraInfoDir = os.path.join(self.output, "extraInfo")
        system("mkdir -p %s" %extraInfoDir)
        system("chmod ug+xrw %s" %extraInfoDir)

        #Merge homologXmls of all regions:
        allHomologXml = "%s/%s-all.xml" %(self.output, "genemapHomolog")
        mergeXmls(genemapHomologXmls, allHomologXml)
        
        genemapHomolog = "%s/%s-*.txt" % (self.output, "genemapHomolog")
        allHomolog = "%s/%s-all.txt" % (self.output, "genemapHomolog")
        system("rm -f %s" %allHomolog)
        system("cat %s > %s" %(genemapHomolog, allHomolog))
        
        #geneToChain = "%s/%s" %(extraInfoDir, "gene2chain")
        geneToChain = "%s/%s" %(self.output, "gene2chain")
        
        genemapChainCommand = "genemapChain.py -o %s -c \"%s\" -i \"%s\" > %s" %(extraInfoDir, "cat",\
                               " ".join(genemapChainXmls), geneToChain)
        system("%s" % genemapChainCommand)

        chainMergeHomolog = "%s/%s" %(extraInfoDir, "chainMergeHomolog")
        chainMergeHomologTex = "%s/%s" %(self.output, "chainVsDup.tex")
        #chainMergeHomologTex = chainMergeHomolog + ".tex" 
        missedGenes = "%s/%s" %(extraInfoDir, "missedGenes")
        genemapMergeCommand = "genemapMerge.py -f c -n %s %s %s %s %s > %s" %(self.options.runName, \
                               allHomolog, geneToChain, chainMergeHomolog, chainMergeHomologTex, missedGenes)
        system("%s" % genemapMergeCommand)

        homologCmp = "%s/%s" %(self.output, "homologCmp")
        homologCmpTex = "%s/%s" %(self.output, "homologCmp.tex")
        homologCmpV = "%s/%s" %(extraInfoDir, "homologCmpV")
        cactusVsMultizCommand = "genemapCactusVsMultiz.py -a %s -d %s %s %s %s > %s" %(extraInfoDir + "/perSpcDiff", \
                        self.options.geneDir + "/all.tx", self.options.multiz, allHomologXml, homologCmp, homologCmpV)
        system("%s" % cactusVsMultizCommand)

        makeLatexTabCommand = "genemapMakeLatexTab.py -s \"%s\" -n %s %s %s" \
                               %(self.species, self.options.runName, homologCmp, homologCmpTex)
        system("%s" % makeLatexTabCommand)
	
        #Cleanup now...
        self.setFollowOnTarget(Cleanup( self.output, extraInfoDir ))

class Cleanup(Target):
    def __init__(self, outdir, extraInfoDir):
        Target.__init__(self, time=0.00025)
        self.outdir = outdir
        self.extraInfoDir = extraInfoDir

    def run(self):
        system("mv %s/%s-all.xml %s" %(self.outdir, "genemapHomolog", self.extraInfoDir))
        system("rm -f %s/%s-*.xml" %(self.outdir, "genemapHomolog"))
        
        system("mv %s/%s-all.txt %s" %(self.outdir, "genemapHomolog", self.extraInfoDir))
        system("rm -f %s/%s*.txt" %(self.outdir, "genemapHomolog"))

        system("mv %s/%s-*.xml %s" %(self.outdir, "genemapChain", self.extraInfoDir))
        system("mv %s/%s %s" %(self.outdir, "gene2chain", self.extraInfoDir))


#------------ UTILITES FUNCTIONS --------------------------#
def check_dir(path):
    """Check if directories on the path, and create them if not."""
    if not os.path.exists(path):
        os.makedirs(path)

def getList(file):
    f = open(file, 'r')
    line = f.readline().strip()
    list = line.split()
    f.close()
    return list

def getFirstLine(file):
    f = open(file, 'r')
    firstline = f.readline().strip()
    f.close()
    return firstline

def getDatabaseString(databaseName, dbTableName):
    dbStr = "'<st_kv_database_conf type=\"mysql\"><mysql host=\"kolossus-10\" port=\"0\" user=\"cactus\" password=\"cactus\" database_name=\"%s\"  table_name=\"%s\" /></st_kv_database_conf>'" %(databaseName, dbTableName)
    return dbStr

def mergeXmls(infiles, outfile):
    outfh = open(outfile, "w")
    num = 0
    for file in infiles:
        infh = open(file, "r")
        firstline = infh.readline()
        if num == 0:
            outfh.write(firstline) 
        for line in infh.readlines():
            if len(line) > 2 and line[0:2] == "</":
                if num == len(infiles) -1:
                    outfh.write(line)
            else:
                outfh.write(line)
        num += 1
        infh.close()
    outfh.close()
    return


#------------ MAIN --------------------#
def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)

    parser.add_option("-c", "--cactusConfig", dest="config", help="Path to config file of cactus. Default = cactus default config.",\
                                           default="default")
    parser.add_option("-i", "--sequenceDir",dest= "seqDir", help= "Directory where the input sequences are.")
    parser.add_option("-s", "--species", dest="species", help="File containing list of input species, space separated")
    parser.add_option("-t", "--tree", dest="tree", help="File containing newick tree")
    parser.add_option("-r", "--regions", dest="regions", help="File containing a list of regions (space separated) to analyze.")
    parser.add_option("-d", "--database", dest="db", help="Mysql database name to write cactus disk to.")
    parser.add_option("-b", "--databaseTableBase", dest="dbTabBase", help="Base name of the tables to write cactus disk to. E.g, 'primDefault', cactusDisk will be written to table 'primDefaultRegionName'. Default=''", default="")
    parser.add_option("-g", "--geneDir", dest="geneDir", help="Location of the refseq directory.")
    parser.add_option("-m", "--multiz", dest="multiz", help="Multiz consEvid file")
    parser.add_option("-e", "--refSpecies", dest="refSpecies", help="The reference species. E.g: hg18")
    parser.add_option("-o", "--output", dest="outdir", help="Output directory. Default= '.'", default=".")
    parser.add_option("-n", "--runName", dest="runName", help="Name of the run (used to name latex tables)", default="")

    options, args = parser.parse_args()

    if options.seqDir == None:
        raise RuntimeError("No sequenceDir (location of input sequences) given")
    if options.species == None:
        raise RuntimeError("No species given")
    if options.tree == None:
        raise RuntimeError("No Newick tree given")
    if options.regions == None:
        raise RuntimeError("No regions given")
    if options.db == None:
        raise RuntimeError("No database given")
    if options.geneDir == None:
        raise RuntimeError("No geneDir (location of refseqs) given")
    if options.multiz == None:
        raise RuntimeError("No multiz consEvid file given")
    if options.refSpecies == None:
        raise RuntimeError("No reference species given")


    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == "__main__":
    from cactus.utilities.genemap.cactus_runGenemap import *    
    main()


