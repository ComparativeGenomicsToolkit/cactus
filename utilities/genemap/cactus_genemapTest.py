#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys, os, re, time

from sonLib.bioio import system

class TestCase(unittest.TestCase):
    runsDir = "/hive/users/nknguyen/mhc/data/ucsc/runs"

    testruns = {"HAR1_primates":"hg19.chr20.63025520.61636691.212254.1", 
                "HCG27-MICB-hc":"hg19.chr6.171115067.31165033.315710.1", 
		"HCG27-MICB-hcmr":"hg19.chr6.171115067.31165033.315710.1", 
		"egfr":"hg19.chr7"}

    refGeneFile = "refGene.bed"

    def test_cactus_geneMap(self):
        for run in self.testruns:
           currdir = os.getcwd()
	   speciesName = self.testruns[run]
	   dir = os.path.join(self.runsDir, run)
	   cactusDisk = os.path.join(dir, "cactusDisk")
           #if not os.path.exists(cactusDisk):
           #    system("cd %s" % dir)
           #    system("make runCactus")
           #    system("cd %s" % currdir)

           databaseString = "'<st_kv_database_conf type=\"tokyo_cabinet\"><tokyo_cabinet database_dir=\"%s\" /></st_kv_database_conf>'" % (cactusDisk)

	   refGeneFilePath = os.path.join(dir, self.refGeneFile)
           #if not os.path.exists(refGeneFilePath):
           #    system("cd %s" % dir)
           #    system("make getRefseq")
           #    system("cd %s" % currdir)

	   #geneMap = os.path.join(dir, "%s_geneMap" % speciesName)
	   geneMap = os.path.join(dir, "geneMap")
	   geneMapXml = os.path.join(dir, "%s.xml" % geneMap)
	   geneMapTex = os.path.join(dir, "%s.tex" % geneMap)

	   system("rm -f %s*" % geneMap)
	   system("cactus_geneMap -c %s -o \"%s\" -s \"%s\" -g \"%s\"" % (databaseString, geneMapXml, speciesName, refGeneFilePath))
	   system("geneMap.py %s %s" %(geneMapXml, geneMapTex))   
	   #system("latex %s" % geneMap)
	   #system("latex %s" % geneMap)
	   #system("latex %s" % geneMap)
           #system("chmod ug+xrw %s*" % geneMap)
	   #system("dvips -Ppdf %s.dvi" % geneMap)
	   #system("ps2pdf %s.ps" %geneMap)
	   #system("rm %s.aux %s.dvi %s.log %s.ps" % (geneMap, geneMap, geneMap, geneMap))

        
def main():
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
