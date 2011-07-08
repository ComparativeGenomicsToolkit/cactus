#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Test the progressive workflow.  Not very thorough and does not use the giant existing cactus 

test framework.  It does use the sonlib data so SON_TRACE_DATASETS is required. 

"""
import unittest

import os
import sys
import random
import math
import copy
import filecmp

from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import TempFileTree
from sonLib.bioio import getRandomAlphaNumericString
from sonLib.bioio import printBinaryTree
from sonLib.bioio import newickTreeParser
from sonLib.bioio import system


class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        self.host = "localhost"
        self.port = "2645"
        self.normalTokyoDB = self.tempDir + "/NORMALTOKYODB"
        self.progressiveTokyoDB = self.tempDir + "/PROGRESSIVETOKYODB"
        self.progressiveKyotoDB = self.tempDir + "/PROGRESSIVETOKYDB"
        self.jt = self.tempDir + "/JT"
        self.exp = self.tempDir + "/EXP.xml"
        self.ref = self.tempDir + "/REF.fa"
        self.tree = "(((HUMAN:0.006969,CHIMP:0.009727):0.025291,BABOON:0.044568):0.11,(MOUSE:0.072818,RAT:0.081244):0.260342);"
        self.mrtree = "(MOUSE:0.072818, RAT:0.081244):0.260342;"
        unittest.TestCase.setUp(self)
            
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        for tempFile in self.tempFiles:
            os.remove(tempFile)
    
    def getBlanchetteSequencePath(self, name, number = "00"):
        seqPath = os.getenv("SON_TRACE_DATASETS") + "/"
        seqPath += "blanchettesSimulation" + "/" + number + ".job/" + name
        return seqPath 
    
    def getTokyoDBString(self, progressive=False):
        dbString = "<st_kv_database_conf type=\"tokyo_cabinet\"> <tokyo_cabinet database_dir="
        if progressive:
            dbString += "\""  + self.progressiveTokyoDB + "\""
        else:
            dbString += "\""  + self.normalTokyoDB + "\""
        dbString += " /> </st_kv_database_conf>"
        return dbString
    
    def getKyotoDBString(self):
        dbstring = "<st_kv_database_conf type=\"kyoto_tycoon\"> <kyoto_tycoon database_dir="
        dbString += "\"" + self.progressiveKyotoDB + "\""
        dbstring += " host=" + "\"" + self.host + "/"
        dbString += " port=" + "\"" + self.port + "/"
        dbString += " /> </st_kv_database_conf>"
        return dbString
        
    def getExperiment(self, sequences, tree, progressive=False, kyoto=False):
        dbstring = ""
        if kyoto:
            dbstring = self.getKyotoDBString()
        else:
            dbstring = self.getTokyoDBString(progressive)    
        expString = "<cactus_workflow_experiment config=\"default\" sequences="
        expString += "\"" + sequences + "\" species_tree=\"" + tree +"\"> <cactus_disk> " + dbstring
        expString += "</cactus_disk> </cactus_workflow_experiment>"
        return expString
    
    def getCmdLine(self, progressive):
        name = ""
        if progressive:
            name = "progressive"
        else:
            name = "workflow"
        cmdString = "cactus_" + name + ".py"
        cmdString += " --experiment " + self.exp
        cmdString += " --setupAndBuildAlignments --buildReference"
        cmdString += " --jobTree " + self.jt
        return cmdString
    
    def runExperiment(self, sequences, tree, progressive, kyoto):
        expString = self.getExperiment(sequences, tree, progressive, kyoto)        
        expFile = open(self.exp, "w")
        expFile.write(expString)
        expFile.close()
        os.system("rm -rf " + self.jt)
        assert os.system(self.getCmdLine(progressive)) == 0
    
    def getReference(self):
        cmdLine = "cactus_getReferenceSeq --flowerName 0"
        cmdLine += " --cactusDisk " + "'" + self.getTokyoDBString(False) + "'"
        cmdLine += " --referenceEventString reference"
        cmdLine += " --outputFile " + self.ref
        assert os.system(cmdLine) == 0
    
    def testProgressiveTokyo(self):
        sequences = self.getBlanchetteSequencePath("MOUSE") + " " + self.getBlanchetteSequencePath("RAT")
        os.system("rm -rf " + self.normalTokyoDB)
        self.runExperiment(sequences, self.mrtree, False, False)
        self.getReference()
         
        os.system("rm -rf " + self.progressiveTokyoDB)
        sequences = " ".join(map(lambda x: self.getBlanchetteSequencePath(x),
                           ["HUMAN", "CHIMP", "BABOON", "MOUSE", "RAT"]))
        self.runExperiment(sequences, self.tree, True, False)
    
        return filecmp.cmp(self.ref, self.progressiveTokyoDB + "/Anc2/Anc2_reference.fa")
 
 # too lazy to write code to automatically set up server and verify file structure in automated test right now   
    '''def testProgressiveKyoto(self):
        sequences = self.getBlanchetteSequencePath("MOUSE") + " " + self.getBlanchetteSequencePath("RAT")
        os.system("rm -rf " + self.normalTokyoDB)
        self.runExperiment(sequences, self.mrtree, False, False)
        self.getReference()
        
        os.system("rm -rf " + self.progressiveKyotoDB)
        sequences = " ".join(map(lambda x: self.getBlanchetteSequencePath(x),
                           ["HUMAN", "CHIMP", "BABOON", "MOUSE", "RAT"]))
        self.runExperiment(sequences, self.tree, True, True)
    
        return filecmp.cmp(self.ref, self.progressiveKyotoDB + "/Anc2/Anc2_reference.fa")'''
        
    
if __name__ == '__main__':
    unittest.main()
        
    