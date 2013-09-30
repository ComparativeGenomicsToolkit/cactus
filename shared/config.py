#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import xml.etree.ElementTree as ET
import os
import sys

from sonLib.bioio import system
from sonLib.bioio import getRandomAlphaNumericString

class CactusWorkflowExperiment:
    """Object used for generating cactus workflow experiment config files,
    using the inputs to generate valid config strings and files.
    """
    def __init__(self, sequences, newickTreeString, 
                 outgroupEvents=None, outputDir=None, databaseName=None, 
                 databaseConf=None, configFile=None, 
                 halFile=None, fastaFile=None,
                 constraints=None, progressive=False, 
                 outputSequences=None):
        self.experiment = ET.Element("cactus_workflow_experiment")
        if databaseName == None:
            self.databaseName = "cactusDisk_%s_%i" % (getRandomAlphaNumericString(), os.getpid()) #Needs to be unique
        else:
            self.databaseName = databaseName
        self.databaseFile = None
        #Do the database first
        database = ET.SubElement(self.experiment, "cactus_disk")
        self.mysql = 0
        self.kyotoTycoon = 0
        if databaseConf != None:
            checkDatabaseConf(databaseConf)
            databaseConf = ET.fromstring(ET.tostring(databaseConf))
            checkDatabaseConf(databaseConf) #This is just a redundant check
            database.append(databaseConf)
            if databaseConf.attrib["type"] == "tokyo_cabinet":
                tokyoCabinet = databaseConf.find("tokyo_cabinet")
                tokyoCabinet.attrib["database_dir"] = os.path.join(tokyoCabinet.attrib["database_dir"], self.databaseName)
                self.databaseFile = tokyoCabinet.attrib["database_dir"]
                assert not os.path.exists(self.databaseFile)
            else:
                self.kyotoTycoon = 1
                assert databaseConf.attrib["type"] == "kyoto_tycoon"
                kyotoTycoon = databaseConf.find("kyoto_tycoon")
                kyotoTycoon.attrib["database_dir"] = os.path.join(kyotoTycoon.attrib["database_dir"], self.databaseName)
                self.databaseFile = kyotoTycoon.attrib["database_dir"]
                assert not os.path.exists(self.databaseFile)
        else:
            databaseConf = ET.SubElement(database, "st_kv_database_conf")
            databaseConf.attrib["type"] = "tokyo_cabinet"
            tokyoCabinet = ET.SubElement(databaseConf, "tokyo_cabinet")
            assert outputDir != None
            assert os.path.exists(outputDir) #Check it exists
            self.databaseFile = os.path.join(outputDir, self.databaseName)
            assert not os.path.exists(self.databaseFile)
            tokyoCabinet.attrib["database_dir"] = self.databaseFile
        #Now add in the user stuff..
        self.experiment.attrib["sequences"] = " ".join(sequences)
        if outputSequences == None:
            assert outputDir != None
            outputSequences = os.path.join(outputDir, "outputSequences")
        self.experiment.attrib["outputSequences"] = outputSequences
        self.experiment.attrib["species_tree"] = newickTreeString
        if outgroupEvents != None:
            self.experiment.attrib["outgroup_events"] = outgroupEvents
        if progressive == True:
            self.experiment.attrib["config"] = "defaultProgressive"
        else:
            self.experiment.attrib["config"] = "default"
        if configFile != None:
            self.experiment.attrib["config"] = configFile
        if halFile != None or mafFile != None or fastaFile != None:
            halElem = ET.SubElement(self.experiment, "hal")
            if halFile != None:
                halElem.attrib["halPath"] = halFile
            if mafFile != None:
                halElem.attrib["mafPath"] = mafFile
            if fastaFile != None:
                halElem.attrib["fastaPath"] = fastaFile
        #Constraints
        if constraints != None:
            if not os.path.exists(constraints):
                raise RuntimeError("Constraints file does not appear to exist: %s" % constraints)
            self.experiment.attrib["constraints"] = constraints
