#!/usr/bin/env python3
import logging
import os
import xml.etree.ElementTree as ET

from sonLib.bioio import absSymPath

from .seqFile import SeqFile
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import cactusRootPath

from cactus.progressive.cactus_createMultiCactusProject import runCreateMultiCactusProject

log = logging.getLogger(__name__)

# Wrap up the cactus_progressive interface:
# - intialize the working directory
# - create Experiment file from seqfile and options
# - create Config file from options
# - run cactus_createMultiCactusProject
# - now ready to launch cactus progressive
class ProjectWrapper:
    alignmentDirName = 'progressiveAlignment'
    def __init__(self, options):
        self.options = options
        self.seqFile = SeqFile(options.seqFile)
        self.workingDir = options.cactusDir
        self.configWrapper = None
        self.expWrapper = None
        self.processConfig()
        self.processExperiment()

    def processConfig(self):
        # read in the default right out of cactus
        if self.options.configFile is not None:
            configPath = self.options.configFile
        else:
            dir = cactusRootPath()
            configPath = os.path.join(dir,
                                      "cactus_progressive_config.xml")
        log.info("Using config from path %s." % configPath)
        configXml = ET.parse(configPath).getroot()
        self.configWrapper = ConfigWrapper(configXml)
        # here we can go through the options and apply some to the config
        self.configWrapper.setBuildHal(True)
        self.configWrapper.setBuildFasta(True)

    def processExperiment(self):
        expXml = self.seqFile.toXMLElement()
        #create the cactus disk
        cdElem = ET.SubElement(expXml, "cactus_disk")
        database = self.options.database
        assert database == "kyoto_tycoon" or database == "tokyo_cabinet"
        confElem = ET.SubElement(cdElem, "st_kv_database_conf")
        confElem.attrib["type"] = database
        ET.SubElement(confElem, database)
        self.expWrapper = ExperimentWrapper(expXml)
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)

    def writeXml(self):
        assert os.path.isdir(self.workingDir)
        configPath = absSymPath(
            os.path.join(self.workingDir, "config.xml"))
        expPath = absSymPath(
            os.path.join(self.workingDir, "expTemplate.xml"))
        self.expWrapper.setConfigPath(configPath)
        self.configWrapper.writeXML(configPath)
        self.expWrapper.writeXML(expPath)

        projPath = os.path.join(self.workingDir,
                                ProjectWrapper.alignmentDirName)
        if len(self.seqFile.outgroups) == 0:
            # No outgroups specified, assume the default outgroup set
            outgroups = None
        else:
            outgroups = self.seqFile.outgroups
        runCreateMultiCactusProject(expPath, projPath,
                                    outgroupNames=outgroups,
                                    root=self.options.root)
