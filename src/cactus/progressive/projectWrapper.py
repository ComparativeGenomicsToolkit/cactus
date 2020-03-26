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
# - initialize the working directory
# - create Experiment file from seqfile and options
# - create Config file from options
# - run cactus_createMultiCactusProject
# - now ready to launch cactus progressive
class ProjectWrapper:
    alignmentDirName = 'progressiveAlignment'
    def __init__(self, options, configPath, ignoreSeqPaths=[]):
        self.options = options
        self.seqFile = SeqFile(options.seqFile)
        self.workingDir = options.cactusDir
        self.configWrapper = None
        self.configPath = configPath
        self.expWrapper = None
        self.processConfig()
        self.processExperiment(ignoreSeqPaths)

    def processConfig(self):
        log.info("Using config from path %s." % self.configPath)
        configXml = ET.parse(self.configPath).getroot()
        self.configWrapper = ConfigWrapper(configXml)
        # here we can go through the options and apply some to the config
        self.configWrapper.setBuildHal(True)
        self.configWrapper.setBuildFasta(True)

    def processExperiment(self, ignoreSeqPaths):
        expXml = self.seqFile.toXMLElement(ignoreSeqPaths)
        #create the cactus disk
        cdElem = ET.SubElement(expXml, "cactus_disk")
        database = self.options.database
        assert database == "kyoto_tycoon"
        confElem = ET.SubElement(cdElem, "st_kv_database_conf")
        confElem.attrib["type"] = database
        ET.SubElement(confElem, database)
        self.expWrapper = ExperimentWrapper(expXml)
        self.expWrapper.setConfigPath(self.configPath)
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
