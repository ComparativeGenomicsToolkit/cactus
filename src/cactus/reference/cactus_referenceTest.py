#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
import argparse
import xml.etree.ElementTree as ET

from toil.common import Toil
from toil.job import Job

from sonLib.bioio import getLogLevelString
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory

from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import silentOnSuccess

from cactus.shared.common import cactusRootPath
from cactus.shared.common import cactus_call

from cactus.pipeline.cactus_workflow import CactusReferencePhase, CactusCheckpointJob, CactusWorkflowArguments

class CactusReferenceTestCheckpoint(CactusCheckpointJob):
    def run(self, fileStore):
        return self.runPhaseWithPrimaryDB(CactusReferencePhase).rv()

class TestCase(unittest.TestCase):
    def setUp(self):
        self.dbDump = os.path.join(cactusRootPath(), "testdata", "blanchette-caf.dbdump")
        self.experimentXML = os.path.join(cactusRootPath(), "testdata", "blanchette_experiment.xml")

    @unittest.skip('')
    @silentOnSuccess
    def testCactus_Blanchette_Greedy(self):
        self.runUsingSavedDB("greedy")

    @unittest.skip('')
    @silentOnSuccess
    def testCactus_Blanchette_Blossum(self):
        self.runUsingSavedDB("blossom5")

    def testCactus_Blanchette_MaxCardinality(self):
        self.runUsingSavedDB("maxCardinality")

    @unittest.skip('')
    @silentOnSuccess
    def testCactus_Blanchette_MaxWeight(self):
        self.runUsingSavedDB("maxCardinality")

    def testCuTest(self):
        cactus_call(parameters=["referenceTests", getLogLevelString()])

    def runUsingSavedDB(self, matchingAlgorithm):
        configFile = getConfigFile(matchingAlgorithm)
        configNode = ET.parse(configFile).getroot()
        fakeOptions = argparse.Namespace()
        fakeOptions.buildHal = True
        fakeOptions.buildFasta = True
        fakeOptions.buildReference = True
        fakeOptions.buildAvgs = False
        fakeOptions.dumpDBsToURL = False
        fakeOptions.intermediateResultsUrl = None
        cactusWorkflowArguments = CactusWorkflowArguments(fakeOptions, experimentFile=self.experimentXML,
                                                          configNode=configNode, seqIDMap=None)
        cactusWorkflowArguments.totalSequenceSize = 1
        job = CactusReferenceTestCheckpoint(ktServerDump=self.dbDump, cactusWorkflowArguments=cactusWorkflowArguments, phaseName="reference")
        toilOpts = Job.Runner.getDefaultOptions(os.path.join(getTempDirectory(), "jobStore"))
        with Toil(toilOpts) as toil:
            experiment = toil.start(job)
            import sys
            print ET.ElementTree(experiment.xmlRoot).write(sys.stdout)
            toil.exportFile(experiment.getReferenceID(), "file:///home/joel/reference.fa")
        os.remove(configFile)

def generateBlanchetteDump(resultsUrl):
    configFile = getConfigFile('greedy')
    runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                 buildReference=True,
                                 configFile=configFile,
                                 intermediateResultsUrl=resultsUrl)
    os.remove(configFile)

def getConfigFile(matchingAlgorithm="greedy", referenceEventName="root"):
    tempConfigFile = getTempFile(rootDir="./", suffix=".xml")
    config = ET.parse(os.path.join(cactusRootPath(), "cactus_progressive_config.xml")).getroot()
    #Set the matching algorithm
    config.find("reference").attrib["matching_algorithm"] = matchingAlgorithm

    config.find("reference").attrib["reference"] = referenceEventName

    #Now print the file..
    fileHandle = open(tempConfigFile, 'w')
    ET.ElementTree(config).write(fileHandle)
    fileHandle.close()
    return os.path.abspath(tempConfigFile)

if __name__ == '__main__':
    unittest.main()
