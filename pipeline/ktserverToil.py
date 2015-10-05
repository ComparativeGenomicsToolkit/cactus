#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
import random
import math
import xml.etree.ElementTree as ET

from sonLib.bioio import getTempFile
from toil.job import Job
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.pipeline.ktserverControl import runKtserver
from cactus.pipeline.ktserverControl import blockUntilKtserverIsRunnning
from cactus.pipeline.ktserverControl import killKtServer
from cactus.pipeline.ktserverControl import getKtServerReport
from cactus.pipeline.ktserverControl import ktServerAlreadyRunning
from cactus.pipeline.ktserverControl import getHostName
from cactus.pipeline.ktserverControl import getLogPath




class ChildWithKtServer(Job):
    """This job launches ktserver as a Toil service, and runs
    a child Cactus workflow job (either a CactusPhasesJob or a 
    CactusRecursionJob). The ktserver will be available for the child
    job and all of its successors, but will be terminated before the follow-on
    jobs of rootJob are run."""
    def __init__(self, rootJob, newChild, isSecondary):
        Job.__init__(self)
        self.rootJob = rootJob
        self.newChild = newChild
        self.isSecondary = isSecondary
    def run(self, fileStore):
        from cactus.pipeline.cactus_workflow import CactusPhasesJob
        from cactus.pipeline.cactus_workflow import CactusRecursionJob
        dbConfString = self.addService(KtServerService(self.rootJob, self.newChild, self.isSecondary))

        #Tell the child job what port and hostname to use for connecting
        #to the database.
        if isinstance(self.newChild, CactusPhasesJob):
            self.newChild.cactusWorkflowArguments.cactusDiskDatabaseString = dbConfString
        elif isinstance(self.newChild, CactusRecursionJob):
            self.newChild.phaseNode.attrib["secondaryDatabaseString"] = dbConfString
        self.addChild(self.newChild)

class KtServerService(Job.Service):

    def __init__(self, rootJob, newChild, isSecondary):
        Job.Service.__init__(self)
        self.rootJob = rootJob
        self.newChild = newChild
        self.isSecondary = isSecondary
        self.blockTimestep = 10
        self.blockTimeout = sys.maxint
        self.killSwitchPath = None
        self.process = None



    def start(self):
        from cactus.pipeline.cactus_workflow import CactusPhasesJob
        from cactus.pipeline.cactus_workflow import CactusRecursionJob

        if self.isSecondary == False:
            assert isinstance(self.newChild, CactusPhasesJob)
            wfArgs = self.newChild.cactusWorkflowArguments
            self.dbElem = ExperimentWrapper(wfArgs.experimentNode)
            experiment = self.dbElem
        else:
            assert isinstance(self.newChild, CactusRecursionJob)
            dbString = self.newChild.getOptionalPhaseAttrib("secondaryDatabaseString")
            assert dbString is not None
            confXML = ET.fromstring(dbString)
            self.dbElem = DbElemWrapper(confXML)

        if not os.path.exists(self.dbElem.getDbDir()):
            os.mkdir(self.dbElem.getDbDir())
        self.killSwitchPath = getTempFile(suffix="_kill.txt",
                                          rootDir=self.dbElem.getDbDir())
        killSwitchFile = open(self.killSwitchPath, "w")
        killSwitchFile.write("init")
        killSwitchFile.close()


        self.process = runKtserver(self.dbElem, self.killSwitchPath)
        assert self.dbElem.getDbHost() != None
        
        if self.isSecondary == False:
            experiment.writeXML(wfArgs.experimentFile)
            wfArgs.cactusDiskDatabaseString = self.dbElem.getConfString()
        else:
            self.newChild.phaseNode.attrib[
                "secondaryDatabaseString"] = self.dbElem.getConfString()
            # added on as a hack to get this into the experiment.xml
            etPath = self.newChild.phaseNode.attrib[
                "experimentPath"]
            experiment = ExperimentWrapper(ET.parse(etPath).getroot())
            experiment.setSecondaryDBElem(self.dbElem)
            experiment.writeXML(etPath)            
        blockUntilKtserverIsRunnning(self.dbElem, self.killSwitchPath, self.blockTimeout, self.blockTimestep)
        return self.dbElem.getConfString()
        

    def stop(self):
        if self.process:
            self.process.kill()
        logPath = getLogPath(self.dbElem)
        if os.path.exists(logPath):
            os.remove(logPath)
        if self.killSwitchPath:
            os.remove(self.killSwitchPath)
    
