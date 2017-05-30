#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import sys
from toil.job import Job
from cactus.pipeline.ktserverControl import runKtserver, blockUntilKtserverIsRunning, stopKtserver, \
    blockUntilKtserverIsFinished

class KtServerService(Job.Service):
    def __init__(self, dbElem, isSecondary, existingSnapshotID=None,
                 memory=None, cores=None, disk=None):
        Job.Service.__init__(self, memory=memory, cores=cores, disk=disk, preemptable=False)
        self.dbElem = dbElem
        self.isSecondary = isSecondary
        self.existingSnapshotID = existingSnapshotID
        self.blockTimestep = 10
        self.blockTimeout = sys.maxint

    def start(self, job):
        snapshotExportID = job.fileStore.jobStore.getEmptyFileStoreID()
        self.dbElem, self.logPath = runKtserver(self.dbElem, fileStore=job.fileStore,
                                                existingSnapshotID=self.existingSnapshotID,
                                                snapshotExportID=snapshotExportID)
        assert self.dbElem.getDbHost() != None
        blockUntilKtserverIsRunning(self.dbElem, self.logPath, self.blockTimeout, self.blockTimestep)
        return self.dbElem.getConfString(), snapshotExportID

    def stop(self, job):
        stopKtserver(self.dbElem)
        blockUntilKtserverIsFinished(self.logPath, timeout=600)

    def check(self):
        return True
    
