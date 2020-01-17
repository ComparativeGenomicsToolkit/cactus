#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import stat
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
        self.failed = False
        self.process = None

    def start(self, job):
        snapshotExportID = job.fileStore.jobStore.getEmptyFileStoreID()
        # We need to run this garbage in case we are on a file-based
        # jobStore with caching enabled. The caching jobStore sets
        # this empty file to be unwritable for some reason. Since we
        # need to write something to it, obviously that won't do.
        path = job.fileStore.readGlobalFile(snapshotExportID)
        os.chmod(path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        self.process, self.dbElem, self.logPath = runKtserver(self.dbElem, fileStore=job.fileStore,
                                                              existingSnapshotID=self.existingSnapshotID,
                                                              snapshotExportID=snapshotExportID)
        assert self.dbElem.getDbHost() != None
        blockUntilKtserverIsRunning(self.logPath)
        self.check()
        return self.dbElem.getConfString(), snapshotExportID

    def stop(self, job):
        self.check()
        try:
            stopKtserver(self.dbElem)
        except:
            # Server is probably already terminated
            pass
        if not self.failed:
            blockUntilKtserverIsFinished(self.logPath, timeout=1200)

    def check(self):
        if self.process.exceptionMsg.empty():
            return True
        else:
            self.failed = True
            msg = self.process.exceptionMsg.get()
            raise RuntimeError(msg)
