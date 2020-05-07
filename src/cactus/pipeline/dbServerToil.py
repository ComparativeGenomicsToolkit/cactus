#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

import os
import stat
from toil.job import Job
from cactus.pipeline.ktserverControl import KtServer
from cactus.pipeline.redisServerControl import RedisServer

class DbServerService(Job.Service):
    def __init__(self, dbElem, isSecondary, existingSnapshotID=None,
                 memory=None, cores=None, disk=None):
        Job.Service.__init__(self, memory=memory, cores=cores, disk=disk, preemptable=False)
        self.dbElem = dbElem
        self.isSecondary = isSecondary
        self.existingSnapshotID = existingSnapshotID
        self.failed = False
        self.process = None
        self.dbServer = None

    def start(self, job):
        snapshotExportID = job.fileStore.jobStore.getEmptyFileStoreID()
        # We need to run this garbage in case we are on a file-based
        # jobStore with caching enabled. The caching jobStore sets
        # this empty file to be unwritable for some reason. Since we
        # need to write something to it, obviously that won't do.
        path = job.fileStore.readGlobalFile(snapshotExportID)
        os.chmod(path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        self.dbServer = getDbServer(self.dbElem, fileStore=job.fileStore,
                                    existingSnapshotID=self.existingSnapshotID,
                                    snapshotExportID=snapshotExportID)
        self.process, self.dbElem, self.logPath = self.dbServer.runServer()
        assert self.dbElem.getDbHost() != None
        self.dbServer.blockUntilServerIsRunning()
        self.check()
        return self.dbElem.getConfString(), snapshotExportID

    def stop(self, job):
        self.check()
        try:
            self.dbServer.stopServer()
        except:
            # Server is probably already terminated
            pass
        if not self.failed:
            self.dbServer.blockUntilServerIsFinished()

    def check(self):
        if self.process.exceptionMsg.empty():
            return True
        else:
            self.failed = True
            msg = self.process.exceptionMsg.get()
            raise RuntimeError(msg)


def getDbServer(dbElem, fileStore=None, existingSnapshotID=None, snapshotExportID=None):
    """
    Get the correct object that handles the database server based on database type
    Each database has a specific class with common functionalities
    """
    if dbElem.getDbType() == "kyoto_tycoon":
        return KtServer(dbElem, fileStore, existingSnapshotID, snapshotExportID)
    elif dbElem.getDbType() == "redis":
        return RedisServer(dbElem, fileStore, existingSnapshotID, snapshotExportID)
    raise RuntimeError("The database type, %s, is not supported" % dbElem.getDbType())