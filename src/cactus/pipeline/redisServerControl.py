#!/usr/bin/env python3
"""
Functions to launch and manage Redis servers.
"""

import os
import random
import signal
import sys
import traceback
from multiprocessing import Process, Queue
from time import sleep

from toil.lib.bioio import logger
from cactus.shared.common import cactus_call
from cactus.pipeline.dbServerCommon import getHostName, findOccupiedPorts


MAX_REDIS_PORT = 65535

# The name of the snapshot that Redis outputs.
REDIS_SNAPSHOT_NAME = "dump.rdb"


class RedisServer:
    def __init__(self, dbElem, fileStore=None, existingSnapshotID=None, snapshotExportID=None):
        self.dbElem = dbElem
        self.logPath = None
        self.fileStore = fileStore
        self.existingSnapshotID = existingSnapshotID
        self.snapshotExportID = snapshotExportID
        self.databaseDir = None

    def runServer(self):
        """
        Run a redis-server. This function launches a separate python process that manages the server.

        Writing to the special key "TERMINATE" signals this thread to safely shut
        down the DB and save the results. After finishing, the data will
        eventually be written to snapshotFile.

        Returns a tuple containing an updated version of the database config dbElem and the
        path to the log file.
        """

        self.databaseDir = self.fileStore.getLocalTempDir()
        # log file can be saved in a subdirectory of where the snapshot is being saved
        self.logPath = os.path.join(self.databaseDir, "redis.log")
        open(self.logPath, 'a').close()

        self.dbElem.setDbHost(getHostName())
        # Find a suitable port to run on.
        try:
            occupiedPorts = findOccupiedPorts()
            unoccupiedPorts = set(range(1025, MAX_REDIS_PORT)) - occupiedPorts
            port = random.choice(list(unoccupiedPorts))
        except:
            logger.warning("Can't find which ports are occupied--likely netstat is not installed."
                           " Choosing a random port to start the DB on, good luck!")
            port = random.randint(1025, MAX_REDIS_PORT)
        self.dbElem.setDbPort(port)
        process = RedisServerProcess(self)
        process.daemon = True
        process.start()

        if not self.blockUntilServerIsRunning():
            try:
                with open(self.logPath) as f:
                    log = f.read()
            except:
                log = ''
            raise RuntimeError("Unable to launch redis-server in time. Log: %s" % log)

        return process, self.dbElem, self.logPath

    def blockUntilServerIsRunning(self, createTimeout=1800):
        """Check status until it's successful, an error is found, or we timeout.

        Returns True if the redis-server is now running, False if something went wrong."""
        success = False
        for i in range(createTimeout):
            if self.isServerFailed():
                logger.critical('Error starting Redis server.')
                success = False
                break
            if self.isServerRunning():
                logger.info('Redis server running.')
                success = True
                break
            sleep(1)
        return success

    def blockUntilServerIsFinished(self, timeout=1800, timeStep=10):
        """Wait for the redis-server log to indicate that it shut down properly.

        Returns True if the server shut down, False if the timeout expired."""
        for i in range(0, timeout, timeStep):
            with open(self.logPath) as f:
                log = f.read()
                if 'ready to exit' in log:
                    return True
            sleep(timeStep)
        raise RuntimeError("Timeout reached while waiting for redis server.")

    def isServerRunning(self):
        """Check if the server started running."""
        success = False
        with open(self.logPath) as f:
            for line in f:
                if line.lower().find("accept connections") >= 0:
                    success = True
        return success

    def isServerFailed(self):
        """Does the server log contain an error?"""
        isFailed = False
        with open(self.logPath) as f:
            for line in f:
                if line.lower().find("error") >= 0:
                    isFailed = True
                    break
        return isFailed

    def getTuningOptions(self):
        """Get the appropriate redis-server tuning parameters (bucket size, etc.)"""
        # these are some hardcoded defaults.  should think about moving to config
        # TODO: check if every necessary tuning option is added (maybe to be merged with getServerOptions())
        tuningOptions = "--maxmemory 50Gb"
        # override default redis-server settings if they are present in the
        # experiment xml file.
        if self.dbElem.getDbTuningOptions() is not None:
            tuningOptions = self.dbElem.getDbTuningOptions()
        if self.dbElem.getDbCreateTuningOptions() is not None:
            tuningOptions = self.dbElem.getDbCreateTuningOptions()
        return tuningOptions

    def getServerOptions(self):
        # these are some hardcoded defaults.  should think about moving to config
        # TODO: check if every necessary option is added
        serverOptions = "--timeout 0 --databases 1 --protected-mode no"
        if self.dbElem.getDbServerOptions() is not None:
            serverOptions = self.dbElem.getDbServerOptions()
        return serverOptions

    def getServerCommand(self, snapshotDir):
        """Get a redis-server command line with the proper options (in popen-type list format)."""
        serverOptions = self.getServerOptions()
        tuning = self.getTuningOptions()
        cmd = ["redis-server", "--port", str(self.dbElem.getDbPort())]
        cmd += serverOptions.split()
        # Configure background snapshots, but set the interval between
        # snapshots to ~ 10 days so it'll never trigger. We are only
        # interested in the snapshot that the DB creates on termination.
        cmd += ["--save", "1000000", "1"]
        cmd += ["--dir", snapshotDir, "--dbfilename", REDIS_SNAPSHOT_NAME]
        cmd += ["--logfile", 'redis.log']
        cmd += tuning.split()
        return cmd

    def getRemoteParams(self):
        """Get parameters to supply to redis-cli to connect to the right DB."""
        # TODO: find why redis-cli raise error when called with an IP as the host
        #host = self.dbElem.getDbHost() or 'localhost'
        host = 'localhost'
        return ['-p', str(self.dbElem.getDbPort()), '-h', host]

    def stopServer(self):
        """Attempt to send the terminate signal to a redis-server."""
        cactus_call(parameters=['redis-cli'] + self.getRemoteParams() + ['set', 'TERMINATE', '1'])


class RedisServerProcess(Process):
    """Independent process that babysits the redis-server process.

    Waits for the TERMINATE flag to be set, then kills the DB and
    copies the final snapshot to snapshotExportID.
    """
    exceptionMsg = Queue()

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        super(RedisServerProcess, self).__init__()

    def run(self):
        """Run the tryRun method, signaling the main thread if an exception occurs."""
        try:
            self.tryRun(*self.args, **self.kwargs)
        except BaseException:
            self.exceptionMsg.put("".join(traceback.format_exception(*sys.exc_info())))
            raise

    def tryRun(self, redisServer):
        snapshotDir = redisServer.databaseDir
        snapshotPath = os.path.join(snapshotDir, REDIS_SNAPSHOT_NAME)
        if redisServer.existingSnapshotID is not None:
            # Extract the existing snapshot to the snapshot
            # directory so it will be automatically loaded
            redisServer.fileStore.readGlobalFile(redisServer.existingSnapshotID, userPath=snapshotPath)
        process = cactus_call(server=True, shell=False,
                              parameters=redisServer.getServerCommand(snapshotDir),
                              port=redisServer.dbElem.getDbPort())

        redisServer.blockUntilServerIsRunning()
        if redisServer.existingSnapshotID is not None:
            # Clear the termination flag from the snapshot
            cactus_call(parameters=["redis-cli"] + redisServer.getRemoteParams() + ["del", "TERMINATE"])
        while True:
            # Check for the termination signal
            terminateFlag = cactus_call(parameters=["redis-cli"] + redisServer.getRemoteParams() + ["get", "TERMINATE"],
                                        swallowStdErr=True, check_output=True)
            if terminateFlag.strip() != '1':
                # No terminate signal sent yet
                pass
            else:
                # Terminate signal received
                break
            # Check that the DB is still alive
            if process.poll() is not None or redisServer.isServerFailed():
                with open(redisServer.logPath) as f:
                    raise RuntimeError("redis server failed. Log: %s" % f.read())
            sleep(60)
        process.send_signal(signal.SIGINT)
        process.wait()
        redisServer.blockUntilServerIsFinished()
        if redisServer.snapshotExportID is not None:
            if not os.path.exists(snapshotPath):
                with open(redisServer.logPath) as f:
                    raise RuntimeError("redis-server did not leave a snapshot on termination,"
                                       " but a snapshot was requested. Log: %s" % f.read())
            # Export the snapshot file to the file store
            redisServer.fileStore.jobStore.updateFile(redisServer.snapshotExportID, snapshotPath)