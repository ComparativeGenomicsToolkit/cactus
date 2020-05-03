#!/usr/bin/env python3
"""
Functions to launch and manage KyotoTycoon servers.
"""

import os
import random
import signal
import sys
import traceback
from glob import glob
from multiprocessing import Process, Queue
from time import sleep

from toil.lib.bioio import logger
from cactus.shared.common import cactus_call
from cactus.pipeline.dbServerCommon import getHostName, findOccupiedPorts

# For some reason ktserver believes there are only 32768 TCP ports.
MAX_KTSERVER_PORT = 32767

# The name of the snapshot that KT outputs.
KTSERVER_SNAPSHOT_NAME = "00000000.ktss"

def runKtserver(dbElem, fileStore, existingSnapshotID=None, snapshotExportID=None):
    """
    Run a KTServer. This function launches a separate python process that manages the server.

    Writing to the special key "TERMINATE" signals this thread to safely shut
    down the DB and save the results. After finishing, the data will
    eventually be written to snapshotFile.

    Returns a tuple containing an updated version of the database config dbElem and the
    path to the log file.
    """
    logPath = fileStore.getLocalTempFile()
    dbElem.setDbHost(getHostName())

    # Find a suitable port to run on.
    try:
        occupiedPorts = findOccupiedPorts()
        unoccupiedPorts = set(range(1025,MAX_KTSERVER_PORT)) - occupiedPorts
        port = random.choice(list(unoccupiedPorts))
    except:
        logger.warning("Can't find which ports are occupied--likely netstat is not installed."
                       " Choosing a random port to start the DB on, good luck!")
        port = random.randint(1025,MAX_KTSERVER_PORT)
    dbElem.setDbPort(port)

    process = ServerProcess(dbElem, logPath, fileStore, existingSnapshotID, snapshotExportID)
    process.daemon = True
    process.start()

    if not blockUntilKtserverIsRunning(logPath):
        try:
            with open(logPath) as f:
                log = f.read()
        except:
            log = ''
        raise RuntimeError("Unable to launch ktserver in time. Log: %s" % log)

    return process, dbElem, logPath

class ServerProcess(Process):
    """Independent process that babysits the ktserver process.

    Waits for the TERMINATE flag to be set, then kills the DB and
    copies the final snapshot to snapshotExportID.
    """
    exceptionMsg = Queue()

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        super(ServerProcess, self).__init__()

    def run(self):
        """Run the tryRun method, signaling the main thread if an exception occurs."""
        try:
            self.tryRun(*self.args, **self.kwargs)
        except BaseException:
            self.exceptionMsg.put("".join(traceback.format_exception(*sys.exc_info())))
            raise

    def tryRun(self, dbElem, logPath, fileStore, existingSnapshotID=None, snapshotExportID=None):
        snapshotDir = os.path.join(fileStore.getLocalTempDir(), 'snapshot')
        os.mkdir(snapshotDir)
        snapshotPath = os.path.join(snapshotDir, KTSERVER_SNAPSHOT_NAME)
        if existingSnapshotID is not None:
            # Extract the existing snapshot to the snapshot
            # directory so it will be automatically loaded
            fileStore.readGlobalFile(existingSnapshotID, userPath=snapshotPath)
        process = cactus_call(server=True, shell=False,
                              parameters=getKtserverCommand(dbElem, logPath, snapshotDir),
                              port=dbElem.getDbPort())

        blockUntilKtserverIsRunning(logPath)
        if existingSnapshotID is not None:
            # Clear the termination flag from the snapshot
            cactus_call(parameters=["ktremotemgr", "remove"] + getRemoteParams(dbElem) + ["TERMINATE"])

        while True:
            # Check for the termination signal
            try:
                cactus_call(parameters=["ktremotemgr", "get"] + getRemoteParams(dbElem) + ["TERMINATE"],
                            swallowStdErr=True)
            except:
                # No terminate signal sent yet
                pass
            else:
                # Terminate signal received
                break
            # Check that the DB is still alive
            if process.poll() is not None or isKtServerFailed(logPath):
                with open(logPath) as f:
                    raise RuntimeError("KTServer failed. Log: %s" % f.read())
            sleep(60)
        process.send_signal(signal.SIGINT)
        process.wait()
        blockUntilKtserverIsFinished(logPath)
        if snapshotExportID is not None:
            if not os.path.exists(snapshotPath):
                with open(logPath) as f:
                    raise RuntimeError("KTServer did not leave a snapshot on termination,"
                                       " but a snapshot was requested. Log: %s" % f.read())
            if len(glob(os.path.join(snapshotDir, "*.ktss"))) != 1:
                # More than one snapshot file. It's not clear what
                # conditions trigger this--if any--but we
                # don't support it right now.
                with open(logPath) as f:
                    raise RuntimeError("KTServer left more than one snapshot. Log: %s" % f.read())

            # Export the snapshot file to the file store
            fileStore.jobStore.updateFile(snapshotExportID, snapshotPath)

def blockUntilKtserverIsRunning(logPath, createTimeout=1800):
    """Check status until it's successful, an error is found, or we timeout.

    Returns True if the ktserver is now running, False if something went wrong."""
    success = False
    for i in range(createTimeout):
        if isKtServerFailed(logPath):
            logger.critical('Error starting ktserver.')
            success = False
            break
        if isKtServerRunning(logPath):
            logger.info('Ktserver running.')
            success = True
            break
        sleep(1)
    return success

def blockUntilKtserverIsFinished(logPath, timeout=1800,
                                 timeStep=10):
    """Wait for the ktserver log to indicate that it shut down properly.

    Returns True if the server shut down, False if the timeout expired."""
    for i in range(0, timeout, timeStep):
        with open(logPath) as f:
            log = f.read()
            if '[FINISH]' in log:
                return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for ktserver.")

def isKtServerRunning(logPath):
    """Check if the server started running."""
    success = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("listening") >= 0:
                success = True
    return success

def isKtServerFailed(logPath):
    """Does the server log contain an error?"""
    isFailed = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("error") >= 0:
                isFailed = True
                break
    return isFailed

def getKtTuningOptions(dbElem):
    """Get the appropriate KTServer tuning parameters (bucket size, etc.)"""
    # these are some hardcoded defaults.  should think about moving to config
    tuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
    # override default ktserver settings if they are present in the
    # experiment xml file.
    if dbElem.getDbTuningOptions() is not None:
        tuningOptions = dbElem.getDbTuningOptions()
    if dbElem.getDbCreateTuningOptions() is not None:
        tuningOptions = dbElem.getDbCreateTuningOptions()
    return tuningOptions

def getKtServerOptions(dbElem):
    # these are some hardcoded defaults.  should think about moving to config
    serverOptions = "-ls -tout 200000 -th 64"
    if dbElem.getDbServerOptions() is not None:
        serverOptions = dbElem.getDbServerOptions()
    return serverOptions

def getKtserverCommand(dbElem, logPath, snapshotDir):
    """Get a ktserver command line with the proper options (in popen-type list format)."""
    serverOptions = getKtServerOptions(dbElem)
    tuning = getKtTuningOptions(dbElem)
    cmd = ["ktserver", "-port", str(dbElem.getDbPort())]
    cmd += serverOptions.split()
    # Configure background snapshots, but set the interval between
    # snapshots to ~ 10 days so it'll never trigger. We are only
    # interested in the snapshot that the DB creates on termination.
    cmd += ["-bgs", snapshotDir, "-bgsc", "lzo", "-bgsi", "1000000"]
    cmd += ["-log", logPath]
    cmd += [":" + tuning]
    return cmd

def getRemoteParams(dbElem):
    """Get parameters to supply to ktremotemgr to connect to the right DB."""
    host = dbElem.getDbHost() or 'localhost'
    return ['-port', str(dbElem.getDbPort()),
            '-host', host]

def stopKtserver(dbElem):
    """Attempt to send the terminate signal to a ktserver."""
    cactus_call(parameters=['ktremotemgr', 'set'] + getRemoteParams(dbElem) + ['TERMINATE', '1'])
