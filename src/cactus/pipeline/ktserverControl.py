#!/usr/bin/env python
"""
Functions to launch and manage KyotoTycoon servers.
"""

import os
import random
import socket
import shutil
import signal
from zipfile import ZipFile
from time import sleep
from multiprocessing import Process

from toil.lib.bioio import logger
from cactus.shared.common import cactus_call

# For some reason ktserver believes there are only 32768 TCP ports.
MAX_KTSERVER_PORT = 32767

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
        unoccupiedPorts = set(xrange(1025,MAX_KTSERVER_PORT)) - occupiedPorts
        port = random.choice(list(unoccupiedPorts))
    except:
        logger.warning("Can't find which ports are occupied--likely netstat is not installed."
                       " Choosing a random port to start the DB on, good luck!")
        port = random.randint(1025,MAX_KTSERVER_PORT)
    dbElem.setDbPort(port)

    process = Process(target=serverProcess,
                      args=(dbElem, logPath, fileStore, existingSnapshotID, snapshotExportID))
    process.daemon = True
    process.start()

    if not blockUntilKtserverIsRunning(logPath):
        raise RuntimeError("Unable to launch ktserver in time.")

    return dbElem, logPath

def serverProcess(dbElem, logPath, fileStore, existingSnapshotID=None, snapshotExportID=None):
    """Independent process that babysits the ktserver process.

    Waits for the TERMINATE flag to be set, then kills the DB and
    copies the final snapshot to snapshotExportID.
    """
    snapshotDir = os.path.join(fileStore.getLocalTempDir(), 'snapshot')
    os.mkdir(snapshotDir)
    if existingSnapshotID is not None:
        # Extract the existing snapshot (a zip file) to the snapshot
        # directory so it will be automatically loaded
        snapshot = fileStore.readGlobalFile(existingSnapshotID)
        zip = ZipFile(snapshot)
        zip.extractall(snapshotDir)
        zip.close()
    process = cactus_call(server=True, shell=False,
                          parameters=__getKtserverCommand(dbElem, logPath, snapshotDir))

    blockUntilKtserverIsRunning(logPath)
    if existingSnapshotID is not None:
        # Clear the termination flag from the snapshot
        cactus_call(parameters=["ktremotemgr", "remove"] + getRemoteParams(dbElem) + ["TERMINATE"])

    while True:
        try:
            cactus_call(parameters=["ktremotemgr", "get"] + getRemoteParams(dbElem) + ["TERMINATE"])
        except:
            # No terminate signal sent yet
            pass
        else:
            # Terminate signal received
            break
        sleep(60)
    process.send_signal(signal.SIGINT)
    process.wait()
    blockUntilKtserverIsFinished(logPath)
    if snapshotExportID is not None:
        # Export the snapshot file(s) to the file store
        zipBasePath = fileStore.getLocalTempFile()
        zipPath = shutil.make_archive(zipBasePath, 'zip', root_dir=snapshotDir)
        fileStore.jobStore.updateFile(snapshotExportID, zipPath)

def blockUntilKtserverIsRunning(logPath, createTimeout=1800):
    """Check status until it's successful, an error is found, or we timeout.

    Returns True if the ktserver is now running, False if something went wrong."""
    success = False
    for i in xrange(createTimeout):
        if __isKtServerFailed(logPath):
            logger.critical('Error starting ktserver.')
            success = False
            break
        if __isKtServerRunning(logPath):
            logger.info('Ktserver running.')
            success = True
            break
        sleep(1)
    return success

def blockUntilKtserverIsFinished(logPath, timeout=1800,
                                 timeStep=10):
    """Wait for the ktserver log to indicate that it shut down properly.

    Returns True if the server shut down, False if the timeout expired."""
    for i in xrange(0, timeout, timeStep):
        with open(logPath) as f:
            log = f.read()
            if '[FINISH]' in log:
                return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for ktserver" %
                       logPath)

def __isKtServerRunning(logPath):
    """Check if the server started running."""
    success = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("listening") >= 0:
                success = True
    return success

def __isKtServerFailed(logPath):
    """Does the server log contain an error?"""
    isFailed = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("error") >= 0:
                isFailed = True
                break
    return isFailed

def __getKtTuningOptions(dbElem):
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

def __getKtServerOptions(dbElem):
    # these are some hardcoded defaults.  should think about moving to config
    serverOptions = "-ls -tout 200000 -th 64"
    if dbElem.getDbServerOptions() is not None:
        serverOptions = dbElem.getDbServerOptions()
    return serverOptions

def __getKtserverCommand(dbElem, logPath, snapshotDir):
    """Get a ktserver command line with the proper options (in popen-type list format)."""
    serverOptions = __getKtServerOptions(dbElem)
    tuning = __getKtTuningOptions(dbElem)
    cmd = ["ktserver", "-host", dbElem.getDbHost(), "-port", dbElem.getDbPort()]
    cmd += serverOptions.split()
    cmd += ["-bgs", snapshotDir, "-bgsc", "lzo"]
    cmd += ["-log", logPath]
    cmd += [":" + tuning]
    return cmd

def getRemoteParams(dbElem):
    """Get parameters to supply to ktremotemgr to connect to the right DB."""
    return ['-port', str(dbElem.getDbPort()),
            '-host', dbElem.getDbHost() or 'localhost']

def stopKtserver(dbElem):
    """Attempt to send the terminate signal to a ktserver."""
    try:
        cactus_call(parameters=['ktremotemgr', 'set'] + getRemoteParams(dbElem) + ['TERMINATE', '1'])
    except:
        # The server is likely already down.
        pass

###############################################################################
# Hostnames of swarm nodes (ex kkr18u57.local) are not visible from
# other swarm nodes (dropping the .local does work).  So we try to
# work with ip addresses whenever possible. Current recipe (to review):
# if IP address does not start with 127. return IP address
# otherwise return primary hostname
###############################################################################
def getHostName():
    hostName = socket.gethostname()
    hostIp = socket.gethostbyname(hostName)
    if hostIp.find("127.") != 0:
        return hostIp
    return hostName

def findOccupiedPorts():
    """Attempt to find all currently taken TCP ports.

    Returns a set of ints, representing taken ports."""
    netstatOutput = cactus_call(parameters=["netstat", "-tuplen"], check_output=True)
    ports = set()
    for line in netstatOutput.split("\n"):
        fields = line.split()
        if len(fields) != 9:
            # Header or other garbage line
            continue
        port = int(fields[3].split(':')[-1])
        ports.add(port)
    logger.debug('Detected ports in use: %s' % repr(ports))
    return ports
