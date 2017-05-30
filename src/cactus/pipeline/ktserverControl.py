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

def runKtserver(dbElem, fileStore, createTimeout=30, loadTimeout=10000,
                existingSnapshotID=None, snapshotExportID=None):
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

    if not __validateKtserver(dbElem, logPath, createTimeout, loadTimeout):
        raise RuntimeError("Unable to launch ktserver.")

    return dbElem, logPath

def serverProcess(dbElem, logPath, fileStore, existingSnapshotID=None, snapshotExportID=None):
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


    blockUntilKtserverIsRunning(dbElem, logPath)
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

###############################################################################
# Check status until it's successful, an error is found, or we timeout
###############################################################################
def __validateKtserver(dbElem, logPath, createTimeout, loadTimeout):
    success = False
    for i in xrange(createTimeout):
        if __isKtServerFailed(logPath):
            success = False
            break
        if __isKtServerRunning(dbElem, logPath):
            success = True
            break
        if __isKtServerReorganizing(logPath):
            raiseTimeout = True
            for j in xrange(loadTimeout):
                if __isKtServerReorganizing(logPath) is False:
                    raiseTimeout = False
                    break
                sleep(1)
            if raiseTimeout is True:
                raise RuntimeError("Reorganization wait timeout failed.")
        sleep(1)
    return success

###############################################################################
# Wait until a ktserver is running
# Note that this function will update dbElem with the currnet host/port
# information of the server
###############################################################################
def blockUntilKtserverIsRunning(dbElem, logPath, timeout=518400,
                                timeStep=10):
    for i in xrange(0, timeout, timeStep):
        if os.path.isfile(logPath) and __isKtServerRunning(dbElem, logPath):
            return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for ktserver" %
                       logPath)

def blockUntilKtserverIsFinished(logPath, timeout=518400,
                                 timeStep=10):
    for i in xrange(0, timeout, timeStep):
        with open(logPath) as f:
            log = f.read()
            if '[FINISH]' in log:
                return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for ktserver" %
                       logPath)

###############################################################################
# Test if a server is running by looking at the log
# if the log looks okay, verify by pinging the server
# note that some information is duplicated across the log and killswitch
# path.  if there are any inconsistencies we raise an exception.
# Note that this function will update dbElem with the currnet host/port
# information of the server
###############################################################################
def __isKtServerRunning(dbElem, logPath):
    success = False
    serverPortFromLog = None
    if os.path.isfile(logPath):
        try:
            outFile = open(logPath, "r")
        except:
            sleep(1)
            outFile = open(logPath, "r")
        for line in outFile.readlines():
            if line.lower().find("listening") >= 0:
                success = True
            if serverPortFromLog is None and line.find("expr=") >= 0:
                try:
                    hostPort = line[line.find("expr="):].split()[0]                    
                    serverPortFromLog = int(hostPort[hostPort.find(":")+1:])
                    if serverPortFromLog < 0:
                        serverPortFromLog = None
                except:
                    serverPortFromLog = None
            if line.lower().find("error") >= 0:
                success = False
                break
    if success is False:
        return False

    if serverPortFromLog != dbElem.getDbPort():
        raise RuntimeError("Ktserver launched on unexpected port %s (expected %s)" % (
            str(serverPortFromLog), str(dbElem.getDbPort()),
            logPath))

    return True

###############################################################################
# Test if the server is reorganizing.  don't know what this means
# except that it can really add to the opening time
###############################################################################
def __isKtServerReorganizing(logPath):
    success = False
    if os.path.exists(logPath):                    
        outFile = open(logPath, "r")
        for line in outFile.readlines():
            if line.lower().find("listening") >= 0:
                success = False
                break
            if line.lower().find("error") >= 0:
                success = False
                break
            if line.lower().find("reorganizing") >= 0 or\
                   line.find("applying a snapshot") >= 0:
                success = True
    return success

def __isKtServerFailed(logPath):
    """Does the server log contain an error?"""
    isFailed = False
    if os.path.exists(logPath):                    
        outFile = open(logPath, "r")
        for line in outFile.readlines():
            if line.lower().find("listening") >= 0:
                isFailed = False
                break
            if line.lower().find("error") >= 0:
                isFailed = True
                break
            if line.lower().find("reorganizing") >= 0 or\
                   line.find("applying a snapshot") >= 0:
                isFailed = False
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
    """Attempt to send the terminate singal to a ktserver."""
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
    logger.debug(netstatOutput)
    for line in netstatOutput.split("\n"):
        fields = line.split()
        if len(fields) != 9:
            # Header or other garbage line
            continue
        port = int(fields[3].split(':')[-1])
        ports.add(port)
    logger.debug('Detected ports in use: %s' % repr(ports))
    return ports
