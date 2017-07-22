#!/usr/bin/env python
"""
Functions to launch and manage Redis servers.
"""

import os
import random
import socket
import signal
import sys
import traceback
from glob import glob
import redis
from multiprocessing import Process, Queue
from time import sleep
from toil.lib.bioio import logger
from cactus.shared.common import cactus_call
from os import listdir
from os.path import isfile, join


_MAX_SERVER_PORT = 65535

# The name of the snapshot that database outputs.
SERVER_SNAPSHOT_NAME = "dump.rdb"

def runserver(dbElem, fileStore, existingSnapshotID=None, snapshotExportID=None):
    """
    Run a DB. This function launches a separate python process that manages the server.

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
        unoccupiedPorts = set(xrange(1025,MAX_SERVER_PORT)) - occupiedPorts
        port = random.choice(list(unoccupiedPorts))
    except:
        logger.warning("Can't find which ports are occupied--likely netstat is not installed."
                       " Choosing a random port to start the DB on, good luck!")
        port = random.randint(1025,MAX_SERVER_PORT)
    dbElem.setDbPort(port)

    process = ServerProcess(dbElem, logPath, fileStore, existingSnapshotID, snapshotExportID)
    process.daemon = True
    process.start()

    if not blockUntilserverIsRunning(dbElem):
        raise RuntimeError("Unable to launch server in time.")

    return process, dbElem, logPath

class ServerProcess(Process):
    """Independent process that babysits the server process.

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
        snapshotPath = os.path.join(snapshotDir, SERVER_SNAPSHOT_NAME)
        if existingSnapshotID is not None:
            # Extract the existing snapshot to the snapshot
            # directory so it will be automatically loaded
            fileStore.readGlobalFile(existingSnapshotID, userPath=snapshotPath)
        
        process = cactus_call(server=True, shell=False,parameters=getserverCommand(dbElem, logPath, snapshotDir))
        blockUntilserverIsRunning(dbElem)

        # the format for debugging must be in redis-cli -p <Port>  MONITOR | head -n <number of lines to read> > <path to file> you cannot take out n or call assumes 0.
        cactus_call(server=True, shell=False, parameters=['redis-cli'] + getRemoteParams(dbElem) + ['MONITOR'], outfile=getFileName(), check_result=True)
        
        if existingSnapshotID is not None:
            # Clear the termination flag from the snapshot
            cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["del", "TERMINATE"])

        cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["get","TERMINATE"])
        while True:
            status=cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["get","TERMINATE"], check_output=True)
            if "1\n" in status:
                break
            else:
                continue
            sleep(1)

        cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["save"])
        cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["shutdown"])
        process.wait()
        if snapshotExportID is not None:
            if not os.path.exists(snapshotPath):
                raise RuntimeError("Redis did not leave a snapshot on termination but a snapshot was requested.")
            if len(glob(os.path.join(snapshotDir, "*.rdb"))) != 1:
                # More than one snapshot file. It's not clear what
                # conditions trigger this--if any--but we
                # don't support it right now.
                raise RuntimeError("Redis left more than one snapshot.")

            # Export the snapshot file to the file store
            fileStore.jobStore.updateFile(snapshotExportID, snapshotPath)


def getFileName():  
    onlyfiles = [f for f in listdir("/home/akul/snapshot/") if isfile(join("/home/akul/snapshot/", f))]
    max=0
    for files in onlyfiles:
       stat= files[7:8]
       if int(stat) > max:
           max = int(stat)
    filename="/home/akul/snapshot/outfile%d.txt" %(max+1)
    print filename
    return filename

def blockUntilserverIsRunning(dbElem, createTimeout=180):
    """Check status until it's successful, an error is found, or we timeout.

    Returns True if the server is now running, False if something went wrong."""
    success = False
    for i in xrange(createTimeout):
        running=cactus_call(shell= False, parameters= ["redis-cli -p", str(dbElem.getDbPort()) , " ping"], check_result=True)
        if (running==0):
           success=True
           break
        sleep(1)
    return success

def blockUntilserverIsFinished(logPath,dbElem, timeout=180,
                                 timeStep=10):
    """Wait for the server to indicate that it shut down properly.

    Returns True if the server shut down, False if the timeout expired."""
    for i in xrange(0, timeout, timeStep):
        running=cactus_call(shell= False, parameters= ["redis-cli -p", str(dbElem.getDbPort()) , " ping"], check_result=True)
        if (running!=0):
            return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for server.")

def isServerRunning(logPath):
    """Check if the server started running."""
    success = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("listening") >= 0:
                success = True
    return success

def isServerFailed(logPath):
    """Does the server log contain an error?"""
    isFailed = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("error") >= 0:
                isFailed = True
                break
    return isFailed

def getTuningOptions(dbElem):
    """Get the appropriate Server tuning parameters (bucket size, etc.)"""
    # these are some hardcoded defaults.  should think about moving to config
    tuningOptions = "#opts=ls#bnum=30m#msiz=50g#opts=p"
    # override default server settings if they are present in the
    # experiment xml file. 
    if dbElem.getDbTuningOptions() is not None:
        tuningOptions = dbElem.getDbTuningOptions()
    if dbElem.getDbCreateTuningOptions() is not None:
        tuningOptions = dbElem.getDbCreateTuningOptions()
    return tuningOptions

def getServerOptions(dbElem):
    # these are some hardcoded defaults.  should think about moving to config
    serverOptions = "-ls -tout 200000 -th 64"
    if dbElem.getDbServerOptions() is not None:
        serverOptions = dbElem.getDbServerOptions()
    return serverOptions


def getRemoteParams(dbElem):
    """Get parameters to supply to connect to the right DB."""
    return ['-p', str(dbElem.getDbPort())]

def stopserver(dbElem):
    """Attempt to send the terminate signal to a server."""
    try:
        cactus_call(parameters=['redis-cli'] + getRemoteParams(dbElem) + ['set','TERMINATE', '1'])
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

