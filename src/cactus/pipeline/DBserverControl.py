#!/usr/bin/env python
"""
Functions to launch and manage Redis DBservers.
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


_MAX_DBSERVER_PORT = 65535

# The name of the snapshot that database outputs.
DBSERVER_SNAPSHOT_NAME = "dump.rdb"

def runDBserver(dbElem, fileStore, existingSnapshotID=None, snapshotExportID=None):
    """
    Run a DB. This function launches a separate python process that manages the DBserver.

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
        unoccupiedPorts = set(xrange(1025,MAX_DBSERVER_PORT)) - occupiedPorts
        port = random.choice(list(unoccupiedPorts))
    except:
        logger.warning("Can't find which ports are occupied--likely netstat is not installed."
                       " Choosing a random port to start the DB on, good luck!")
        port = random.randint(1025,_MAX_DBSERVER_PORT)
    dbElem.setDbPort(port)

    process = DBServerProcess(dbElem, logPath, fileStore, existingSnapshotID, snapshotExportID)
    process.daemon = True
    process.start()

    if not blockUntilDBserverIsRunning(dbElem):
        raise RuntimeError("Unable to launch DBserver in time.")

    return process, dbElem, logPath

class DBServerProcess(Process):
    """Independent process that babysits the DBserver process.

    Waits for the TERMINATE flag to be set, then kills the DB and
    copies the final snapshot to snapshotExportID.
    """
    exceptionMsg = Queue()

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        super(DBServerProcess, self).__init__()

    def run(self):
        """Run the tryRun method, signaling the main thread if an exception occurs."""
        try:
            self.tryRun(*self.args, **self.kwargs)
        except BaseException:
            self.exceptionMsg.put("".join(traceback.format_exception(*sys.exc_info())))
            raise

    def tryRun(self, dbElem, logPath, fileStore, existingSnapshotID=None, snapshotExportID=None)
        while True:
            cactus_call(server=True, shell=False, parameters=['redis-cli'] + getRemoteParams(dbElem) + ['MONITOR'], outfile=getFileName(),
                    check_result=True)
        
            if existingSnapshotID is not None:
                # Clear the termination and save flag from the snapshot
                cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["del", "TERMINATE"])
                cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["del", "SAVE"])


            cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["get","TERMINATE"])
            status = ""
            while True:
                status=cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["get","TERMINATE"], check_output=True)
                running= cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["get","SAVE"], check_output=True)
                if "1\n" in status or running::
                    break
                else:
                    continue
                sleep(1)

            if "1\n" in status:
                break

            #saves the database
            waitUntilSnapshotTaken(dbElem)

            #updates the database save file
            fileStore.jobStore.updateFile(snapshotExportID, snapshotPath)
   
        #saves the database and shuts it down
        cactus_call(parameters=["redis-cli"] + getRemoteParams(dbElem) + ["shutdown","save"])
        process.wait()

def waitUntilSnapshotTaken(dbElem, createTimeout=180):
    """Checks to see that the "Save" command works and snapshot was taken
    """
    success= False
    for i in xrange(createTimeout):
        running=cactus_call(shell= False, parameters= ["redis-cli -p", str(dbElem.getDbPort()) , " save"], check_output=True)
        if "ok\n" in running:
            success=True
            break
    return success

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

def blockUntilDBserverIsRunning(dbElem, createTimeout=180):
    """Check status until it's successful, an error is found, or we timeout.

    Returns True if the DBserver is now running, False if something went wrong."""
    success = False
    for i in xrange(createTimeout):
        running=cactus_call(shell= False, parameters= ["redis-cli -p", str(dbElem.getDbPort()) , " ping"], check_result=True)
        if (running==0):
            loading=cactus_call(shell= False, parameters= ["redis-cli -p", str(dbElem.getDbPort()) , " ping"], check_output=True)
            if 'LOADING' not in loading:
                success=True
                break
        sleep(1)
    return success

def blockUntilDBserverIsFinished(logPath,dbElem, timeout=180,
                                 timeStep=10):
    """Wait for the DBserver to indicate that it shut down properly.

    Returns True if the DBserver shut down, False if the timeout expired."""
    for i in xrange(0, timeout, timeStep):
        running=cactus_call(shell= False, parameters= ["redis-cli -p", str(dbElem.getDbPort()) , " ping"], check_result=True)
        if (running!=0):
            return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for DBserver.")

def isDBServerRunning(logPath):
    """Check if the DBserver started running."""
    success = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("listening") >= 0:
                success = True
    return success

def isDBServerFailed(logPath):
    """Does the DBserver log contain an error?"""
    isFailed = False
    with open(logPath) as f:
        for line in f:
            if line.lower().find("error") >= 0:
                isFailed = True
                break
    return isFailed

def getTuningOptions(dbElem):
    """Get the appropriate DBServer tuning parameters (bucket size, etc.)"""
    # these are some hardcoded defaults.  should think about moving to config
    tuningOptions = "#opts=ls#bnum=30m#msiz=50g#opts=p"
    # override default DBserver settings if they are present in the
    # experiment xml file. 
    if dbElem.getDbTuningOptions() is not None:
        tuningOptions = dbElem.getDbTuningOptions()
    if dbElem.getDbCreateTuningOptions() is not None:
        tuningOptions = dbElem.getDbCreateTuningOptions()
    return tuningOptions

def getDBserverCommand(dbElem, logPath, snapshotPath):
    """Get a ktserver command line with the proper options (in popen-type list format)."""
    cmd = ["redis-server", "--port",str(dbElem.getDbPort()),"--save \"\" --dir", snapshotPath]
    return cmd

def getDBServerOptions(dbElem):
    # these are some hardcoded defaults.  should think about moving to config
    DBserverOptions = "-ls -tout 200000 -th 64"
    if dbElem.getDbDBServerOptions() is not None:
        DBserverOptions = dbElem.getDbDBServerOptions()
    return DBserverOptions


def getRemoteParams(dbElem):
    """Get parameters to supply to connect to the right DB."""
    return ['-p', str(dbElem.getDbPort())]

def stopDBserver(dbElem):
    """Attempt to send the terminate signal to a DBserver."""
    try:
        cactus_call(parameters=['redis-cli'] + getRemoteParams(dbElem) + ['set','TERMINATE', '1'])
    except:
        # The DBserver is likely already down.
        pass

def saveDBserver(dbElem):
    """Attempt to send the terminate signal to a DBserver."""
    try:
        cactus_call(parameters=['redis-cli'] + getRemoteParams(dbElem) + ['set','SAVE', '1'])
    except:
        # The DBserver is likely already down.
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

