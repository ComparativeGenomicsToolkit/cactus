#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Some functions to launch, wait for, kill, and test ktservers
"""

import os
import sys
import random
import math
import subprocess
import signal
import psutil
from time import sleep
from optparse import OptionParser
import threading
import xml.etree.ElementTree as ET
from cactus.progressive.experimentWrapper import DbElemWrapper
from cactus.progressive.experimentWrapper import ExperimentWrapper

###############################################################################
# run a server until killSwitchPath gets deleted
# throws an exception if we were unable to launch the server
###############################################################################
def runKtserver(dbElem, killSwitchPath, maxPortsToTry=100, readOnly = False,
                createTimeout=30, loadTimeout=10000, killTimeout=sys.maxint,
                killPingInterval=5):
    if not os.path.isfile(killSwitchPath):
        raise RuntimeError("Kill switch file not found, can't " +
                           "launch without it %s" % killSwitchPath)
    dbPathExists = False
    if dbElem.getDbInMemory() == False:
        if os.path.splitext(dbElem.getDbName())[1] == ".kch":
            raise RuntimeError("Expected path to end in .kch" %
                               dbElem.getDbName())
        dbPathExists = os.path.exists(dbPath)

    logPath = getLogPath(dbElem)
    basePort = dbElem.getDbPort()
    success = False
    process = None
    procWaiter = None
    try:        
        for port in xrange(basePort, basePort + maxPortsToTry):
            dbElem.setDbPort(port)
            if os.path.exists(logPath):
                os.remove(logPath)               
            cmd = __getKtserverCommand(dbElem, dbPathExists, readOnly)
            process = subprocess.Popen(cmd, shell=True, 
                                       stdout=subprocess.PIPE,
                                       stderr=sys.stderr, bufsize=-1)
            procWaiter = ProcessWaiter(process)
            procWaiter.start()
            success = __validateKtserver(process, dbElem, createTimeout, loadTimeout)
            if success is True:
                break
            else:
                process.kill()
                process = None

        if success is False:
            raise RuntimeError("Unable to launch ktserver.  "+
                               "Server log is: %s" % logPath)
        assert process is not None
        for i in xrange(0, killTimeout, killPingInterval):
            if process.returncode is not None:
                raise RuntimeError("ktserver finished before it could be " +
                                   "killed with the file. " +
                                   "  Server log is: %s" % logPath)
            if not os.path.isfile(killSwitchPath):
                process.terminate()
                return True
            sleep(killPingInterval)

        process.terminate()
        raise RuntimeError("Kill timeout %d reached." % killTimeout)
    
    except Exception as e:
        # if we don't kill the spawned process, the waiter thread will keep
        # this process alive which we don't want in the case of an error
        if procWaiter is not None and procWaiter.is_alive():
            process.terminate()
        raise e

###############################################################################
# Check status until it's successful, an error is found, or we timeout
###############################################################################
def __validateKtserver(process, dbElem, createTimeout, loadTimeout):
    success = False
    for i in xrange(createTimeout):
        if process.returncode is not None:
            break
        if isKtServerFailed(dbElem) or isKtServerOnTakenPort(dbElem):
            success = False
            break
        if isKtServerRunning(dbElem):
            success = True
            break
        if isKtServerReorganizing(dbElem):
            raiseTimeout = True
            for j in xrange(loadTimeout):
                if isKtServerReorganizing(dbElem) is False:
                    raiseTimeout = False
                    break
                sleep(1)
            if raiseTimeout is True:
                raise RuntimeError("Reorganization wait timeout failed. " +
                                   "Server log is: %s" % logPath)
        sleep(1)
    return success

###############################################################################
# Kill a server by deleting the given kill switch file.  Check that it's
# no longer running using the timeout.  If it's being saved to disk, it
# may take a while for serialization to complete.  This method isn't going to
# wait.
###############################################################################
def killKtServer(dbElem, killSwitchPath, killTimeout=10):
    if not os.path.isfile(logPath):
        raise RuntimeError("Can't kill server because file" +
                           " not found %s" % killSwitchPath)
    os.remove(logPath)
    success = False
    for i in xrange(killTimeout):
        if isKtServerRunning(dbElem):
            sleep(1)
        else:
            success = True
    if not success:
        raise RuntimeError("Failed to kill server. " +
                           "Server log is %s" % getLogPath(dbElem))
    return True

###############################################################################
# Test if a server is running by looking at the log
# if the log looks okay, verify by pinging the server
###############################################################################
def isKtServerRunning(dbElem):
    logPath = getLogPath(dbElem)
    success = False
    port = None
    if os.path.exists(logPath):                    
        outFile = open(logPath, "r")
        for line in outFile.readlines():
            if line.lower().find("listening") >= 0:
                success = True
            if port is None and line.find("expr=") >= 0:
                try:
                    hostPort = line[line.find("expr="):].split()[0]                    
                    port = int(hostPort[hostPort.find(":")+1:])
                    if port < 0:
                        port = None
                except:
                    port = None
            if line.lower().find("error") >= 0:
                success = False
                break
    if success is True and port is not None:
        success = subprocess.call(['ktremotemgr', 'report',
                                   '-port', str(port),
                                   '-host', dbElem.getDbHost()],
                                  shell=False, bufsize=-1,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE) == 0        
    return success is True and port is not None

###############################################################################
# Test if the server is reorganizing.  don't know what this means
# except that it can really add to the opening time
###############################################################################
def isKtServerReorganizing(dbElem):
    logPath = getLogPath(dbElem)
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

###############################################################################
# Test if the server log has an error.
###############################################################################
def isKtServerFailed(dbElem):
    logPath = getLogPath(dbElem)
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

###############################################################################
# Check if a server has the same port as another server.  The only way to
# do this is by checking running processes (in this case with psutil). 
# Ktserver can sometimes catch these errors, but often doesn't.  Once you
# have two servers on the same port running all bets are off.
###############################################################################
def isKtServerOnTakenPort(dbElem):
    if isKtServerReorganizing(dbElem) or isKtServerRunning(dbElem):
        logPath = getLogPath(dbElem)
        pidList = scrapePids([logPath])
        if len(pidList) > 1:
            raise RuntimeError("Ktserver already found running with log %s" %
                               logPath)
        pidList = scrapePids(['port %d' % dbElem.getDbPort()])
        if len(pidList) > 1:
            return True
    return False        

###############################################################################
# All access to the server log should use this method
###############################################################################
def getLogPath(dbElem):
    return os.path.join(dbElem.getDbDir(), "ktout.log")

###############################################################################
# Get the database tuning options.  They can be different depending
# one whether or not we are creating a new database
###############################################################################
def __getKtTuningOptions(dbElem, exists = False, readOnly = False):
    # these are some hardcoded defaults.  should think about moving to config
    createTuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
    readTuningOptions = "#opts=ls#ktopts=p"
    # override default ktserver settings if they are present in the
    # epxeriment xml file. 
    if dbElem.getDbServerOptions() is not None:
        serverOptions = dbElem.getDbServerOptions()
    if dbElem.getDbTuningOptions() is not None:
        createTuningOptions = dbElem.getDbTuningOptions()
        readTuningOptions = dbElem.getDbTuningOptions()
    if dbElem.getDbCreateTuningOptions() is not None:
        createTuningOptions = dbElem.getDbCreateTuningOptions()
    if dbElem.getDbReadTuningOptions() is not None:
        readTuningOptions = dbElem.getDbReadTuningOptions()
    tuning = createTuningOptions
    if exists or readOnly:
        tuning = self.readTuningOptions
    return tuning

###############################################################################
# Get the ktserver options
###############################################################################
def __getKtServerOptions(dbElem):
    # these are some hardcoded defaults.  should think about moving to config
    serverOptions = "-ls -tout 200000 -th 64"
    if dbElem.getDbServerOptions() is not None:
        serverOptions = dbElem.getDbServerOptions()
    return serverOptions

###############################################################################
# Construct the ktserver command line from the xml database element.
###############################################################################
def __getKtserverCommand(dbElem, exists = False, readOnly = False):
    logPath = getLogPath(dbElem)
    port = dbElem.getDbPort()
    serverOptions = __getKtServerOptions(dbElem)
    tuning = __getKtTuningOptions(dbElem, exists, readOnly)
    cmd = "ktserver -log %s -port %d %s" % (logPath, port, serverOptions)
    if readOnly is True and dbElem.getDbSnapshot() == False:
        cmd += " -ord -onr"
    if dbElem.getDbHost() is not None:
        cmd += " -host %s" % dbElem.getDbHost()
    if dbElem.getDbSnapshot() == True:
        cmd += " -bgs %s -bgsi 100000000" % dbElem.getDbDir()
    if dbElem.getDbInMemory() == False:
        cmd += " %s" % os.path.join(dbElem.getDbDir(), dbElem.getDbName())
    else:
        cmd += " :"
    cmd += tuning
    return cmd

###############################################################################
# By having a (non-daemon)thread waiting on the ktserver process at all times,
# we hope to guarantee that the server proc gets aborted if
# anything happens in the parent
###############################################################################
class ProcessWaiter(threading.Thread):
    def __init__(self, process):
        threading.Thread.__init__(self)
        self.__process = process
        self.daemon = False
    def run (self):
        self.__process.wait()

###############################################################################
# Use psutil to find all the running ktserver processes on the current
# machine.
###############################################################################
def __scrapePids(keywordList = []):
    pidList = []
    for proc in psutil.process_iter():
        try:
            if proc.cmdline[0].find("ktserver") >= 0:
                cmd = ' '.join(proc.cmdline)
                foundTerms = True
                for word in keywordList:
                    if cmd.find(word) < 0:
                        foundTerms = False
                        break
                if foundTerms:
                    pidList.append(proc.pid)
        except Exception:
            # probably somebody else's process we can't access
            pass
    return pidList

def main():
    try:
        usage = "usage: %prog <dbElem file> <killFile>"
        description = "Open ktserver of a given dbElem until <killFile> erased"
        parser = OptionParser(usage=usage, description=description)
        
        options, args = parser.parse_args()
        
        if len(args) != 2:
            parser.print_help()
            raise RuntimeError("Wrong number of arguments")

        exp = ExperimentWrapper(ET.parse(args[0]).getroot())

        runKtserver(exp, args[1])
        return 0    
    except RuntimeError, e:
        print "ERROR " + str(e)
        return 1
    
if __name__ == '__main__':    
    main()
    
