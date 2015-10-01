#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Some functions to launch, wait for, kill, and test ktservers.  They are
intended to be used as follows:
- Controlling process creates a unique file on a globally accessible disk
  to govern the lifespan of the server.
  
- Subsequent process calls runKtserver() to launch a ktserver on whatever
  machine it is run from.  This server will run until the unique file
  is erased or until the process running runKtserver() is terminated.

- An independent process (which depends on the server existing) calls
  blockUntilKtserverIsRunnning() to wait until runKtserver() has
  been successfully launched (keep in mind that we do not know a
  prirori in which order runKtserver() and blockUntil..() actually
  get executed.

- After blockUntilKtserverIsRunnning() returns, subsequent proceseses can
  access the server.  Once these are done, the server can be killed
  using killKtServer().
"""

import os
import sys
import random
import math
import subprocess
import signal
import psutil
import socket
from sonLib.bioio import logger
from time import sleep
from optparse import OptionParser
import threading
import xml.etree.ElementTree as ET
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.experimentWrapper import ExperimentWrapper

###############################################################################
# run a server until killSwitchPath gets deleted
# throws an exception if we were unable to launch the server
###############################################################################
def runKtserver(dbElem, killSwitchPath, maxPortsToTry=100, readOnly = False,
                createTimeout=30, loadTimeout=10000, killTimeout=518400,
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
    logDir = os.path.dirname(logPath)
    if os.path.exists(logDir) is False:
        os.makedirs(logDir)
    assert os.path.isdir(logDir)
    basePort = dbElem.getDbPort()
    host = dbElem.getDbHost()
    if host is not None:
        dbElem.setDbHost(host)
    else:
        dbElem.setDbHost(getHostName())

    success = False
    process = None
    procWaiter = None
    try:        
        for port in xrange(basePort, basePort + maxPortsToTry):
            dbElem.setDbPort(port)
            if os.path.exists(logPath):
                os.remove(logPath)
            if __isKtServerOnTakenPort(dbElem, killSwitchPath, pretest=True):
                continue
            cmd = __getKtserverCommand(dbElem, dbPathExists, readOnly)
            process = subprocess.Popen(cmd.split(), shell=False, 
                                       stdout=subprocess.PIPE,
                                       stderr=sys.stderr, bufsize=-1)
            procWaiter = ProcessWaiter(process)
            procWaiter.start()
            __writeStatusToSwitchFile(dbElem, process.pid, killSwitchPath)
            success = __validateKtserver(process, dbElem, killSwitchPath,
                                         createTimeout, loadTimeout)
            if success is True:
                break
            else:
                if process.returncode is None:
                    process.kill()
                process = None

        if success is False:
            raise RuntimeError("Unable to launch ktserver.  "+
                               "Server log is: %s" % logPath)
            
    except Exception as e:
        # make an attempt to alert the world of the launch failure by
        # writing some -1's to the switch file
        switchFile = open(killSwitchPath, "w")
        switchFile.write("-1\n-1\n-1\n")
        switchFile.close()

        # if we don't kill the spawned process, the waiter thread will keep
        # this process alive which we don't want in the case of an error
        if process is not None and process.returncode is None:
            process.terminate()
        #assert procWaiter is None or procWaiter.is_alive() is False
        
        raise e
    return process

###############################################################################
# Check status until it's successful, an error is found, or we timeout
###############################################################################
def __validateKtserver(process, dbElem, killSwitchPath,
                       createTimeout, loadTimeout):
    success = False
    for i in xrange(createTimeout):
        if process.returncode is not None:
            break
        if __isKtServerFailed(dbElem) or __isKtServerOnTakenPort(dbElem,
                                                                 killSwitchPath):
            success = False
            break
        if __isKtServerRunning(dbElem, killSwitchPath):
            success = True
            break
        if __isKtServerReorganizing(dbElem):
            raiseTimeout = True
            for j in xrange(loadTimeout):
                if __isKtServerReorganizing(dbElem) is False:
                    raiseTimeout = False
                    break
                sleep(1)
            if raiseTimeout is True:
                raise RuntimeError("Reorganization wait timeout failed. " +
                                   "Server log is: %s" % logPath)
        sleep(1)
    return success

###############################################################################
# We use the kill-switch file to store some vital information about the
# kterver that's not always obvious to scrape from the log file
###############################################################################
def __writeStatusToSwitchFile(dbElem, serverPid, killSwitchPath):
    try:
        switchFile = open(killSwitchPath, "w")
        switchFile.write("%s\n%d\n%d\n" % (dbElem.getDbHost(), 
                                           int(dbElem.getDbPort()),
                                           int(serverPid)))
        switchFile.close()
        return True
    except:
        raise RuntimeError("Unexpected error updating killswitch file " %
                          killSwitchPath)
        
###############################################################################
# We use the kill-switch file to store some vital information about the
# kterver that's not always obvious to scrape from the log file
# Returns false if the killswitch file is empty, true of the information was
# able to be read. 
###############################################################################
def __readStatusFromSwitchFile(dbElem, serverPidAsList, killSwitchPath):
    try:
        assert isinstance(serverPidAsList, list)
        assert len(serverPidAsList) == 0
        if not os.path.isfile(killSwitchPath):
            return False
        switchFile = open(killSwitchPath, "r")
        host = switchFile.readline().strip()
        port = switchFile.readline().strip()
        serverPid = switchFile.readline().strip()
        if host == '' or port == '' or serverPid == '':
            return False
        port = int(port)
        serverPid = int(serverPid)    
    except:
        raise RuntimeError("Unexpected error reading killswitch file %s" %
                          killSwitchPath)

    if port < 0 or serverPid < 0:
        # try not to leave the switch file in an error state
        switchFile = open(killSwitchPath, "w")
        switchFile.write("")
        switchFile.close()
        raise RuntimeError("Ktserver polling detected fatal error")
    
    serverPidAsList.append(serverPid)    
    dbElem.setDbPort(port)
    dbElem.setDbHost(host)
    switchFile.close()
    return True
                                    

###############################################################################
# Wait until a ktserver is running
# Note that this function will update dbElem with the currnet host/port
# information of the server
###############################################################################
def blockUntilKtserverIsRunnning(dbElem, killSwitchPath, timeout=518400,
                                 timeStep=10):
    logPath = getLogPath(dbElem)
    for i in xrange(0, timeout, timeStep):
        if os.path.isfile(logPath) and __isKtServerRunning(dbElem,
                                                           killSwitchPath):
            return True
        sleep(timeStep)
    raise RuntimeError("Timeout reached while waiting for ktserver" %
                       logPath)

###############################################################################
# Kill a server by deleting the given kill switch file.  Check that it's
# no longer running using the timeout.  If it's being saved to disk, it
# may take a while for serialization to complete.  This method isn't going to
# wait.
# Note that this function will update dbElem with the currnet host/port
# information of the server
###############################################################################
def killKtServer(dbElem, killSwitchPath, killTimeout=10000):
    if not os.path.isfile(killSwitchPath):
        raise RuntimeError("Can't kill server because file" +
                           " not found %s" % killSwitchPath)
    logPath = getLogPath(dbElem)
    isRunning =  __isKtServerRunning(dbElem, killSwitchPath)
    os.remove(killSwitchPath)
    logPath = getLogPath(dbElem)
    if not isRunning:
        raise RuntimeError("Can't find running server to kill %s" % logPath)

    success = False
    for i in xrange(killTimeout):
        try:
            if pingKtServer(dbElem) or len(__scrapePids([logPath])) > 0:
                logger.critical("Waiting for ktserver to die, but still running with logPath: %s, pingKtServer returned %s, dB port: %s, __scrapePids returned: %s" % (logPath, pingKtServer(dbElem), dbElem.getDbPort(), __scrapePids([logPath])))
                sleep(1)
            else:
                success = True
                break
        except RuntimeError:
            logger.critical("Got runtime error while trying to kill ktserver, putting it down to bad luck and carrying on")
    if not success:
        raise RuntimeError("Failed to kill server within timeout. " +
                           "Server log is %s" % logPath)
    if os.path.exists(logPath):
        os.remove(logPath)
    return True

###############################################################################
# Test if a server is running by looking at the log
# if the log looks okay, verify by pinging the server
# note that some information is duplicated across the log and killswitch
# path.  if there are any inconsistencies we raise an exception.
# Note that this function will update dbElem with the currnet host/port
# information of the server
###############################################################################
def __isKtServerRunning(dbElem, killSwitchPath):
    logPath = getLogPath(dbElem)
    success = False
    serverPidAsList = []
    serverPidFromLog = None
    serverPortFromLog = None
    if os.path.isfile(logPath):
        try:
            outFile = open(logPath, "r")
        except:
            sleep(1)
            outFile = open(logPath, "r")
        for line in outFile.readlines():
            if line.lower().find("listening") >= 0:
                success = __readStatusFromSwitchFile(dbElem, serverPidAsList,
                                                     killSwitchPath)
            if serverPortFromLog is None and line.find("expr=") >= 0:
                try:
                    hostPort = line[line.find("expr="):].split()[0]                    
                    serverPortFromLog = int(hostPort[hostPort.find(":")+1:])
                    if serverPortFromLog < 0:
                        serverPortFromLog = None
                except:
                    serverPortFromLog = None
            if serverPidFromLog is None and line.find("pid=") >=0:
                try:
                    pidToken = line[line.find("pid=") + 4:].split()[0]
                    serverPidFromLog = int(pidToken)
                except:
                    serverPidFromLog = None
            if line.lower().find("error") >= 0:
                success = False
                break
    if success is False:
        return False

    if serverPidFromLog != serverPidAsList[0]:
        raise RuntimeError("Pid %s != %s (former from %s, lastter %s)" % (
            str(serverPidFromLog), str(serverPidAsList[0]),
            logPath, killSwitchPath))
    if serverPortFromLog != dbElem.getDbPort():
        raise RuntimeError("Port %s != %s (former from %s, lastter %s)" % (
            str(serverPortFromLog), str(dbElem.getDbPort()),
            logPath, killSwitchPath))
    
    return pingKtServer(dbElem)

###############################################################################
# Query a running server
###############################################################################
def pingKtServer(dbElem):
    canPing = subprocess.call(['ping', '-c', '1', dbElem.getDbHost()],
                              shell=False, bufsize=-1,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
    if canPing != 0:
        raise RuntimeError("Unable to ping ktserver host %s from %s" % (
            dbElem.getDbHost(), getHostName()))

    return subprocess.call(['ktremotemgr', 'report',
                            '-port', str(dbElem.getDbPort()),
                            '-host', dbElem.getDbHost()],
                           shell=False, bufsize=-1,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) == 0

###############################################################################
# Get a report on a running server
###############################################################################
def getKtServerReport(dbElem):
    assert pingKtServer(dbElem) is True
    process = subprocess.Popen(['ktremotemgr', 'report',
                                '-port', str(dbElem.getDbPort()),
                                '-host', dbElem.getDbHost()],
                               shell=False, bufsize=-1,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    report, nothing = process.communicate()

    process = subprocess.Popen(['ls', '-lh',
                                dbElem.getDbDir()],
                               shell=False, bufsize=-1,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    dirList, nothing = process.communicate()

    return "Report for %s:%d:\n%s\nContents of %s:\n%s\n" % (
        dbElem.getDbHost(), dbElem.getDbPort(), report,
        dbElem.getDbDir(), dirList)
    
###############################################################################
# Test if the server is reorganizing.  don't know what this means
# except that it can really add to the opening time
###############################################################################
def __isKtServerReorganizing(dbElem):
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
def __isKtServerFailed(dbElem):
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

def ktServerAlreadyRunning(dbElem):
    logPath = getLogPath(dbElem)
    pidList = __scrapePids([logPath])
    if len(pidList) > 0:
        return True
    return False
###############################################################################
# Check if a server has the same port as another server.  The only way to
# do this is by checking running processes (in this case with psutil). 
# Ktserver can sometimes catch these errors, but often doesn't.  Once you
# have two servers on the same port running all bets are off.
###############################################################################
def __isKtServerOnTakenPort(dbElem, killSwitchPath, pretest = False):
    if __isKtServerReorganizing(dbElem) or __isKtServerRunning(dbElem,
                                                               killSwitchPath):
        logPath = getLogPath(dbElem)
        pidList = __scrapePids([logPath])
        if pretest is False:
            thresh = 1
        else:
            thresh = 0
        if len(pidList) > thresh:
            raise RuntimeError("Ktserver already found running with log %s" %
                               logPath)
        pidList = __scrapePids(['port %d' % dbElem.getDbPort()])
        if len(pidList) > thresh:
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
    
