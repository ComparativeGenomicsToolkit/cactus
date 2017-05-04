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
  blockUntilKtserverIsRunning() to wait until runKtserver() has
  been successfully launched (keep in mind that we do not know a
  prirori in which order runKtserver() and blockUntil..() actually
  get executed.

- After blockUntilKtserverIsRunning() returns, subsequent proceseses can
  access the server.  Once these are done, the server can be killed
  using killKtServer().
"""

import os
import subprocess
import psutil
import socket
import shutil
import math
import glob
from multiprocessing import Process
from time import sleep

from toil.lib.bioio import logger
import threading
from sonLib.bioio import getTempFile, getTempDirectory, system, popenCatch

from cactus.shared.common import cactus_call, pullCactusImage

###############################################################################
# run a server until killSwitchPath gets deleted
# throws an exception if we were unable to launch the server
###############################################################################
def runKtserver(dbElem, killSwitchPath, maxPortsToTry=100, readOnly = False,
                createTimeout=30, loadTimeout=10000, killTimeout=518400,
                killPingInterval=5, fileStore = None):
    if not os.path.isfile(killSwitchPath):
        raise RuntimeError("Kill switch file not found, can't " +
                           "launch without it %s" % killSwitchPath)

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
    assert os.path.exists(dbElem.getDbDir())

    # Ensure the cactus image is pulled here, because we have a 10
    # second limit on the ktserver calls, which probably won't be
    # enough to download the image
    pullCactusImage()
    success = False
    process = None
    procWaiter = None
    for port in xrange(basePort, basePort + maxPortsToTry):
        dbElem.setDbPort(port)
        logger.info("Trying port: %i" % port)
        if os.path.exists(logPath):
            os.remove(logPath)

        if __isKtServerOnTakenPort(dbElem, killSwitchPath, pretest=True):
            logger.info("Ktserver already on port %i" % port)
            continue
        process = cactus_call(server=True, shell=False, work_dir=os.path.dirname(logPath),
                              parameters=['ktserver', '-log', logPath,
                                          '-port', dbElem.getDbPort(),
                                          __getKtServerOptions(dbElem)])

        procWaiter = ProcessWaiter(process)
        procWaiter.start()
        __writeStatusToSwitchFile(dbElem, process.pid, killSwitchPath)
        validated = __validateKtserver(process, dbElem, killSwitchPath,
                                     createTimeout, loadTimeout)
        logger.info("Was able to validate ktserver: %r" % validated)
        success = validated
        if success is True:
            break
        else:
            if process.returncode is None:
                process.kill()
            process = None

    if success is False:
        raise RuntimeError("Unable to launch ktserver.  "+
                           "Server log is: %s" % logPath)
            
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
        if __isKtServerFailed(dbElem):
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
                raise RuntimeError("Reorganization wait timeout failed.")
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
def blockUntilKtserverIsRunning(dbElem, killSwitchPath, timeout=518400,
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

    #if serverPidFromLog != serverPidAsList[0]:
    #    raise RuntimeError("Pid %s != %s (former from %s, latter %s)" % (
    #        str(serverPidFromLog), str(serverPidAsList[0]),
    #        logPath, killSwitchPath))
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

    return cactus_call(check_result=True,
                       parameters=['ktremotemgr', 'report'] + getRemoteParams(dbElem)) == 0

###############################################################################
# Get a report on a running server
###############################################################################
def getKtServerReport(dbElem):
    assert pingKtServer(dbElem) is True
    process = subprocess.Popen(['ktremotemgr', 'report'] + getRemoteParams(dbElem),
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
    work_dir = os.path.dirname(os.path.abspath(logPath))
    cmd = "docker run --interactive -v %s:/data --log-driver=none -p %d:%d --net=host quay.io/comparative-genomics-toolkit/ktserver:%s -log %s -port %d %s" % (work_dir, port, port, cactus_commit, os.path.basename(logPath), port, serverOptions)
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

def getRemoteParams(dbElem):
    """Get parameters to supply to ktremotemgr to connect to the right DB."""
    return ['-port', str(dbElem.getDbPort()),
            '-host', dbElem.getDbHost() or 'localhost']

def chunkFileByLines(path, pieces):
    """Chunk up the file by lines into N pieces."""
    numLines = int(popenCatch("wc -l %s" % path).split()[0])
    linesPerFile = int(math.ceil(float(numLines)/pieces))
    tempDir = getTempDirectory()
    system("split -l %s %s %s/part" % (linesPerFile, path, tempDir))
    return tempDir

def dumpKtServer(dbElem, path, numProcesses=10):
    """Dump the KT server to 'path' as a gzipped TSV file."""
    temp = getTempFile()
    cactus_call(outfile=temp, parameters=['ktremotemgr', 'list', '-px'] + getRemoteParams(dbElem))
    # The ktremotemgr hex output includes spaces between each byte that it
    # can't handle when reading hex input
    system("sed -i 's/ //g' %s" % temp)
    tempDir = chunkFileByLines(temp, numProcesses)
    os.remove(temp)
    # Launch several processes.
    processes = []
    files = []
    for subfile in glob.glob(os.path.join(tempDir, '*')):
        outfile = getTempFile()
        files.append(outfile)
        process = Process(target=lambda: cactus_call(infile=subfile, outfile=outfile,
                                                     parameters=['xargs', '-n', '50', 'ktremotemgr',
                                                                 'getbulk', '-sx', '-px']
                                                                + getRemoteParams(dbElem)))
        process.start()
        processes.append(process)

    map(lambda x: x.join(), processes)
    if any([p.exitcode != 0 for p in processes]):
        raise RuntimeError("Saving the DB failed")
    system("rm -fr %s" % tempDir)

    # Again, have to remove the spaces from the ktremotemgr output so it can parse it properly
    system("sed 's/ //g' %s > %s" % (" ".join(files), path))
    system("gzip %s" % path)
    shutil.move(path + '.gz', path)

def clearKtServer(dbElem):
    """Remove all data in the KT server."""
    cactus_call(parameters=['ktremotemgr', 'clear'] + getRemoteParams(dbElem))

def isGzipped(path):
    """Check for GZIP magic number in the file at path "path"."""
    gzipped = False
    with open(path) as f:
        if f.read(2) == '\x1f\x8b':
            gzipped = True
    return gzipped

def restoreKtServer(dbElem, path, numProcesses=30):
    """Load a KT server with data from 'path' (a gzipped TSV file)."""
    gzipped = isGzipped(path)
    if gzipped:
        temp = getTempFile()
        os.remove(temp)
        shutil.copy(path, temp + '.gz')
        system("gzip -d %s" % temp + '.gz')
        path = temp
    tempDir = chunkFileByLines(path, numProcesses)
    # Launch several processes.
    processes = []
    for subfile in glob.glob(os.path.join(tempDir, '*')):
        call = ['ktremotemgr', 'import', '-sx'] + getRemoteParams(dbElem) + [subfile]
        if os.environ.get("CACTUS_DOCKER_MODE") != "0":
            # Silence the ktremotemgr process, the amount of output
            # chokes Docker
            call += ['cactus-redirect', '/dev/null']
        process = Process(target=lambda: cactus_call(parameters=call))
        process.start()
        processes.append(process)

    map(lambda x: x.join(), processes)
    if any([p.exitcode != 0 for p in processes]):
        raise RuntimeError("Loading the DB failed")
    system("rm -fr %s" % tempDir)
    if gzipped:
        os.remove(temp)

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
