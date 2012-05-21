#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Help with spawning ktservers.  Something more sophisticated
may well be necessary down the road (multiple hostnames?)... 
dependent on psutil package
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
import xml.etree.ElementTree as ET
from sonLib.bioio import spawnDaemon
from cactus.progressive.experimentWrapper import ExperimentWrapper

class KtserverLauncher:
    def __init__(self):
        self.maxRunningServers = 20
        self.waitFromMax = 30
        self.rangeSize = 100
        self.listenWaitIntervals = 1000
        self.listenWait = 1
        self.checkWaitIntervals = 30
        self.checkWait = 1
        self.killWaitIntervals = 1000
        self.killWait = 10
        self.createTuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
        # it seems that using msiz to open an existing db can 
        # cause an error sometimes, especially on kolossus
        self.readTuningOptions = "#opts=ls#ktopts=p"
        self.serverOptions = "-ls -tout 200000 -th 64"
    
    # use psutil to find all the running ktserver processes
    # todo: check user matches current user (though exceptions
    # seem to be working for now)
    def scrapePids(self, keywordList = []):
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
                        
    # determine if there is a ktserver process running
    # on a given port.  if there are too many, wait a bit
    # before aborting because maybe one is in process of
    # closing
    def isServerOnPort(self, port):
        serverPids = self.scrapePids(['port %d' % port])
        numServers = len(serverPids)
        for i in xrange(self.checkWaitIntervals):
            if numServers == 0 or numServers == 1:
                return numServers == 1
            else:
                sleep(self.checkWait)
                serverPids = self.scrapePids(['port %d' % port])
                numServers = len(serverPids)
                
        raise RuntimeError("%d servers found open on port %s" % (numServers, 
                                                                 port))
    
    # hang until fewer than self.maxRunningServers are running
    def waitOnTotalNumberOfServers(self):
        block = True
        while block:
            numServers = len(self.scrapePids())
            if numServers < self.maxRunningServers:
                block = False
            else:
                sleep(self.waitFromMax)
    
    # a valid server is a) running and b) listening 
    # we check for listening by looking into the output file 
    # return pid if running or -1 if fail       
    def validateServer(self, dbPath, outPath, port):     
        # sometimes a server can be slow to start up
        # so we poll the file where stdout was piped for a
        # few seconds waiting to hear that it's "listening"
        success = None
        for i in range(0, self.listenWaitIntervals):
            sleep(self.listenWait)
            if success == None:
                if os.path.exists(outPath):                    
                    outFile = open(outPath, "r")
                    for line in outFile.readlines():
                        if line.find("listening") >= 0:
                            success = True
                            break
                        if line.find("ERROR") >= 0 or line.find("error") >= 0:
                            success = False
                            break
                    outFile.close()
            else:
                break
        
        if success == True:
            pids = self.scrapePids(['port %d' % port, dbPath])
            if len(pids) == 1:
                return int(pids[0])
            else:
                raise RuntimeError("failed to open ktserver on %s" % dbPath)
                
        else:
            for i in range(0, self.checkWaitIntervals):
                sleep(self.checkWait)
                pids = self.scrapePids(['port %d' % port, dbPath])
                if len(pids) == 0:
                    return -1               
            raise RuntimeError("failed ktserver stuck on %s" % dbPath)
        
    def ktserverCmd(self, experiment, outputPath, port, exists):
        tuning = self.createTuningOptions
        if exists:
            tuning = self.readTuningOptions
        cmd = "ktserver -log %s -port %d %s" % (outputPath, port, self.serverOptions)
        if experiment.getDbHost() is not None:
            cmd += " -host %s" % experiment.getDbHost()
        if experiment.getDbSnapshot() == True:
            cmd += " -bgs %s -bgsi 100000000" % experiment.getDbDir()
        if experiment.getDbInMemory() == False:
            cmd += " %s%s" % (os.path.join(experiment.getDbDir(), 
                                           experiment.getDbName()), tuning)
        return cmd
    
    # launch the ktserver as a new daemon process.  
    def spawnServer(self, experiment):        
        dbPath = os.path.join(experiment.getDbDir(), experiment.getDbName())
        if (len(self.scrapePids([dbPath])) != 0):
            raise RuntimeError("ktserver already running on %s" % dbPath)
        
        # override default ktserver settings if they are present in the
        # epxeriment xml file. 
        if experiment.getDbServerOptions() is not None:
            self.serverOptions = experiment.getDbServerOptions()
        if experiment.getDbTuningOptions() is not None:
            self.createTuningOptions = experiment.getDbTuningOptions()
            self.readTuningOptions = experiment.getDbTuningOptions()
        if experiment.getDbCreateTuningOptions() is not None:
            self.createTuningOptions = experiment.getDbCreateTuningOptions()
        if experiment.getDbReadTuningOptions() is not None:
            self.readTuningOptions = experiment.getDbReadTuningOptions()
        
        self.waitOnTotalNumberOfServers()        
        basePort = experiment.getDbPort()
        outputPath = os.path.join(experiment.getDbDir(), "ktout.log")
        
        dbPathExists = False
        if experiment.getDbInMemory() == False:
            assert os.path.splitext(experiment.getDbName())[1] == ".kch"
            dbPathExists = os.path.exists(dbPath)
            
        for port in range(basePort, basePort + self.rangeSize):
            if self.isServerOnPort(port) == False:
                if os.path.exists(outputPath):
                    os.remove(outputPath)
                # shouldn't be necessary but try to reduce the occurrence of
                # freak concurrency issues by taking a quick nap
                sleep(random.uniform(0, 1))
                spawnDaemon(self.ktserverCmd(experiment, outputPath, port, 
                                             dbPathExists))                                             
                pid = self.validateServer(dbPath, outputPath, port)
                if pid >= 0:
                    experiment.setDbPort(port)
                    break                
        
        if pid < 0:
            raise RuntimeError("Failed to launch ktserver for %s" % dbPath)
        
        if (len(self.scrapePids([dbPath])) > 1):
            raise RuntimeError("Multiple ktservers running on %s" % dbPath)
                 
    def killServer(self, experiment):
        dbPath = os.path.join(experiment.getDbDir(), experiment.getDbName())
        port = experiment.getDbPort()
        pids = self.scrapePids(['port %d' % port, dbPath])
        if len(pids) == 0:
            raise RuntimeError("Can't find ktserver to kill for %s on port %d" \
                               % (dbPath, port))
        elif len(pids) > 1:
            raise RuntimeError("Multiple ktservers running on %s" % dbPath)
        else:
            pid = pids[0]
            assert pid is not None
            assert pid > -1
        os.kill(pid, signal.SIGTERM)
        
        # wait until it's really dead (can take a while for larger dbs since 
        # they have to get serialized at this point
        for i in xrange(self.killWaitIntervals):
            serverPids = self.scrapePids(['port %d' % port, dbPath])
            numServers = len(serverPids)
            if numServers == 0:
                return
            else:
                sleep(self.killWait)
        raise RuntimeError("Failed to kill ktserver on port %d and path %s" % (port, dbPath))
            
        
def main():
    try:
        usage = "usage: %prog <experiment file>"
        description = "Open ktserver of a given experiment"
        parser = OptionParser(usage=usage, description=description)
        
        options, args = parser.parse_args()
        
        if len(args) != 1:
            parser.print_help()
            raise RuntimeError("Wrong number of arguments")
        
        exp = ExperimentWrapper(ET.parse(args[0]).getroot())
        kts = KtserverLauncher()
        kts.spawnServer(exp)
        print exp.getDiskDatabaseString()
        return 0    
    except RuntimeError, e:
        print "ERROR " + str(e)
        return 1
    
if __name__ == '__main__':    
    main()
    
