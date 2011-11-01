#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Help with spawning ktservers.  Something more sophisticated
may well be necessary down the road (multiple hostnames?)... 
for now, everything is based on os.system and grepping ps x
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
        self.maxRunningServers = 10
        self.waitFromMax = 30
        self.rangeSize = 100
        self.listenWaitIntervals = 100
        self.listenWait = 1
        self.tuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
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
                        
    # use ps to determine if there is a ktserver process running
    # on a given port
    def isServerOnPort(self, port):
        serverPids = self.scrapePids(['port %d' % port])
        numServers = len(serverPids)
        assert numServers == 0 or numServers == 1
        return numServers == 1
    
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
        success = False
        for i in range(0, self.listenWaitIntervals):
            sleep(self.listenWait)
            if success == False:
                if os.path.exists(outPath):                    
                    outFile = open(outPath, "r")
                    for line in outFile.readlines():
                        if line.find("listening") >= 0:
                            success = True
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
            for i in range(0, self.listenWaitIntervals):
                sleep(self.listenWait)
                pids = self.scrapePids(['port %d' % port, dbPath])
                if len(pids) == 0:
                    return -1               
            raise RuntimeError("failed ktserver stuck on %s" % dbPath)
    
    def ktserverCmd(self, dbPath, outputPath, port, exists):
        tuning = self.tuningOptions
        if exists:
            tuning = self.readTuningOptions
        return "ktserver -log %s -port %d %s %s%s" % (outputPath, port, self.serverOptions, 
                                              dbPath, tuning)
    
    # launch the ktserver as a new process. The process is 
    # orphaned using disown so that it can live beyond (and doesn't
    # deadlock) its parent jobTree job's process.  
    def spawnServer(self, experiment):
        dbPath = os.path.join(experiment.getDbDir(), experiment.getDbName())
        assert os.path.splitext(experiment.getDbName())[1] == ".kch"
        if (len(self.scrapePids([dbPath])) != 0):
            raise RuntimeError("ktserver already running on %s" % dbPath)
        
        # shouldn't be necessary but try to reduce the occurrence of
        # freak concurrency issues by taking a quick nap
        sleep(random.uniform(0, 1))
        
        self.waitOnTotalNumberOfServers()        
        basePort = experiment.getDbPort()
        outputPath = "%s.log" % dbPath
        if os.path.exists(outputPath):
            os.remove(outputPath)
        dbPathExists = os.path.exists(dbPath)
            
        for port in range(basePort, basePort + self.rangeSize):
            if self.isServerOnPort(port) == False:
                spawnDaemon(self.ktserverCmd(dbPath, outputPath, port, dbPathExists))
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
        os.kill(pid, signal.SIGKILL)
        
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
    