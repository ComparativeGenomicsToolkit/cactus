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
from time import sleep
from optparse import OptionParser
import xml.etree.ElementTree as ET

from cactus.progressive.experimentWrapper import ExperimentWrapper

class KtserverLauncher:
    def __init__(self):
        self.pid = -1  
        self.maxRunningServers = 8
        self.waitFromMax = 30
        self.rangeSize = 100
        self.listenWaitIntervals = 100
        self.listenWait = 1
        self.tuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
        # it seems that using msiz to open an existing db can 
        # cause an error sometimes, especially on kolossus
        self.readTuningOptions = "#opts=ls#ktopts=p"
        self.serverOptions = "-ls -tout 200000 -th 64"
    
    # return list of pids for ktservers with given keywords
    # in command line.  pids are stored in strings, and
    # empty list return if nothing found.
    def scrapePids(self, keywordList = []):
        gstring = "| grep ktserver"
        for keyword in keywordList:
            gstring += "| grep \"%s\"" % keyword
        gstring += "| grep -v grep"
        cmdline = "ps x %s | awk \'{print $1}\' | xargs" % gstring
        outPipe = subprocess.Popen(cmdline, shell=True, 
                                   stdout=subprocess.PIPE).stdout
        pString = outPipe.read()
        return pString.split()
    
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
        pids = self.scrapePids(['port %d' % port, dbPath])
        assert len(pids) < 2
        if len(pids) == 0:
            # server not running
            return -1
        pid = int(pids[0])
        
        # sometimes a server can be slow to start up
        # so we poll the file where stdout was piped for a
        # few seconds waiting to hear that it's "listening"
        assert os.path.exists(outPath)
        success = False
        for i in range(0, self.listenWaitIntervals):
            if success == False:
                sleep(self.listenWait)
                outFile = open(outPath, "r")
                for line in outFile.readlines():
                    if line.find("listening") >= 0:
                        success = True
                        break
                outFile.close()
        if success == True:
            return pid
        else:
            self.killServer()
            return -1
    
    def ktserverCmd(self, dbPath, port, exists):
        tuning = self.tuningOptions
        if exists:
            tuning = self.readTuningOptions
        return "ktserver -port %d %s %s%s" % (port, self.serverOptions, 
                                              dbPath, tuning)
    
    # launch the ktserver as a new process. The process is 
    # orphaned using disown so that it can live beyond (and doesn't
    # deadlock) its parent jobTree job's process. the pid is
    # remembered in self.pid and the port is written back to the
    # input experiment object         
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
        outputPath = "%s.stdout" % dbPath
        if os.path.exists(outputPath):
            os.remove(outputPath)
        dbPathExists = os.path.exists(dbPath)
        
        for port in range(basePort, basePort + self.rangeSize):
            if self.isServerOnPort(port) == False:
                os.system("nohup %s > %s 2>&1 < /dev/null & disown" % 
                          (self.ktserverCmd(dbPath, port, dbPathExists),
                           outputPath))
                self.pid = self.validateServer(dbPath, outputPath, port)
            if self.pid >= 0:
                experiment.setDbPort(port)
                break
        
        if self.pid < 0:
            raise RuntimeError("Failed to launch ktserver for %d" % dbPath)
 
    def killServer(self):
        assert self.pid is not None
        assert self.pid > -1
        os.kill(self.pid, signal.SIGKILL)
        
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
        print "PID: %d" % kts.pid
        return 0    
    except RuntimeError, e:
        print "ERROR " + str(e)
        return 1

if __name__ == '__main__':    
    main()
    