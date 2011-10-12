#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Help with spawning ktservers.  Something more sophisticated
may well be necessary down the road... for now, everything is 
based on subprocess and ps x

"""
import unittest

import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import random
import math
import copy
import filecmp
import subprocess
from time import sleep
from optparse import OptionParser


from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.multiCactusProject import MultiCactusProject

class KtserverLauncher:
    def __init__(self):
        self.proc = None 
        self.max = 8
        self.waitFromMax = 10
        self.rangeSize = 100
        self.rangeWait = 2
        self.waitToKill = 2
        self.listenWaitIntervals = 100
        self.listenWait = 1
        self.tuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
        # it seems that using msiz to open an existing db can 
        # cause an error sometimes, especially on kolossus
        self.readTuningOptions = "#opts=ls#ktopts=p"
        self.serverOptions = "-ls -tout 200000 -th 64"
    
    def countRunningServers(self, keyword = None):
        gstring = ""
        if keyword:
            gstring = "| grep \"%s\"" % keyword
        cmdline = "ps x | grep ktserver %s | grep -v grep | wc -l" % gstring
        outPipe = subprocess.Popen(cmdline, shell=True, 
                                   stdout=subprocess.PIPE).stdout
        pString = outPipe.read()
        numServers = int(pString)
        assert numServers >= 0
        return numServers
    
    def isServerOnPort(self, port):
        numServers = self.countRunningServers('port %d' % port)
        assert numServers == 0 or numServers == 1
        return numServers == 1
    
    # if there are already too many servers running, wait...
    # otherwise, cycle through ports, beginning at the one 
    # specified in the experiment parameter, until a 
    # server is successfully launched.  experiment
    # is updated with the new port (and i guess should
    # be saved
    def spawnServer(self, experiment):
        assert self.countRunningServers(experiment.getDbDir()) == 0
        assert os.path.splitext(experiment.getDbName())[1] == ".kch"

        exists = os.path.isfile(os.path.join(experiment.getDbDir(), 
                                             experiment.getDbName()))
        block = True
        while block:
            sleep(random.uniform(0, self.rangeWait))
            numServers = self.countRunningServers()
            if numServers < self.max:
                block = False
            else:
                sleep(self.waitFromMax)
        basePort = experiment.getDbPort()
        
        # hangs unless i use an actual file !?
        listenFilePath = os.path.join(experiment.getDbDir(), "listen.txt")
        if os.path.isfile(listenFilePath):
            os.remove(listenFilePath)
        listenFile = open(listenFilePath, "w+")
        
        for port in range(basePort, basePort + self.rangeSize):
            sleep(random.uniform(0, self.rangeWait))
            if not self.isServerOnPort(port):
                experiment.setDbPort(port)
                cmdLine = self.ktserverCmd(experiment, exists)
                self.proc = subprocess.Popen(cmdLine.split(), shell=False,
                                             stdout=listenFile)
                sleep(self.rangeWait)
                if self.proc.poll() is None:
                    break
        assert self.proc.poll() is None
        
        success = False
        for i in range(0, self.listenWaitIntervals):
            if not success:
                listenRead = open(listenFilePath, "r")
                for line in listenRead.readlines():
                    if line.find("listening") >= 0:
                        success = True
                        break
                listenRead.close()
                sleep(self.listenWait)
        os.remove(listenFilePath)
        assert success
        
                        
    def ktserverCmd(self, experiment, exists):
        tuning = self.tuningOptions
        if exists:
            tuning = self.readTuningOptions
        return "ktserver -port %s %s %s/%s%s" % (experiment.getDbPort(),
                                                  self.serverOptions,
                                                  experiment.getDbDir(),
                                                  experiment.getDbName(),
                                                  tuning)
    
    def killServer(self):
        assert self.proc is not None
        assert self.proc.poll() is None
        sleep(self.waitToKill)
        self.proc.kill()

def main():
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

if __name__ == '__main__':    
    main()
    