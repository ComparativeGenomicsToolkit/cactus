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
        self.tuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
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

        block = True
        while block:
            sleep(random.uniform(0, self.rangeWait))
            numServers = self.countRunningServers()
            if numServers < self.max:
                block = False
            else:
                sleep(self.waitFromMax)
        basePort = experiment.getDbPort()

        for port in range(basePort, basePort + self.rangeSize):
            sleep(random.uniform(0, self.rangeWait))
            if not self.isServerOnPort(port):
                experiment.setDbPort(port)
                cmdLine = self.ktserverCmd(experiment)
                self.proc = subprocess.Popen(cmdLine.split(), shell=False)
                sleep(self.rangeWait)
                if self.proc.poll() is None:
                    break
        assert self.proc.poll() is None
                        
    def ktserverCmd(self, experiment):
         return "ktserver -port %s %s %s/%s%s" % (experiment.getDbPort(),
                                                  self.serverOptions,
                                                  experiment.getDbDir(),
                                                  experiment.getDbName(),
                                                  self.tuningOptions)
    
    def killServer(self):
        assert self.proc is not None
        assert self.proc.poll() is None
        self.proc.kill()

    