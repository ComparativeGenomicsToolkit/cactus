#!/usr/bin/env python

#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Bunch of functions used for progressive alignment.  In particular, they provide logic to 

decompose a phylogenetic tree into clades of a given maximum size.  Functions to create and 

manage the progressive directory structure are all here as well. 

"""

import os
import xml.etree.ElementTree as ET
import math
from optparse import OptionParser
from collections import deque
import random
import copy
import xml.etree.ElementTree as ET
from tempfile import mkstemp
from shutil import move
from os import remove, close
import subprocess
from time import sleep

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from sonLib.bioio import getLogLevelString
from sonLib.tree import getBinaryTreeNodes

from jobTree.src.bioio import getLogLevelString
from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions

from cactus.progressive.progressiveSplitUtils import getCladeLeaves

# we just hardcode these for now.  it's possible that they may
# need to be dynamically set (msiz in particular) or at least specified
# in the xml in the future.
tuningOptions = "#opts=ls#bnum=30m#msiz=50g#ktopts=p"
serverOptions = "-ls -tout 200000 -th 64"
portRangeSize = 100

def isKyotoTycoon(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    return dbConfElem.attrib["type"] == "kyoto_tycoon"

def getPort(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    return dbTypeElem.attrib["port"]

def setPort(options, port):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    dbTypeElem.attrib["port"] = str(port)
    options.cactusDiskDatabaseString = ET.tostring(dbConfElem)

def getHost(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    return dbTypeElem.attrib["host"]

def getDatabaseDir(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    return dbTypeElem.attrib["database_dir"]

def getDatabaseName(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    return dbTypeElem.attrib["database_name"]

def testServerExists(host, port):
    r = os.system("ktremotemgr report -port " + str(port) + " -host " + host +
                  "> /dev/null 2>&1")
    return r == 0

# Create a mapping between node IDs and port numbers
def creatKTPortMap(tree, options):
    basePort = int(getPort(options))
    dfStack = [tree]
    lookup = dict()
    while dfStack:
        node = dfStack.pop(len(dfStack)-1)
        if node is not None and node.internal:
            lookup[node.iD] = basePort
            basePort += 1
            dfStack.extend(getCladeLeaves(node, options))
    return lookup

def ktServerCommandLine(cladeOptions):
    return "ktserver -port %s %s %s/%s%s" % (getPort(cladeOptions),
                                                            serverOptions,
                                                            getDatabaseDir(cladeOptions),
                                                            getDatabaseName(cladeOptions),
                                                            tuningOptions)
                                                 
    return "ktserver -port %s %s \"%s/%s%s\"" % (getPort(cladeOptions),
                                                 serverOptions,
                                                 getDatabaseDir(cladeOptions),
                                                 getDatabaseName(cladeOptions),
                                                 tuningOptions)

# leave this in the working directory to help kill trailing servers
def updateKillScript(cladeOptions):
    scriptPath = getDatabaseDir(cladeOptions) + "/../../killKtServers.sh"
    os.system("echo kill -KILL %i >> %s" % (cladeOptions.serverProcess.pid, scriptPath))
    os.system("chmod u+x %s" % scriptPath)
    
# Spawn a server process on the next available port in range   
# Note: cladeOptions are updated with the successful port
def spawnLocalKtserver(node, cladeOptions):
    assert cladeOptions.autoKtserver
    assert isKyotoTycoon(cladeOptions)
    
    basePort = cladeOptions.portMap[node.iD]
    procHandle = None
    host = getHost(cladeOptions)
    
    if not os.path.isdir(getDatabaseDir(cladeOptions)):
        os.makedirs(getDatabaseDir(cladeOptions))
                
    for port in range(basePort, basePort + portRangeSize):
        
        if testServerExists(host, port) == False:
            setPort(cladeOptions, str(port))          
            procHandle = subprocess.Popen(ktServerCommandLine(cladeOptions).split(), shell=False)
            sleep(2)
            
            if procHandle.poll() is None:
                logger.info("Spawned ktserver with command: %s" % ktServerCommandLine(cladeOptions))      
                break
            else:
                logger.info("ktserver failed on port %i... retrying" % port)
                sleep(2)
                
    if procHandle is None:
        logger.critical("Could not open server in port range [%i, %i]" % (basePort, 
                                                                          basePort + portRangeSize-1))  

    # make sure the server is still running
    assert procHandle.poll() is None
    cladeOptions.serverProcess = procHandle
    updateKillScript(cladeOptions)

# kill the server
def killLocalKtserver(cladeOptions):
    assert cladeOptions.autoKtserver
    assert isKyotoTycoon(cladeOptions)
    assert testServerExists(getHost(cladeOptions), getPort(cladeOptions))
    assert cladeOptions.serverProcess.poll() is None
    cladeOptions.serverProcess.kill()
    logger.info("Killed ktserver started by command: %s" % ktServerCommandLine(cladeOptions))
