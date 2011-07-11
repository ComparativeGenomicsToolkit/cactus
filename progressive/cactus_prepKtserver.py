#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Kyto Tycoon requires a server to be running before cactus alignment is begun.  The full path of 

each database must be passed to the server when it is launched (via ktserver).  The directory of

each database must exist beforehand.  This script parses a cactus experiment xml to determine 

the names and paths of each database (one per clade) required to run cactus_progressive with those

options.  A ktserver command is printed to stdout that can be used to manually initalise the server,

and all required database directories are created if they do not already exist.   

"""

import os
from optparse import OptionParser
import xml.etree.ElementTree as ET

from sonLib.bioio import newickTreeParser

from jobTree.src.bioio import getLogLevelString
from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions

from cactus.progressive.progressiveSplitUtils import nameUnlabeledInternalNodes
from cactus.progressive.progressiveSplitUtils import getCladeLeaves
from cactus.progressive.progressiveSplitUtils import getCladeDatabaseDir
from cactus.progressive.progressiveSplitUtils import getCladeDatabaseName

def recurseClades(root, options, cladeNodes):
    if root.left is not None and root.right is not None:
        dbDir = getCladeDatabaseDir(root, options)
        if not os.access(dbDir, os.F_OK):
            os.makedirs(dbDir)        
        cladeNodes.append(root)
        leaves = getCladeLeaves(root, options)
        for leaf in leaves:
            recurseClades(leaf, options, cladeNodes)
    elif root.left is not None or root.right is not None:
        assert False

def getDbPaths(cladeNodes, options):
    paths = ""
    for node in cladeNodes:
        paths += "\"" + getCladeDatabaseDir(node, options) + "/" + getCladeDatabaseName(node, options)
        paths += "#opts=ls#bnum=30m#msiz=50g#ktopts=p" + "\"" + " " 
        
    return paths

def getDbDirs(cladeNodes, options):
    paths = ""
    for node in cladeNodes:
        paths += getCladeDatabaseDir(node, options) + " "
    return paths

def getPort(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    return dbTypeElem.attrib["port"]

def getHost(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    return dbTypeElem.attrib["host"]

def testServerExists(options):
    r = os.system("ktremotemgr report -port " + getPort(options) + " -host " + getHost(options) +
                  "> /dev/null 2>&1")
    return r == 0
    
def main():
    parser = OptionParser()
    parser.add_option("--experiment", dest="experimentFile", 
                       help="The file containing a link to the experiment parameters")
    
    parser.add_option("--cladeSize", dest="cladeSize", type="int", 
                      help="Max number of sequences to align at a time", default=2)
    
    options, args = parser.parse_args()
    #setLoggingFromOptions(options)
    
    options.experimentFile = ET.parse(options.experimentFile).getroot()
    #Get the database string
    options.cactusDiskDatabaseString = ET.tostring(options.experimentFile.find("cactus_disk").find("st_kv_database_conf"))
    #Get the species tree
    options.speciesTree = options.experimentFile.attrib["species_tree"]
   
    
    root = newickTreeParser(options.speciesTree)
    nameUnlabeledInternalNodes(root)
    
    clades = []
    recurseClades(root, options, clades)
    print ""
    if testServerExists(options):
        print "WARNING: ktserver already running on specified port and host\n"
    print "ktserver" + " -port " + getPort(options) + " -host " + getHost(options) \
    + " -ls -tout 200000 -th 64 " + getDbPaths(clades, options) + "\n"
    

if __name__ == '__main__':    
    #from cactusTools.progressive.prepKtserverDbs import *
    main()