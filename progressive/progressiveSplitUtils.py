#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
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

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from sonLib.bioio import getLogLevelString
from sonLib.tree import getBinaryTreeNodes

from jobTree.src.bioio import getLogLevelString
from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions


# Extract a subtree with given root and leaves
def copySubtree(root, leaves):
    def cpyChildren(nodeCpy, node, leaves):
        # Note: slow list search won't scale to huge clades
        if node not in leaves:
            if node.left is not None:
                nodeCpy.left = copy.deepcopy(node.left)
                cpyChildren(nodeCpy.left, node.left, leaves)
            else:
                nodeCpy.left = None
            if node.right is not None:
                nodeCpy.right = copy.deepcopy(node.right)
                cpyChildren(nodeCpy.right, node.right, leaves)
            else:
                nodeCpy.right = None
        else:
            nodeCpy.left = None
            nodeCpy.right = None
            nodeCpy.internal = False
    
    rootCpy = copy.deepcopy(root)
    cpyChildren(rootCpy, root, leaves)
    return rootCpy

# Get the list nodes corresponding to the leaves of a clade
def getCladeLeaves(root, options):
    # Get children of a set of nodes
    def allChildren(nodes):
        children = []
        for i in nodes:
            if i.left is not None:
                children.append(i.left)
            if i.right is not None:    
                children.append(i.right)
            elif i.left is None:
                children.append(i)
        return children
    
    curLevel = []
    nextLevel = allChildren([root])
    while (len(nextLevel) <= options.cladeSize and len(nextLevel) > len(curLevel)):
        curLevel = nextLevel
        nextLevel = allChildren(curLevel)
        
    return curLevel

# Assign name to unlabeled internal nodes
def nameUnlabeledInternalNodes(root, prefix = "Anc", startIdx = 0):
    bfQueue = [root]
    count = startIdx
    while bfQueue:
        node = bfQueue.pop(0)
        if node is not None:
            if node.left is not None or node.right is not None:
                if node.iD is None:
                    node.iD = "Anc" + str(count)
                    bfQueue.append(node.left)
                    bfQueue.append(node.right)
                count += 1
            # tack on a dot so mafJoin will work
            # will have to do a review down the road for multi-chromosome cases 
            # down the road...    
            # if os.path.splitext(node.iD)[1] == '':
            #    node.iD += ".0"

def createProgWorkDir(options, root):
    workDir = getWorkingDir(options)
    if not os.access(workDir, os.F_OK):
        os.makedirs(workDir)
    assert os.access(workDir, os.W_OK)
    open(getSpeciesTreePath(options), "w").write(printBinaryTree(root, True))
    
# the following methods describe the progressive directory structure, all relative
# to the database_dir specified in the original experiment file
def getWorkingDir(options):
    dbConfElem = options.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    baseDbDir = dbTypeElem.attrib["database_dir"]
    return baseDbDir

def getCladeWorkingDir(root, options):
    return getWorkingDir(options) + "/" + root.iD

def getCladeReferenceName(root):
    name, ext = os.path.splitext(root.iD)
    return name + "_reference" + ext

def getCladeReferencePath(root, options):
    return getCladeWorkingDir(root, options) + "/" + getCladeReferenceName(root) + ".fa"

def getCladeExperimentPath(root, options):
    return getCladeWorkingDir(root, options) + "/" + root.iD + "_experiment.xml"

def getCladeDatabaseDir(root, options):
    return getCladeWorkingDir(root, options) + "/" + root.iD + "_DB"

def getCladeDatabaseName(root, options):
    return "data_" + root.iD + ".kch"

def getCladeMAFPath(root, options, isJoined):
    path = getCladeWorkingDir(root, options) + "/" + root.iD
    if not isJoined:
        path += "_local"
    return path + ".maf"

def getCladeMAFJoinTempPath(root, leaf, options):
    return getCladeWorkingDir(root, options) + "/" + root.iD + "_" + leaf.iD + ".maf"

def getCladeConfigPath(root, options):
    return getCladeWorkingDir(root, options) + "/" + root.iD + "_config.xml"

def getSpeciesTreePath(options):
    return getWorkingDir(options) + "/species.tree"
    
# Create a mapping between node IDs and sequence files
def createSeqeunceMap(tree, options, sequences):
    dfStack = [tree]
    nameIterator = iter(sequences)
    lookup = dict()
    while dfStack:
        node = dfStack.pop(len(dfStack)-1)
        if node is not None:
            if node.left is None and node.right is None:
                lookup[node.iD] = nameIterator.next()
            else:
                lookup[node.iD] = getCladeReferencePath(node, options)
                dfStack.append(node.right)
                dfStack.append(node.left)
    
    l = []
    getBinaryTreeNodes(tree, l)    
    assert len(lookup) == len(l)
    return lookup

def getCladeSequences(leaves, options):
    return map(lambda x: options.lookup[x.iD], leaves)

# Create new options object for the subproblem, specifically overriding
# the original species tree and diskdatabasestring, experimentfile
def createCladeOptions(root, leaves, options):  
    cladeName = root.iD
    assert cladeName is not None
    
    #copy the options and the clade tree from input into some new options
    cladeOptions = copy.deepcopy(options)
    cladeTree = copySubtree(root, leaves)
    cladeOptions.speciesTree = printBinaryTree(cladeTree, True)
    
    #modifiy sequences attribute in new XML
    assert options.lookup is not None
    sequences = map(lambda x: options.lookup[x.iD], leaves)
    sequences = reduce(lambda x,y: x + " " + y, sequences)
    
    cladeOptions.experimentFile.attrib["sequences"] = sequences
    cladeOptions.experimentFile.attrib["species_tree"] = cladeOptions.speciesTree    
     
    #modify database_dir attribute in new XML
    dbConfElem = cladeOptions.experimentFile.find("cactus_disk").find("st_kv_database_conf")
    dbTypeElem = dbConfElem.find(dbConfElem.attrib["type"])
    dbTypeElem.attrib["database_dir"] = getCladeDatabaseDir(root, options)
    dbTypeElem.attrib["database_name"] = getCladeDatabaseName(root, options)
    cladeOptions.cactusDiskDatabaseString = ET.tostring(dbConfElem)
    
    #modify the config attribute in new XML
    cladeOptions.experimentFile.attrib["config"]  = getCladeConfigPath(root, options)
    cladeOptions.config.find("reference").attrib["reference"] = getCladeReferenceName(root)     
    
    #add field for server process for auto ktserver spawning
    cladeOptions.serverProcess = None
    
    return cladeOptions 
    
# Create the directory and xml for cactus on clade   
def createCladeFileStructure(root, leaves, options, cladeOptions):
    cladeName = root.iD
    #Create the work directory
    workDir = getCladeWorkingDir(root, options)
    if not os.access(workDir, os.F_OK):
        os.makedirs(workDir)
    assert os.access(workDir, os.W_OK)
    dbDir = getCladeDatabaseDir(root, options)
    if not os.access(dbDir, os.F_OK):
        os.makedirs(dbDir)
    assert os.access(dbDir, os.W_OK)
    
    #Copy over the config.xml file, making sure to update the reference attribute
    ET.ElementTree(cladeOptions.config).write(getCladeConfigPath(root, options))
    
    #Write out a local experimentfile xml file
    ET.ElementTree(cladeOptions.experimentFile).write(getCladeExperimentPath(root, options))
    
def getReferenceSeqOptions(root, options, cladeOptions):
    optString = "--logLevel " + getLogLevelString()
    optString += " --cactusDisk " + "'" + cladeOptions.cactusDiskDatabaseString + "'"
    optString += " --flowerName 0"
    optString += " --referenceEventString " + getCladeReferenceName(root)
    optString += " --outputFile " + getCladeReferencePath(root, options)
    return optString

# Make a string of options for running cactus_MAFGenerator on the given clade
def getMAFGeneratorOptions(root, options, cladeOptions):
    optString = "--logLevel " + getLogLevelString()
    optString += " --cactusDisk " + "'" + cladeOptions.cactusDiskDatabaseString + "'"
    optString += " --flowerName 0"
    optString += " --outputFile " + getCladeMAFPath(root, options, False)
    # reference string should be picked up automatically
    return optString
    
def getMAFJoinOptions(root, leaf, options, cladeOptions, rootIsTreeMAF):
    def fn(leaf):
        leaves = getCladeLeaves(leaf, options)
        for i in leaves:
            if i.left is not None or i.right is not None:
                return True
        return False
    leafIsTreeMAF = fn(leaf)
    optString = " -maxBlkWidth=10000 -maxInputBlkWidth=1000"
    if not rootIsTreeMAF:
        optString += " -treelessRoot1=" + getCladeReferenceName(root) 
    if not leafIsTreeMAF:
        optString += " -treelessRoot2=" + getCladeReferenceName(leaf)
    optString += " " + getCladeReferenceName(leaf)
    
    optString += " " + getCladeMAFPath(root, options, rootIsTreeMAF)   
    optString += " " + getCladeMAFPath(leaf, options, leafIsTreeMAF)
    optString += " " + getCladeMAFJoinTempPath(root, leaf, options)
    return optString
    
#cheap temporary hack copy pasted from itnerent
def addDotsToMAF(root, options):
    file = getCladeMAFPath(root, options, False)
     #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file)
    for line in old_file:
        linelist = line.split()
        if len(linelist) > 0:
            if linelist[0] == "s":
                line = line.replace(linelist[1], linelist[1] + ".1", 1)
        new_file.write(line)
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file)
    #Move new file
    move(abs_path, file)
               
    
    