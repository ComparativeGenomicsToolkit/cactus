#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Remove outgroup from a maf file. 

"""

import os
import sys

from sonLib.nxtree import NXTree
from sonLib.nxnewick import NXNewick
from cactus.progressive.multiCactusTree import MultiCactusTree

def removeOutgroupLeafFromTree(tree, event):
    for leaf in tree.getLeaves():
        if tree.getName(leaf).find(event) == 0:
            tree.nxDg.remove_node(leaf)
            return
    assert False

# cactus_createMulti.... guarantees one event name is never
# a prefix of another EXCEPT in the case of the _self. 
# so any sequence with the given outgroup name as a prefix
# except those with outgroup_self as the prefix belong
# to the outgroup and need to be filtered
def eventInLine(name, line):
    tokens = line.split()
    if line[0] == 's' and len(line) >= 2:
        word = tokens[1]
        if word.find(name) == 0 and \
        word.find(name + MultiCactusTree.self_suffix) != 0:
            return True
    return False
    
# quick and dirty friday afternoon! 
def removeOutgroupFromMaf(mafPath, events):
    if event is None or len(events) == 0:
        return
    # only support 1 outgroup for now
    assert len(events) == 1
    event = events[0].split('.')[0]
    tempPath = "%s.og" % mafPath
    inFile = open(mafPath)
    outFile = open(tempPath, 'w')
    num = 1
    newick = NXNewick() 
    for line in inFile:
        if num == 2:
            treeString = line.split()[2]
            tree = MultiCactusTree(newick.parseString(treeString))
            removeOutgroupLeafFromTree(tree, event)
            line = line.replace(treeString, newick.writeString(tree))            
        if line[0] == '#' or not eventInLine(event, line):
            outFile.write(line)
        num += 1
    inFile.close()
    outFile.close()
    os.rename(tempPath, mafPath)
