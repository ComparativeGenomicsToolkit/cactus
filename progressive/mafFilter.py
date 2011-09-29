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

def removeOutgroupLeafFromTree(tree, event):
    for leaf in tree.getLeaves():
        if tree.getName(leaf) == event:
            tree.nxDg.remove_node(leaf)
            return
    assert False
    
# quick and dirty friday afternoon! 
def removeOutgroupFromMaf(mafPath, event):
    if event is None or event == "":
        return
    tempPath = "%s.og" % mafPath
    inFile = open(mafPath)
    outFile = open(tempPath, 'w')
    num = 1
    newick = NXNewick() 
    for line in inFile:
        if num == 2:
            treeString = line.split()[2]
            tree = newick.parseString(treeString)
            removeOutgroupLeafFromTree(tree, event)
            line = line.replace(treeString, newick.writeString(tree))            
        if line[0] == '#' or line.find(event) < 0:
            outFile.write(line)
        num += 1
    inFile.close()
    outFile.close()
    os.rename(tempPath, mafPath)