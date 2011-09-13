#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Remove outgroup from a maf file. 

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

from optparse import OptionParser

from sonLib.bioio import printBinaryTree
from sonLib.bioio import newickTreeParser
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper

def removeOutgroupLeafFromTree(tree, event):
    if tree:
        if tree.left and tree.left.iD == event and not tree.left.internal:
            tree.left = None            
        elif tree.right and tree.right.iD == event and not tree.right.internal:
            tree.right = None
        else:
            removeOutgroupLeafFromTree(tree.left, event)
            removeOutgroupLeafFromTree(tree.right, event)
    
# quick and dirty friday afternoon! 
def removeOutgroupFromMaf(mafPath, event):
    if event is None or event == "":
        return
    tempPath = "%s.og" % mafPath
    inFile = open(mafPath)
    outFile = open(tempPath, 'w')
    num = 1 
    for line in inFile:
        if num == 2:
            treeString = line.split()[2]
            tree = newickTreeParser(treeString, reportUnaryNodes=True)
            removeOutgroupLeafFromTree(tree, event)
            line = line.replace(treeString, printBinaryTree(tree, True))            
        if line[0] == '#' or line.find(event) < 0:
            outFile.write(line)
        num += 1
    inFile.close()
    outFile.close()
    os.rename(tempPath, mafPath)