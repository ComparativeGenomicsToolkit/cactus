#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Compute outgroups for subtrees. Currently uses simple heuristic:
for a leaf: use the nearest other leaf
otherwise node x: find nearest node at height lower than x that
is not in a subtree of x.  

"""

import os
import xml.etree.ElementTree as ET
import sys
import math
import copy

from optparse import OptionParser

from sonLib.bioio import printBinaryTree
from sonLib.bioio import newickTreeParser
from sonLib.tree import binaryTree_depthFirstNumbers
from sonLib.tree import getDistanceMatrix

class OutgroupFinder:
    def __init__(self, mcTree):
        self.mcTree = mcTree
        self.ogMap = dict()
    
    def computeDistances(self):
        binaryTree_depthFirstNumbers(self.mcTree.tree)
        self.dm = getDistanceMatrix(self.mcTree.tree)
        # throw in some topological sort information 
        def fn(node, idx):
            if node:
                node.traversalID.sort = idx
                fn(node.left, idx + 1)
                fn(node.right, idx + 1)
        fn(self.mcTree.tree, 0) 
    
    def computeHeights(self):
        self.heights = dict()
        def fn(node):
            height = 0
            if node.internal:
               height = 1 + max(fn(node.left), fn(node.right))
            self.heights[node] = height
            return height
        fn(self.mcTree.tree)
    
    def computeBelow(self):
        def computeBelowRecursive(node, aboveStack):
            if node:
                for aboveNode in aboveStack:
                    self.below.add((node, aboveNode))
                aboveStack.append(node)
                computeBelowRecursive(node.left, aboveStack)
                computeBelowRecursive(node.right, aboveStack)
                aboveStack.pop()
        self.below = set()
        computeBelowRecursive(self.mcTree.tree, [])
        
    # really naive inefficient placeholder function
    # assigns outgroup as closest node that is at least one 
    # level lower in the tree (or the same for leaf).
    def findAllNearestBelow(self):
        self.computeHeights()
        self.computeDistances()
        self.computeBelow()
        byHeight = dict()
        
        # sort nodes by height
        for name, node in self.mcTree.subtreeRoots.items():
            h = self.heights[node]
            if h in byHeight:
                byHeight[h].append(node)
            else:
                byHeight[h] = [node]
        byHeight[0] = []
        for node, height in self.heights.items():
            if height == 0:
                byHeight[0].append(node)
        
        for name, node in self.mcTree.subtreeRoots.items():
            h = self.heights[node]
            maxOh = max(0, h-1)
            minD = sys.maxint
            outgroup = None
            id = node.traversalID.mid
            for oh in range(maxOh, -1, -1):                
                for otherNode in byHeight[oh]:
                    if otherNode != node and \
                    (otherNode, node) not in self.below:
                        d = self.dm[(otherNode.traversalID.mid, id)]                     
                        if d < minD:
                            minD = d
                            outgroup = otherNode
            assert name not in self.ogMap
            if outgroup:
                distance = self.dm[(outgroup.traversalID.mid, id)]
                self.ogMap[name] = (outgroup.iD, distance)
            else:
                assert node == self.mcTree.tree 
        
        
    
        