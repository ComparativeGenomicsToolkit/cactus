#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Wrap a tree and add some simple partitioning functionality.  note
that nameUnlabeledInternalNodes() and computeSubtreeRoots() need to be
called (in that order) for anything to work...  

"""

import os
import xml.etree.ElementTree as ET
import sys
import math
import copy

from optparse import OptionParser

from sonLib.nxtree import NXTree
from sonLib.nxnewick import NXNewick

class MultiCactusTree(NXTree):
    self_suffix = "_self"
    def __init__(self, tree = None, subtreeSize = 2):
        if isinstance(tree, NXTree):
            NXTree.__init__(self, tree.nxDg)
        else:   
            NXTree.__init__(self, tree)
        # ids of all subtree roots for fast checking
        self.subtreeRoots = set()
        # map of names to node ids
        self.nameToId = dict()
        for node in self.breadthFirstTraversal():
            if self.hasName(node):
                self.nameToId[self.getName(node)] = node
        # size a subtree (in number of leaves)
        self.subtreeSize = subtreeSize
        
    # fill in unlabeled node ids with a breadth-first
    # traversal numbering from the root
    def nameUnlabeledInternalNodes(self, prefix = "Anc", startIdx = 0):
        count = startIdx
        numInternal = 0
        width = 0
        for node in self.breadthFirstTraversal():
            if self.isLeaf(node) is False:
                numInternal += 1
        if numInternal > 0:
            width = int(math.log10(numInternal)) + 1
        for node in self.breadthFirstTraversal():
            if not self.isLeaf(node) and not self.hasName(node):
                self.setName(node, "%s%s" % (prefix, str(count).zfill(width)))
                count += 1
            self.nameToId[self.getName(node)] = node
    
    # identify roots of subclades in the tree and 
    # add them to the self.claderoots dicitonary
    def computeSubtreeRoots(self, root = None):
        if root is None:
            root = self.rootId
            self.subtreeRoots = set()
        assert root not in self.subtreeRoots
        self.subtreeRoots.add(root)
        leaves = self.getSubtreeLeaves(root)    
        for subtreeLeaf in leaves:
            if not self.isLeaf(subtreeLeaf):
                self.computeSubtreeRoots(subtreeLeaf)
        
    # blindly read in the roots from given list of names 
    def assignSubtreeRootNames(self, rootNames):
        self.subtreeRoots = set()
        for node in self.breadthFirstTraversal():
            if self.getName(node) in rootNames:
                self.subtreeRoots.add(node)
            self.nameToId[self.getName(node)] = node
                
    def getSubtreeRootNames(self):
        return [self.getName(x) for x in self.subtreeRoots]
    
    # generate eall nodes beneath (and including) given
    # root
    def traverseSubtree(self, root, node):
        yield node
        if node == root or node not in self.subtreeRoots:
            for child in self.getChildren(node):
                for i in self.traverseSubtree(root, child):
                    yield i
            
    # copy a subtree rooted at node with given name
    def extractSubTree(self, name):
        root = self.nameToId[name]
        subtree = [i for i in self.traverseSubtree(root, root)]
        cpy = self.nxDg.subgraph(subtree).copy()
        mcCpy = MultiCactusTree(cpy, 2)
        mcCpy.assignSubtreeRootNames(self.getSubtreeRootNames())
        return mcCpy
        
    # find the root of the subtree containing the given node
    # as leaf (slowly.. for nwo) (returns event name, not node id)
    def getSubtreeRoot(self, name):
        assert len(self.nameToId) > 0
        node = self.nameToId[name]
        if node == self.rootId:
            return name
        parent = self.getParent(node)
        while parent is not None:
            if parent in self.subtreeRoots:
                return self.getName(parent)
            parent = self.getParent(parent)
        return None
        
    # find the leaves of af subtree, subject to 
    # 1) number of leaves maximal but less than self.subtreeSize
    # 2) if a node is returned, its sibling must me as well
    def getSubtreeLeaves(self, node):
        curLevel = []
        nextLevel = self.getChildren(node)
        while True:
            curLevel = nextLevel
            curLevel.sort()
            nextLevel = []
            for node in curLevel:
                if self.isLeaf(node):
                    nextLevel.append(node)
                else:
                    nextLevel += self.getChildren(node)
            nextLevel.sort()
            if len(nextLevel) > self.subtreeSize or curLevel == nextLevel:
                break
        return curLevel
    
    # safe id to insert is current max + 1
    def getNextIndex(self):
        return sorted([i for i in self.breadthFirstTraversal()])[-1] + 1
    
    # insert a new node above a specified node in the tree
    def insertAbove(self, node, newNode, newName = "", newWeight= None):
        parent = self.getParent(node)
        if parent is not None:
            oldWeight = self.getWeight(parent, node)
            self.nxDg.remove_edge(parent, node)
            self.nxDg.add_edge(parent, newNode)
            if oldWeight is not None:
                self.setWeight(parent, newNode, oldWeight)
        else:
            assert node == self.rootId
            self.rootId = newNode
        self.nxDg.add_node(newNode)
        self.setName(newNode, newName)
        self.nxDg.add_edge(newNode, node)
        if newWeight is not None:
            self.setWeight(newNode, node, newWeight)
        if len(newName) > 0:
            self.nameToId[newName] = newNode
            
    # insert a node with id (name_self) directly above 
    # every node in the tree
    # should be run after subtreeroots are computed (otherwise
    # won't work
    def addSelfEdges(self):
        nextIndex = self.getNextIndex()
        traversal = [i for i in self.breadthFirstTraversal()]
        for node in traversal:
            if (node in self.subtreeRoots or self.isLeaf(node)) and\
            node != self.rootId:
                newNode = nextIndex
                nextIndex += 1
                parent = self.getParent(node)
                weight = None
                if parent is not None:
                    weight = self.getWeight(parent, node)
                    assert weight is not None
                assert self.self_suffix not in self.getName(node)
                newName = self.getName(node) + self.self_suffix
                self.insertAbove(node, newNode, newName, weight)
                self.subtreeRoots.add(newNode) 
    
    # tack an outgroup onto the root
    # if root is a leaf, we make a new root above. 
    def addOutgroup(self, ogName, distance):
        assert ogName not in self.nameToId
        if self.isLeaf(self.rootId):
            newNode = self.getNextIndex()
            self.insertAbove(self.rootId, newNode, "", distance / 2)
            distance = distance / 2
        newNode = self.getNextIndex()
        self.nxDg.add_edge(self.rootId, newNode )
        self.setName(newNode, ogName)
        self.nameToId[ogName] = newNode
        self.setWeight(self.rootId, newNode, distance)

    # make an attempt at keeping nameToId up to date
    def setName(self, node, name):
        super(MultiCactusTree, self).setName(node, name)
        self.nameToId[name] = node
    
    # map from event name to networkx node ID
    def getNodeId(self, name):
        return self.nameToId[name]
