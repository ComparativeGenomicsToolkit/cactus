#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Wrap a tree and add some simple partitioning functionality 

"""

import os
import xml.etree.ElementTree as ET
import sys
import math
import copy

from optparse import OptionParser

from sonLib.bioio import printBinaryTree
from sonLib.bioio import newickTreeParser


class MultiCactusTree:
    def __init__(self, tree, subtreeSize = 2):
        self.tree = tree
        self.subtreeSize = subtreeSize
        # maps node id to node for all subtree roots
        self.subtreeRoots = dict()
        # maps node to the subtree root that contains it
        self.nearestRoots = dict()
    
    # fill in unlabeled node ids with a breadth-first
    # traversal numbering from the root
    def nameUnlabeledInternalNodes(self, prefix = "Anc", startIdx = 0):
        bfQueue = [self.tree]
        count = startIdx
        while bfQueue:
            node = bfQueue.pop(0)
            if node is not None and node.iD is None:
                if node.internal:
                    node.iD = prefix + str(count)
                    bfQueue.append(node.left)
                    bfQueue.append(node.right)
                count += 1
    
    # identify roots of subclades in the tree and 
    # add them to the self.claderoots dicitonary
    def computeSubtreeRoots(self):
        def computeSubtreeRootsRecursive(node):
            assert node.iD is not None
            self.subtreeRoots[node.iD] = node
            leaves = self.getSubtreeLeaves(node)
            for subtreeLeaf in leaves:
                if subtreeLeaf.internal:
                    computeSubtreeRootsRecursive(subtreeLeaf)
        self.subtreeRoots = dict()
        computeSubtreeRootsRecursive(self.tree)
    
    # map a node to the closest subtree root that contains it
    def computeNearestRoots(self):
        def computeNearestRootsRecursive(node, root):
            if node:
                self.nerestRoots[node] = root
                if node in self.subtreeRoots:
                    root = node
                computeNearestRootsRecursive(node.left, root)
                computeNearesetRootsRecursive(node.right, root)
        self.nearestRoots = dict()
        computeSubtreeRootsRecursive(self.tree, None)
        
    # blindly read in the roots from a given map or set 
    def assignSubtreeRoots(self, roots):
        def assignSubtreeRootsRecursive(roots, node):
            if node is None:
                return
            if node == self.tree:
                self.subtreeRoots = dict()
            if node.iD and node.iD in roots:
                self.subtreeRoots[node.iD] = node
            assignSubtreeRootsRecursive(roots, node.left)
            assignSubtreeRootsRecursive(roots, node.right)
        assignSubtreeRootsRecursive(roots, self.tree)
          
    # copy a subtree rooted at node with given name
    def extractSubTree(self, name):
        def copyRecursive(root, node):
            if node is None:
                return None
            cpy = copy.deepcopy(node)
            if (node != root and node.iD in self.subtreeRoots) \
                or not node.internal:
                cpy.left = None
                cpy.right = None
                cpy.internal = False
            else:
                cpy.left = copyRecursive(root, node.left)
                cpy.right = copyRecursive(root, node.right)
            return cpy
        root = self.subtreeRoots[name]
        cpy = copyRecursive(root, root)
        assert cpy is not None
        return cpy
        
    # find the root of the subtree containing the given node
    # as leaf (slowly.. for nwo)
    def getSubtreeRoot(self, node):
        def getSubtreeRootRecursive(node, cur):
            assert node != self.tree
            if cur is None:
                return None
            leaves = getSubtreeLeaves(cur)
            for leaf in leaves:
                if leaf.iD == node.iD:
                    return cur
                par = getSubtreeRootRecursive(node, leaf)
                if par is not None:
                    return par
            return None
        return getSubtreeRootRecursive(node, self.tree)
        
    # find the leaves of af subtree, subject to 
    # 1) number of leaves maximal but less than self.subtreeSize
    # 2) if a node is returned, its sibling must me as well
    def getSubtreeLeaves(self, node):
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
        nextLevel = allChildren([node])
        while (len(nextLevel) <= self.subtreeSize and len(nextLevel) > len(curLevel)):
            curLevel = nextLevel
            nextLevel = allChildren(curLevel)    
        return curLevel
    
    # check if given tree is a subtree of self.tree        
    def checkSubtree(self, tree):
        assert tree.iD is not None
        subtree = self.extractSubTree(tree.iD)
        
        def treeComp(node1, node2, outgroup):
            if node1 is None and node2 is None:
                return True
            elif node1 is not None and node2 is not None:
                return node1.iD == node2.iD and\
                     treeComp(node1.left, node2.left) and\
                     treeComp(node1.right, node2.right)    
            else:
                return False
        
        assert treeComp(subtree, tree, outgroup)    


    