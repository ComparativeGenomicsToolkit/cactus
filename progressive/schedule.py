#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Schedule jobs based on the dependencies between their associated
events.  These dependencies are based 1, on the phylogeny and 2,
on outgroups.  The input is a DAG, the output is a TREE 

"""

import os
import xml.etree.ElementTree as ET
import sys
import math
import copy
import networkx as NX

from optparse import OptionParser

from sonLib.bioio import printBinaryTree
from sonLib.bioio import newickTreeParser
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper

class Schedule:
    def __init__(self):
        self.inGraph = None
        self.depTree = None
        self.heights = None
        self.rootName = None
    
    # read the experiments, compute the dependency dag
    def compute(self, mcProject):
        self.inGraph = NX.DiGraph()
        for name, expPath in mcProject.expMap.items():
            exp = ExperimentWrapper(ET.parse(expPath).getroot())
            self.loadBinaryTree(exp.getTree(), name) 
        self.checkInGraph()
        self.rootName = mcProject.mcTree.tree.iD
        self.depTree = NX.DiGraph(self.inGraph)
        self.heights = dict()
        self.computeHeights(self.rootName)
        self.makeTree()
        self.checkDepTree()
        self.inGraph = None
        self.heights = None
        
     # for a given event name, get the names of all the 
    # events that are directly dependent on it in the 
    # schedule
    def deps(self, name):
        assert name in self.depTree
        
        edges = self.depTree.out_edges(name)
        depList = []
        for edge in edges: 
            depList.append(edge[1])
        return depList
    
    # write to graphviz dot file
    def writeToFile(self, path):
        NX.drawing.nx_agraph.write_dot(self.depTree, path)
    
    # read from graphviz dot file
    def readFromFile(self, path):
        self.depTree = NX.read_edgelist(path)
        self.checkDepTree()
        
    # PRIVATE BELOW:
    
    # load only the root and leaves
    # if there's an outgroup, handle it and flag its dependency
    # with info='outgroup'
    def loadBinaryTree(self, tree, rootName):
        def getLeaves(node):
            leaves = []
            if not node.internal:
                assert node.left is None and node.right is None
                leaves.append(node)
            elif node.left:
                leftLeaves = getLeaves(node.left)
                for i in leftLeaves:
                    leaves.append(i)
            if node.right:
                rightLeaves = getLeaves(node.right)
                for i in rightLeaves:
                    leaves.append(i)
            return leaves
        
        # if there's an outgroup, the event is one node below the 
        # tree root (don't think we need)
        root = tree
        og = None
        if root.iD != rootName:
            assert tree.left.iD == rootName or tree.right.iD == rootName
            if tree.left.iD == rootName:
                root = tree.left
                og = tree.right
            else:
                root = tree.right
                og = tree.left    
        assert og is None or not og.internal
        
        leaves = getLeaves(root)
        for leaf in leaves:
            assert len(leaf.iD) > 0
            self.inGraph.add_edge(root.iD, leaf.iD)
        if og:
            self.inGraph.add_edge(root.iD, og.iD, info='outgroup')
    
    # compute the children below a node in the phylogenetic tree
    # note we explicitly skip outgroups
    def phyloChildren(self, name):
        edges = self.inGraph.out_edges(name)
        children = []
        for edge in edges:
            if 'info' not in self.inGraph.adj[edge[0]][edge[1]] or\
             self.inGraph.adj[edge[0]][edge[1]]['info'] != 'outgroup':
                children.append(edge[1])
        return children
    
    # repeatedly break undirected cycles, moving them "up" the tree.  
    # claim: this should converge
    # note: can be sped up by keeping a priority queue of nodes with 
    # indegree > 1
    def makeTree(self):
        N = self.depTree.size() + 1
    
        for i in range(0, N):
            self.dirty = False
            self.breakUndirectedCycles(self.rootName)
            if not self.dirty:
                break
            else:
                assert i < N-1
        
    # jobtree doesn't support two nodes sharing the same dependency
    # so we make sure to remove such cases by recursively pushing
    # up undirected cycles starting at the leaves
    def breakUndirectedCycles(self, rootName):
        children = self.phyloChildren(rootName)
        children = sorted(children, key = lambda x : -self.heights[x])
        for i in children:
            self.breakUndirectedCycles(i)
        inEdges = self.depTree.in_edges([rootName])
        if len(inEdges) > 0:
            parents = []
            for inEdge in inEdges:
                parents.append(inEdge[0])
            parents = sorted(parents, key = lambda x: self.heights[x])            
            keeper = parents[0]
            for inEdge in inEdges:
                if inEdge[0] != keeper:
                    print "remove %s %s" % (inEdge[0], inEdge[1])
                    self.depTree.remove_edge(inEdge[0], inEdge[1])
                    if not self.depTree.has_edge(inEdge[0], keeper):
                        print "add %s %s" % (inEdge[0], keeper)
                        self.depTree.add_edge(inEdge[0], keeper)
                        if self.depTree.in_degree(keeper) > 1:
                            self.dirty = True
        
        
    # fill up the map htable with the height
    # of each node in the phylogenetic tree (not considering
    # outgroup edges
    def computeHeights(self, rootName):
        #leaves = self.phyloChildren(rootName)
        leaves = self.deps(rootName)
        if len(leaves) == 0:
            self.heights[rootName] = 0
        else:
            leafHeights = []
            for i in leaves:
                self.computeHeights(i)
                leafHeights.append(self.heights[i])
            self.heights[rootName] = max(leafHeights) + 1
    
    def checkInGraph(self):
        assert NX.is_directed_acyclic_graph(self.inGraph)
        
    def checkDepTree(self):
        def checkTree(node):
            assert self.depTree.in_degree(node) <= 1
            for i in self.deps(node):
                checkTree(i)
        checkTree(self.rootName)
        assert NX.is_directed_acyclic_graph(self.depTree)
        
def main():
    usage = "usage: %prog <project> <output graphviz .dot file>"
    description = "TEST: create schedule from project file"
    parser = OptionParser(usage=usage, description=description)
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    proj = MultiCactusProject()
    proj.readXML(args[0])
    schedule = Schedule()
    schedule.compute(proj)
    schedule.writeToFile(args[1])    

if __name__ == '__main__':    
    main()