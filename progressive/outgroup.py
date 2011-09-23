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
import networkx as NX

from optparse import OptionParser

from sonLib.bioio import printBinaryTree
from sonLib.tree import binaryTree_depthFirstNumbers
from sonLib.tree import getDistanceMatrix
from sonLib.tree import BinaryTree

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree

class GreedyOutgroup:
    def __init__(self):
        self.dag = None
        self.dm = None
        self.dmDirected = None
        self.rootName = None
        self.ogMap = None
        
    # add edges from sonlib tree to self.dag
    # compute self.dm: an undirected distance matrix
    def importTree(self, mcTree):
        def importNode(node):        
            if node and node.left:
                w = node.left.distance
                self.dag.add_edge(node.iD, node.left.iD, weight=w)
                importNode(node.left)
            if node and node.right:
                w = node.right.distance
                self.dag.add_edge(node.iD, node.right.iD, weight=w)
                importNode(node.right)
        self.dag = NX.DiGraph()
        importNode(mcTree.tree)
        self.rootName = mcTree.tree.iD
        self.stripNonEvents(self.rootName, mcTree.subtreeRoots)
        self.dmDirected = NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(self.dag)
        graph = NX.Graph(self.dag)
        self.dm = NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(graph)
 
    # get rid of any node that's not an event
    def stripNonEvents(self, name, subtreeRoots):
        children = []
        for outEdge in self.dag.out_edges(name):
            children.append(outEdge[1])
        parentCount = len(self.dag.in_edges(name))
        assert parentCount <= 1
        parent = None
        if parentCount == 1:
            parent = self.dag.in_edges[0][0]
        if name not in subtreeRoots and len(children) >= 1:
            assert parentCount == 1
            self.dag.remove_node(name)
            for child in children:
                self.dag.add_edge(parent, child)
                self.stripNonEvents(child, subtreeRoots)
    
    # are source and sink son same path to to root? if so,
    # they shouldn't be outgroups of each other.             
    def onSamePath(self, source, sink):
        if source in self.dmDirected:
            if sink in self.dmDirected[source]:
                return True
        if sink in self.dmDirected:
            if source in self.dmDirected[sink]:
                return True 
        return False
    
    # greedily assign closest possible valid outgroup
    # all outgroups are stored in self.ogMap
    # edges between leaves ARE NOT kept in the dag         
    def greedy(self):
        orderedPairs = []
        for source, sinks in self.dm.items():
            for sink, dist in sinks.items():
                if source != self.rootName and sink != self.rootName:
                    orderedPairs.append((dist, (source, sink)))
        orderedPairs.sort(key = lambda x: x[0])
        finished = set()
        self.ogMap = dict()
        
        for candidate in orderedPairs:
            source = candidate[1][0]
            sink = candidate[1][1]
            dist = candidate[0]
            
            # skip leaves (as sources)
            if len(self.dag.out_edges(source)) == 0:
                finished.add(source)
    
            if source not in finished and \
            not self.onSamePath(source, sink):
                self.dag.add_edge(source, sink, weight=dist, info='outgroup')
                if NX.is_directed_acyclic_graph(self.dag):
                    finished.add(source)
                    self.ogMap[source] = (sink, dist)                    
                else:
                    self.dag.remove_edge(source, sink)
    
             

def main():
    usage = "usage: %prog <project> <output graphviz .dot file>"
    description = "TEST: draw the outgroup DAG"
    parser = OptionParser(usage=usage, description=description)
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    proj = MultiCactusProject()
    proj.readXML(args[0])
    outgroup = GreedyOutgroup()
    outgroup.importTree(proj.mcTree)
    outgroup.greedy()
    NX.drawing.nx_agraph.write_dot(outgroup.dag, args[1])
    return 0

if __name__ == '__main__':    
    main()