#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Compute outgroups for subtrees by greedily finding
the nearest valid candidate.  has option to only assign
leaves as outgroups

"""

import os
import xml.etree.ElementTree as ET
import sys
import math
import copy
import networkx as NX

from optparse import OptionParser

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree

class GreedyOutgroup:
    def __init__(self):
        self.dag = None
        self.dm = None
        self.dmDirected = None
        self.root = None
        self.ogMap = None
        self.mcTree = None
        
    # add edges from sonlib tree to self.dag
    # compute self.dm: an undirected distance matrix
    def importTree(self, mcTree):
        self.mcTree = mcTree
        self.dag = mcTree.nxDg.copy()
        self.root = mcTree.rootId
        self.stripNonEvents(self.root, mcTree.subtreeRoots)
        self.dmDirected = NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(self.dag)
        graph = NX.Graph(self.dag)
        self.dm = NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(graph)
 
    # get rid of any node that's not an event
    def stripNonEvents(self, id, subtreeRoots):
        children = []
        for outEdge in self.dag.out_edges(id):
            children.append(outEdge[1])
        parentCount = len(self.dag.in_edges(id))
        assert parentCount <= 1
        parent = None
        if parentCount == 1:
            parent = self.dag.in_edges[0][0]
        if id not in subtreeRoots and len(children) >= 1:
            assert parentCount == 1
            self.dag.remove_node(id)
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
    
    # fill up a dictionary of node id -> height in tree where
    # leaves have height = 0
    def heightTable(self, node, htable):
        children = [x[1] for x in self.dag.out_edges(node)]
        if len(children) == 0:
            htable[node] = 0
        else:            
            for i in children:
                self.heightTable(i, htable)
            htable[node] = max([htable[i] for i in children]) + 1

    # check the candidate using the set and and fraction
    def inCandidateSet(self, node, candidateChildFrac):
        if self.candidateMap is None or len(self.candidateMap) == 0:
            return True
        if self.mcTree.getName(node) in self.candidateMap:
            return self.candidateMap[self.mcTree.getName(node)]
        children = self.mcTree.breadthFirstTraversal(node)
        leaves = []
        for child in children:
            if self.mcTree.isLeaf(child):
                leaves.append(child)
        candidateLeaves = 0
        for leaf in leaves:
            if self.mcTree.getName(leaf) in self.candidateMap:
                candidateLeaves += 1
        if len(leaves) == 0:
            self.candidateMap[self.mcTree.getName(node)] = False
            return False
        frac = float(candidateLeaves) / float(len(leaves))
        if frac >= candidateChildFrac:
            self.candidateMap[self.mcTree.getName(node)] = True
            return True
        self.candidateMap[self.mcTree.getName(node)] = False
        return False
                                    
    # greedily assign closest possible valid outgroup
    # all outgroups are stored in self.ogMap
    # edges between leaves ARE NOT kept in the dag    
    # the threshold parameter specifies how much parallelism can
    # be sacrificed by the selection of an outgroup
    # threshold = None : just greedy with now constraints
    # threshold = 0 : depth of schedule guaranteed to be unaffected by outgroup
    # threshold = k : depth increases by at most k per outgroup
    # candidateSet : names of valid outgroup genomes. (all if is None)
    # candidateChildFrac : min fraction of children of ancestor in
    # candidatSet in order for the ancestor to be an outrgoup candidate
    # if > 1, then only members of the candidate set and none of their
    # ancestors are chosen
    def greedy(self, threshold = None, candidateSet = None,
               candidateChildFrac = 2.):
        orderedPairs = []
        for source, sinks in self.dm.items():
            for sink, dist in sinks.items():
                if source != self.root and sink != self.root:
                    orderedPairs.append((dist, (source, sink)))
        orderedPairs.sort(key = lambda x: x[0])
        finished = set()
        self.ogMap = dict()
        self.candidateMap = dict()
        if candidateSet is not None:
            assert isinstance(candidateSet, set)
            for candidate in candidateSet:
                self.candidateMap[candidate] = True
        
        htable = dict()
        self.heightTable(self.root, htable)
        
        for candidate in orderedPairs:
            source = candidate[1][0]
            sink = candidate[1][1]
            dist = candidate[0]
            
            # skip leaves (as sources)
            if len(self.dag.out_edges(source)) == 0:
                finished.add(source)
                
            # skip nodes that aren't in the candidate set (if specified)
            # or don't have enough candidate children
            if not self.inCandidateSet(sink, candidateChildFrac):
                continue
            
            # canditate pair exceeds given threshold, so we skip
            if threshold is not None and \
            htable[sink] - htable[source] + 1 > threshold:
                continue
    
            if source not in finished and \
            not self.onSamePath(source, sink):
                self.dag.add_edge(source, sink, weight=dist, info='outgroup')
                if NX.is_directed_acyclic_graph(self.dag):
                    finished.add(source)
                    sourceName = self.mcTree.getName(source)
                    sinkName = self.mcTree.getName(sink)
                    self.ogMap[sourceName] = (sinkName, dist)
                    htable[source] = max(htable[source], htable[sink] + 1)                 
                else:
                    self.dag.remove_edge(source, sink)

def main():
    usage = "usage: %prog <project> <output graphviz .dot file>"
    description = "TEST: draw the outgroup DAG"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--justLeaves", dest="justLeaves", action="store_true", 
                      default = False, help="Assign only leaves as outgroups")
    parser.add_option("--threshold", dest="threshold", type='int',
                      default = None, help="greedy threshold")
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    proj = MultiCactusProject()
    proj.readXML(args[0])
    outgroup = GreedyOutgroup()
    outgroup.importTree(proj.mcTree)
    if options.justLeaves:
        candidates = set([proj.mcTree.getName(x)
                          for x in proj.mcTree.getLeaves()])
    else:
        candidates = None
    outgroup.greedy(candidates, threshold=options.threshold,
                    candidateChildFrac=1.1)
    NX.drawing.nx_agraph.write_dot(outgroup.dag, args[1])
    return 0

if __name__ == '__main__':    
    main()
