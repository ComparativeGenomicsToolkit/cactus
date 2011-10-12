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


from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.experimentWrapper import ExperimentWrapper

class Schedule:
    def __init__(self):
        # input graph
        self.inGraph = None
        # output tree
        self.depTree = None
    
    # read the experiments, compute the dependency dag
    def loadProject(self, mcProject):
        self.inGraph = NX.DiGraph()
        for name, expPath in mcProject.expMap.items():
            exp = ExperimentWrapper(ET.parse(expPath).getroot())
            tree = exp.getTree()
            for leaf in tree.getLeaves():
                self.inGraph.add_edge(name, tree.getName(leaf))
        assert NX.is_directed_acyclic_graph(self.inGraph)
    
    # break all the cycles in reverse topological order
    def compute(self):
        self.depTree = self.inGraph.copy()
        tsort = NX.topological_sort(self.depTree)
        tsort.reverse()
        for node in tsort:
            parents = [i[0] for i in self.depTree.in_edges(node)]
            if len(parents) > 1:
                chosenParentIndex = min([tsort.index(i) for i in parents])
                chosenParent = tsort[chosenParentIndex]
                for parent in parents:
                    if parent != chosenParent:
                        self.depTree.remove_edge(parent, node)
                        if chosenParent not in self.deps(parent):
                            self.depTree.add_edge(parent, chosenParent)
            assert len(self.depTree.in_edges(node)) < 2
        
        for node in tsort:
            assert len(self.depTree.in_edges(node)) < 2
        assert NX.is_directed_acyclic_graph(self.depTree)
 
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
    schedule.loadProject(proj)
    schedule.compute()
    schedule.writeToFile(args[1])    

if __name__ == '__main__':    
    main()