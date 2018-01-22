#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import networkx as NX
import random

from cactus.progressive.schedule import Schedule

from sonLib.nxtreeTest import randomTreeSet

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def testRunsOnRandom(self):
        for tree in randomTreeSet():
            if tree.size() < 120:
                dag = self.__addDagEdges(tree)
                sched = Schedule()
                sched.inGraph = dag
                sched.compute()
                
    def __addDagEdges(self, tree):
        count = tree.size() / random.randrange(1,10)
        tsort = NX.topological_sort(tree)
        for i in range(0, count):
            source = random.randrange(0, len(tsort)-1)            
            sink = random.randrange(source+1, len(tsort))
            sourceNode = tsort[source]
            sinkNode = tsort[sink]
            if (sourceNode, sinkNode) not in tree.out_edges():
                assert sourceNode != sinkNode
                tree.add_edge(sourceNode, sinkNode)
        assert NX.is_directed_acyclic_graph(tree)
        return tree
            
            
def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
