#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import os
import sys
import copy
import xml.etree.ElementTree as ET
import networkx as NX
import random
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.progressive.schedule import Schedule

from sonLib.nxnewick import NXNewick
from sonLib.nxtreeTest import randomTreeSet

class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

    @TestStatus.shortLength
    def testRunsOnRandom(self):
        for tree in randomTreeSet():
            if tree.size() < 120:
                dag = self.__addDagEdges(tree)
                sched = Schedule()
                sched.inGraph = dag
                sched.maxParallel = 2
                sched.compute()

    def __addDagEdges(self, tree):
        count = tree.size() // random.randrange(1,10)
        tsort = list(NX.topological_sort(tree))
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
