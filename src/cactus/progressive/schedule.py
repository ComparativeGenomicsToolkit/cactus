#!/usr/bin/env python3

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


from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import cactusRootPath
from cactus.progressive.seqFile import SeqFile

from sonLib.nxnewick import NXNewick

class Schedule:
    def __init__(self):
        # input graph
        self.inGraph = None
        # output tree
        self.depTree = None
        # maximum number of nodes that can be computed in parallel
        # we will change the topology such that no more than this
        # number of internal nodes can ever be dependency-free at
        # the same time
        self.maxParallelSubtrees = None

    # read the tree and outgroup map, compute the dependency dag
    def loadProject(self, mc_tree, og_map, max_parallel_subtrees):
        self.inGraph = NX.DiGraph()
        globTree = mc_tree
        self.maxParallelSubtrees = max_parallel_subtrees
        leafEvents = [globTree.getName(i) for i in globTree.getLeaves()]

        # add dependency edge from root to children (including outgroups)
        for node in globTree.getSubtreeRoots():
            name = globTree.getName(node)
            ogs = []
            if name in og_map:
                for og in og_map[name]:
                    ogs.append(og)
            igs = [globTree.getName(child) for child in globTree.getChildren(node)]

            for dep in igs + ogs:
                if not mc_tree.isLeaf(mc_tree.getNodeId(dep)):
                    self.inGraph.add_edge(name, dep)

        assert NX.is_directed_acyclic_graph(self.inGraph)

    # break all the cycles in reverse topological order
    def compute(self):
        def getName(pref, id):
            while '%s%d' % (pref, id) in self.depTree:
                id += 1
            return '%s%d' % (pref, id)
        self.depTree = self.inGraph.copy()
        self.transitveReduction(self.depTree)
        tsort = list(reversed(list(NX.topological_sort(self.depTree))))
        nextId = 0
        for node in tsort:
            #strip out cycles within parents
            parents = set([i[0] for i in self.depTree.in_edges(node)])
            if len(parents) > 1:
                for parent in parents:
                    for pedge in self.depTree.out_edges(parent):
                        if pedge[1] in parents:
                            self.depTree.remove_edge(parent, node)
                            break

            parents = [i[0] for i in self.depTree.in_edges(node)]
            if len(parents) > 1:
                # insert new "virtual" node as parent to all parents
                vpNode = getName("Virtual", nextId)
                self.depTree.add_node(vpNode, virtual='1')
                # insert new "virtual follow-on" node as follow up to above node
                vfoNode = getName("VirtualF", nextId)
                self.depTree.add_node(vfoNode, virtual='1')
                self.depTree.add_edge(vpNode, vfoNode, followOn='1')

                nextId += 1
                for parent in parents:
                    # insert vpNode above all the parents
                    for inEdge in list(self.depTree.in_edges(parent)):
                        if not self.isFollowOn(inEdge[0], inEdge[1]):
                            self.depTree.add_edge(inEdge[0], vpNode)
                            self.depTree.remove_edge(inEdge[0], inEdge[1])
                    # add parent as child of the followOn
                    self.depTree.add_edge(vfoNode, parent)
                    # remove node from parent
                    self.depTree.remove_edge(parent, node)
                    # add parent's childrent to vpNode
                    for child in list(self.depTree.successors(parent)):
                        if not self.isFollowOn(parent, child):
                            self.depTree.remove_edge(parent, child)
                            self.depTree.add_edge(vpNode, child)
                            if self.depTree.has_edge(child, vpNode):
                                self.depTree.remove_edge(child, vpNode)
                #add node to vpNode
                self.depTree.add_edge(vpNode, node)
                #lazy clean-up of inserted transitive edges
                #will need to not do this brute-force for bigger trees
                self.transitveReduction(self.depTree)
                assert(len(self.depTree.in_edges(node)) == 1)
                # process virtual node next
                tsort.insert(tsort.index(node) + 1, vpNode)

            assert len(self.depTree.in_edges(node)) < 2

        for node in tsort:
            assert len(self.depTree.in_edges(node)) < 2
        assert NX.is_directed_acyclic_graph(self.depTree)
        self.enforceMaxParallel()

    # friday afternoon!
    def transitveReduction(self, digraph):
        paths = dict(NX.all_pairs_shortest_path(digraph))
        def hasPath(node1, node2):
            if node1 == node2:
                return False
            return (node1 in paths and node2 in paths[node1] and
                    len(paths[node1][node2]) > 0)
        for x in digraph.nodes():
            for y in digraph.nodes():
                for z in digraph.nodes():
                    if (x != z and digraph.has_edge(x, z) and hasPath(x, y) and
                        hasPath(y, z)):
                        digraph.remove_edge(x, z)


    # add dependencies to ensure that more than self.maxParallelSubtrees
    # different jobs can never be scheduled at the same time
    def enforceMaxParallel(self):
        if self.maxParallelSubtrees is None or len(self.depTree.nodes()) <= 2:
            return
        assert self.maxParallelSubtrees > 0
        tree = self.depTree.copy()
        # remove followOn edges
        for edge in list(tree.edges()):
            if self.isFollowOn(edge[0], edge[1]):
                tree.remove_edge(edge[0], edge[1])
        # get roots
        roots = []
        for node in tree.nodes():
            if len(tree.in_edges(node)) == 0:
                roots.append(node)
        # for each root, independently reduce the number of leaves
        for root in roots:
            leaves = []
            for node in NX.bfs_tree(tree, root):
                if (len(tree.out_edges(node)) == 0 and
                    not self.isVirtual(node)):
                    leaves.append(node)
            numLeaves = len(leaves)
            # sort leaves by chain length (lowest at end)
            leafChains = [self.getChainParent(tree, x) for x in leaves]
            leafChains.sort()
            leaves = [x[1] for x in leafChains]
            leaves.reverse()
            while numLeaves > self.maxParallelSubtrees:
                # try to insert the last leaf of the list above
                # its sibling
                leaf = leaves[-1]
                chainParent = self.getChainParent(tree, leaf)[2]
                parents = [x for x in tree.predecessors(chainParent)]
                assert len(parents) == 1
                parent = parents[0]
                children = [x for x in tree.successors(parent)]
                for child in children:
                    if child != chainParent and chainParent != root:
                        tree.add_edge(leaf, child)
                        self.depTree.add_edge(leaf, child)
                        tree.remove_edge(parent, child)
                        self.depTree.remove_edge(parent, child)
                        numLeaves -= 1
                        break
                del leaves[-1]
            assert numLeaves <= self.maxParallelSubtrees

    # walk up a path (not bifurcations) from a leaf as far as we can
    def getChainParent(self, tree, node):
        chainLen = 0
        chainParent = node
        while True:
            parents = list(tree.predecessors(chainParent))
            if len(parents) == 1 and len(list(tree.successors(parents[0]))) == 1:
                chainParent = parents[0]
                chainLen += 1
            else:
                break
        return (chainLen, node, chainParent)

    # for a given event name, get the names of all the
    # events that are directly dependent on it in the
    # schedule
    def deps(self, name):
        assert name in self.depTree
        edges = self.depTree.out_edges(name)
        depList = []
        for edge in edges:
            if not self.isFollowOn(edge[0], edge[1]):
                depList.append(edge[1])
        return depList

    # get the follow on node if it exists
    def followOn(self, name):
        assert name in self.depTree
        edges = self.depTree.out_edges(name)
        depList = []
        for edge in edges:
            if self.isFollowOn(edge[0], edge[1]):
                depList.append(edge[1])
        assert len(depList) < 2
        if len(depList) == 1:
            return depList[0]
        return None

    # test if an edge is a follow on
    def isFollowOn(self, parent, child):
        edge = self.depTree[parent][child]
        return 'followOn' in edge and str(edge['followOn']) == '1'

    # test if node corresponds to genome event or was added to break a virtual
    def isVirtual(self, name):
        node = self.depTree.nodes[name]
        return 'virtual' in node and str(node['virtual']) == '1'

    # write to graphviz dot file
    def writeToFile(self, path):
        NX.drawing.nx_agraph.write_dot(self.depTree, path)

    # read from graphviz dot file
    def readFromFile(self, path):
        self.depTree = NX.read_edgelist(path)
        self.checkDepTree()


def main():
    usage = "usage: %prog <seqfile> <output graphviz .dot file>"
    description = "TEST: create schedule from project file"
    parser = OptionParser(usage=usage, description=description)
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))    

    options, args = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    # load up the seqfile and figure out the outgroups and schedule
    # can't use progressive_decomposition.py for circular deps...
    config_node = ET.parse(options.configFile).getroot()
    config_wrapper = ConfigWrapper(config_node)
    config_wrapper.substituteAllPredefinedConstantsWithLiterals()
    seq_file = SeqFile(args[0])
    mc_tree = MultiCactusTree(seq_file.tree)
    mc_tree.nameUnlabeledInternalNodes(config_wrapper.getDefaultInternalNodePrefix())
    mc_tree.computeSubtreeRoots()
    
    # todo:
    # move to own module so can use progressive_decomposition.py
    assert(False)

if __name__ == '__main__':
    main()
