#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Wrap a tree and add some simple partitioning functionality.  note
that nameUnlabeledInternalNodes() and computeSubtreeRoots() need to be
called (in that order) for anything to work...

"""

import math

from sonLib.nxtree import NXTree

import networkx as nx
from networkx.algorithms.shortest_paths.weighted import dijkstra_path

class MultiCactusTree(NXTree):
    self_suffix = "_self"
    def __init__(self, tree = None):
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
    # add them to the self.subtreeRoots set
    def computeSubtreeRoots(self):
        self.subtreeRoots = set(node for node in self.breadthFirstTraversal() if not self.isLeaf(node))

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

    # Extracts a tree spanning the nodes with the given names.
    def extractSpanningTree(self, nodes):
        assert len(nodes) > 1
        nodeIds = [self.nameToId[name] for name in nodes]
        paths = [dijkstra_path(self.nxDg.to_undirected(), source=nodeIds[0], target=x) for x in nodeIds[1:]]
        nodesToInclude = set()
        for path in paths:
            for node in path:
                nodesToInclude.add(node)

        cpy = self.nxDg.subgraph(nodesToInclude).copy()
        # Get rid of nodes that have only 1 children
        graphWasModified = True
        while graphWasModified:
            graphWasModified = False
            for node in nx.nodes(cpy):
                if cpy.out_degree(node) == 1 and cpy.in_degree(node) == 1:
                    if node not in nodeIds:
                        # This is a spurious node in the species tree,
                        # we can and should remove
                        childEdge = list(cpy.out_edges(node, data=True))[0]
                        parentEdge = list(cpy.in_edges(node, data=True))[0]
                        child = childEdge[1]
                        childDist = childEdge[2]['weight']
                        parent = parentEdge[0]
                        assert parent != node
                        parentDist = parentEdge[2]['weight']
                        cpy.remove_node(node)
                        cpy.add_edge(parent, child, weight = childDist + parentDist)
                        graphWasModified = True
                        break

        mcCpy = MultiCactusTree(cpy)
        mcCpy.nameUnlabeledInternalNodes(prefix="thisPrefixShouldNeverAppear")
        mcCpy.computeSubtreeRoots()
        return mcCpy

    # get names of children below this node
    def getChildNames(self, name):
        id = self.nameToId[name]
        subtree = [i for i in self.traverseSubtree(id, id)]
        # remove the root from the set of children
        subtree.remove(id)
        names = [self.getName(i) for i in subtree]
        return names

    # copy a subtree rooted at node with given name
    def extractSubTree(self, name):
        root = self.nameToId[name]
        subtree = [i for i in self.traverseSubtree(root, root)]
        cpy = self.nxDg.subgraph(subtree).copy()
        mcCpy = MultiCactusTree(cpy)
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
        # de-activating assert because new outgroup munging in
        # cactus_createMultiCactusProject will temporarility duplicate
        # note in tree (is normal)
        #assert ogName not in self.nameToId
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
