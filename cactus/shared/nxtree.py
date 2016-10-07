#!/usr/bin/env python

#Copyright (C) 2006-2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""A more general (ie arbitrary degree) tree to replace the BinaryTree
in sonLib.tree.py.  Implemented as a lightweight wrapper over a 
networkx digraph.  nodes are accessed using unique ids.    

for now, there is no interface (beyond directly working the 
networkx graph or reading from newick) to modify the structure of the tree.
update -- added function to remove leaves (and their parents if deg. 2)
"""
import sys
import os
import networkx as NX
from optparse import OptionParser

class NXTree(object):    
    def __init__(self, nxDg = None):
        self.nxDg = nxDg
        if self.nxDg is None:
            self.nxDg = NX.DiGraph()
        self.isTree()    
        self.rootId = None
        self.computeRootId()
        
    def isTree(self):
        assert NX.is_directed_acyclic_graph(self.nxDg)
        for node in self.nxDg.nodes():
            assert len(self.nxDg.in_edges(node)) < 2
            
    def computeRootId(self):
        self.rootId = None
        for node in self.nxDg.nodes():
            if len(self.nxDg.in_edges(node)) == 0:
                self.rootId = node
                break
    
    def loadNetworkXTree(self, nxDg):
        self.nxDg = nxDg
        self.isTree()
        self.computeRootId()            
            
    def getChildren(self, id):
        assert id in self.nxDg
        return sorted([x[1] for x in self.nxDg.out_edges(id)])
    
    def isLeaf(self, id):
        return len(self.getChildren(id)) == 0
    
    def getLeaves(self):
        leaves = []
        for i in self.breadthFirstTraversal():
            if self.isLeaf(i):
                leaves.append(i)
        return leaves
    
    def hasParent(self, id):
        return self.getParent(id) is not None

    def getParent(self, id):
        assert id in self.nxDg
        edges = self.nxDg.in_edges(id)
        assert len(edges) < 2
        if len(edges) == 0:
            return None
        else:
            return edges[0][0]
    
    def getName(self, id):
        assert id in self.nxDg
        node = self.nxDg.node[id]
        if 'name' in node:
            return node['name']
        return ""
    
    def setName(self, id, name):
        assert id in self.nxDg
        self.nxDg.node[id]['name'] = name
    
    def hasName(self, id):
        assert id in self.nxDg
        node = self.nxDg.node[id]
        return 'name' in node
        
    def getWeight(self, parentId, childId, defaultValue=None):
        assert parentId in self.nxDg and childId in self.nxDg
        edge = self.nxDg[parentId][childId]
        if 'weight' in edge:
            return edge['weight']
        return defaultValue
    
    def setWeight(self, parentId, childId, weight):
        assert parentId in self.nxDg and childId in self.nxDg
        self.nxDg[parentId][childId]['weight'] = float(weight)
    
    def getRootId(self):
        return self.rootId
    
    def getRootName(self):
        return self.getName(self.getRootId())
        
    def preOrderTraversal(self, root = None):
        if root is None:
            root = self.rootId
        yield root 
        for child in self.getChildren(root):
            for node in self.preOrderTraversal(child):
                yield node
          
    def postOrderTraversal(self, root = None):
        if root is None:
            root = self.rootId
        for child in self.getChildren(root):
            for node in self.postOrderTraversal(child):
                yield node
        yield root
        
    def breadthFirstTraversal(self, root = None):
        if root is None:
            root = self.rootId
        bfQueue = [root]
        while len(bfQueue) > 0:
            node = bfQueue.pop(0)
            yield node
            for child in self.getChildren(node):
                bfQueue.append(child)

    def removeDegree2Vertex(self, nodeId):
        if self.hasParent(nodeId) is True:
            parentId = self.getParent(nodeId)
            topWeight = self.getWeight(parentId, nodeId)
            children = self.getChildren(nodeId)
            if len(children) == 1:
                childId = children[0]
                botWeight = self.getWeight(nodeId, childId)
                if topWeight is None and botWeight is None:
                    newWeight = None
                elif topWeight is None:
                    newWeight = botWeight
                elif botWeight is None:
                    newWeight = topWeight
                else:                    
                    newWeight = topWeight + botWeight
                self.nxDg.remove_node(nodeId)
                self.nxDg.add_edge(parentId, childId)
                self.nxDg[parentId][childId]['weight'] = float(newWeight)

    def removeLeaf(self, leafId):
        assert self.isLeaf(leafId)
        parentId = self.getParent(leafId)
        self.nxDg.remove_node(leafId)
        if len(self.getChildren(parentId)) == 1 and \
           self.hasParent(parentId) is True:
            self.removeDegree2Vertex(parentId)

    def removeEdge(self, parentId, childId):
        self.nxDg.remove_edge(parentId, childId)
        
    def reroot(self, newRootId):
        flipEdges = []
        node = newRootId        
        while self.hasParent(node) is True:
            parent = self.getParent(node)
            flipEdges.append((parent, node))
            node = parent
        for edge in flipEdges:
            self.nxDg.add_edge(edge[1], edge[0])
            self.nxDg[edge[1]][edge[0]] = self.nxDg[edge[0]][edge[1]]
            self.nxDg.remove_edge(edge[0], edge[1])
        self.rootId = newRootId

        

            
                
                
        
        
    
    
            
