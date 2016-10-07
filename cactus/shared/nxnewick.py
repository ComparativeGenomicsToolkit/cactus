#!/usr/bin/env python

#Copyright (C) 2006-2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""read and write newick trees to and from networkx graphs (as wrapped by nxtree). 
"""
import sys
import os
import re
import math
import random
from string import whitespace as ws
from cactus.shared.misc import close
import bioio
import networkx as NX
from optparse import OptionParser
from cactus.shared.nxtree import NXTree


class NXNewick(object):
    def __init__(self, nxTree = None):
        self.nxTree = None
        self.bracketMatch = None
        self.inString = None
        self.nextId = 0
        self.outString = None
        
    def parseFile(self, path):
        inFile = open(path)
        self.parseString(inFile.read())
        inFile.close()
        return self.nxTree
    
    def parseString(self, newickString, addImpliedRoots = True):
        self.nxTree = NXTree()
        self.inString = self.__filterWhitespace(newickString)
        self.__createBracketTable()
        self.nextId = 0
        assert self.inString[-1] == ';'
        self.__addNode(0, len(self.inString)-1, None, addImpliedRoots)
        self.nxTree.isTree()
        return self.nxTree
    
    def writeString(self, nxTree = None):
        if nxTree:
            self.nxTree = nxTree
        self.outString = ""
        self.__writeNode(self.nxTree.getRootId(), None)
        self.outString += ";"
        return self.outString
    
    def writeFile(self, path, nxTree = None):
        outFile = open(path, "w")
        outFile.write(self.writeString(nxTree))
        outFile.write("\n")
        outFile.close()
        return None
 
    #### PRIVATE WRITING FUNCTIONS ####
    def __writeNode(self, node, parent = None):
        children = self.nxTree.getChildren(node)
        if len(children) > 0:
            self.outString += "("
            for child in children:
                self.__writeNode(child, node)
                if child != children[-1]:
                    self.outString += ","
        if len(children) > 0:
            self.outString += ")"
        name = self.nxTree.getName(node)
        
        if len(name) > 0:
            containsSpace = True in [c1 in name for c1 in ws]
            if containsSpace:
                self.outString += "\""
            self.outString += name
            if containsSpace:
                self.outString += "\""
        if parent is not None:
            weight = self.nxTree.getWeight(parent, node, defaultValue=None)
            if weight is not None:
                self.outString += ":%s" % str(weight)      
        
    #### PRIVATE READING FUNCTIONS ####       
    
    def __filterWhitespace(self, newickString):
        filteredString = ""
        inQuote = False
        for c in newickString:
            if c == "\'" or c == "\"":
                inQuote = not inQuote
            elif inQuote or c not in ws:
                filteredString += c
        return filteredString
            
    def __createBracketTable(self):
        bracketStack = []
        self.bracketMatch = dict()
        index = 0
        for c in self.inString:
            if c == '(':
                bracketStack.append(index)
            elif c == ')':
                leftIndex = bracketStack.pop()
                self.bracketMatch[leftIndex] = index
            index += 1
        assert len(bracketStack) == 0
    
    def __childRanges(self, start, length):
        ranges = []
        currentStart = start
        i = currentStart
        while i < start + length + 1:
            if self.inString[i] == ',' or i == start + length:
                ranges.append((currentStart, i - currentStart))
                currentStart = i + 1
            if i in self.bracketMatch:
                i = self.bracketMatch[i]
            i += 1
        return ranges
    
    def __parseName(self, nameString):
        if nameString == ';':
            return ('','')
        tokens = nameString.split(':')
        assert len(tokens) == 1 or len(tokens) == 2
        name = tokens[0]
        weight = ''
        if len(tokens) == 2:
            weight = tokens[1]
        return (name, weight)
        
    def __addNode(self, start, length, parent = None, addImpliedRoots = True):
        # parse the children (..,...,..)
        children = []
        if self.inString[start] == '(':
            assert start in self.bracketMatch
            chLength = self.bracketMatch[start] - start - 1
            children = self.__childRanges(start+1, chLength)      
            start = self.bracketMatch[start] + 1
            length -= (chLength + 2)
            
        # prase the name abc:123
        name, weight = self.__parseName(self.inString[start:start+length])
        id = self.nextId
        self.nextId += 1
        self.nxTree.nxDg.add_node(id)
        if len(name) > 0:
            self.nxTree.nxDg.node[id]['name'] = name
       
        #update the graph
        if parent is not None:
            self.nxTree.nxDg.add_edge(parent, id)
            if len(weight) > 0:
                self.nxTree.nxDg[parent][id]['weight'] = float(weight)
       
        #update the root (implied roots are added as a new node)
        if self.nxTree.getRootId() is None:
            assert parent is None
            root = id
            if len(weight) > 0:
                if addImpliedRoots:
                    root = self.nextId
                    self.nextId += 1
                    self.nxTree.nxDg.add_edge(root, id)
                    self.nxTree.setWeight(root,id, weight)
            self.nxTree.rootId = root
        
        # recurse on children
        for child in children:
            self.__addNode(child[0], child[1], id)

def main():
    usage = "usage: %prog <tree file> <output graphviz .dot file> <output tree file>"
    description = "TEST: convert newick tree to graphviz tree"
    parser = OptionParser(usage=usage, description=description)
    options, args = parser.parse_args()
    
    if len(args) != 3:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    parser = NXNewick()
    parser.parseFile(args[0])
    NX.drawing.nx_agraph.write_dot(parser.nxTree.nxDg, args[1])
    parser.writeFile(args[2])
    print "PRE"
    for i in parser.nxTree.preOrderTraversal():
        print ("%d %s" % (i, parser.nxTree.getName(i)))
    print "POST"
    for i in parser.nxTree.postOrderTraversal():
        print ("%d %s" % (i, parser.nxTree.getName(i)))
    print "BFS"
    for i in parser.nxTree.breadthFirstTraversal():
        print ("%d %s" % (i, parser.nxTree.getName(i)))
    return 0

if __name__ == '__main__':    
    main()
