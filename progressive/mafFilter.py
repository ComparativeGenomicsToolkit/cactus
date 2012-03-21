#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Remove outgroup from a maf file.  
If specified in the config, remove internal nodes too 

"""

import os
import sys
import xml.etree.ElementTree as ET

from sonLib.nxtree import NXTree
from sonLib.nxnewick import NXNewick
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.configWrapper import ConfigWrapper
from cactus.progressive.experimentWrapper import ExperimentWrapper

# remove an event from the tree
# resulting tree won't be binary
def removeNode(tree, event):
    os.system("echo \'%s xxx %s\' > ~/blin.txt" % (event, NXNewick().writeString(tree) ))
    id = None
    for node in tree.breadthFirstTraversal():
        if tree.getName(node) == event:
            id = node
    if id is not None:
        assert id != tree.getRootId()
        parent = tree.getParent(id)
        assert parent is not None
        for child in tree.getChildren(id):
            tree.nxDg.remove_edge(id, child)
            tree.nxDg.add_edge(parent, child)
        tree.nxDg.remove_node(id)

# cactus_createMulti.... guarantees one event name is never
# a prefix of another EXCEPT in the case of the _self. 
# so any sequence with the given outgroup name as a prefix
# except those with outgroup_self as the prefix belong
# to the outgroup and need to be filtered
def eventsInLine(names, line):
    tokens = line.split()
    if line[0] == 's' and len(line) >= 2:
        word = tokens[1]
        for name in names:
            if word.find(name) == 0 and \
            word.find(name + MultiCactusTree.self_suffix) != 0:
                return True
    return False
    
# quick and dirty friday afternoon! 
def removeEventsFromMaf(mafPath, events):
    if events is None or len(events) == 0:
        return
    for i in range(len(events)):
        events[i] = events[i].split('.')[0]        
    tempPath = "%s.og" % mafPath
    inFile = open(mafPath)
    outFile = open(tempPath, 'w')
    num = 1
    newick = NXNewick() 
    for line in inFile:
        if num == 2:
            treeString = line.split()[2]
            tree = MultiCactusTree(newick.parseString(treeString))
            for event in events:
                removeNode(tree, event)
            line = line.replace(treeString, newick.writeString(tree))            
        if line[0] == '#' or not eventsInLine(events, line):
            outFile.write(line)
        num += 1
    inFile.close()
    outFile.close()
    os.rename(tempPath, mafPath)

def mafFilterOutgroup(experiment):
    events = experiment.getOutgroupEvents()
    removeEventsFromMaf(experiment.getMAFPath(), events)
    
def mafFilterInternalNodes(event, options, project, experiment):
    config = ConfigWrapper(ET.parse(experiment.getConfig()).getroot())
    filterInternal = config.getFilterMafInternalEvents()
    if filterInternal == True:
        events = experiment.getOutgroupEvents()
        tree = project.mcTree
        rootId = tree.nameToId[event]
        for node in tree.breadthFirstTraversal(root=rootId):
            if node != rootId and tree.isLeaf(node) == False:
                events.append(tree.getName(node))
        if len(events) > 0:
            removeEventsFromMaf(experiment.getMAFPath(), events)
    