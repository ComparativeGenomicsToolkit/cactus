#!/usr/bin/env python

"""Example script to build trees for a reconstruction problem.
"""

import random
import logging
from sonLib.bioio import logger
import xml.etree.ElementTree as ET

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions

from sonLib.tree import makeRandomBinaryTree

def main():
    ##########################################
    #Construct the arguments.
    ##########################################  
    
    parser = getBasicOptionParser("usage: %prog [options]", "%prog 0.1")
    
    parser.add_option("--absolutePathPrefix", dest="absolutePathPrefix", 
                      help="The path to the root of the reconstruction tree problem",
                      default="None")
    
    parser.add_option("--reconstructionProblem", dest="reconstructionProblem", 
                      help="File containing the reconstruction problem",
                      default="None")
    
    parser.add_option("--uniqueNamePrefix", dest="uniqueNamePrefix", 
                      help="An alpha-numeric prefix which, when appended with any alpha-numeric characters is guaranteed to produce a unique name",
                      default=None)
    
    options, args = parseBasicOptions(parser)
    assert options.uniqueNamePrefix != None
    logger.setLevel(logging.DEBUG)
    
    logger.info("Parsed arguments")
    assert len(args) == 0
    
    absoluteReconstructionProblem = "%s/%s" % (options.absolutePathPrefix, options.reconstructionProblem)
    
    logger.info("input file: %s" % absoluteReconstructionProblem)
    
    ##########################################
    #Temp files
    ##########################################

    #It might not be obvious, but the parser automatically parses the dir for the temp files.
    logger.info("The temp dir root: %s" % options.tempDirRoot)
    
    ##########################################
    #Parse the reconstruction problem
    ##########################################
    
    reconstructionProblemTag = ET.parse(absoluteReconstructionProblem).getroot()
    
    logger.info("Parsed the reconstruction problem")
    
    ##########################################
    #Stuff to make random trees
    ##########################################
    
    ###You actually need to construct chains, get alignments for each atom/chain
    #and then make trees.
    
    eventTree = reconstructionProblemTag.find("event_tree").find("phylogeny").find("clade")
    atoms = reconstructionProblemTag.find("atoms")
    
    logger.info("Got the event tree and atoms")
    
    #Make a set of the leaf events, to check stuff
    def leafEvents(eventNode, leafEventsSet):
        if len(eventNode.findall("clade")) > 0:
            for childEventNode in eventNode.findall("clade"):
                leafEvents(childEventNode, leafEventsSet)
        else:
            leafEventsSet.add(eventNode.attrib["event"])
    leafEventsSet = set()
    leafEvents(eventTree, leafEventsSet)
    
    #Make a map of the contigs and their events
    contigsToEventsMap = {}
    sequences = reconstructionProblemTag.find("sequences")
    for sequence in sequences:
        contigsToEventsMap[sequence.attrib["contig"]] = sequence.attrib["event"]
        assert sequence.attrib["event"] in leafEventsSet
    
    logger.info("Made a map of contigs and events")
    
    #Build a map of the all the child event nodes for each event node.
    def makeAllEventsMap(eventNode, eventsMap):
        eventNodes = [ eventNode.attrib["event"] ]
        for childEventNode in eventNode.findall("clade"):
            eventNodes += makeAllEventsMap(childEventNode, eventsMap) 
        eventsMap[eventNode] = eventNodes
        return eventNodes
    eventsMap = {}
    makeAllEventsMap(eventTree, eventsMap)
    
    #Builds a map of parent events
    def makeParentEventMap(eventNode, eventMap):
        for childEventNode in list(eventNode.findall("clade")):
            eventMap[childEventNode] = eventNode
            makeParentEventMap(childEventNode, eventMap)
    parentEventMap = {}
    makeParentEventMap(eventTree, parentEventMap)
    
    extraEventIndex = [0]
    for atom in atoms:
        logger.info("Building a random tree for atom: %s" % atom.text)
        
        atomTree = ET.SubElement(atom, "phylogeny")
        atomTreeRoot = ET.SubElement(atomTree, "clade") #We deal with the root seperately
        atomTreeRoot.attrib["event"] = eventTree.attrib["event"]
        
        #Get a random tree of the appropriate size, we have to iterate
        #a bit with random trees to get a random tree that has appropriate size.
        atomInstances = list(atom.find("atom_instances").findall("atom_instance"))
        binaryTree = makeRandomBinaryTree()
        while True:
            def fn(binaryTree):
                #Leaf number
                if binaryTree.internal == True:
                    return fn(binaryTree.left) + fn(binaryTree.right)
                return 1
            if fn(binaryTree) == len(atomInstances):
                break
            binaryTree = makeRandomBinaryTree()
        
        logger.info("Built a random tree for the atom instances")
        
        def getCommonAncestorEvent(parentEventNode, event1, event2):
            """Gets the most recent common ancestor of two event nodes in the eventtree.
            """
            assert event1 in eventsMap[parentEventNode]
            assert event2 in eventsMap[parentEventNode]
            for childEventNode in parentEventNode.findall("clade"):
                events = eventsMap[childEventNode]
                if event1 in events and event2 in events:
                    return getCommonAncestorEvent(childEventNode, event1, event2)
            return parentEventNode
        
        def fn2(binaryTree):
            """Function walks bottom up the atom tree labelling it with atom instances (at the leaves),
            and events.
            """
            if binaryTree.internal:
                eventNode = getCommonAncestorEvent(eventTree, fn2(binaryTree.left), fn2(binaryTree.right))
                if binaryTree.left.event == eventNode.attrib["event"] or \
                binaryTree.right.event == eventNode.attrib["event"]:
                    #We must event an ancestral event to this point
                    parentEventNode = parentEventMap[eventNode]
                    if random.random() > 0.95 or \
                    parentEventNode.attrib["event"] == eventTree.attrib["event"] or \
                    len(parentEventNode.findall("clade")) != 1: #This last checks prevents dups coalescening during 'speciation' event.
                        extraNode = ET.SubElement(parentEventNode, "clade")
                        parentEventNode.remove(eventNode)
                        extraNode.append(eventNode)
                        extraNode.attrib["event"] = "%s%i" % (options.uniqueNamePrefix, extraEventIndex[0]) 
                        extraEventIndex[0] += 1
                        eventNode = extraNode
                        #Rejig the event map, given the extra event
                        eventsMap.clear()
                        makeAllEventsMap(eventTree, eventsMap)
                         #Reset the parent event map
                        parentEventMap.clear()
                        makeParentEventMap(eventTree, parentEventMap)
                    else:
                        eventNode = parentEventNode
                binaryTree.event = eventNode.attrib["event"]
            else:
                atomInstance = atomInstances.pop()
                binaryTree.atomInstance = atomInstance
                binaryTree.event = contigsToEventsMap[atomInstance.find("coordinates").attrib["contig"]]
            return binaryTree.event
            
        fn2(binaryTree)
        assert len(atomInstances) == 0
        logger.info("Labelled the tree with the atom instances")
        
        #Now add to the XML
        def fn4(atomNode, binaryTree):
            atomNode = ET.SubElement(atomNode, "clade")
            atomNode.attrib["event"] = binaryTree.event
            if binaryTree.internal:
                assert atomNode.attrib["event"] not in leafEventsSet
                fn4(atomNode, binaryTree.left)
                fn4(atomNode, binaryTree.right)
            else:
                assert atomNode.attrib["event"] in leafEventsSet
                atomNode.attrib["instance"] = binaryTree.atomInstance.text
        fn4(atomTreeRoot, binaryTree)
        logger.info("Added the atom tree to the XML")
    
    logger.info("Added trees to the atoms")
    
    ##########################################
    #Write out the updated xml
    ##########################################
    
    fileHandle = open(absoluteReconstructionProblem, 'w')
    tree = ET.ElementTree(reconstructionProblemTag)
    tree.write(fileHandle)
    fileHandle.close()
    
    logger.info("Written out the updated reconstruction problem")
            
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
