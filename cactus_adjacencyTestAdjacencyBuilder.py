#!/usr/bin/env python

"""Example script to build random reconstructions.
"""

import random
import logging
from sonLib.bioio import logger
import xml.etree.ElementTree as ET

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions

def main():
    ##########################################
    #Construct the arguments.
    ##########################################  
    
    parser = getBasicOptionParser("usage: %prog [options]", "%prog 0.1")
    
    parser.add_option("--absolutePathPrefix", dest="absolutePathPrefix", 
                      help="The path to the root of the reconstruction tree problem")
    
    parser.add_option("--reconstructionProblem", dest="reconstructionProblem", 
                      help="File containing the reconstruction problem")
    
    parser.add_option("--uniqueNamePrefix", dest="uniqueNamePrefix", 
                      help="An alpha-numeric prefix which, when appended with any alpha-numeric characters is guaranteed to produce a unique name")

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
    #(Reverse evolution) Recursively transform the local break point graph into a progressively more ancestral form.
    ##########################################
    
    ###
    # A note on naming.
    # I call the nodes in the event tree 'events'
    # I call the nodes in cap trees 'cap nodes'
    # Where a string identifies an object I call it a 'name'
    # These naming conventions for the objects in the reconstruction should be adhered to.
    ####
    
    ####Global variables (stuff to be kept up to date by all the functions
    
    #Map of non-root event and cap nodes to there parent event/cap nodes.
    parentMap = {}
    #Retrieve a node by its instance name
    instanceNamesToCapNodes = {}
    #Retrieve an event by its instance name
    eventNamesToEvents = {}
    #Get the cap instance nodes for an event
    capNodesForEvents = {}
    #The name of the root event
    rootEventName = reconstructionProblemTag.find("event_tree").find("phylogeny").find("clade").attrib["event"]
    ####End global variables
    
    def getCapList():
        #Gets all caps in a reconstruction problem
        capList = {}
        for cap in reconstructionProblemTag.find("caps"):
            capList[cap.text] = cap
        for atom in reconstructionProblemTag.find("atoms"):
            for cap in atom.findall("cap"):
                capList[cap.text] = cap
        return capList
    
    #####Global variable Setup functions
    #Put cap trees in list;
    def setupGlobalVariables():
        capList = getCapList()
              
        #Hash for instance names to nodes
        # and a map of parent pointers for the cap trees
        for cap in capList.values():
            def fn(capNode):
                instanceNamesToCapNodes[capNode.attrib["instance"]] = capNode
                for childCapNode in capNode.findall("clade"):
                    parentMap[childCapNode] = capNode
                    fn(childCapNode)
            fn(cap.find("phylogeny").find("clade"))
            
        #Hash for event names to events
        # and a map of parent pointers for the event tree.
        def fn2(event):
            eventNamesToEvents[event.attrib["event"]] = event
            for childEvent in event.findall("clade"):
                parentMap[childEvent] = event
                fn2(childEvent)
        fn2(reconstructionProblemTag.find("event_tree").find("phylogeny").find("clade"))
    
    setupGlobalVariables()
    
    def getNodesForEvents():
        #Hash cap nodes in the cap trees to the names of events.
        capList = getCapList()
        def fn(event):
            capNodesForEvents[event.attrib["event"]] = set()
            for childEvent in event.findall("clade"):
                fn(childEvent)
        fn(reconstructionProblemTag.find("event_tree").find("phylogeny").find("clade"))
        for cap in capList.values():
            def fn2(capNode):
                capNodesForEvents[capNode.attrib["event"]].add(capNode)
                for childCapNode in capNode.findall("clade"):
                    fn2(childCapNode)
            fn2(cap.find("phylogeny").find("clade"))
    getNodesForEvents()
    
    logger.info("Setup the global variables")
    #####End global variable setup functions
     
    ####Functions used to build the adjacencies and operations 
    
    def linkCapNodes(capNode1, capNode2):
        assert capNode1.attrib["event"] == capNode2.attrib["event"]
        if capNode1.attrib.has_key("adjacency"):
            assert capNode1.attrib["adjacency"] == capNode2.attrib["instance"]
            assert capNode2.attrib["adjacency"] == capNode1.attrib["instance"]
        else:
            assert not capNode1.attrib.has_key("adjacency")
        capNode1.attrib["adjacency"] = capNode2.attrib["instance"]
        capNode2.attrib["adjacency"] = capNode1.attrib["instance"]
        
    def linkCapNodesSecondarily(capNode1, capNode2):
        assert capNode1.attrib["event"] == capNode2.attrib["event"]
        assert not capNode1.attrib.has_key("adjacency_2")
        assert not capNode2.attrib.has_key("adjacency_2")
        capNode1.attrib["adjacency_2"] = capNode2.attrib["instance"]
        capNode2.attrib["adjacency_2"] = capNode1.attrib["instance"]
        
    def isCycle(component):
        for i in xrange(1, len(component)):
            assert parentMap[component[i-1][1]] == parentMap[component[i][0]]
        #Returns true if the component is a cycle
        return parentMap[component[0][0]] == parentMap[component[-1][1]] and \
        parentMap[component[0][0]].attrib["event"] == parentMap[component[0][1]].attrib["event"] #This ensures the end parent node is part of the cycle of nodes
    
    def mostAncestralEvent(eventNames):
        #Gets the most ancestral event of a set of events
        def fn(eventName1, eventName2): #Sorts by position in tree using traversal to root.
            if eventName1 == eventName2: #Case of equal
                return 0
            event1 = eventNamesToEvents[eventName1]
            event2 = eventNamesToEvents[eventName2]
            while event1 != event2 and event1.attrib["event"] != rootEventName:
                event1 = parentMap[event1]
            if event1 == event2:
                return 1
            return -1
        eventNames = list(eventNames)
        eventNames.sort(fn)
        return eventNames[0]
    
    makeChildCapNode_Index = [0]
    def makeChildCapNode(eventName, parentCapNode):
        #Makes a child cap node hanging from the given parent cap node
        newCapNode = ET.SubElement(parentCapNode, "clade")
        assert mostAncestralEvent([ parentCapNode.attrib["event"], eventName ]) == parentCapNode.attrib["event"]
        assert parentCapNode.attrib["event"] != eventName
        newCapNode.attrib["event"] = eventName
        newCapNodeName = "%s.%s%i" % (parentCapNode.attrib["instance"].split(".")[0], options.uniqueNamePrefix, makeChildCapNode_Index[0])
        makeChildCapNode_Index[0] += 1
        newCapNode.attrib["instance"] = newCapNodeName
        newCapNode.attrib["augmented"] = "true"
        #Update the global variables
        parentMap[newCapNode] = parentCapNode
        instanceNamesToCapNodes[newCapNodeName] = newCapNode
        if not capNodesForEvents.has_key(eventName):
            capNodesForEvents[eventName] = set()
        capNodesForEvents[eventName].add(newCapNode)
        return newCapNode
        
    def addAncestralEventToCap(childCapNode, newEventName):
        #Adds a unary node to the cap immediately above the current one.
        if childCapNode.attrib["event"] != newEventName:
            assert mostAncestralEvent([ childCapNode.attrib["event"], newEventName ]) == newEventName
            parentCapNode = parentMap[childCapNode]
            if parentCapNode.attrib["event"] != newEventName:
                newCapNode = makeChildCapNode(newEventName, parentCapNode)
                parentCapNode.remove(childCapNode)
                newCapNode.insert(0, childCapNode)
                #Update the global variables
                parentMap[childCapNode] = newCapNode
                if not capNodesForEvents.has_key(newEventName):
                    capNodesForEvents[newEventName] = set()
                capNodesForEvents[newEventName].add(newCapNode)
                return newCapNode
            return parentCapNode
        return childCapNode
    
    createExtraStubCap_Index = [0]
    def createExtraStubCap(derivedEventName):
        #Creates an extra  stub cap containing the given event at its leaf. 
        capName = "$%s%i" % (options.uniqueNamePrefix, createExtraStubCap_Index[0])
        createExtraStubCap_Index[0] += 1
        
        assert derivedEventName != rootEventName
        assert mostAncestralEvent([ derivedEventName, rootEventName ]) == rootEventName
        
        cap = ET.SubElement(reconstructionProblemTag.find("caps"), "cap")
        cap.text = capName
        #ap.attrib['is_stub'] = 'true'
        
        capTree = ET.SubElement(cap, "phylogeny")
        capRootNode = ET.SubElement(capTree, "clade")
        #Add in the event name to identify with the species tree.
        capRootNode.attrib["event"] = rootEventName
        capRootNode.attrib["instance"] = "%s.1" % capName
        
        capLeafNode = ET.SubElement(capRootNode, "clade")
        #Add in the event name to identify with the species tree.
        capLeafNode.attrib["event"] = derivedEventName
        capLeafNode.attrib["instance"] = "%s.0" % capName
        
        #Update the global variables
        parentMap[capLeafNode] = capRootNode
        instanceNamesToCapNodes[capRootNode.attrib["instance"]] = capRootNode
        instanceNamesToCapNodes[capLeafNode.attrib["instance"]] = capLeafNode
        capNodesForEvents[rootEventName].add(capRootNode)
        capNodesForEvents[derivedEventName].add(capLeafNode)
        
        return capLeafNode
    
    addEventToEventTree_Index = [0]
    def addEventToEventTree(derivedEventName):
        #Create a new event in the event tree between the ancestral and derived events
        derivedEvent = eventNamesToEvents[derivedEventName]
        ancestralEvent = parentMap[derivedEvent]
        
        newEventName = "%s%i" % (options.uniqueNamePrefix, addEventToEventTree_Index[0])
        addEventToEventTree_Index[0] += 1
        newEvent = ET.SubElement(ancestralEvent, "clade")
        newEvent.attrib["event"] = newEventName
        ancestralEvent.remove(derivedEvent)
        newEvent.insert(0, derivedEvent)
        
        #Update the global variables
        parentMap[derivedEvent] = newEvent
        parentMap[newEvent] = ancestralEvent
        eventNamesToEvents[newEventName] = newEvent
        return newEvent
    
    operationIndex = [0]    
    def createOperation(configuration1, configuration2, configuration2IsAncestral):
        #Creates an operation structure for the operation.
        operation = ET.SubElement(reconstructionProblemTag.find("operations"), "operation")
        
        for capNode in configuration1[1:]:
            assert capNode.attrib["event"] == configuration1[0].attrib["event"]
        
        for capNode in configuration2[1:]:
            assert capNode.attrib["event"] == configuration2[0].attrib["event"]    
        
        #Give the operation an id.
        operation.attrib["operation_index"] = str(operationIndex[0])
        if configuration2IsAncestral:
            l = configuration1
        else:
            l = configuration1 + configuration2
        for capNode in l:
            capNode.attrib["operation_index"] = str(operationIndex[0])
        
        operationIndex[0] += 1
        
        #Write out the adjacencies for the derived configuration
        def fn(capNodes, configuration):
            configuration.attrib["event"] = capNodes[0].attrib["event"]
            seen = set()
            adjacencies = ET.SubElement(configuration, "adjacency_pairs")
            for capNode in capNodes:
                capNodeName = capNode.attrib["instance"]
                otherCapName = capNode.attrib["adjacency"]
                if (otherCapName, capNodeName) not in seen:
                    assert (capNodeName, otherCapName) not in seen
                    seen.add((capNodeName, otherCapName))
                    adjacency_pair = ET.SubElement(adjacencies, "adjacency_pair")
                    adjacency_pair.attrib["left"] = capNodeName
                    adjacency_pair.attrib["right"] = otherCapName
            return adjacencies
        
        configuration1Tag = ET.SubElement(operation, "configuration")
        fn(configuration1, configuration1Tag)
        
        configuration2Tag = ET.SubElement(operation, "configuration")
        fn(configuration2, configuration2Tag)
    
    def createPolarisedOperation(childCapNodes):
        createOperation(childCapNodes, [ parentMap[childCapNode] for childCapNode in childCapNodes ], True)
        
    def buildComponent(event, childCapNode, component, capNodesForEvent):
        #Builds the component associated the the operation
        adjacentChildCapNode = instanceNamesToCapNodes[childCapNode.attrib["adjacency"]]
        parentAdjacentChildCapNode = parentMap[adjacentChildCapNode]
        eventName = event.attrib["event"]
        
        assert childCapNode.attrib["event"] == adjacentChildCapNode.attrib["event"]
        assert mostAncestralEvent([ parentAdjacentChildCapNode.attrib["event"], eventName ]) == \
        parentAdjacentChildCapNode.attrib["event"]
        assert mostAncestralEvent([ adjacentChildCapNode.attrib["event"], eventName ]) == eventName
        assert (adjacentChildCapNode, childCapNode) not in component
        
        if (childCapNode, adjacentChildCapNode) not in component:
            component.append((childCapNode, adjacentChildCapNode))
            if parentAdjacentChildCapNode in capNodesForEvent:
                for nextChildCapNode in parentAdjacentChildCapNode.findall("clade"):
                    if nextChildCapNode != adjacentChildCapNode:
                        buildComponent(event, nextChildCapNode, \
                                       component, capNodesForEvent)
    
    def buildComponents(event, capNodesForEvent):
        #Build the local break point graph components for a cap node.
        seen = set()
        components = []
        for capNode in capNodesForEvent:
            if capNode not in seen:
                childCapNodes = list(capNode.findall("clade"))
                assert len(childCapNodes) == 2 or len(childCapNodes) == 1
                component = []
                buildComponent(event, childCapNodes[0], component, capNodesForEvent)
                if len(childCapNodes) == 2 and not isCycle(component): #Hash other child branch (not the end of a snake)
                    l = []
                    buildComponent(event, childCapNodes[1], l, capNodesForEvent)
                    assert not isCycle(l)
                    l.reverse()
                    assert not isCycle(component)
                    component = [ (capNode2, capNode1) for (capNode1, capNode2) in l ] + component
                    assert not isCycle(component)
                components.append(component)
                for childCap1, childCap2 in component:
                    seen.add(parentMap[childCap1])
                    seen.add(parentMap[childCap2])
        return components
    
    def addAdjacenciesForTrivialComponent(event, component):
        #Links the two vertices in a trivial component, that are adjacent in both children.
        assert len(component) == 1 or (isCycle(component) and len(component) == 2)
        #Get the nodes
        childCapNode1, childCapNode2 = component[0]
        parentCapNode1 = addAncestralEventToCap(childCapNode1, event.attrib["event"])
        parentCapNode2 = addAncestralEventToCap(childCapNode2, event.attrib["event"]) 
        linkCapNodes(parentCapNode1, parentCapNode2)
    
    def breakUnevenCycle(component):
        #Breaks a cycle of adjacencies of uneven degree by introducing an extra operation.
        
        assert isCycle(component) and len(component) % 2 == 1 #is odd numbered cycle
        #Get the nodes for the first link in the component chain
        childCapNode1, childCapNode2 = component[0]
        
        assert parentMap[childCapNode1].attrib["event"] == parentMap[childCapNode2].attrib["event"]
        assert childCapNode1.attrib["adjacency"] == childCapNode2.attrib["instance"]
        assert childCapNode2.attrib["adjacency"] == childCapNode1.attrib["instance"]
        assert childCapNode1.attrib["event"] == childCapNode2.attrib["event"]
        
        derivedEventName = childCapNode1.attrib["event"]
        
        #Create two caps whose youngest event is the derived event of the extra op and 
        #whose intermediate is the ancestral event
        childExtraCapNode1 = createExtraStubCap(derivedEventName) #Returns the child of the new stub cap
        childExtraCapNode2 = createExtraStubCap(derivedEventName)
        
        #Link the two caps in the derived state.
        linkCapNodes(childExtraCapNode1, childExtraCapNode2)
        
        operation = [ childCapNode1, childCapNode2, childExtraCapNode1, childExtraCapNode2 ]
         
        #Invent new event
        newEvent = addEventToEventTree(derivedEventName)
        
        #Add ancestral event nodes to the chosen node caps.
        childCapNode1 = addAncestralEventToCap(childCapNode1, newEvent.attrib["event"])
        childCapNode2 = addAncestralEventToCap(childCapNode2, newEvent.attrib["event"])
        childExtraCapNode1 = addAncestralEventToCap(childExtraCapNode1, newEvent.attrib["event"])
        childExtraCapNode2 = addAncestralEventToCap(childExtraCapNode2, newEvent.attrib["event"])
        
        #Link the ancestral events
        linkCapNodes(childExtraCapNode1, childCapNode1)
        linkCapNodes(childExtraCapNode2, childCapNode2)
        
        #Add the ops to the op list
        createPolarisedOperation(operation)
        
        #Remove the first link from the component and replace with the two new edges. And Done.
        return [ (childExtraCapNode2, childCapNode2) ] + \
        component[1:] + [ (childCapNode1, childExtraCapNode1) ]
        
    def extendEvenSnake(event, component):
        #Adds a stub to the end of a snake with an even number of edges, so that, by choosing the odd edges to break,
        #the operations does not create any null adjacencies.
        assert (not isCycle(component)) and len(component) % 2 == 0
        #Find the leaf event for the new stub
        l = []
        for capNode1, capNode2 in component[1::2]:
            assert capNode1.attrib["event"] == capNode2.attrib["event"]
            #l += [ capNode1.attrib["event"], capNode2.attrib["event"] ]
            l.append(capNode1.attrib["event"])
        derivedEventName = mostAncestralEvent(set(l))
        #Make a new cap
        extraCapNode = createExtraStubCap(derivedEventName)
        #Make a new branch on the child branch to connect to.
        childCapNode = component[0][0]
        assert childCapNode.attrib["event"] != event.attrib["event"]
        assert parentMap[childCapNode].attrib["event"] != event.attrib["event"]
        addAncestralEventToCap(childCapNode, event.attrib["event"])
        assert parentMap[childCapNode].attrib["event"] == event.attrib["event"]
        otherChildCapNode = makeChildCapNode(derivedEventName, parentMap[childCapNode])
        #Link the new cap node and the new stub
        linkCapNodes(extraCapNode, otherChildCapNode)
        #Add the cap link to the component
        return [ (extraCapNode, otherChildCapNode) ] + component
    
    def scheduleOperation(event, component):
        oddComponent = component[1::2]
        evenComponent = component[::2]
        
        def fn(componentConfiguration):
            #Find the switch point range for the chosen operation 
            l = []
            for capNode1, capNode2 in componentConfiguration:
                #l += [ capNode1.attrib["event"], capNode2.attrib["event"] ]
                assert capNode1.attrib["event"] == capNode2.attrib["event"]
                l.append(capNode1.attrib["event"])
            derivedEventName = mostAncestralEvent(set(l))
            
            #Link the derived nodes
            l = []
            for childCapNode1, childCapNode2 in componentConfiguration:
                #Make child events for the switch range
                childCapNode1 = addAncestralEventToCap(childCapNode1, derivedEventName)
                childCapNode2 = addAncestralEventToCap(childCapNode2, derivedEventName)
                #Link children
                linkCapNodes(childCapNode1, childCapNode2)
                l += [ (childCapNode1, childCapNode2) ]
            return l
        
        oddComponent = fn(oddComponent)
        configuration1 = []
        for childCapNode1, childCapNode2 in oddComponent:
            configuration1 += [ childCapNode1, childCapNode2 ]
        
        if event.attrib["event"] == rootEventName: #Make ambiguuity at the root
            evenComponent = fn(evenComponent)
            configuration2 = []
            for childCapNode1, childCapNode2 in evenComponent:
                configuration2 += [ childCapNode1, childCapNode2 ]
        
        for childCapNode1, childCapNode2 in evenComponent:
            parentCapNode1 = addAncestralEventToCap(childCapNode1, event.attrib["event"])
            parentCapNode2 = addAncestralEventToCap(childCapNode2, event.attrib["event"])
            linkCapNodes(parentCapNode1, parentCapNode2)
        
        if event.attrib["event"] == rootEventName: #Make ambiguuity at the root
            for childCapNode1, childCapNode2 in oddComponent:
                parentCapNode1 = addAncestralEventToCap(childCapNode1, event.attrib["event"])
                parentCapNode2 = addAncestralEventToCap(childCapNode2, event.attrib["event"])
                linkCapNodesSecondarily(parentCapNode1, parentCapNode2)
            createOperation(configuration1, configuration2, False)
        else:
            createPolarisedOperation(configuration1)
    
    #For each event in the event tree (in DFS order), ancestral graph to children. (recursive function)
    def traverseEventTree(event):
        if len(list(event.findall("clade"))) > 0: #Not a leaf event
            #Do the children first
            for childEvent in event.findall("clade"):
                traverseEventTree(childEvent)
        
            #For each unique component do history
            capNodesForEvent = capNodesForEvents[event.attrib["event"]]
            for component in buildComponents(event, capNodesForEvent):
                #Deal with components that do not involve an operation (trivial)
                if (len(component) == 1 and not isCycle(component)) or \
                (isCycle(component) and len(component) == 2):
                    addAdjacenciesForTrivialComponent(event, component)
                    continue
                
                #Handle uneven cycles by adding extra op, this also deals with self loops by creating (with cost of an op)
                if isCycle(component) and len(component) % 2 == 1: #is odd numbered cycle
                    component = breakUnevenCycle(component)
                    
                #Handle even chains by adding a stub cap to one end
                if not isCycle(component) and len(component) % 2 == 0:
                    component = extendEvenSnake(event, component)
                
                ####Now schedule the main operation
                scheduleOperation(event, component)
    
    traverseEventTree(reconstructionProblemTag.find("event_tree").\
                      find("phylogeny").find("clade"))
    
    logger.info("Finished doing the adjacency reconstruction")
    
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
