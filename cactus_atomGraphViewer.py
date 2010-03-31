#!/usr/bin/python

"""Creates a graphviz format graph document of a reconstruction problem as a master break point graph; the graph can be visualised with neato.
"""

import sys
import os
import re
import math
import xml.etree.ElementTree as ET

from sonLib.bioio import logger
from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions

from sonLib.bioio import setupGraphFile
from sonLib.bioio import finishGraphFile
from sonLib.bioio import addNodeToGraph
from sonLib.bioio import addEdgeToGraph

getNodeName_Integer = [0]
def getNodeName():
    """Gets a name for a node.
    """
    nodeName = "n%in" % getNodeName_Integer[0]
    getNodeName_Integer[0] += 1
    return nodeName

getColour_Index = [0]
def getColour():
    """Returns a valid colour.
    """
    getColour_Index[0] += 1
    colours = [ "red", "blue", "green", "yellow", "cyan", "magenta", "orange", "purple", "brown", "grey80" ]
    return colours[getColour_Index[0] % len(colours)]

def addCapNodeAdjacencyEdges(capNode, nodeNames, capColours, graphFileHandle, includeInternalAdjacencies):
    """Adds adjacencies between the cap nodes.
    """
    capNodeName = nodeNames[capNode.text]
    if capNode.text in capColours:
        capColour = capColours[capNode.text]
    else:
        capColour = "red"
    def fn(capInstanceNode):
        adjacentCap = capInstanceNode.attrib["adjacency"].split(".")[0]
        if adjacentCap in nodeNames:
            adjacentCapNodeName = nodeNames[adjacentCap]
            if includeInternalAdjacencies or len(capInstanceNode.findall("clade")) == 0:
                if capNodeName < adjacentCapNodeName:
                    addEdgeToGraph(capNodeName, adjacentCapNodeName, graphFileHandle, 
                                   colour=capColour, length=2.5, weight=100, dir="none")
        #Add the child instance adjacencies in
        for childCapInstanceNode in capInstanceNode.findall("clade"):
            fn(childCapInstanceNode)
    fn(capNode.find("phylogeny").find("clade"))

def linkAtomCaps(atomNode, nodeNames, graphFileHandle, includeInternalAdjacencies):
    """Links the atom ends with sequence edges.
    """
    capNode1, capNode2 = atomNode.findall("cap")
    for atomInstance in atomNode.find("atom_instances").findall("atom_instance"): 
        if includeInternalAdjacencies or len(atomInstance.findall("coordinates")) == 1:
            addEdgeToGraph(nodeNames[capNode1.text], nodeNames[capNode2.text], graphFileHandle, 
                           colour="black", length=2.5, weight=100, dir="none")

def makeReconstructionGraph(absPathPrefix, reconstructionProblem, graphFileName, includeCaps=False, includeInternalAdjacencies=False):
    """Processes a reconstruction problem.
    """
    logger.info("Processing the reconstruction problem: %s" % reconstructionProblem)
    
    graphFileHandle = open(graphFileName, 'w')
    setupGraphFile(graphFileHandle)
    
    #Load the reconstruction problem
    reconNode = ET.parse(os.path.join(absPathPrefix, reconstructionProblem)).getroot()
    nodeNames = {}
    capColours = {}
    #Get the cap nodes.
    if includeCaps:
        capNodes = reconNode.find("caps").findall("cap")
    else:
        capNodes = []
    for atomNode in reconNode.find("atoms").findall("atom"):
        capNodes += atomNode.findall("cap")
    
    #Make the colours for each adjacency component
    for adjacencyComponent in reconNode.find("adjacency_components").findall("adjacency_component"):
        colour = getColour()
        for capNodeName in adjacencyComponent.text.split():
            capColours[capNodeName] = colour
    
    #Get the cap nodes.
    for capNode in capNodes:
        nodeNames[capNode.text] = getNodeName()
        addNodeToGraph(nodeNames[capNode.text], graphFileHandle, label=capNode.text)
    
    logger.info("Named all the caps nodes.")
    
    #Add the adjacency edges between the caps.
    for capNode in capNodes:
        addCapNodeAdjacencyEdges(capNode, nodeNames, capColours, graphFileHandle, includeInternalAdjacencies)
    
    logger.info("Added the adjacency edges.")
    
    #Link the atom nodes
    for atomNode in reconNode.find("atoms").findall("atom"):
        linkAtomCaps(atomNode, nodeNames, graphFileHandle, includeInternalAdjacencies)
    
    logger.info("Added the sequence edges.")
    
    finishGraphFile(graphFileHandle)
    graphFileHandle.close()

def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = getBasicOptionParser("usage: %prog [options]", "%prog 0.1")
    
    parser.add_option("--absolutePathPrefix", dest="absPathPrefix", help="Path to the base directory of reconstruction tree") 
    parser.add_option("--reconstructionProblem", dest="reconstructionProblem", help="The relative path to the file containing the reconstruction problem.") 
    parser.add_option("--graphFile", dest="graphFileName", help="The file to write the output graph in") 
    parser.add_option("--includeCaps", dest="includeCaps", default=False, action="store_true", help="Include the caps in the atom graph") 
    parser.add_option("--includeInternalAdjacencies", dest="includeInternalAdjacencies", default=False, action="store_true", help="Include the internal adjacencies in the atom graph") 
    
    options, args = parseBasicOptions(parser)

    logger.info("Parsed arguments")
    
    makeReconstructionGraph(options.absPathPrefix, options.reconstructionProblem, options.graphFileName, options.includeCaps, options.includeInternalAdjacencies)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
