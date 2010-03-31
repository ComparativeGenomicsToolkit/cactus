#!/usr/bin/python

"""Creates a graphviz format graph document of a reconstruction tree; the graph can be visualised with dot.
"""

import sys
import os
import re
import math
import xml.etree.ElementTree as ET
import random

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
    colours = [ "red", "blue", "green", "yellow", "cyan", "magenta", "orange", "purple", "brown", "black", "grey80" ]
    return colours[getColour_Index[0] % len(colours)]

def proportionOfTotalSequenceLength(reconNode):
    """Returns the total sequence length the reconstruction problem covers.
    """
    totalPossibleSequenceLength = 0.0
    totalSequenceLength = 0.0
    for sequence in reconNode.find("sequences").findall("sequence"):
        totalPossibleSequenceLength += float(sequence.attrib["length"])
    for string in reconNode.find("strings").findall("string"):
        totalSequenceLength += float(string.attrib["length"])
    return totalSequenceLength/totalPossibleSequenceLength

def calculateAtomStats(reconNode):
    """Calculates the number, average length and depth of a reconstruction problem,
    and returns a tuple of the results.
    """
    atomNumber = float(len(reconNode.find("atoms").findall("atom")))
    atomTotalLength = 0.0
    atomInstanceNumber = 0.0
    for atom in reconNode.find("atoms").findall("atom"):
        atomTotalLength += float(atom.attrib["length"])
        atomInstanceNumber += len(atom.find("atom_instances").findall("atom_instance"))
    if atomNumber > 0:
        return atomNumber, atomTotalLength/atomNumber, atomInstanceNumber/atomNumber
    return 0.0, 0.0, 0.0

def addNodeToGraph2(nodeName, graphFileHandle, scalingFactor, label):
    """Adds a node to the graph.
    """
    height = 1
    width = 0.1
    if scalingFactor >= 1:
        height = 4 * math.sqrt(scalingFactor)
        width = 0.25 * math.sqrt(scalingFactor)
    addNodeToGraph(nodeName, graphFileHandle, label=label, width=width, height=height, shape="box", colour="green", fontsize=14)

def processReconstructionProblem(absPathPrefix, reconstructionProblem, parentNodeName, parentEdgeColour, graphFileHandle, nodesProportionalTo):
    """Processes a reconstruction tree.
    """
    logger.info("Processing the reconstruction problem: %s" % reconstructionProblem)
    #Load the reconstruction problem
    reconNode = ET.parse(os.path.join(absPathPrefix, reconstructionProblem)).getroot()
    #Make the name of the node
    nodeName = getNodeName()
    atomNumber, averageAtomLength, averageAtomDepth = calculateAtomStats(reconNode)
    propOfTotalSequenceLength = proportionOfTotalSequenceLength(reconNode)
    if nodesProportionalTo == "atoms":
        scalingFactor = atomNumber
    elif nodesProportionalTo == "sequences":
        scalingFactor = propOfTotalSequenceLength * 200
    else:
        raise RuntimeError("Nodes to proportional to, unrecognised option: %s" % nodesProportionalTo)
    label = "%.1f %.1f %.1f %.3f" % (atomNumber, averageAtomLength, averageAtomDepth, propOfTotalSequenceLength)
    addNodeToGraph2(nodeName, graphFileHandle, scalingFactor, label)
    #Write the edge to the parent
    if parentNodeName != None:
        addEdgeToGraph(parentNodeName, nodeName, graphFileHandle, colour=parentEdgeColour, length=10.5, weight=1.0, dir="forward")
    edgeColour = getColour()
    #Now iterate through children
    for adjacencyComponent in reconNode.find("adjacency_components").findall("adjacency_component"):
        childReconstructionProblem = adjacencyComponent.attrib["child_file"]
        processReconstructionProblem(absPathPrefix, childReconstructionProblem, nodeName, 
                                     edgeColour, graphFileHandle, nodesProportionalTo)

def makeReconstructionGraph(absPathPrefix, reconstructionProblem, graphFileName, nodesProportionalTo):
    """Function to create a graph viz dot file representation of a reconstruction problem.
    """
    graphFileHandle = open(graphFileName, 'w')
    setupGraphFile(graphFileHandle)
    processReconstructionProblem(absPathPrefix, reconstructionProblem, None, None, graphFileHandle, nodesProportionalTo)
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
    parser.add_option("--nodesProportionalTo", dest="nodesProportionalTo", default="atoms", help="Size nodes are proportional to: either 'atoms'/'sequences'") 
    
    options, args = parseBasicOptions(parser)

    logger.info("Parsed arguments")
    
    makeReconstructionGraph(options.absPathPrefix, options.reconstructionProblem, options.graphFileName, options.nodesProportionalTo)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
