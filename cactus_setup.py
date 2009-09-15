#!/usr/bin/env python

"""Script for creating the top level reconstruction tags and first problem node for a 
recursive reconstruction.
"""

import os
import sys
from sonLib.bioio import logger
import xml.etree.ElementTree as ET

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import system
from sonLib.bioio import getTempFile

#Tree stuff
from sonLib.bioio import newickTreeParser

def main():
    ##########################################
    #Construct the arguments.
    ##########################################  
    
    parser = getBasicOptionParser("usage: %prog [options] [sequences in fasta format]", "%prog 0.1")
    
    parser.add_option("--speciesTree", dest="speciesTree", 
                      help="Newick tree string of species tree with which to associate the caps")
    
    parser.add_option("--reconstructionTree", dest="reconstructionTree", 
                      help="Directory to contain the reconstruction tree \n\
(setup will create this direction and a file in this directory called reconstructionProblem.xml)")
    
    parser.add_option("--uniqueNamePrefix", dest="uniqueNamePrefix", 
                      help="An alpha-numeric prefix which, when appended with any alpha-numeric characters is guaranteed to produce a unique name",
                      default=None)
    
    options, args = parseBasicOptions(parser)
    logger.info("Parsed arguments")
   
    ##########################################
    #Create the XML for the top level reconstruction problem
    ##########################################
    
    root = ET.Element("reconstruction_problem")
    
    assert options.speciesTree != None
    assert options.reconstructionTree != None
    if os.path.isdir(options.reconstructionTree):
        system("rm -rf %s" % options.reconstructionTree)
    os.mkdir(options.reconstructionTree)
    
    logger.info("Built the reconstruction directory")
    
    ##########################################
    #The parent structure to hold the map of fasta headers and sequence/contig names + event names.
    ##########################################
    
    nameMap = ET.Element("name_map")
    eventMap = ET.SubElement(nameMap, "event_map")
    fastaMap = ET.SubElement(nameMap, "fasta_map")
    
    eventNames = set() 
    def addEventToEventMap(originalName, reconName):
        assert originalName not in eventNames #We have seen a duplicated input event name in the input tree.
        eventNames.add(originalName)
        ET.SubElement(eventMap, "event", 
                  { "original_name":originalName, "event":reconName })
    
    fastaHeaders = set() 
    def addFastaToFastaMap(fastaHeader, contigName, eventName, seqFile):
        seqFile = os.path.split(seqFile)[1] #Use only file name for seq file
        assert fastaHeader not in fastaHeaders #Check we can identify the fasta header.
        fastaHeaders.add(fastaHeader)
        ET.SubElement(fastaMap, "fasta", 
                          { "header":fastaHeader, "contig":contigName, "event":eventName, "sequence_file":seqFile })
      
    ##########################################
    #Function for getting unique names.
    ##########################################

    nameIndex = [0]
    def getUniqueName():
        """Gets a unique name for the element
        """
        name = options.uniqueNamePrefix + str(nameIndex[0])
        nameIndex[0] += 1
        return name
    
    ##########################################
    #Create the event tree, using the species tree
    ##########################################
    
    newickTree = newickTreeParser(options.speciesTree)
    
    eventTree = ET.SubElement(root, "event_tree")
    tree = ET.SubElement(eventTree, "phylogeny")
    rootBranch = ET.SubElement(tree, "clade")
    #Add in the event name
    rootEventName = getUniqueName() #This is the root event.
    rootBranch.attrib["event"] = rootEventName 
    addEventToEventMap("root", rootEventName)
   
    #Add in the branch length
    rootBranch.attrib["branch_length"] = str(sys.maxint)
    
    eventNameMap = {}
    def fn(newickTree, node):
        """Function to convert newick tree into phylogeny.
        """
        if newickTree != None:
            node = ET.SubElement(node, "clade")
            #Add in the event name
            if newickTree.iD == None:
                node.attrib["event"] = getUniqueName()
            else:
                eventName = getUniqueName()
                node.attrib["event"] = eventName
                #Add to fasta header map
                addEventToEventMap(newickTree.iD, eventName)
                eventNameMap[newickTree.iD] = eventName
            #Add in the branch length
            assert newickTree.distance != None
            node.attrib["branch_length"] = str(float(newickTree.distance))
            #Call recursively
            fn(newickTree.left, node)
            fn(newickTree.right, node)
    
    fn(newickTree, rootBranch)
    
    def fn2(tree):
        """Gets the leaf event names in left to right order.
        """
        if tree.internal:
            return fn2(tree.left) + fn2(tree.right)
        else:
            return (eventNameMap[tree.iD],)
    orderedEventNames = fn2(newickTree)
    
    logger.info("Created the event tree")
    
    ##########################################
    #Create the sequences, caps and strings
    ##########################################
    
    sequences = ET.SubElement(root, "sequences")
    ET.SubElement(root, "atoms")
    caps = ET.SubElement(root, "caps")
    strings = ET.SubElement(root, "strings")
    ET.SubElement(root, "adjacency_components")
    ET.SubElement(root, "chains")
    ET.SubElement(root, "operations")
    
    sequences.attrib['alphabet'] = 'DNA'
    
    logger.info("Constructed the basic elements")
    
    ##########################################
    #Process the fasta files
    ##########################################
    
    fastaProcess_sequenceNumber = [0]
    def fastaProcess(sequenceFile):
        tempFile = getTempFile(suffix=".txt", rootDir=options.tempDirRoot)
        system("cactus_setup_processFasta %s %s %s %i" % (options.reconstructionTree, tempFile, sequenceFile, fastaProcess_sequenceNumber[0]))
        fileHandle = open(tempFile, 'r')
        line = fileHandle.readline()
        while line != '':
            fastaHeader = line[:-1]
            sequenceFile = fileHandle.readline()[:-1]
            assert sequenceFile != ''
            length = int(fileHandle.readline())
            yield fastaHeader, sequenceFile, length
            fastaProcess_sequenceNumber[0] += 1
            line = fileHandle.readline()
        fileHandle.close()
        os.remove(tempFile)
               
    def processSequence(parentSequenceFile, eventName, capCount):
        for fastaHeader, sequenceFile, sequenceLength in fastaProcess(parentSequenceFile):
            #Get a unique a fasta name.
            contigName = getUniqueName()
            #Add to fasta header map
            addFastaToFastaMap(fastaHeader, contigName, eventName, parentSequenceFile)
            
            #Make sequence
            sequence = ET.SubElement(sequences, "sequence")
            sequence.attrib['contig'] = contigName
            sequence.attrib['length'] = str(sequenceLength)
            sequence.attrib['event'] = eventName
            sequence.attrib['sequence_file'] = sequenceFile
            
            #Make base sequence
            string = ET.SubElement(strings, "string")
            string.attrib["contig"] = contigName
            string.attrib["start"] = "0"
            string.attrib["length"] = str(sequenceLength)
            
            other = "$%i.0" % (capCount+1)
            for attrib in [ "left", "right" ]:
                #Make cap name / cap instance name
                capName = "$%i" % capCount
                capCount += 1
                capInstanceName = "%s.0" % capName
                #The cap
                cap = ET.SubElement(caps, "cap")
                cap.text = capName
                
                capTree = ET.SubElement(cap, "phylogeny")
                capRootBranch = ET.SubElement(capTree, "clade")
                #Add in the event name to identify with the species tree.
                capRootBranch.attrib["event"] = rootEventName
                capRootBranch.attrib["instance"] = "%s.1" % capName
                
                capLeafBranch = ET.SubElement(capRootBranch, "clade")
                #Add in the event name to identify with the species tree.
                capLeafBranch.attrib["event"] = eventName
                capLeafBranch.attrib["instance"] = "%s.0" % capName
                capLeafBranch.attrib["adjacency"] = other
                other = capLeafBranch.attrib["instance"]
                
                #The string
                string.attrib[attrib] = capInstanceName
        return capCount
    
    ###Process fasta files.
    capCount = 0
    assert len(args) == len(orderedEventNames) #Check the number of leaves in the species tree is the same as the number of sequence arguments.
    for sequenceThing, eventName in zip(args, orderedEventNames):
        if os.path.isdir(sequenceThing):
            for relSequenceFile in os.listdir(sequenceThing):
                sequenceFile = os.path.join(sequenceThing, relSequenceFile)
                capCount = processSequence(sequenceFile, eventName, capCount)
        else:
            capCount = processSequence(sequenceThing, eventName, capCount)               
    
    logger.info("Processed the sequences")
    
    ##########################################
    #Write the reconstruction_problem
    ##########################################
    
    fileHandle = open(os.path.join(options.reconstructionTree, "reconstructionProblem.xml"), 'w')
    tree = ET.ElementTree(root)
    tree.write(fileHandle)
    fileHandle.close()
    
    logger.info("Written out the reconstruction problem")
    
    ##########################################
    #Write the fasta map
    ##########################################
    
    fileHandle = open(os.path.join(options.reconstructionTree, "fastaMap.xml"), 'w')
    tree = ET.ElementTree(nameMap)
    tree.write(fileHandle)
    fileHandle.close()
    
    logger.info("Written out the fastaMap.xml")
   
    return 0

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
