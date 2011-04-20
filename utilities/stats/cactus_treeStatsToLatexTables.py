#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script for generating latex tables from a cactus tree stats file.
"""

import xml.etree.ElementTree as ET
import os
import math

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import logger

def formatFloat(string, decimals=1):
    f = float(string)
    if f == 2147483647:
        return "NaN"
    return ("%." + str(decimals) + "f") % f

def writeDocumentPreliminaries(fileHandle):
    fileHandle.write("\\documentclass[11pt,oribibl]{llncs}\n")

    fileHandle.write("\\pagenumbering{arabic}\n")
    fileHandle.write("\\pagestyle{plain}\n")
    
    fileHandle.write("\\usepackage{epsfig}\n")
    fileHandle.write("\\usepackage{url}\n")
    fileHandle.write("\\usepackage{rotating}\n")
    
    fileHandle.write("\\usepackage{multirow}\n")
    
    fileHandle.write("\\setlength{\\evensidemargin}{0in}\n")
    fileHandle.write("\\setlength{\\oddsidemargin}{0in}\n")
    fileHandle.write("\\setlength{\\marginparwidth}{1in}\n")
    fileHandle.write("\\setlength{\\textwidth}{6.5in}\n")
    fileHandle.write("\\setlength{\\topmargin}{-0.5in}\n")
    fileHandle.write("\\setlength{\\textheight}{9in}\n")

    fileHandle.write("\\begin{document}\n")
    

def writeDocumentEnd(fileHandle):
    fileHandle.write("\\end{document}\n")

def writePreliminaries(columnNumber, fileHandle):
    fileHandle.write("\\begin{table}\n\\centering\n")
    fileHandle.write("\\begin{tabular}{" + ("|c"*columnNumber) + "|}\n")

def writeEnd(fileHandle, tableLabel, caption):
    fileHandle.write("\hline\n")
    fileHandle.write("\\end{tabular}\n")
    fileHandle.write("\caption{%s}\n" % caption)
    fileHandle.write("\label{%s}\n" % tableLabel)
    fileHandle.write("\end{table}\n")

def writeLine(columnNumber, rowNumber, entries, fileHandle):
    updatedEntries = []
    for name, x1, x2, y1, y2 in entries:
        columnNumber = y2 - y1 + 1
        updatedEntries.append((y1, x1, x2, name, columnNumber, y2 - y1 == 0))
        while y2 - y1 > 0: #Make multiple new entries:
            y1 += 1
            updatedEntries.append((y1, x1, x2, "", columnNumber, y2 - y1 == 0))
    updatedEntries.sort()
    fileHandle.write("\hline\n")
    start = True
    currentRow = 0
    row = []
    for y1, x1, x2, name, columnNumber, cLine in updatedEntries:
        if y1 != currentRow:
            fileHandle.write(" \\\\ %s\n" % " ".join([ "\\cline{%i-%i}" % (x3+1, x4+1) for x3, x4 in row ]))
            currentRow = y1
            row = []
        else:
            if not start:
                fileHandle.write(" & ")
            else:
                start = False
        if cLine:
            row.append((x1, x2))
        fileHandle.write("\multicolumn{%i}{|c|}{\multirow{%i}{*}{%s}}" % (x2-x1+1, columnNumber, name))
    fileHandle.write(" \\\\\n")

def writeRow(entries, fileHandle):
    fileHandle.write("\hline\n%s \\\\\n" % " & ".join(entries))

def writeFlowerTable(stats, fileHandle):
    columnNumber = 15
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Flowers", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Bp Size", 1, 1, 0, 1), 
                                 ("T. Flowers", 2, 4, 0, 0), 
                                 ("All", 2, 2, 1, 1), 
                                 ("Links", 3, 3, 1, 1), 
                                 ("Tangles", 4, 4, 1, 1), 
                                 ("Relative Entropy", 5, 7, 0, 0),
                                 ("P(X)/Z", 5, 5, 1, 1),
                                 ("Q(X)/Z", 6, 6, 1, 1),
                                 ("NRE", 7, 7, 1, 1),
                                 ("Children", 8, 10, 0, 0),
                                 ("Max", 8, 8, 1, 1),
                                 ("Avg.", 9, 9, 1, 1),
                                 ("Med.", 10, 10, 1, 1),
                                 ("Depth", 11, 14, 0, 0),
                                 ("Min", 11, 11, 1, 1),
                                 ("Max", 12, 12, 1, 1),
                                 ("Avg.", 13, 13, 1, 1),
                                 ("Med.", 14, 14, 1, 1)), fileHandle)
    for statNode, regionName in stats:
        #The nodes we use
        flowersNode = statNode.find("flowers")
        childrenNode = flowersNode.find("children")
        tangleChildrenNode = flowersNode.find("tangle_children")
        linkChildrenNode = flowersNode.find("link_children")
        depthNode = flowersNode.find("depths")
        relativeEntropyNode = statNode.find("relative_entropy_stats")
       
        #Flower/group type numbers
        totalFlowers = float(childrenNode.attrib["sum"])
        totalTangles = float(tangleChildrenNode.attrib["sum"])
        totalLinks = float(linkChildrenNode.attrib["sum"])
        
        totalSequenceLength = float(statNode.attrib["total_sequence_length"])
        
        #Relative entropy
        z = totalSequenceLength
        p = math.log(z)/math.log(2)
        nre = float(relativeEntropyNode.attrib["normalised_relative_entropy"])
        q = nre + p
        
        writeRow((regionName, 
                  formatFloat(totalSequenceLength, decimals=0),
                  formatFloat(1 + totalFlowers, decimals=0),
                  formatFloat(totalLinks, decimals=0),
                  formatFloat(totalTangles, decimals=0),
                  formatFloat(str(p), decimals=2), 
                  formatFloat(str(q), decimals=2), 
                  formatFloat(str(nre), decimals=2), 
                  formatFloat(childrenNode.attrib["max"], decimals=0),
                  formatFloat(childrenNode.attrib["avg"], decimals=2),
                  formatFloat(childrenNode.attrib["median"], decimals=0),
                  formatFloat(depthNode.attrib["min"], decimals=0),
                  formatFloat(depthNode.attrib["max"], decimals=0),
                  formatFloat(depthNode.attrib["avg"], decimals=2),
                  formatFloat(depthNode.attrib["median"], decimals=0)), fileHandle)
    writeEnd(fileHandle, "flowers_table", "Statistics on the flowers of the cactus trees. \
    Region: region name. \
    Bp size: total number of basepairs in the input sequences. \
    T. flowers: Total flowers in the cactus tree, either 'all', including all flowers in tree, 'links', including \
    only links groups or 'tangles', including only tangle groups. Note the sum \
    of link and tangle groups is equal to all minus one (for the root node).\
    Norm. relative entropy: Let $N$ be a flower in the set of all flowers $T$ in a cactus tree $X$. Furthermore let $N_0 \ldots N_{n-1}, N_{n}$ \
    denote the ancestral path of nodes from the root flower $N_0$ of the cactus tree to the flower $N_n$. Let $P(X) = \sum_{N_n \in T} |b(N_n)| \
    (log_2(|b(N_n)|) + \sum_{i=0}^{n-1} log_2(|c(N_i)|))$ and $Q(X) = Z log_2(Z)$, where $Z$ is the total number of basepairs in the input sequences, \
    $b(N)$ is the set of basepairs contained in blocks of the flower $N$, $|b(N)|$ is the the size of $b(N)$, $c(N)$ is the set of child flowers (direct descendants) \
    of $N$ and $| c(N) |$ is the size of $c(N)$. The total relative entropy is $P(X) - Q(X)$ and the normalised relative entropy (NRE) is $ ( P(X) - Q(X) ) / Z $. \
    The measure therefore reflects the balance of the tree.\
    Children: The children of a flower are its direct descendants flowers in the subsequent flower layer of the (multi layered) cactus tree. Results given for non-terminal flowers only. \
    Depth: The depth of a flower is the number of nodes (excluding itself) on the path from it to the root node. Results for terminal flowers only. \
    (A leaf flower is a terminal flower in the multi-layered cactus tree)")

def writeBlocksTable(stats, fileHandle):
    columnNumber = 15
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Blocks", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Min. Block Degree", 1, 1, 0, 1),
                                 ("Total", 2, 2, 0, 1), 
                                 ("Per Flower", 3, 5, 0, 0),
                                 ("Max", 3, 3, 1, 1),
                                 ("Avg.", 4, 4, 1, 1),
                                 ("Med.", 5, 5, 1, 1),
                                 ("Length", 6, 8, 0, 0),
                                 ("Max", 6, 6, 1, 1),
                                 ("Avg.", 7, 7, 1, 1),
                                 ("Med.", 8, 8, 1, 1),
                                 ("Degree", 9, 11, 0, 0),
                                 ("Max", 9, 9, 1, 1),
                                 ("Avg.", 10, 10, 1, 1),
                                 ("Med.", 11, 11, 1, 1),
                                 ("Coverage", 12, 14, 0, 0),
                                 ("Max", 12, 12, 1, 1),
                                 ("Avg.", 13, 13, 1, 1),
                                 ("Med.", 14, 14, 1, 1)), fileHandle)
    for statNode, regionName in stats:
        l = []
        l.append((regionName, 0, 0, 0, 1))
        i = 0
        for blocksNode in statNode.findall("blocks"):
            countsNode = blocksNode.find("counts")
            lengthsNode = blocksNode.find("lengths")
            degreesNode = blocksNode.find("leaf_degrees")
            coverageNode = blocksNode.find("leaf_coverage")
            l.append((blocksNode.attrib["minimum_leaf_degree"], 1, 1, i, i))
            j = 2
            for field, decimals in ((countsNode.attrib["sum"], 0),
                                   (countsNode.attrib["max"], 0), 
                                   (countsNode.attrib["avg"], 2),
                                   (countsNode.attrib["median"], 0),
                                   (lengthsNode.attrib["max"], 0),
                                   (lengthsNode.attrib["avg"], 2),
                                   (lengthsNode.attrib["median"], 0),
                                   (degreesNode.attrib["max"], 0),
                                   (degreesNode.attrib["avg"], 2),
                                   (degreesNode.attrib["median"], 0),
                                   (coverageNode.attrib["max"], 0),
                                   (coverageNode.attrib["avg"], 2),
                                   (coverageNode.attrib["median"], 0)):
                l.append((formatFloat(field, decimals=decimals), j, j, i, i))
                j+=1
            i+=1
        writeLine(columnNumber, 2, l, fileHandle)
    writeEnd(fileHandle, "blocks_table", "Statistics on the blocks of the cactus trees. \
    Region: region name. \
    Min. Block Degree: the minimum number of leaf sequences in a block considered this round. \
    Total: total number of blocks in the cactus tree. \
    Per Flower: numbers of blocks in the child chains of each flower, excluding terminal flowers. \
    Length: number of basepairs in a block. \
    Degree: number of leaf sequences in a block. \
    Coverage: block's length * degree.")

def writeChainsTable(stats, fileHandle):
    columnNumber = 15
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Chains", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Type", 1, 1, 0, 1), 
                                
                                ("Total", 2, 2, 0, 1),
                                 
                                 ("Per Flower", 3, 5, 0, 0), 
                                 ("Max", 3, 3, 1, 1),
                                 ("Avg.", 4, 4, 1, 1),
                                 ("Med.", 5, 5, 1, 1),
                                 
                                 ("Link Number", 6, 8, 0, 0), 
                                 ("Max", 6, 6, 1, 1),
                                 ("Avg.", 7, 7, 1, 1),
                                 ("Med.", 8, 8, 1, 1),
                                 
                                 ("Block Bp length", 9, 11, 0, 0), 
                                 ("Max", 9, 9, 1, 1),
                                 ("Avg.", 10, 10, 1, 1),
                                 ("Med.", 11, 11, 1, 1),
                                 
                                 ("Instance Length", 12, 14, 0, 0), 
                                 ("Max", 12, 12, 1, 1),
                                 ("Avg.", 13, 13, 1, 1),
                                 ("Med.", 14, 14, 1, 1)), fileHandle)
    
    for statNode, regionName in stats:
        l = [ (regionName, 0, 0, 0, 1) ]
        i = 0
        for chainsNode, blockType in ((statNode.findall("chains")[0], "all"), (statNode.findall("chains")[1], "$>=2$ B.")):
            l.append((blockType, 1, 1, i, i))
            l.append((formatFloat(chainsNode.find("counts").attrib["sum"], decimals=0), 2, 2, i, i))
            j = 3
            for typeNode in (chainsNode.find("counts"), chainsNode.find("base_block_lengths"), chainsNode.find("link_numbers"), chainsNode.find("avg_instance_base_length")):
                l.append((formatFloat(typeNode.attrib["max"], decimals=0), j, j, i, i))
                l.append((formatFloat(typeNode.attrib["avg"], decimals=2), j+1, j+1, i, i))
                l.append((formatFloat(typeNode.attrib["median"], decimals=0), j+2, j+2, i, i))
                j += 3
            i += 1 
        writeLine(columnNumber, 2, l, fileHandle)
    writeEnd(fileHandle, "chains_table", "Statistics on the chains of the cactus trees. \
    Region: region name. \
    Type: categories of chains, either `all', which includes all chains or `$>=2$ B.', \
    which includes only chains containing a minimum of two blocks. \
    Total: total number of chains in the cactus tree. \
    Per Flower: numbers of child chains in each non-terminal flower. \
    Link Number: number of links in chain. \
    Block Bp length: number of basepairs in blocks of chain. \
    Instance length: average number of basepairs in an instance of the chain, including both its blocks and intervening links.")

def writeTerminalGroupsTable(stats, fileHandle):
    columnNumber = 4
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Terminal Groups", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("End number", 1, 2, 0, 0), 
                                ("All", 1, 1, 1, 1), 
                                ("No Free Stubs", 2, 2, 1, 1), 
                                 ("End degrees ", 3, 3, 0, 1),), fileHandle)
    for statNode, regionName in stats:
        netNode = statNode.find("nets")
        assert(netNode != None)
        l = [ (regionName, 0, 0, 0, 0) ]
        l.append((formatFloat(netNode.find("total_end_numbers_per_terminal_group").attrib["avg"], decimals=2), 1, 1, 0, 0))
        l.append((formatFloat(netNode.find("total_non_free_stub_end_numbers_per_terminal_group").attrib["avg"], decimals=2), 2, 2, 0, 0))
        l.append((formatFloat(netNode.find("end_degrees_per_terminal_group").attrib["avg"], decimals=2), 3, 3, 0, 0))
        #l.append((formatFloat(netNode.find("total_groups_per_net").attrib["avg"], decimals=2), 4, 4, 0, 0))
        writeLine(4, 1, l, fileHandle)
                     
    writeEnd(fileHandle, "nets_table", "Statistics on the terminal groups of the cactus tree. \
    Region: region name. \
    End number: ends per terminal group. \
    End degrees: the number of distinct ends an end in a terminal group is adjacent to.\
    Groups/net: the number of terminal groups per net")
    
def writeFacesTable(stats, fileHandle):
    columnNumber = 13
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Faces", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Group", 1, 1, 0, 1), 
                                 ("Per Group", 2, 4, 0, 0), 
                                 ("Max", 2, 2, 1, 1),
                                 ("Avg.", 3, 3, 1, 1),
                                 ("Med.", 4, 4, 1, 1),
                                 ("Cardinality", 5, 7, 0, 0), 
                                 ("Max", 5, 5, 1, 1),
                                 ("Avg.", 6, 6, 1, 1),
                                 ("Med.", 7, 7, 1, 1),
                                 ("Breakpoint reuse", 8, 10, 0, 0), 
                                 ("Max", 8, 8, 1, 1),
                                 ("Avg.", 9, 9, 1, 1),
                                 ("Med.", 10, 10, 1, 1),
                                 ("Prop. Reg.", 11, 11, 0, 1), 
                                 ("Prop. Can.", 12, 12, 0, 1)), fileHandle)
    for statNode, regionName in stats:
        facesNodes = statNode.findall("faces")
        l = [ (regionName, 0, 0, 0, 3) ]
        i = 0
        for groupType in ("all", "tangles", "links"):
            l.append((groupType, 1, 1, i, i))
            l.append((formatFloat(facesNodes[i].find("number_per_group").attrib["max"], decimals=0), 2, 2, i, i))
            l.append((formatFloat(facesNodes[i].find("number_per_group").attrib["avg"], decimals=2), 3, 3, i, i))
            l.append((formatFloat(facesNodes[i].find("number_per_group").attrib["median"], decimals=0), 4, 4, i, i))
            l.append((formatFloat(facesNodes[i].find("cardinality").attrib["max"], decimals=0), 5, 5, i, i))
            l.append((formatFloat(facesNodes[i].find("cardinality").attrib["avg"], decimals=2), 6, 6, i, i))
            l.append((formatFloat(facesNodes[i].find("cardinality").attrib["median"], decimals=0), 7, 7, i, i))
            l.append((formatFloat(facesNodes[i].find("faces_per_face_associated_end").attrib["max"], decimals=0), 8, 8, i, i))
            l.append((formatFloat(facesNodes[i].find("faces_per_face_associated_end").attrib["avg"], decimals=2), 9, 9, i, i))
            l.append((formatFloat(facesNodes[i].find("faces_per_face_associated_end").attrib["median"], decimals=0), 10, 10, i, i))
            l.append((formatFloat(facesNodes[i].find("is_regular").attrib["avg"], decimals=2), 11, 11, i, i))
            l.append((formatFloat(facesNodes[i].find("is_canonical").attrib["avg"], decimals=0), 12, 12, i, i))
            i += 1
        writeLine(columnNumber, 3, l, fileHandle)
    
    writeEnd(fileHandle, "faces_table", "Statistics on the faces of the terminal AVGs. \
    Region: region name. \
    Group: either `links', `tangles' or both (`all') groups. \
    Per Group: numbers of faces in a group. \
    Cardinality: the cardinality of faces. \
    Breakpoint reuse: Let E be the set of ends in AVGs such that each end in E has involved in atleast one\
    nontrivial face. Breakpoint reuse is the number of distinct faces associated with members of E. \
    Prop. Reg.: Proportion of non-trivial faces which are regular.\
    Prop. Can.: Proportion of non-trivial faces which are canononical.")

def writeReferenceTable(stats, fileHandle):
    columnNumber = 10
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Reference", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Type", 1, 1, 0, 1), 
                                
                                 ("P.chr.", 2, 5, 0, 0), 
                                 ("Total", 2, 2, 1, 1),
                                 ("Avg.", 3, 3, 1, 1),
                                 ("Med.", 4, 4, 1, 1),
                                 ("Max", 5, 5, 1, 1),
                                 
                                 ("Contigs", 6, 9, 0, 0), 
                                 ("Total", 6, 6, 1, 1),
                                 ("Avg.", 7, 7, 1, 1),
                                 ("Med.", 8, 8, 1, 1),
                                 ("Max", 9, 9, 1, 1),), fileHandle)
    for statNode, regionName in stats:
        referenceNodes = statNode.findall("reference2")
        l = [ (regionName, 0, 0, 0, len(referenceNodes)-1) ]
        for i in xrange(len(referenceNodes)):
            l.append((referenceNodes[i].attrib["method"], 1, 1, i, i))
            
            l.append((formatFloat(referenceNodes[i].find("top_level_pseudo_chromosome_lengths").attrib["total"], decimals=0), 2, 2, i, i))
            l.append((formatFloat(referenceNodes[i].find("top_level_pseudo_chromosome_lengths").attrib["avg"], decimals=2), 3, 3, i, i))
            l.append((formatFloat(referenceNodes[i].find("top_level_pseudo_chromosome_lengths").attrib["median"], decimals=0), 4, 4, i, i))
            l.append((formatFloat(referenceNodes[i].find("top_level_pseudo_chromosome_lengths").attrib["max"], decimals=0), 5, 5, i, i))
            
            l.append((formatFloat(referenceNodes[i].find("contig_lengths_filtered").attrib["total"], decimals=2), 6, 6, i, i))
            l.append((formatFloat(referenceNodes[i].find("contig_lengths_filtered").attrib["avg"], decimals=0), 7, 7, i, i))
            l.append((formatFloat(referenceNodes[i].find("contig_lengths_filtered").attrib["median"], decimals=0), 8, 8, i, i))
            l.append((formatFloat(referenceNodes[i].find("contig_lengths_filtered").attrib["max"], decimals=2), 9, 9, i, i))
            
        writeLine(columnNumber, len(referenceNodes), l, fileHandle)
    
    writeEnd(fileHandle, "reference_table", "Statistics on the reference lengths of pseudo-chromosomes and contigs.\
    Region: region name. \
    Type: reference algorithm. \
    Columns (abbreviated 'P.chr'.):base lengths of pseudo chromosomes.\
    Contigs: base lengths of contigs.")

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = getBasicOptionParser("usage: %prog [options] treeStatsFiles", "%prog 0.1")
    
    parser.add_option("--outputFile", dest="outputFile", 
                      help="File to put the latex tables in.")

    options, args = parseBasicOptions(parser)
        
    logger.info("Parsed arguments")
    
    assert options.outputFile != None
    
    ##########################################
    #Get the input data etc.
    ##########################################
    
    assert len(args) % 2 == 0
    stats = [ (ET.parse(statsFile).getroot(), regionName) for statsFile, regionName in zip(args[::2], args[1::2]) ] 
    fileHandle = open(options.outputFile, "w")
    
    ##########################################
    #Make the document
    ##########################################
    
    writeDocumentPreliminaries(fileHandle)
    writeFlowerTable(stats, fileHandle)
    writeBlocksTable(stats, fileHandle)
    writeChainsTable(stats, fileHandle)
    writeTerminalGroupsTable(stats, fileHandle)
    writeFacesTable(stats, fileHandle)
    writeReferenceTable(stats, fileHandle)
    writeDocumentEnd(fileHandle)
    
    ##########################################
    #Cleanup
    ##########################################
    
    fileHandle.close()
    
if __name__ == "__main__":
    main()
