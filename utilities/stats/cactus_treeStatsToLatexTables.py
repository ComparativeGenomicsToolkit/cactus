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

def writeNetTable(stats, fileHandle):
    columnNumber = 15
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Nets", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Bp Size", 1, 1, 0, 1), 
                                 ("T. Nets", 2, 4, 0, 0), 
                                 ("All", 2, 2, 1, 1), 
                                 ("Nets", 3, 3, 1, 1), 
                                 ("Chains", 4, 4, 1, 1), 
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
        netsNode = statNode.find("nets")
        totalNets = float(netsNode.attrib["total_net_number"])
        nonTrivialGroups = float(statNode.find("non_trivial_groups").attrib["total"])
        relativeEntropyNode = statNode.find("relative_entropy_stats")
        largestChildNode = statNode.find("largest_child")
        z = float(statNode.attrib["total_sequence_length"])
        p = math.log(z)/math.log(2)
        nre = float(relativeEntropyNode.attrib["normalised_relative_entropy"])
        q = nre + p
        writeRow((regionName, 
                  formatFloat(statNode.attrib["total_sequence_length"], decimals=0),
                  formatFloat(str(totalNets), decimals=0),
                  formatFloat(str(totalNets - nonTrivialGroups - 1), decimals=0),
                  formatFloat(str(nonTrivialGroups), decimals=0),
                  formatFloat(str(p), decimals=2), 
                  formatFloat(str(q), decimals=2), 
                  formatFloat(str(nre), decimals=2), 
                  formatFloat(netsNode.attrib["max_children"], decimals=0),
                  formatFloat(netsNode.attrib["avg_children"], decimals=2),
                  formatFloat(netsNode.attrib["median_children"], decimals=0),
                  formatFloat(netsNode.attrib["min_depth"], decimals=0),
                  formatFloat(netsNode.attrib["max_depth"], decimals=0),
                  formatFloat(netsNode.attrib["avg_depth"], decimals=2),
                  formatFloat(netsNode.attrib["median_depth"], decimals=0)), fileHandle)
    writeEnd(fileHandle, "nets_table", "Statistics on the nets of the cactus trees. \
    Region: region name. \
    Bp size: total number of basepairs in the input sequences. \
    T. nets: Total nets in the cactus tree, either 'all', including all nets in tree, 'nets', including \
    only net contained groups or 'chains', including only chain contained groups.. Note the sum \
    of net and chain contained groups is equal to all minus one (for the root node).\
    Norm. relative entropy: Let $N$ be a net in the set of all nets $T$ in a cactus tree $X$. Furthermore let $N_0 \ldots N_{n-1}, N_{n}$ \
    denote the ancestral path of nodes from the root net $N_0$ of the cactus tree to the net $N_n$. Let $P(X) = \sum_{N_n \in T} |b(N_n)| \
    (log_2(|b(N_n)|) + \sum_{i=0}^{n-1} log_2(|c(N_i)|))$ and $Q(X) = Z log_2(Z)$, where $Z$ is the total number of basepairs in the input sequences, \
    $b(N)$ is the set of basepairs contained in blocks of the net $N$, $|b(N)|$ is the the size of $b(N)$, $c(N)$ is the set of child nets (direct descendants) \
    of $N$ and $| c(N) |$ is the size of $c(N)$. The total relative entropy is $P(X) - Q(X)$ and the normalised relative entropy (NRE) is $ ( P(X) - Q(X) ) / Z $. \
    The measure therefore reflects the balance of the tree.\
    Children: The children of a net are its direct descendants nets in the subsequent net layer of the (multi layered) cactus tree. Results given for non-leaf nets only. \
    Depth: The depth of a net is the number of nodes (excluding itself) on the path from it to the root node. Results for leaf nets only. \
    (A leaf net is net with only chain descendants in the multi-layered cactus tree)")

def writeBlocksTable(stats, fileHandle):
    columnNumber = 14
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Blocks", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                 ("Total", 1, 1, 0, 1), 
                                 ("Per Net", 2, 4, 0, 0),
                                 ("Max", 2, 2, 1, 1),
                                 ("Avg.", 3, 3, 1, 1),
                                 ("Med.", 4, 4, 1, 1),
                                 ("Length", 5, 7, 0, 0),
                                 ("Max", 5, 5, 1, 1),
                                 ("Avg.", 6, 6, 1, 1),
                                 ("Med.", 7, 7, 1, 1),
                                 ("Degree", 8, 10, 0, 0),
                                 ("Max", 8, 8, 1, 1),
                                 ("Avg.", 9, 9, 1, 1),
                                 ("Med.", 10, 10, 1, 1),
                                 ("Coverage", 11, 13, 0, 0),
                                 ("Max", 11, 11, 1, 1),
                                 ("Avg.", 12, 12, 1, 1),
                                 ("Med.", 13, 13, 1, 1)), fileHandle)
    for statNode, regionName in stats:
        blocksNode = statNode.find("blocks")
        writeRow((regionName, formatFloat(blocksNode.attrib["total_number"], decimals=0),
               formatFloat(blocksNode.attrib["max_number_per_net"], decimals=0), 
               formatFloat(blocksNode.attrib["average_number_per_net"], decimals=2),
               formatFloat(blocksNode.attrib["median_number_per_net"], decimals=0),
               formatFloat(blocksNode.attrib["max_length"], decimals=0),
               formatFloat(blocksNode.attrib["average_length"], decimals=2),
               formatFloat(blocksNode.attrib["median_length"], decimals=0),
               formatFloat(blocksNode.attrib["max_degree"], decimals=0),
               formatFloat(blocksNode.attrib["average_degree"], decimals=2),
               formatFloat(blocksNode.attrib["median_degree"], decimals=0),
               formatFloat(blocksNode.attrib["max_coverage"], decimals=0),
               formatFloat(blocksNode.attrib["average_coverage"], decimals=2),
               formatFloat(blocksNode.attrib["median_coverage"], decimals=0)), fileHandle)
    writeEnd(fileHandle, "blocks_table", "Statistics on the blocks of the cactus trees. \
    Region: region name. \
    Total: total number of blocks in the cactus tree. \
    Per Net: numbers of blocks in the child chains of each net. \
    Length: number of basepairs in a block. \
    Degree: number of sequences in a block. \
    Coverage: block's length * degree.")

def writeChainsTable(stats, fileHandle):
    columnNumber = 15
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Chains", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Type", 1, 1, 0, 1), 
                                
                                ("Total", 2, 2, 0, 1),
                                 
                                 ("Per Net", 3, 5, 0, 0), 
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
        chainsNode = statNode.find("chains")
        chains2BNode = statNode.find("chains_with_two_or_more_blocks")
        
        writeLine(columnNumber, 2, ((regionName, 0, 0, 0, 1),
                  ("all", 1, 1, 0, 0),
                  (formatFloat(chainsNode.attrib["total_number"], decimals=0), 2, 2, 0, 0),
                  
                  (formatFloat(chainsNode.attrib["max_number_per_net"], decimals=0), 3, 3, 0, 0), 
                  (formatFloat(chainsNode.attrib["average_number_per_net"], decimals=2), 4, 4, 0, 0),
                  (formatFloat(chainsNode.attrib["median_number_per_net"], decimals=0), 5, 5, 0, 0),
                  
                  (formatFloat(chainsNode.attrib["max_length"], decimals=0), 6, 6, 0, 0),
                  (formatFloat(chainsNode.attrib["average_length"], decimals=2), 7, 7, 0, 0),
                  (formatFloat(chainsNode.attrib["median_length"], decimals=0), 8, 8, 0, 0),
                  
                  (formatFloat(chainsNode.attrib["max_base_length"], decimals=0), 9, 9, 0, 0),
                  (formatFloat(chainsNode.attrib["average_base_length"], decimals=2), 10, 10, 0, 0),
                  (formatFloat(chainsNode.attrib["median_base_length"], decimals=0), 11, 11, 0, 0),
                  
                  (formatFloat(chainsNode.attrib["max_instance_length"], decimals=0), 12, 12, 0, 0),
                  (formatFloat(chainsNode.attrib["average_instance_length"], decimals=2), 13, 13, 0, 0),
                  (formatFloat(chainsNode.attrib["median_instance_length"], decimals=0), 14, 14, 0, 0),
                  
                  ("$>=2$ B.", 1, 1, 1, 1),
                  (formatFloat(chains2BNode.attrib["total_number"], decimals=0), 2, 2, 1, 1),
                  
                  (formatFloat(chains2BNode.attrib["max_number_per_net"], decimals=0), 3, 3, 1, 1), 
                  (formatFloat(chains2BNode.attrib["average_number_per_net"], decimals=2), 4, 4, 1, 1),
                  (formatFloat(chains2BNode.attrib["median_number_per_net"], decimals=0), 5, 5, 1, 1),
                  
                  (formatFloat(chains2BNode.attrib["max_length"], decimals=0), 6, 6, 1, 1),
                  (formatFloat(chains2BNode.attrib["average_length"], decimals=2), 7, 7, 1, 1),
                  (formatFloat(chains2BNode.attrib["median_length"], decimals=0), 8, 8, 1, 1),
                  
                  (formatFloat(chains2BNode.attrib["max_base_length"], decimals=0), 9, 9, 1, 1),
                  (formatFloat(chains2BNode.attrib["average_base_length"], decimals=2), 10, 10, 1, 1),
                  (formatFloat(chains2BNode.attrib["median_base_length"], decimals=0), 11, 11, 1, 1),
                  
                  (formatFloat(chains2BNode.attrib["max_instance_length"], decimals=0), 12, 12, 1, 1),
                  (formatFloat(chains2BNode.attrib["average_instance_length"], decimals=2), 13, 13, 1, 1),
                  (formatFloat(chains2BNode.attrib["median_instance_length"], decimals=0), 14, 14, 1, 1)), 
                  fileHandle)
    
    writeEnd(fileHandle, "chains_table", "Statistics on the chains of the cactus trees. \
    Region: region name. \
    Type: categories of chains, either `all', which includes all chains or `$>=2$ B.', \
    which includes only chains containing a minimum of two blocks. \
    Total: total number of chains in the cactus tree. \
    Per Net: numbers of child chains in each net. \
    Link Number: number of links in chain. \
    Block Bp length: number of basepairs in blocks of chain. \
    Instance length: average number of basepairs in an instance of the chain, including both its blocks and intervening links.")

def writeEndsTable(stats, fileHandle):
    columnNumber = 8
    writePreliminaries(columnNumber, fileHandle)
    writeLine(columnNumber, 1, (("Ends", 0, columnNumber-1, 0, 0),), fileHandle)
    writeLine(columnNumber, 2, (("Region", 0, 0, 0, 1), 
                                ("Type", 1, 1, 0, 1), 
                                 ("Per Group", 2, 4, 0, 0), 
                                 ("Max", 2, 2, 1, 1),
                                 ("Avg.", 3, 3, 1, 1),
                                 ("Med.", 4, 4, 1, 1),
                                 ("Connectivity", 5, 7, 0, 0), 
                                 ("Max", 5, 5, 1, 1),
                                 ("Avg.", 6, 6, 1, 1),
                                 ("Med.", 7, 7, 1, 1)), fileHandle)
    for statNode, regionName in stats:
        endsNode = statNode.find("ends")
        endsNonTrivialNode = statNode.find("ends_non_trivial")
        endsTrivialNode = statNode.find("ends_trivial")
        
        writeLine(columnNumber, 3, ((regionName, 0, 0, 0, 2),
                  ("all", 1, 1, 0, 0),
                  (formatFloat(endsNode.attrib["max_number_per_group"], decimals=0), 2, 2, 0, 0),
                  (formatFloat(endsNode.attrib["average_number_per_group"], decimals=2), 3, 3, 0, 0), 
                  (formatFloat(endsNode.attrib["median_number_per_group"], decimals=0), 4, 4, 0, 0),
                  (formatFloat(endsNode.attrib["max_degree"], decimals=0), 5, 5, 0, 0),
                  (formatFloat(endsNode.attrib["average_degree"], decimals=2), 6, 6, 0, 0),
                  (formatFloat(endsNode.attrib["median_degree"], decimals=0), 7, 7, 0, 0),
                  
                  ("links", 1, 1, 1, 1),
                  (formatFloat(endsTrivialNode.attrib["max_number_per_group"], decimals=0), 2, 2, 1, 1),
                  (formatFloat(endsTrivialNode.attrib["average_number_per_group"], decimals=2), 3, 3, 1, 1), 
                  (formatFloat(endsTrivialNode.attrib["median_number_per_group"], decimals=0), 4, 4, 1, 1),
                  (formatFloat(endsTrivialNode.attrib["max_degree"], decimals=0), 5, 5, 1, 1),
                  (formatFloat(endsTrivialNode.attrib["average_degree"], decimals=2), 6, 6, 1, 1),
                  (formatFloat(endsTrivialNode.attrib["median_degree"], decimals=0), 7, 7, 1, 1),
                  
                  ("tangles", 1, 1, 2, 2),
                  (formatFloat(endsNonTrivialNode.attrib["max_number_per_group"], decimals=0), 2, 2, 2, 2),
                  (formatFloat(endsNonTrivialNode.attrib["average_number_per_group"], decimals=2), 3, 3, 2, 2), 
                  (formatFloat(endsNonTrivialNode.attrib["median_number_per_group"], decimals=0), 4, 4, 2, 2),
                  (formatFloat(endsNonTrivialNode.attrib["max_degree"], decimals=0), 5, 5, 2, 2),
                  (formatFloat(endsNonTrivialNode.attrib["average_degree"], decimals=2), 6, 6, 2, 2),
                  (formatFloat(endsNonTrivialNode.attrib["median_degree"], decimals=0), 7, 7, 2, 2)), fileHandle)
                     
    writeEnd(fileHandle, "ends_table", "Statistics on the ends of the cactus trees. \
    Region: region name. \
    Type: categories of end, either `all', which includes all ends, `links', which includes all ends within a link group or \
    `tangles', which includes ends found in tangle groups. \
    Per Group: numbers of ends in a group. \
    Connectivity: the number of distinct ends an end is adjacent to.")

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
    writeNetTable(stats, fileHandle)
    writeBlocksTable(stats, fileHandle)
    writeChainsTable(stats, fileHandle)
    writeEndsTable(stats, fileHandle)
    writeDocumentEnd(fileHandle)
    
    ##########################################
    #Cleanup
    ##########################################
    
    fileHandle.close()
    
if __name__ == "__main__":
    main()
