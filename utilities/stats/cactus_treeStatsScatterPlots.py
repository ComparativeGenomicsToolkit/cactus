#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/python

"""Script for generating scatter plots from a cactus tree stats file.
"""

import xml.etree.ElementTree as ET
import os

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import system

def plot(xDists, yDists, seriesNames, outputFile, xLabel, yLabel, title):
    tempFile = getTempFile()
    assert len(xDists) == len(yDists)
    assert len(xDists) == len(seriesNames)
    tempFiles = []
    plotCommand = "plot"
    for i, j, series in zip(xDists, yDists, seriesNames):
        tempFile = getTempFile()
        fileHandle = open(tempFile, 'w')
        assert len(i) == len(j)
        for k, l in zip(i, j):
            fileHandle.write("%f\t%f\n" % (k, l))
        fileHandle.close()
        tempFiles.append(tempFile)
        plotCommand += ' "%s" using 1:2 title "%s",' % (tempFile, series)
    if plotCommand[-1] == ',':
        plotCommand = plotCommand[:-1]
    system("""echo '\
    set logscale x; set logscale y; set xlabel "%s"; \
    set ylabel "%s"; %s' | gnuplot""" % 
    (xLabel, yLabel, plotCommand))
    os.remove(tempFile)

def chainScatterPlots(stats):
    linkLengths = []
    baseLengths = []
    instanceLengths = []
    regionNames = []
    for statNode, regionName in stats:
        chainsNode = statNode.find("chains")
        regionNames.append(regionName)
        linkLengths.append([ int(i) for i in chainsNode.find("link_numbers").text.split() ])
        baseLengths.append([ int(i) for i in chainsNode.find("base_block_lengths").text.split() ])
        instanceLengths.append([ int(i) for i in chainsNode.find("avg_instance_base_length").text.split() ])
    
    plot(linkLengths, baseLengths, regionNames, "chains_linkLengths_blockLengths.ps", "Links", "Blocks Combined Basepair Length", "Chains")
    plot(linkLengths, instanceLengths, regionNames, "chains_linkLengths_instanceLengths.ps", "Links", "Avg. Basepair Instance Length", "Chains")
    plot(baseLengths, instanceLengths, regionNames, "chains_blockLengths_instanceLengths.ps", "Blocks Combined Basepair Length", "Avg. Basepair Instance Length", "Chains")

def blockScatterPlots(stats):
    blockDegrees = []
    blockLengths = []
    blockCoverage = []
    regionNames = []
    for statNode, regionName in stats:
        blocksNode = statNode.find("blocks")
        regionNames.append(regionName)
        blockDegrees.append([ int(i) for i in blocksNode.find("leaf_degrees").text.split() ])
        blockLengths.append([ int(i) for i in blocksNode.find("lengths").text.split() ])
        blockCoverage.append([ int(i) for i in blocksNode.find("leaf_coverage").text.split() ])

    plot(blockLengths, blockDegrees, regionNames, "blocks_lengths_degrees.ps", "Lengths", "Degrees", "Blocks")
    plot(blockLengths, blockCoverage, regionNames, "blocks_lengths_coverage.ps", "Lengths", "Coverage", "Blocks")
    plot(blockDegrees, blockCoverage, regionNames, "blocks_degrees_covergae.ps", "Degrees", "Coverage", "Blocks")

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = getBasicOptionParser("usage: %prog [options] treeStatsFiles", "%prog 0.1")

    options, args = parseBasicOptions(parser)
        
    logger.info("Parsed arguments")
    
    ##########################################
    #Get the input data etc.
    ##########################################
    
    assert len(args) % 2 == 0
    stats = [ (ET.parse(statsFile).getroot(), regionName) for statsFile, regionName in zip(args[::2], args[1::2]) ] 
    
    ##########################################
    #Make the scatter plots
    ##########################################
    
    chainScatterPlots(stats)
    blockScatterPlots(stats)
    
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
