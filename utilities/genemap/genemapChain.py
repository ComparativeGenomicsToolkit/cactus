#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""
nknguyen@soe.ucsc.edu
Jan 13 2011: Edit so that the script can combine genemap results from separate cactus runs
Dec 13 2010: 
Revised: Sep 14 2010
June 01 2010
Score how well genes mapped to cactus structure.
Input: xml output file from cactus_geneMap
Output: Latex summary table
"""

import os, re, sys, subprocess
import xml.etree.ElementTree as ET
from optparse import OptionParser
from sonLib.bioio import system

def usage():
    sys.stderr.write("Usage: geneMap.py <config.xml>\n")
    sys.exit(2)

#returns aggregated elements of multiple trees:
def getAggElements(trees, elemTag):
    if len(trees) == 0:
        return None
    root = trees[0].getroot()
    list = root.getiterator(elemTag)
    for i in range(1, len(trees)):
        root = trees[i].getroot()
        list.extend(root.getiterator(elemTag))
    return list

def getOverlap(s1, e1, s2, e2):
    s = s1
    e = e1
    if s < s2:
        s = s2
    if e > e2:
        e = e2
    return e - s

def countChains(root):
    """Return number of different chains in the tree starting from input root
    """
    chains = root.getiterator("chainName")
    visited = []
    count = 0
    for chain in chains:
        if chain.text not in visited:
            count += 1
	visited.append(chain.text)
    return count

def countChains2(root):
    """Return a dictionary of level to different chains at that level
        in the tree starting from input root
    """
    blocks = root.getiterator("block")
    levelToChains = {}
    for block in blocks:
        blockName = block.find('blockName').text
        chain = block.find('chainName').text
        if chain == "NA":
            chain += blockName
        level = block.find('level').text
        if level not in levelToChains:
            levelToChains[level] = [chain]
        else:
            if chain not in levelToChains[level]:
                 levelToChains[level].append(chain)
    return levelToChains

def getBreakLenDist(exon):
    #Return a dictionary of breakLen to number of breaks with those lengths.
    lenToBreakCounts = {}
    levelToChains = countChains2(exon)
    levels = sorted(levelToChains.keys())
    if len(levels) == 0:
        sys.stderr.write("Exon %s, start: %s, end: %s\n" %(exon.get("id"), exon.get("start"), exon.get("end")))
        return -1

    #Calculate total number of bases at lower levels
    start = int(exon.get("start"))
    end = int(exon.get("end"))
    numLowerLevelBases = 0
    topLevel = levels[0]
    level = topLevel
    blocks = exon.findall("block")
    for block in blocks:
        blockStart = int(block.find("start").text)
        blockEnd = int(block.find("end").text)
        level = block.find("level").text
        if level == topLevel:
            if numLowerLevelBases > 0:
                if numLowerLevelBases in lenToBreakCounts:
                    lenToBreakCounts[numLowerLevelBases] +=1
                else:
                    lenToBreakCounts[numLowerLevelBases] = 1
                numLowerLevelBases = 0
        else:
            numLowerLevelBases += getOverlap(start, end, blockStart, blockEnd)
    
    if numLowerLevelBases > 0:
        if numLowerLevelBases in lenToBreakCounts:
            lenToBreakCounts[numLowerLevelBases] +=1
        else:
            lenToBreakCounts[numLowerLevelBases] = 1
    
    return lenToBreakCounts

def getNonTripleBreaks(exon):
    breakLenDist = getBreakLenDist(exon)
    numNtb = 0
    for len in breakLenDist:
        if len%3 != 0:
            numNtb += breakLenDist[len]
    return numNtb

def mult3totalIndels(exon):
    breakLenDist = getBreakLenDist(exon)
    sumBreakLen = 0
    for len in breakLenDist:
        sumBreakLen += len
    
    if sumBreakLen%3 == 0:
        return "true"
    else:
        return "false"

def writeDocumentStart(f):
    f.write("\\documentclass[11pt]{article}\n") 
    f.write("\\usepackage{epsfig}\n")
    f.write("\\usepackage{multirow}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{array}\n")
    f.write("\\usepackage{slashbox}\n")
    f.write("\n")

    f.write("\\newcommand{\\figref}[1]{Figure~\\ref{fig:#1}}\n")
    f.write("\\newcommand{\\tabref}[1]{Table~\\ref{tab:#1}}\n")
    f.write("\n")

    f.write("\\textwidth=6.5in\n")
    f.write("\\textheight=9in\n")
    f.write("\\oddsidemargin=0in\n")
    f.write("\\evensidemargin=0in\n")
    f.write("\\topmargin=0in\n")
    f.write("\\topskip=0in\n")
    f.write("\\headheight=0in\n")
    f.write("\\headsep=0in\n")
    f.write("\n")
    
    f.write("\\begin{document}\n")
    return

def writeDocumentEnd(f):
    f.write("\\end{document}\n")


#============ GET GENE STRUCTURE TABLE =========
def score(gene, f, max):
    l2c = countChains2(gene)
    allLevels = sorted(l2c.keys())
    name = gene.get("name")

    level2Str = {}
    #initializing the hash
    if len(allLevels) <= 0:
        sys.stderr.write("Numlevel = 0!\n")
        return
    level2Str[allLevels[0]] = "\\multirow{%d}{*}{%s} & %s &" %(len(allLevels) + 1, name, allLevels[0])
    for i in range(1, len(allLevels)):
        level2Str[allLevels[i]] = " & %s & " %(allLevels[i])

    nonTripleBreaks = " & NTB & "
    exons = gene.findall("exon")
    for exon in exons:
        levelToChains = countChains2(exon)
        #levels = sorted(levelToChains.keys())
        for level in allLevels:
            if level in levelToChains:
                level2Str[level] += str(len(levelToChains[level]))
            else:
                level2Str[level] += " - "
            level2Str[level] += " & "
        nonTB = getNonTripleBreaks(exon)
        if nonTB == 0:
            nonTripleBreaks += " - & "
        else:
            nonTripleBreaks += str(nonTB) + " & "
    
    for i in range(len(exons), max):
        for level in allLevels:
            level2Str[level] += "  & "
        nonTripleBreaks += "  & "

    for level in allLevels:
        f.write("%s %d\\\\\n" %(level2Str[level], len(l2c[level])))

    f.write("%s\\\\\n" %(nonTripleBreaks))
    f.write("\\hline\n")
    return     

def getkey(elem):
     return elem.get("name")


def getCombinedGeneStructureTab(trees, f, max):
    for tree in trees:
        getGeneStructureTab(tree, f, max)

def getGeneStructureTab(tree, f, max):
    root = tree.getroot()

    genes = root.findall("gene")
    genes[:] = sorted(genes, key=getkey)
    #get score for each gene
    for gene in genes:
        score(gene, f, max)

def getMaxExonNum(tree):
    max = 0
    root = tree.getroot()
    genes = root.findall("gene")
    for gene in genes:
        exonCount = int(gene.get("exonCount"))
        if max < exonCount:
            max = exonCount
    return max

def geneStrucTableHeader(exonNum, f):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    exonC = ""
    exonStr = ""
    for i in range(exonNum):
        exonC += ">{\small}c"
        exonStr += "exon" + str(i) + " & "
    f.write("\\scalebox{0.7}{%\n")
    f.write("\\begin{tabular}{|>{\small}c|>{\small}c|%s|>{\small}c|}\n" % (exonC))
    f.write("\\hline\n")
    f.write("Gene & CLevel & %s TNumDiffChains\\\\\n" %(exonStr))
    f.write("\\hline\n")

def tableCloser(f, captionStr):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\end{center}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{table}\n")

def printGeneStructure(outdir, file, trees):
    file = os.path.join(outdir, file)
    f = open(file, "w")

    maxExonNum = 0
    for tree in trees:
        currMax = getMaxExonNum(tree)
        if maxExonNum < currMax:
            maxExonNum = currMax

    writeDocumentStart(f)

    geneStrucTableHeader(maxExonNum, f)
    getCombinedGeneStructureTab(trees, f, maxExonNum)
    captionStr = "Mapping cactus chain structure to RefSeq genes. CLevel is the chain level. The number corresponding to a CLevel i and an exon j is the number of different chains at level i within exon j. NTB is the number of Non-Triplet Breaks within an exon. TNumDiffChains of CLevel i is the number of different chains with level i across the whole gene."
    tableCloser(f, captionStr)
    
    writeDocumentEnd(f)
    f.close()
#========== END GENE STRUCTURE ========

#============= HISTOGRAM =========
#x-axis: number of chain-level a gene has
#y-axis: number of genes that have x chain-levels
def writeData(outdir, name, header, list):
    dataFile = "data-" + name + ".txt"
    file = os.path.join(outdir, dataFile)
    f = open(file, "w")
    f.write(header + "\n")
    for item in list:
        f.write("%d\n" %(item))
    f.close()

def writeRscript(outdir, name, titles, xlabels, ylabels, numPlots, numHistsPerPlot, outputFormat='pdf'):
    #credit to the script 'cactus_treeStatsPlotter.py' from dearl@soe.ucsc.edu 
    #numHistsPerPlot: Number of HISTOGRAMS (ON TOP OF EACH OTHER) per plot
    #numPlots: number of plots printed out to the pdf (or png)

    #if isLog:
    #    ylabel = "log of " + ylabel
    if numPlots <= 0 or numHistsPerPlot <= 0:
        return
    if outputFormat == 'pdf':
        imgName = "hist-%s.pdf" %(name)
        imgFile = os.path.join(outdir, imgName)
        startImgDeviceStr = "pdf('%s', height=4, width=8)\n" %(imgFile)
    elif outputFormat == 'png':
        imgName = "hist-%s.png" %(name)
        imgFile = os.path.join(outdir, imgName)
        startImgDeviceStr = "png('%s', height=400, width=800)\n" %(imgFile)
    else:
        sys.stderr.write("Unknown image output format: %s\n" %(outputFormat))
        usage()

    dataFiles = []
    for i in range(0, numPlots):
        for j in range(0, numHistsPerPlot):
            dataFiles.append( os.path.join(outdir, "data-" + name + str(i) + str(j) + ".txt") )

    #Start writing script:
    scriptFile = os.path.join(outdir, name + ".R")
    f = open(scriptFile, "w")
    f.write("#This script is automatically generated by geneMap.py\n")
    for i in range(0, numPlots):
        for j in range(0, numHistsPerPlot):
            f.write("table%d%d = read.table('%s', header=TRUE, sep='\\t')\n"\
                     %(i, j, dataFiles[i*numHistsPerPlot + j]))
            f.write("x%d%d=table%d%d[,1]\n" %(i,j, i,j))
            f.write("label%d%d = names(table%d%d)[1]\n" %(i,j, i,j))
    #writing the legend labels
    for i in range(0, numPlots):
        f.write("labels%d = c(label%d0" %(i, i))
        for j in range(1, numHistsPerPlot):
            f.write(", label%d%d" %(i,j))
        f.write(")\n")


    cols = ["blue", "green", "red", "black", "pink", "purple", "orange"]
    if len(cols) < numHistsPerPlot:
        sys.stderr.write("Function 'writeRscript': Not enough color to draw histogram!\n")
        sys.exit(1)
    mycols = cols[0: numHistsPerPlot]
    f.write("mycols = c('%s'" %(mycols[0]))
    if numHistsPerPlot > 1:
        for j in range(1, numHistsPerPlot):
            f.write(", '%s'" %(mycols[j]))
    f.write(")\n")

    f.write(startImgDeviceStr)
    if numPlots == 1:
        numCols = 1
        numRows = 1
    else:
        numCols = 2
        numRows = numPlots/numCols
        if (numPlots % numCols) != 0:
            numRows +=1
    f.write("par(mfrow=c(%d, %d))\n" %(numRows, numCols))
        
    for i in range(0, numPlots):
        f.write("breaks%d = (min(x%d0)-0.5):(max(x%d0)+0.5)\n" %(i, i, i))
        f.write("if ( length(breaks%d) > 36 ) {step = max(x%d0)/20; breaks%d = step*0:20}\n" %(i, i, i))
        f.write("hist(x%d0, breaks=breaks%d, col=mycols[1], main='%s', xlab='%s',ylab='%s')\n" %(i, i, titles[i], xlabels[i], ylabels[i]))

        if numHistsPerPlot > 1:
            for j in range(1, numHistsPerPlot):
                f.write("hist(x%d%d, breaks=breaks%d, col=mycols[%d], add=TRUE)\n" %(i,j, i, j+1))
    if numHistsPerPlot > 1 and numPlots == 1:
        f.write("legend(\"topright\", labels0, col=mycols, lty=rep(1,length(mycols)),lwd=rep(5, length(mycols)))\n")
    
    f.write("dev.off()\n")
    f.close()

def runRscript(outdir, name, numPlots, numHistsPerPlot, cleanUp):
    scriptName = name + ".R"
    scriptFile = os.path.join(outdir, scriptName)
    #command = "R --no-save --no-restore < %s >/dev/null/ 2>&1" %(os.path.join(os.getcwd(), scriptName))
    command = "R --no-save --no-restore < %s > /dev/null" %(scriptFile)
    status = subprocess.Popen(command, shell=True)
    status.wait()
    if status.returncode:
        sys.stderr.write("%s:%s failed to run properly in R. Return code is %d\n" %(sys.argv[0], scriptName, status.returncode))
    
    if cleanUp:
        os.remove(scriptFile)
        for i in range(0, numPlots):
            for j in range(0, numHistsPerPlot):
                dataName = "data-" + name + str(i) + str(j) + ".txt"
                dataFile = os.path.join(outdir, dataName)
                os.remove(dataFile)

#============= END OF HISTOGRAM =================

#======= GET DISTRIBUTION OF EXONS ACROSS DIFFERENT NUMBER OF NonTripletBreaks
def getNtbDist(outdir, file, trees):
    #root = tree.getroot()
    #allExonNtb = []
    brokenExonNtb = []
    
    #exons = root.getiterator("exon")
    exons = getAggElements(trees, "exon")
    for exon in exons:
        ntb = getNonTripleBreaks(exon)
        #allExonNtb.append(ntb)
        levelToChains = countChains2(exon)
        levels = sorted(levelToChains.keys())
        if not (len(levels) == 1 and len(levelToChains[levels[0]]) == 1):
            brokenExonNtb.append(ntb)
    #writeData(file + "0", "ntb", allExonNtb)
    writeData(outdir, file + "00", "ntb", brokenExonNtb)
    titles = ["Distribution of non-triplet breaks"]
    xlabels = ["Number of non-triplet breaks"]
    ylabels = ["Number of exons"]
    numPlots = 1
    numHistsPerPlot = 1
    writeRscript(outdir, file, titles, xlabels, ylabels, numPlots,numHistsPerPlot, outputFormat='pdf')
    runRscript(outdir, file, numPlots, numHistsPerPlot, None)

#======= GET BREAKS LENGTH DISTRIBUTION ================
def getBreakLenHist(outdir, file, trees):
    #root = tree.getroot()
    breakLens = []
    #exons = root.getiterator("exon")
    exons = getAggElements(trees, "exon")
    for exon in exons:
        lenDist = getBreakLenDist(exon)
        if len(lenDist) == 0:
            continue
        for bl in lenDist:#bl = breaklength
            for i in range(0, lenDist[bl]):
                breakLens.append(bl)
            #if bl <= 36:
            #    for i in range(0, lenDist[bl]):
            #        breakLens.append(bl)

    writeData(outdir, file + "00", "bld", breakLens)
    titles = ["Distribution of break-lengths"]
    xlabels = ["Number of bases"]
    ylabels = ["Number of breaks"]
    numPlots = 1
    numHistsPerPlot = 1
    writeRscript(outdir, file, titles, xlabels, ylabels, numPlots,numHistsPerPlot, outputFormat='pdf')
    runRscript(outdir, file, numPlots, numHistsPerPlot, None)

#========== GET conservedGeneDist ==============
def getGeneInfo(gene):
    exonCount = int(gene.get("exonCount"))
    range = int(gene.get("end")) - int(gene.get("start"))
    exons = gene.findall("exon")
    exonicLen = 0
    for exon in exons:
        exonicLen += int(exon.get("end")) - int(exon.get("start"))
    return (exonCount, exonicLen, range) 

def  getConservedGeneDist(outdir, file, trees):
    """Looks at the distributions of genes that have no broken exon.\
    This is to check how much of a bias these 'conserved genes' are toward\
    whether they are with small number of exons, or with\
    small number of exonic bases, or span smaller regions of the genome."""
    #root = tree.getroot()
    #genes = root.findall("gene")
    genes = getAggElements(trees, "gene")

    numExons = [] #list of total number of exons each conserved gene has
    exonicLens = [] #list of total number of exonic bases each conserved gene has
    ranges = [] #list of (geneEnd - geneStart)
    
    totalNumExons = [] #list of total number of exons each gene (conserved and broken) has
    totalExonicLens = [] #list of total number of exonic bases each gene has
    totalRanges = [] #list of (geneEnd - geneStart)
    #Finding all the 'conserved' (no broken exon) genes and their 
    for gene in genes:
        levelToChains = countChains2(gene)
        levels = sorted(levelToChains.keys())
        (numExon, exonicLen, range) = getGeneInfo(gene)
        totalNumExons.append(numExon)
        totalExonicLens.append(exonicLen)
        totalRanges.append(range)
        if len(levels) == 1 and len(levelToChains[levels[0]]) == 1:
            numExons.append(numExon)
            exonicLens.append(exonicLen)
            ranges.append(range)
    
    writeData(outdir, file + "00", "Total exons", totalExonicLens)
    writeData(outdir, file + "01", "Conserved exons", exonicLens)
    writeData(outdir, file + "10", "Total exons", totalRanges)
    writeData(outdir, file + "11", "Conserved exons", ranges)
    writeData(outdir, file + "20", "Total exons", totalNumExons)
    writeData(outdir, file + "21", "Conserved exons", numExons)
    
    titles = ["", "", ""]
    xlabels = ["Number of exonic bases", "Number of genomic bases", "Number of exons"]
    ylabels = ["Number of genes", "Number of genes", "Number of genes"]
    numPlots = 3
    numHistsPerPlot = 2
    writeRscript(outdir, file, titles, xlabels, ylabels, numPlots, numHistsPerPlot, outputFormat='pdf')
    runRscript(outdir, file, numPlots, numHistsPerPlot, None)


#============= NUMBER OF CHAIN-LEVELS HISTOGRAM =========
def getNumChainLevelHist(outdir, file, list1, list2):
    writeData(outdir, file + "00", "Total genes", list1)
    writeData(outdir, file + "01", "Genes with one top chain", list2)
    titles = ["Distribution of the number of chain-levels each gene has"]
    xlabels = ["Number of chain-levels"]
    ylabels = ["Number of genes"]
    numPlots = 1
    numHistsPerPlot = 2
    writeRscript(outdir, file, titles, xlabels, ylabels, numPlots, numHistsPerPlot, outputFormat='pdf')
    runRscript(outdir, file, numPlots, numHistsPerPlot, None)

#============= END NUMBER OF CHAIN-LEVELS HISTOGRAM =========

#============ 
def getFreq(list):
    dict = {}
    for item in list:
        if item in dict:
            dict[item] += 1
        else:
            dict[item] = 1
    return dict

def getFreqTable(outdir, file, list1, list2):
    freqdict1 = getFreq(list1)
    freqdict2 = getFreq(list2)
   
    file = os.path.join(outdir, file)
    f = open(file, "w")
    f.write("NumberOfChainLevel\tgeneCounts\tcleanGeneCounts\n")
    for key in freqdict1:
        f.write("%d\t%d\t" %(key, freqdict1[key]))
        if key in freqdict2:
            freq2 = freqdict2[key]
        else:
            freq2 =0
        f.write("%d\n" %(freq2))
    f.close()
#=============
        
#============ GET CHAININFO ==============
def statsTableHeader(f):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{0.7}{%\n")
    #| | #CL/genes | #CL/exon | #ConservedExons/gene | #BreakExons/gene | #BreakExonsWithNTB/gene
    f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c||>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\hline\n")
    f.write(" & NumCL/gene & NumCL/exon & NumConservedExons/gene &NumBrokenExons/gene & NumBrokenExonsWithNTB/gene\\\\\n")
    f.write("\\hline\n")

def getExonStats(genes):
    conservedEPG = [] #list of number of conserved (1 chain) exon per gene (EPG)
    brokenEPG = []
    ntbEPG = []
    for gene in genes:
        exons = gene.findall("exon")
        conservedExon = 0
        brokenExon = 0
        ntbExon = 0
        for exon in exons:
            levelToChains = countChains2(exon)
            levels = sorted(levelToChains.keys())
            if len(levelToChains) == 1 and len(levelToChains[levels[0]]) == 1:
                conservedExon += 1
            else:
                brokenExon += 1
                if getNonTripleBreaks(exon) > 0:
                    ntbExon +=1
        conservedEPG.append(conservedExon)
        brokenEPG.append(brokenExon)
        ntbEPG.append(ntbExon)
        #conservedEPG.append(conservedExon*100.0/len(exons))
        #brokenEPG.append(brokenExon*100.0/len(exons))
        #ntbEPG.append(ntbExon*100.0/len(exons))
    return (conservedEPG, brokenEPG, ntbEPG)

def getNumChainLevelList(rootlist):
    #sys.stderr.write("Getting the list of number of chain-levels per gene/exon\n")
    #sys.stderr.write("Number of items (genes or exons): %d\n" %len(rootlist))
    numCLs = [] #list of number of chain-level per gene or per exon
    for root in rootlist:#For each gene (or each exon)
        levelToChains = countChains2(root)
        numCLs.append(len(levelToChains))
    return numCLs

def getNumChainLevelList2(rootlist):
    #return a list of number of chain-Level for each gene/exon
    #that has the lowest-chain-level that only has 1 chain
    numCLs = [] #list of number of chain-level per gene or per exon
    for root in rootlist:#For each gene (or each exon)
        levelToChains = countChains2(root)
        levels = sorted(levelToChains.keys())
        if len(levels) <= 0:
            sys.stderr.write("Gene or Exon without any chain-level!\n")
            sys.exit(2)
        lowestLevel = levels[0]
        lowestChains = levelToChains[lowestLevel]
        if len(lowestChains) == 1:    
            numCLs.append(len(levels))
    return numCLs

def getMaxMinAvg(list):
    if len(list) <= 0:
        sys.stderr.write("Gene or Exon without any chain!\n")
        sys.exit(2)

    avg = sum(list)*1.0/len(list)
    return (max(list), min(list), avg)

def printStats(outdir, file, trees):
    file2 = os.path.join(outdir, file + ".tex")
    f = open(file2, "w")
    writeDocumentStart(f)
    statsTableHeader(f)
    
    #root = tree.getroot()

    #Number of chain-levels/gene
    #genes = root.findall("gene")
    genes = getAggElements(trees, "gene")
    geneNumCLs = getNumChainLevelList(genes)
    (maxGeneCL, minGeneCL, avgGeneCL) = getMaxMinAvg(geneNumCLs)

    #Number of chain-levels/exon
    #exons = root.getiterator("exon")
    exons = getAggElements(trees, "exon")
    eNumCLs = getNumChainLevelList(exons)
    (maxExonCL, minExonCL, avgExonCL) = getMaxMinAvg(eNumCLs)
    
    #Number of conserved, broken, nonTripletBreak exons per gene
    (conservedEPG, brokenEPG, ntbEPG) = getExonStats(genes)
    (maxConservedEPG, minConservedEPG, avgConservedEPG) = getMaxMinAvg(conservedEPG)
    (maxBrokenEPG, minBrokenEPG, avgBrokenEPG) = getMaxMinAvg(brokenEPG)
    (maxNtbEPG, minNtbEPG, avgNtbEPG) = getMaxMinAvg(ntbEPG)
   
    #Print out the latex rows
    f.write("Max & %d & %d & %d & %d & %d\\\\\n" %(maxGeneCL, maxExonCL,\
             maxConservedEPG, maxBrokenEPG, maxNtbEPG)) 
    f.write("\\hline\n")
    f.write("Min & %d & %d & %d & %d & %d\\\\\n" %(minGeneCL, minExonCL,\
             minConservedEPG, minBrokenEPG, minNtbEPG)) 
    f.write("\\hline\n")
    f.write("Average & %.2f & %.2f & %.2f & %.2f & %.2f\\\\\n" %(avgGeneCL, avgExonCL,\
             avgConservedEPG, avgBrokenEPG, avgNtbEPG)) 
    f.write("\\hline\n")
     
    captionStr = ""
    tableCloser(f, captionStr)
    writeDocumentEnd(f)
    f.close()

    #Get Number-Chain-level distribution
    cleanGeneNumCLs = getNumChainLevelList2(genes)
    getNumChainLevelHist(outdir, file, geneNumCLs, cleanGeneNumCLs)
    freqFile = "clFreq-" + file
    getFreqTable(outdir, freqFile, geneNumCLs, cleanGeneNumCLs)

#================== PRINT CATEGORIES of GENES =============
def categoriesTableHeader(f):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\hline\n")
    f.write("\\multicolumn{2}{|c|}{\\multirow{2}{*}{\\backslashbox{NumLevels}{NumTopChains}}} & \\multirow{2}{*}{1} & \\multicolumn{2}{|c|}{$>=2$}\\\\\n")
    f.write("\\cline{4-5}\n")
    f.write("\\multicolumn{2}{|c|}{} &  & No MulTCExon & With MulTCExon\\\\\n")
    f.write("\\hline\n")


def hasMulTopChainsExon(gene):
    exons = gene.findall("exon")
    for exon in exons:
        levelToChains = countChains2(exon)
        levels = sorted(levelToChains.keys())
        if len(levels) <= 0:
            sys.stderr.write("Exon with no chain!\n")
            sys.exit(1)
        topChains = levelToChains[levels[0]]
        if len(topChains) > 1:
            return "true"
    return "false"


def hasNTB(gene):
    exons = gene.findall("exon")
    for exon in exons:
        ntb = getNonTripleBreaks(exon)
        if ntb > 0:
            return "true"
    return "false"

def hasNTB2(gene):
    exons = gene.findall("exon")
    for exon in exons:
        if mult3totalIndels(exon) == "false":
            return "true"
    return "false"

def allExonsHaveSameTopChain(gene):
    #return true if all exons have the same top-level chain. Otherwise return false
    exons = gene.findall("exon")
    toplevel = -1

    for exon in exons:
        levelToChains = countChains2(exon)
        levels = sorted(levelToChains.keys())
        if toplevel == -1:
            toplevel = levels[0]
        elif toplevel != levels[0]:
            return "false"
    return "true"

def getCategory(gene):
    #Catagory is represented as the binary ABCDE, based on 4 aspects (each X
    #corresponds to one aspect: 
    #1/ Number of chain-level: A = 0 if has 1 level, otherwise A = 1
    #2/ Number of top-level chains. B = 0 if has 1 top-chain, otherwise B = 1
    #3/ If the gene has NO exon that has more than one top-chains, C = 0, else C = 1
    #4/ If the gene has NO non-triplet break, D = 0, else D = 1
    #5/ If some exon(s) stay in lower level compared with other exons, E = 1, else E = 0

    levelToChains = countChains2(gene)
    levels = sorted(levelToChains.keys())
    if len(levels) <= 0:
        sys.stderr.write("Gene with no chain!\n")
        sys.exit(1)

    topChains = levelToChains[levels[0]]
    hasNtb = hasNTB2(gene)
    #hasNtb = hasNTB(gene)
    hasMulTCexon = hasMulTopChainsExon(gene)
    sameExonTopLevel = allExonsHaveSameTopChain(gene)
    
    if len(levels) == 1:
        A = "0"
    else: 
        A = "1"

    if len(topChains) == 1:
        B = "0"
    else:
        B = "1"
   
    if hasMulTCexon == "true":
        C = "1"
    else:
        C = "0"
    
    if hasNtb == "true":
        D = "1"
    else:
        D = "0"

    if sameExonTopLevel == "true":
        E = "0"
    else:
        E = "1"

    
    cat = [A , B , C ,D ,E]
    return cat

def printCategories(outdir, file, trees):
    file = os.path.join(outdir, file + ".tex")
    f = open(file, "w")
    writeDocumentStart(f)
    categoriesTableHeader(f)
    
    #root = tree.getroot()
    #genes = root.findall("gene")
    genes = getAggElements(trees, "gene")

    #key=category, val = count
    counts = {"00000":0, "10000":0, "10010":0,\
              "100*1":0, "*10**":0, "*11**":0}

    #Get the category distribution
    sys.stdout.write("Category\tGene\n")
    for gene in genes:
        cat = getCategory(gene)
        catStr = "".join(cat)
        if catStr == "00000":
            counts[catStr] +=1
        elif catStr == "10000":
            counts[catStr] +=1
        elif catStr == "10010":
            counts[catStr] +=1
        elif cat[0] == 1 and cat[1] == 0 and cat[2] == 0 and cat[4] ==1:
            counts["100*1"] +=1
        elif cat[1] == 1 and cat[2] == 0:
            counts["*10**"] +=1
        elif cat[1] == 1 and cat[2] == 1:
            counts["*11**"] +=1
        sys.stdout.write("%s\t%s\n" %(gene.get("name"), catStr))
        #counts[catg] += 1 
    
    #sys.stdout.write("\nFREQ\n")
    #for key in counts:
    #    sys.stdout.write("%s\t%d\n" %(key, counts[key]))
    
    #Print the latex table
    #f.write("\\multicolumn{2}{|c|}{1} & %d & %d & %d\\\\\n" \
    #         %(counts["0000"],counts["0100"], counts["0110"]))
    #f.write("\\hline\n")
    #f.write("\\multirow{2}{*}{$>=2$} & no NTB & %d & %d & %d\\\\\n" \
    #         %(counts["1000"],counts["1100"], counts["1110"]))
    #f.write("\\cline{2-5}\n")
    #f.write("& with NTB & %d & %d & %d\\\\\n" \
    #         %(counts["1001"],counts["1101"], counts["1111"]))
    #f.write("\\hline\n")

    #captionStr = ""
    #tableCloser(f, captionStr)
    writeDocumentEnd(f)
    f.close()


def main(): 
    #if len(sys.argv) < 3:
    #    sys.stderr.write("Too few input arguments.\n")
    #    usage()
    #
    #usg = "usage: %prog <list of inputing xmls> [options] > geneCategories"
    usg = "usage: %prog [options] > geneCategories"
    parser = OptionParser(usage=usg)
    parser.add_option("-a", "--geneStructure", dest="geneStructure",\
           help="If specified, outputs to specified file a general table\
           which describes how each exon of each gene maps to cactus\n")
    parser.add_option("-b", "--stats", dest="stats", help="Print to specified\
                      file a reduced summary (latex table)\n")
    parser.add_option("-c", "--categories", dest="categories", help="print to the\
           specified file latex table that provides count of genes in different\
           categories\n")
    parser.add_option("-d", "--ntb", dest="ntb", help="Generates to specified file\
                      non-triplet breaks histogram\n")
    parser.add_option("-e", "--conservedGeneDist", dest="conservedGeneDist", \
                       help="Generates histograms of distribution of genes that perfectly mapped\
                       across number of exons, exonic bases, gene range, to specified file\n")
    parser.add_option("-f", "--breakLenDist", dest="breakLenDist", help="Generates to specified file\
                      distribution of breaks' length\n")
    parser.add_option("-o", "--outputDir", dest="outdir", help="Output directory\n")
    parser.add_option("-i", "--inputFiles", dest="inFiles", help="(Required argument) List of inputing genemap.xml files. E.g \"run1.xml run2.xml\"\n")
    (options, args) = parser.parse_args()

    #Parsing the input XMLs
    #tree = ET.parse(sys.argv[1])
    if options.inFiles == None:
        sys.stderr.write("Option --inputFiles is required. Please list the list of inputing XMLs\n")
        sys.exit(2)

    infiles = options.inFiles.split()
    trees = []
    for file in infiles:
        trees.append(ET.parse(file))

    if options.outdir:
        outdir = options.outdir
        system("rm -rf %s" % outdir)
        os.mkdir(outdir)
    else:
        outdir = os.getcwd()
    if options.geneStructure:
        printGeneStructure(outdir, options.geneStructure, trees)
    if options.stats:
        printStats(outdir, options.stats, trees)
    if options.categories:
        printCategories(outdir, options.categories, trees)
    if options.ntb:
        getNtbDist(outdir, options.ntb, trees)
    if options.conservedGeneDist:
        getConservedGeneDist(outdir, options.conservedGeneDist, trees)
    if options.breakLenDist:
        getBreakLenHist(outdir, options.breakLenDist, trees)

if __name__ == "__main__":
    main()
