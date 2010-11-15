#!/usr/bin/env python

"""
Revised: Sep 14 2010
June 01 2010: nknguyen@soe.ucsc.edu
Score how well genes mapped to cactus structure.
Input: xml output file from cactus_geneMap
Scoring scheme:
    geneScore = PercentCoverageOfLargestChain - Penalty
    Penalty = 

    (OLD)geneScore = [Sum(exonScore for all exons) - interExonChainCost*(numberOfDiffChains -1)]/(number of exons)
    exonScore = (coveredScore*numberOfCoveredBases - Penalty)/exonLength
    Penalty = breakCost*(numberNA + numberChains -1)
"""

import os, re, sys
import xml.etree.ElementTree as ET

def usage():
    sys.stderr.write("Usage: geneMap.py <config.xml>\n")
    sys.exit(2)

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

def getNonTripleBreaks(exon):
    levelToChains = countChains2(exon)
    #f.write("\tExon %s\n" %(exon.get("id")))
    levels = sorted(levelToChains.keys())
    if len(levels) == 0:
        sys.stderr.write("Exon %s, start: %s, end: %s\n" %(exon.get("id"), exon.get("start"), exon.get("end")))
        return -1
    #for level in levels:
    #    numChains = len(levelToChains[level])
    #    f.write("\t\t%s-chain: %d\n" %(level, numChains))

    #Calculate total number of bases at lower levels
    start = int(exon.get("start"))
    end = int(exon.get("end"))
    numLowerLevelBases = 0
    topLevel = levels[0]
    level = topLevel
    blocks = exon.findall("block")
    numBreaksNot3Mul = 0
    #if len(levels) >1:
    #    f.write("\t\t")
    for block in blocks:
        blockStart = int(block.find("start").text)
        blockEnd = int(block.find("end").text)
        level = block.find("level").text
        if level == topLevel:
            if numLowerLevelBases > 0:
                #f.write("%d\t" %(numLowerLevelBases)) 
                if numLowerLevelBases%3 != 0:
                    numBreaksNot3Mul += 1
                numLowerLevelBases = 0
        else:
            numLowerLevelBases += getOverlap(start, end, blockStart, blockEnd)
    if numLowerLevelBases > 0:
        #f.write("%d\t" %(numLowerLevelBases))
        if numLowerLevelBases%3 != 0:
            numBreaksNot3Mul += 1
    #if len(levels) >1:
        #f.write("\n")
    #if numBreaksNot3Mul > 0:
        #f.write("\t\tNonMul3 breaks: %d\n" %(numBreaksNot3Mul))
        
    return numBreaksNot3Mul

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

def getScores(tree, f, max):
    #where scoreList = [coveredScore, openGapCost, gapCost, intraExonChainCost, interExonChainCost]
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

def writeDocumentStart(f):
    f.write("\\documentclass[11pt]{article}\n") 
    f.write("\\usepackage{epsfig}\n")
    f.write("\\usepackage{multirow}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{array}\n")
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

def tableHeader(exonNum, f):
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

def tableCloser(f):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\end{center}\n")
    f.write("\\caption{}\n")
    f.write("\\end{table}\n")


def main(): 
    if len(sys.argv) < 3:
        sys.stderr.write("Too few input arguments.\n")
        usage()
    #
    tree = ET.parse(sys.argv[1])
    f = open(sys.argv[2], "w")
    maxExonNum = getMaxExonNum(tree)
    writeDocumentStart(f)
    tableHeader(maxExonNum, f)
    getScores(tree, f, maxExonNum)
    tableCloser(f)
    writeDocumentEnd(f)
    f.close()

if __name__ == "__main__":
    main()
