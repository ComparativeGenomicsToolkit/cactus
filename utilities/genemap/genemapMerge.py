#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser

def readHomologs(file, filter):
    f = open(file, "r")
    NAcount = 0
    missedCount = 0
    totalCount = 0
    
    gene2homolog = {}
    for line in f.readlines():
        totalCount += 1
        items = line.strip().split('\t')
        if len(items) == 2 and items[1] == 'NA':
            sys.stdout.write("%s\tNA\n" %(items[0]))
            NAcount += 1
            continue
        elif len(items) != 4 and len(items) != 5:
            totalCount -= 1
            sys.stderr.write("wrong format, line *%s*\n" %(line))
        else:
            if items[0] in gene2homolog:
                sys.stderr.write("some repeated gene: %s, please check...\n" %(items[0]))
            else:
                if (filter == 'a' and len(items) == 4) or (filter== 'b' and len(items) == 5) or (filter =='c'):
                    gene2homolog[items[0]] =  '\t'.join([items[1], items[2], items[3]])
            
            if len(items) == 5 and items[4] == "missed":
                missedCount += 1
                sys.stdout.write("%s\tmissed\n" %(items[0]))
    
    sys.stdout.write("\nTotal number of genes: %d\n" %(totalCount))
    sys.stdout.write("Number of NA genes: %d\n" %(NAcount))
    #sys.stdout.write("Number of genes that are absent in at least one species: %d\n" %(missedCount))
    f.close()
    return gene2homolog

def readChains(file):
    f = open(file, "r")
    gene2chain = {}
    for line in f.readlines():
        items = line.strip().split('\t')
        if len(items) != 2:
            sys.stderr.write("wrong format, line *%s*\n" %(line))
        if items[0] in gene2chain:
            sys.stderr.write("some repeated gene: %s, please check...\n" %(items[0]))
        else:
            gene2chain[items[0]] = items[1]
    f.close()
    return gene2chain

def mergeResults(gene2homolog, gene2chain, outFile):
    f = open(outFile, "w")
    cat2counts = {"nono":{"000":0, "001":0, "011": 0, "100":0, "101":0, "110":0, "111":0}, \
                  "noyes":{"000":0, "001":0, "011":0, "100":0, "101":0, "110":0, "111":0}, \
                  "yesno":{"000":0, "001":0, "011":0, "100":0, "101":0, "110":0, "111":0}, \
                  "yesyes":{"000":0, "001":0, "011":0,  "100":0, "101":0, "110":0, "111":0}}

    for gene in gene2homolog:
        if gene not in gene2chain:
            sys.stderr.write("gene %s does not have chain info\n" %(gene))
        else:
            f.write("%s\t%s\t%s\n" %(gene, gene2chain[gene], gene2homolog[gene]))
            dups = gene2homolog[gene].split('\t')
            dupcat = dups[1] + dups[2]
            chaincat = gene2chain[gene]
            chaincatStr = chaincat[1] + chaincat[2] + chaincat[4]
            cat2counts[dupcat][chaincatStr] += 1

    #for gene in gene2chain:
    #    if gene not in gene2homolog:
    #        sys.stderr.write("gene %s does not have homolog info\n" %(gene))
    f.close()
    return cat2counts

def getTotal(tab, dupcat):
    total = 0
    for k in tab[dupcat]:
        total += tab[dupcat][k]
    return total

def getTotal2(tab, chaincat):
    total = 0
    for dupcat in tab:
        total += tab[dupcat][chaincat]
    return total

def getTotal3(tab):
    total = 0
    for dupcat in tab:
        for chaincat in tab[dupcat]:
            total += tab[dupcat][chaincat]
    return total


def printCategoryTable(tab, f, name):
    tableHeader(f)
    f.write(" 1TC-1TL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["000"], tab["noyes"]["000"], tab["yesno"]["000"], tab["yesyes"]["000"], getTotal2(tab, "000") ))
    f.write("\\hline\n")
    f.write(" 1TC-MTL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["001"], tab["noyes"]["001"], tab["yesno"]["001"], tab["yesyes"]["001"], getTotal2(tab, "001") ))
    f.write("\\hline\n")
    f.write(" 1TC-MTCE-MTL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["011"], tab["noyes"]["011"], tab["yesno"]["011"], tab["yesyes"]["011"], getTotal2(tab, "011") ))
    f.write("\\hline\n")
    f.write(" MTC-1TL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["100"], tab["noyes"]["100"], tab["yesno"]["100"], tab["yesyes"]["100"], getTotal2(tab, "100") ))
    f.write("\\hline\n")
    f.write(" MTC-MTL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["101"], tab["noyes"]["101"], tab["yesno"]["101"], tab["yesyes"]["101"], getTotal2(tab, "101") ))
    f.write("\\hline\n")
    f.write(" MTC-MTCE-1TL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["110"], tab["noyes"]["110"], tab["yesno"]["110"], tab["yesyes"]["110"], getTotal2(tab, "110") ))
    f.write("\\hline\n")
    f.write(" MTC-MTCE-MTL & %d & %d & %d & %d & %d\\\\\n" %( tab["nono"]["111"], tab["noyes"]["111"], tab["yesno"]["111"], tab["yesyes"]["111"], getTotal2(tab, "111") ))
    f.write("\\hline\n")
    f.write("Total & %d & %d & %d & %d & %d\\\\\n" %(getTotal(tab, "nono"), getTotal(tab, "noyes"), getTotal(tab, "yesno"), getTotal(tab, "yesyes"), getTotal3(tab)))
    f.write("\\hline\n")

    capStr = "%s, assessment of cactus chain structure. The rows are categories of cactus chain structure of the genes. 1TC means there is one top chain across the gene, and MTC indicates there are two or more top chains across the gene. 1TL means all exons of the gene has the same top chain level. On the opposite, MTL means the exons of the gene have different top chain levels (Multiple Top Levels). MTCE indicates that the gene contains at least one exon with multiple top chains. The columns are categories of the genes by the locations of duplications happening across the CDS region. None means there is no duplication infered. Intronic means there are duplications in the intronic region of the gene. Exonic means there are duplications in the exonic region of the gene. Intronic and Exonic indicate that there are duplications in both exonic and intronic region of the gene." %(name)
    tableCloser(f, capStr)

def tableHeader(f):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\hline\n")
    f.write("Chain Structure & None & Intronic & Exonic & Intronic and Exonic & Total\\\\\n") 
    f.write("\\hline\n")

def tableCloser(f, captionStr):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\end{center}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{table}\n")

def writeDocumentStart(f):
    f.write("\\documentclass[11pt]{article}\n")
    f.write("\\usepackage{epsfig}\n")
    f.write("\\usepackage{multirow}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{array}\n")
    f.write("\\usepackage{slashbox}\n")
    f.write("\\usepackage{color}\n")
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



def main():
    usg = "usage: %prog homologInfoFile chainInfoFile outputFile latexFile"
    parser = OptionParser(usage=usg)
    parser.add_option("-f", "--filter", dest="filter", help="'a' to include only genes that present in all species.\
                       'b' to only include genes that present in at least one other species, but not present in all.\
                       'c' to include both genes in 'a' and 'b'\n", default="c")
    parser.add_option("-n", "--name", dest="name", help="Name of the run. Eg 'Primates'", default="c")
    (options, args) = parser.parse_args()
    gene2homolog = readHomologs(args[0], options.filter)
    gene2chain = readChains(args[1])
    tab = mergeResults(gene2homolog, gene2chain, args[2])

    f = open(args[3], "w")
    writeDocumentStart(f)
    printCategoryTable(tab, f, options.name)
    writeDocumentEnd(f)
    f.close()

if __name__ == "__main__":
    main()



