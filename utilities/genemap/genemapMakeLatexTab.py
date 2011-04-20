#!/usr/bin/env python

import os, re, sys, subprocess
from optparse import OptionParser

def cactusVSmultizTab(f, tab, tabName, altColor):
    keys = ["CC", "NC", "CN", "NN"]
    for k in keys:
        if k not in tab:
            tab[k] = [0, 0.00]
    totalRowYes = tab["CC"][0] + tab["CN"][0]
    totalRowNo = tab["NC"][0] + tab["NN"][0]
    totalRowYesPercent = tab["CC"][1] + tab["CN"][1]
    totalRowNoPercent = 100.0 - totalRowYesPercent

    totalColYes = tab["CC"][0] + tab["NC"][0]
    totalColNo = tab["CN"][0] + tab["NN"][0]
    totalColYesPercent = tab["CC"][1] + tab["NC"][1]
    totalColNoPercent = 100.0 - totalColYesPercent

    total = totalColYes + totalColNo
    totalPercent = 100.00
   
    if not altColor:
        f.write("\\multirow{3}{*}{%s} & \\cellcolor[gray]{0.9} Conserved & \\cellcolor[gray]{0.9} \\textbf{%d} & \\cellcolor[gray]{0.9} \\textbf{%.2f}\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% &\\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%%\\\\\n" % (tabName, tab["CC"][0], tab["CC"][1], tab["CN"][0], tab["CN"][1], totalRowYes, totalRowYesPercent) )
        #f.write("\\cline{2-8}\n")
        f.write("& Non-conserved & %d & %.2f\\%% & \\textbf{%d} & \\textbf{%.2f}\\%% & %d & %.2f\\%%\\\\\n" % (tab["NC"][0], tab["NC"][1], tab["NN"][0], tab["NN"][1], totalRowNo, totalRowNoPercent) )
        #f.write("\\cline{2-8}\n")
        f.write("& \\cellcolor[gray]{0.9} Total & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%%\\\\\n" % (totalColYes, totalColYesPercent, totalColNo, totalColNoPercent, total, totalPercent) )
    
    else:
        f.write("\\multirow{3}{*}{%s} & Conserved & \\textbf{%d} & \\textbf{%.2f}\\%% & %d &  %.2f\\%% & %d & %.2f\\%%\\\\\n" % (tabName, tab["CC"][0], tab["CC"][1], tab["CN"][0], tab["CN"][1], totalRowYes, totalRowYesPercent) )
        #f.write("\\cline{2-8}\n")
        f.write("& \\cellcolor[gray]{0.9} Non-conserved & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} \\textbf{%d} & \\cellcolor[gray]{0.9} \\textbf{%.2f}\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%%\\\\\n" % (tab["NC"][0], tab["NC"][1], tab["NN"][0], tab["NN"][1], totalRowNo, totalRowNoPercent) )
        #f.write("\\cline{2-8}\n")
        f.write("& Total & %d & %.2f\\%% & %d & %.2f\\%% & %d & %.2f\\%%\\\\\n" % (totalColYes, totalColYesPercent, totalColNo, totalColNoPercent, total, totalPercent) )

    f.write("\\hline\n")
    f.write("\n")


def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")
    
def cactusVSmultizTableHeader(f, title):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n") 
    f.write("\\scalebox{1}{%\n")
    #f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\begin{tabular}{cccccccc}\n")
    f.write("\\multicolumn{8}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("\\multicolumn{2}{c}{\\multirow{2}{*}{\\backslashbox{Cactus}{Multiz}}} & \\multicolumn{2}{c}{Conserved} & \\multicolumn{2}{c}{Non-conserved} & \\multicolumn{2}{c}{Total} \\\\\n")
    #f.write("\\multicolumn{2}{c}{\\multirow{2}{*}{}} & \\multicolumn{2}{c}{Conserved} & \\multicolumn{2}{c}{Non-conserved} & \\multicolumn{2}{c}{Total} \\\\\n")
    f.write("\\cline{3-8}\n")
    f.write("\\multicolumn{2}{c}{} & Count & Percentage & Count & Percentage & Count & Percentage\\\\\n")
    f.write("\\hline\n")

def speciesTable(f, exonTab, geneTab, species):
    #for spc in geneTab:
    spcList = species.split(" ")
    altColor = -1
    for spc in spcList:
        if spc not in exonTab:
            continue
        elist = exonTab[spc]
        glist = geneTab[spc]
        if altColor == 1:
            f.write("%s & %s/%s & %s\\%% & %s/%s & %s\\%% & %s/%s & %s\\%% & %s/%s & %s\\%%\\\\\n" %\
                 (spc, glist[0], glist[1], glist[2], glist[3], glist[4], glist[5], \
                      elist[0], elist[1], elist[2], elist[3], elist[4], elist[5]))
        else:
            f.write("\\cellcolor[gray]{0.9} %s &\\cellcolor[gray]{0.9} %s/%s & \\cellcolor[gray]{0.9} %s\\%% & \\cellcolor[gray]{0.9} %s/%s & \\cellcolor[gray]{0.9} %s\\%% & \\cellcolor[gray]{0.9} %s/%s &\\cellcolor[gray]{0.9} %s\\%% & \\cellcolor[gray]{0.9} %s/%s & \\cellcolor[gray]{0.9} %s\\%%\\\\\n" %\
                 (spc, glist[0], glist[1], glist[2], glist[3], glist[4], glist[5], \
                      elist[0], elist[1], elist[2], elist[3], elist[4], elist[5]))

        altColor *= (-1)
    
    f.write("\\hline\n")

def speciesTableHeader(f, title):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n") 
    f.write("\\scalebox{0.9}{%\n")
    #f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\begin{tabular}{ccccccccc}\n")
    f.write("\\multicolumn{9}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("\\multirow{3}{*}{Species} & \\multicolumn{4}{c}{Genes} & \\multicolumn{4}{c}{Exons} \\\\\n")
    f.write("\\cline{2-9}\n")
    f.write("& \\multicolumn{2}{c}{Cactus} & \\multicolumn{2}{c}{Multiz} & \\multicolumn{2}{c}{Cactus} & \\multicolumn{2}{c}{Multiz}\\\\\n")
    f.write("\\cline{2-9}\n")
    f.write("& Count & Percentage & Count & Percentage & Count & Percentage & Count & Percentage\\\\\n")
    f.write("\\hline\n")

#======== supplement table ==========
def cactusVSmultizFullTableHeader(f, title):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n") 
    f.write("\\scalebox{0.85}{%\n")
    #f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\begin{tabular}{cccccccccc}\n")
    f.write("\\multicolumn{10}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("\\multicolumn{2}{c}{ \\multirow{2}{*}{\\backslashbox{Cactus}{Multiz}} } & \\multicolumn{2}{c}{Conserved} & \\multicolumn{2}{c}{Non-conserved} & \\multicolumn{2}{c}{Missing-data} & \\multicolumn{2}{c}{Total} \\\\\n")
    f.write("\\cline{3-10}\n")
    f.write("\\multicolumn{2}{c}{} & Count & Percentage & Count & Percentage & Count & Percentage & Count & Percentage\\\\\n")
    f.write("\\hline\n")


def cactusVSmultizFullTab(f, tab, tabName, altColor):
    keys = ["CC", "NC", "CN", "NN", "MDC", "CMD", "MDN", "NMD", "MDMD", "NAN", "NNA", "NAC", "CNA", "NANA", "NAMD", "MDNA"]
    for k in keys:
        if k not in tab:
            tab[k] = [0, 0.00]

    #Count NA as non-conserved:
    rowC = [0, 0.00]
    rowN = [0, 0.00]
    rowMD = [0, 0.00]
    colC = [0, 0.00]
    colN = [0, 0.00]
    colMD = [0, 0.00]
    total = [0, 0.00]
    for i in range(2):
        tab["CN"][i] = tab["CN"][i] + tab["CNA"][i]
        tab["NC"][i] = tab["NC"][i] + tab["NAC"][i]
        tab["NN"][i] = tab["NN"][i] + tab["NNA"][i] + tab["NAN"][i] + tab["NANA"][i]
        tab["MDN"][i] = tab["MDN"][i] + tab["MDNA"][i]
        tab["NMD"][i] = tab["NMD"][i] + tab["NAMD"][i]
        
        rowC[i] = tab["CC"][i] + tab["CN"][i] + tab["CMD"][i]
        rowN[i] = tab["NC"][i] + tab["NN"][i] + tab["NMD"][i]
        rowMD[i] = tab["MDC"][i] + tab["MDN"][i] + tab["MDMD"][i]
        
        colC[i] = tab["CC"][i] + tab["NC"][i] + tab["MDC"][i]
        colN[i] = tab["CN"][i] + tab["NN"][i] + tab["MDN"][i]
        colMD[i] = tab["CMD"][i] + tab["NMD"][i] + tab["MDMD"][i]


    total[0] = colC[0] + colN[0] + colMD[0]
    total[1] = 100.00
    rowMD[1] = total[1] - (rowC[1] + rowN[1])
    colMD[1] = total[1] - (colC[1] + colN[1])
    
    #cellColorStr = "\cellcolor[gray]{0.9}"
    if not altColor:
        f.write("\\multirow{4}{*}{%s} & \\cellcolor[gray]{0.9} Conserved & \\cellcolor[gray]{0.9} \\textbf{%d} & \\cellcolor[gray]{0.9} \\textbf{%.2f}\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% &\\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% \\\\\n" % (tabName, tab["CC"][0], tab["CC"][1], tab["CN"][0], tab["CN"][1], tab["CMD"][0], tab["CMD"][1], rowC[0], rowC[1]) )
        f.write("& Non-conserved & %d & %.2f\\%% & \\textbf{%d} & \\textbf{%.2f}\\%% & %d & %.2f\\%% & %d & %.2f\\%% \\\\\n" % (tab["NC"][0], tab["NC"][1], tab["NN"][0], tab["NN"][1], tab["NMD"][0], tab["NMD"][1], rowN[0], rowN[1]) )
        f.write("& \\cellcolor[gray]{0.9} Missing-data & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} \\textbf{%d} & \\cellcolor[gray]{0.9} \\textbf{%.2f} \\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f \\%% \\\\\n" % (tab["MDC"][0], tab["MDC"][1], tab["MDN"][0], tab["MDN"][1], tab["MDMD"][0], tab["MDMD"][1], rowMD[0], rowMD[1]) )
        f.write("& Total & %d & %.2f\\%% & %d & %.2f\\%% & %d & %.2f\\%% & \\textbf{%d} & \\textbf{%.2f}\\%% \\\\\n" % (colC[0], colC[1], colN[0], colN[1], colMD[0], colMD[1], total[0], total[1]) )
    
    else:
        f.write("\\multirow{4}{*}{%s} & Conserved & \\textbf{%d} & \\textbf{%.2f}\\%% & %d & %.2f\\%% & %d & %.2f\\%% & %d & %.2f\\%% \\\\\n" % (tabName, tab["CC"][0], tab["CC"][1], tab["CN"][0], tab["CN"][1], tab["CMD"][0], tab["CMD"][1], rowC[0], rowC[1]) )
        f.write("& \\cellcolor[gray]{0.9} Non-conserved & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} \\textbf{%d} & \\cellcolor[gray]{0.9} \\textbf{%.2f}\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% \\\\\n" % (tab["NC"][0], tab["NC"][1], tab["NN"][0], tab["NN"][1], tab["NMD"][0], tab["NMD"][1], rowN[0], rowN[1]) )
        f.write("& Missing-data &  %d &  %.2f\\%% & %d & %.2f\\%% & \\textbf{%d} & \\textbf{%.2f} \\%% &  %d & %.2f \\%% \\\\\n" % (tab["MDC"][0], tab["MDC"][1], tab["MDN"][0], tab["MDN"][1], tab["MDMD"][0], tab["MDMD"][1], rowMD[0], rowMD[1]) )
        f.write("& \\cellcolor[gray]{0.9} Total & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} %d & \\cellcolor[gray]{0.9} %.2f\\%% & \\cellcolor[gray]{0.9} \\textbf{%d} & \\cellcolor[gray]{0.9} \\textbf{%.2f}\\%% \\\\\n" % (colC[0], colC[1], colN[0], colN[1], colMD[0], colMD[1], total[0], total[1]) )


    f.write("\\hline\n")
    f.write("\n")



#======== end of supplement table =====

def writeDocumentStart(f):
    f.write("\\documentclass[11pt]{article}\n")
    f.write("\\usepackage{epsfig}\n")
    f.write("\\usepackage{multirow}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{array}\n")
    f.write("\\usepackage{slashbox}\n")
    f.write("\\usepackage{color}\n")
    f.write("\\usepackage[table]{xcolor}\n")
    f.write("\\usepackage{rotating}\n")
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

def readTabs(file):
    f = open(file, "r")
    exonTab = {}
    geneTab = {}
    spcExonTab = {}
    spcGeneTab = {}

    exonTabSup = {}
    geneTabSup = {}

    for line in f.readlines():
        if re.search('#', line):#skip comment line
            continue
        items = line.strip().split('\t')
        if items[0] == "exonTab":
            exonTab[items[1]] = [int(items[2]), float(items[3])]
        if items[0] == "geneTab":
            geneTab[items[1]] = [int(items[2]), float(items[3])]
        if items[0] == "spcExonTab":
            #spcExonTab[items[1]] = [int(items[2]), int(items[3]), float(items[4]), int(items[5]), int(items[6]), float(items[7])]
            spcExonTab[items[1]] = items[2:]
            #sys.stderr.write("Length: %d\n" %(len(spcExonTab[items[1]])))
        if items[0] == "spcGeneTab":
            spcGeneTab[items[1]] = items[2:]
            #sys.stderr.write("Length: %d\n" %(len(spcGeneTab[items[1]])))

        #Supplement tables:
        if items[0] == "exonTabSupplement":
            exonTabSup[items[1]] = [int(items[2]), float(items[3])]
        if items[0] == "geneTabSupplement":
            geneTabSup[items[1]] = [int(items[2]), float(items[3])]


    f.close()
    return (exonTab, geneTab, spcExonTab, spcGeneTab, exonTabSup, geneTabSup)

def main():
    usg = "usage: %prog [options] resultsFile outputFile"
    parser = OptionParser(usage=usg)
    parser.add_option("-n", "--name", dest = "name", help="Name of the run\n", default="")
    parser.add_option("-s", "--species", dest = "species", help="Specify the order of the species rows\n")
    (options, args) = parser.parse_args()

    (exonTab, geneTab, spcExonTab, spcGeneTab, exonTabSup, geneTabSup) = readTabs(args[0])
    f = open(args[1], "w")
    
    writeDocumentStart(f)
    title = "%s genic and exonic cross-species conservation" %options.name
    cactusVSmultizTableHeader(f, title)
    cactusVSmultizTab(f, geneTab, "Gene", False)
    f.write("\\hline\n")
    cactusVSmultizTab(f, exonTab, "Exon", True)
    captionStr = "Comparisons of genic and exonic cross-species conservation between Cactus and Multiz using the %s dataset. The `Conserved' and `Non-conserved' columns indicate conservation status inferred from Multiz. Similarly, the `Conserved' and `Non-conserved' rows show conservation status inferred from Cactus. There are totally four categories: Cactus conserved and Multiz conserved, Cactus conserved and Multiz non-conserved, Cactus non-conserved and Multiz conserved, and Cactus non-conserved and Multiz non-conserved. The first section, `Gene', shows the gene statistics, and the second section, `Exon', show the exon statistics. `Count' is the number of genes or exons that fall into each category. `Percentage' is the relative metric of `Count', which shows the percentage of the total genes or exons included in the analysis that fall into each category. The `Total' column shows the total counts (or percentages) of the row (Cactus) statistics, while the `Total' row shows the total counts (or percentages) of the column (Multiz) statistics." %(options.name)
    label = "encodeCrossSpeciesCons%s" %options.name
    tableCloser(f, captionStr, label)
   
    
    title = "%s conserved genes and exons per species" %options.name
    speciesTableHeader(f, title)
    speciesTable(f, spcExonTab, spcGeneTab, options.species)
    captionStr = "Comparisons of per-species conserved genes and exons between Cactus and Multiz using the %s dataset. The rows are different species. `Cactus' or `Multiz' indicates the alignment used to infer the conserved status. `Count': the numerator is the number of the reference (human) genes or exons that are inferred as conserved in each species, and the denominator is the total genes or exons of each species that were included in the analysis. `Percentage' is the corresponding relative metric of `Count'." %(options.name)
    label = "encodeSpeciesCons%s" %options.name
    tableCloser(f, captionStr, label)
    
    #========== supplement tables ============
    #title = "Cross-species conservation of genes and exons in %s, Cactus versus Multiz" %(options.name)
    #cactusVSmultizFullTableHeader(f, title)
    #cactusVSmultizFullTab(f, geneTabSup, "Gene", False)
    #f.write("\\hline\n")
    #cactusVSmultizFullTab(f, exonTabSup, "Exon", True)
    #captionStr = "Comparisons of cross-species conservation of genes and exons between Cactus and Multiz. `Conserved': exon/gene must be conserved in all species. `Non-conserved': exon/gene is non-conserved or absent in at least one of the species. `Missing-data': exon/gene is conserved in some species, but contains missing data or ambiguous bases in the rest."
    #label = "encodeCrossSpeciesCons%s" % options.name
    #tableCloser(f, captionStr, label)


    writeDocumentEnd(f)
    f.close()

if __name__ == "__main__":
    main()


