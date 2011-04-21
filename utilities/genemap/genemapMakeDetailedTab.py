#!/usr/bin/env python

import sys, re, os
from optparse import OptionParser

def writeDocumentStart(f):                                                                    
    f.write("\\documentclass[11pt]{article}\n")                                               
    f.write("\\usepackage{epsfig}\n")                                                         
    f.write("\\usepackage{multirow}\n")                                                       
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{tabularx}\n")
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
                                           
def tableHeader(f, title, mode):
    f.write("\n")
    #f.write("\\begin{sidewaystable}\n")
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{0.75}{%\n")
    if mode == "gene":
        #f.write("\\begin{tabularx}{\\textwidth}{XXXXX}\n")
        #f.write("\\begin{tabular}{llll>{\\small}l}\n")
        f.write("\\begin{tabular}{rrrrr}\n")
        f.write("\multicolumn{5}{c}{%s} \\\\\n" %title)
    elif mode == "exon":
        #f.write("\\begin{tabular}{lllll>{\\small}l}\n")
        f.write("\\begin{tabular}{rrrrrr}\n")
        #f.write("\\begin{tabularx}{\\textwidth}{XXXXXX}\n")
        f.write("\multicolumn{6}{c}{%s} \\\\\n" %title)
    
    f.write("\\hline\n")
    f.write("\\hline\n")
    
    if mode == "gene":
        f.write("Gene & Accession & Cactus & Multiz & Comment \\\\\n")
    elif mode == "exon":
        f.write("Exon Id & Gene & Accession & Cactus & Multiz & Comment \\\\\n")
    
    f.write("\\hline\n")

def tableCloser(f, caption, label):
    #f.write("\\end{tabularx}\n")
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %caption)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    #f.write("\\end{sidewaystable}\n")
    f.write("\\end{table}\n")

def table(f, infile, mode):
    alt = -1
    infh = open(infile, "r")
    for line in infh.readlines():
        if re.search('#', line):
            continue
        line = line.strip().replace('_', '\_')
        items = line.split('\t')

        #hack
        if (mode == "exon" and len(items) == 8) or (mode == "gene" and len(items) == 7):
            lastitem = items[len(items) -1].split()
            commentStr = ' '.join(lastitem[1:])
            items[len(items) - 1] = "\t".join([lastitem[0], commentStr])

        line = '\t'.join(items[3:])
        #if (mode == "gene" and len(items) != 7) or (mode == "exon" and len(items) != 8):
        #    continue
               
        if alt == 1:
            line = line.replace('\t', ' & ')
            #f.write("%s & \\\\\n" % line)
            f.write("%s \\\\\n" % line)
        else:
            line = line.replace('\t', ' & \\cellcolor[gray]{0.9} ')
            #f.write("\\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} \\\\\n" % line)
            f.write("\\cellcolor[gray]{0.9} %s \\\\\n" % line)
        
        alt *= -1
    
    f.write("\\hline\n\n")
    infh.close()
    return

def writeGeneTable(f, options):
    title = "%s" % options.name
    tableHeader(f, title, "gene")
    table(f, options.input, "gene")
    caption = "%s" % options.name
    label = "detailedGeneTab%s" % options.name
    tableCloser(f, caption, label)

def writeExonTable(f, options):
    title = "%s" % options.name
    tableHeader(f, title, "exon")
    table(f, options.input, "exon")
    caption = "%s" % options.name
    label = "detailedGeneTab%s" % options.name
    tableCloser(f, caption, label)

def main():
    usg = "usage: %prog [options] inputFiles outputFile.tex"
    parser = OptionParser(usage=usg)
    parser.add_option("-n", "--name", dest="name", help="Name of the run\n", default="")
    parser.add_option("-o", "--output", dest="output", help="Name of the output file\n", default="out.tex")
    parser.add_option("-i", "--input", dest="input", help="Input file\n")
    parser.add_option("-m", "--mode", dest="mode", help="Specify 'gene' to print the geneTab, and 'exon' to print the exonTab\n")
    (options, args) = parser.parse_args()
    
    if not options.input:
        sys.stderr.write("Need input files\n")
        exit(2)
    f = open(options.output, "w")
    writeDocumentStart(f)

    if options.mode == "gene":
        writeGeneTable(f, options)
    elif options.mode == "exon":
        writeExonTable(f, options)

    writeDocumentEnd(f)
    f.close()

if __name__ == "__main__":
    main()
