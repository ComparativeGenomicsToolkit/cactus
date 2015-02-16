#!/usr/bin/env python

## USE LASTZ TO SOFTMASK REPEATS OF A GIVEN FASTA SEQUENCE FILE.  

############################################################
##  NOTE:   NEW VERSION OF LASTZ REQUIRED (MINIMUM 1.02.40)
##          lastz/tools MUST BE IN SYSTEM PATH
############################################################

import os
import re
import sys
from optparse import OptionParser
import tempfile
import shutil
import random
from sonLib.bioio import system
from sonLib.bioio import catFiles

def testExec(exeName):
    for dir in os.getenv("PATH").split(':'):                                           
        if (os.path.exists(os.path.join(dir, exeName))):
            return True
    return False

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <query> <target> <output>\n\n" + \
            "    <query>:  fasta sequence to search for repeats\n" + \
            "    <target>: fasta sequence to mask\n" + \
            "    <output>: softmasked version of <target>\n\n" + \
            "Example: %prog genome.fa chunk.fa chunk.masked.fa\n\n" 
    description = "softrepeat mask a fasta file using lastz.\n" + \
                    "NOTE: NEW VERSION OF LASTZ REQUIRED (MINIMUM 1.02.40)\n" + \
                    "lastz tools/ dir MUST BE IN SYSTEM PATH"
    parser = OptionParser(usage=usage, description=description)

    #output stuff
    parser.add_option("--fragment", dest="fragment", type="int",
                     help="The size of chunks passed to lastz (must be at least twice as big as overlap)",
                     default=200)
    
    parser.add_option("--minPeriod", dest="period", type="int",
                     help="minimum number of occurrences of a sequence for it to be masked",
                     default=10)
    
    parser.add_option("--lastzOpts", dest="lastzOptions", type="string",
                      help="lastz options for repeat identification",
                      default="")
    
    parser.add_option("--lastzCmd", dest="lastzCmd", type="string",
                      help="lastz executable",
                      default="cPecanLastz")
    
    parser.add_option("--unmaskInput", dest="unmaskInput", action="store_true",
                      help="Makes any previous masking of the input sequence invisible to the repeat masking process",
                      default=False)
    
    parser.add_option("--unmaskOutput", dest="unmaskOutput", action="store_true",
                      help="Discards any previous masking from the output sequence, uses just the masking discovered by lastz",
                      default=False)
    
    parser.add_option("--proportionSampled", dest="proportionSampled", type="float",
                     help="The amount of the genome that is being sampled for masking, used to adjust the minPeriod parameter according to sampling",
                     default="1.0")
    
    parser.add_option("--tempDir", dest="tempDir",
                     help="Location in which to place to temporary files",
                     default=os.getcwd())
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        return 1
    
    queryFile = args[0]
    outputFile = args[1]
    targetFiles = sys.stdin.readline().split() #Read them from stdin
    
    assert os.path.isfile(queryFile)
    assert len(targetFiles) >= 1
    for targetFile in targetFiles:
        assert os.path.isfile(targetFile)
    assert options.fragment > 1
    
    #Adjust the period parameter using the amount of genome sampled
    options.period = max(1, round(options.proportionSampled * options.period))
    
    # make sure fragment size is even so they can overlap by exactly one half. 
    if options.fragment % 2:
        options.fragment += 1
    
    # make temporary working directory in same path as output
    tempDir = tempfile.mkdtemp(dir=options.tempDir)
    maskInfoFile = os.path.join(tempDir, "maskFile.dat")
    targetFile = os.path.join(tempDir, "target.fa")

    try:
        #Make temporary target file, if more than one file
        catFiles(targetFiles, targetFile)
        
        # chop up input fasta file into into fragments of specified size.  fragments overlap by 
        # half their length. 
        fragCmdLine = 'cat ' + queryFile + ' | cactus_fasta_fragments.py ' + '--fragment=' + \
                        str(options.fragment) + ' --step=' + str(options.fragment / 2) + " --origin=zero "
        
        # lastz each fragment against the entire input sequence.  Each time a fragment aligns to a base
        # in the sequence, that base's match count is incremented.  
        # the plus three for the period parameter is a fudge to ensure sufficient alignments are found
        lastZSequenceHandling  = '[multiple][nameparse=darkspace] /dev/stdin[nameparse=darkspace] '
        if options.unmaskInput:
            lastZSequenceHandling  = '[multiple,unmask][nameparse=darkspace] /dev/stdin[unmask][nameparse=darkspace] '
        lastzCmdLine = options.lastzCmd + ' ' + targetFile + \
        lastZSequenceHandling + options.lastzOptions + \
        (' --querydepth=keep,nowarn:%i --format=general:name1,zstart1,end1,name2,zstart2+,end2+ --markend ' % \
         (options.period+3))

        #This runs Bob's covered intervals program, which combins the lastz alignment info into intervals of the query.
        coveredIntervalsCmdLine = "cactus_covered_intervals --queryoffsets --origin=one M=%s > %s" % (int(options.period*2), maskInfoFile)

        system(fragCmdLine + ' | ' + lastzCmdLine + ' | ' + coveredIntervalsCmdLine)

        #open(maskInfoFile, "w").close()
        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
        unmaskString = ""
        if options.unmaskOutput:
            unmaskString = "--unmask"
        softMaskCmdLine = 'cat ' + queryFile + (' | cactus_fasta_softmask_intervals.py --origin=one %s ' % unmaskString) + \
            maskInfoFile +  ' > ' + outputFile
    
        system(softMaskCmdLine)
    
    except Exception, e:
        # delete the temporary files
        shutil.rmtree(tempDir)
        raise e

    # delete the temporary files
    shutil.rmtree(tempDir)
        
if __name__ == '__main__':
    exit(main())
