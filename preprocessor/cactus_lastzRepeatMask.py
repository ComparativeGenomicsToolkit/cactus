#!/usr/bin/env python

## USE LASTZ TO SOFTMASK REPEATS OF A GIVEN FASTA SEQUENCE FILE.  

############################################################
##  NOTE:   NEW VERSION OF LASTZ REQUIRED (MINIMUM 1.02.40)
##          lastz/tools MUST BE IN SYSTEM PATH
############################################################

import os
import re
from optparse import OptionParser
import tempfile
import shutil
import random
from sonLib.bioio import system

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
                     default=100)
    
    parser.add_option("--minPeriod", dest="period", type="int",
                     help="minimum number of occurrences of a sequence for it to be masked",
                     default="12")
    
    parser.add_option("--lastzOpts", dest="lastzOptions", type="string",
                      help="lastz options for repeat identification",
                      default="")
    
    parser.add_option("--lastzCmd", dest="lastzCmd", type="string",
                      help="lastz executable",
                      default="cactus_lastz")
    
    options, args = parser.parse_args()
    
    if len(args) != 3:
        parser.print_help()
        return 1
    
    queryFile = args[0]
    outputFile = args[1]
    targetFiles = args[2:]
    
    assert os.path.isfile(queryFile)
    assert os.path.isfile(targetFile)
    assert options.fragment > 1
    
    if testExec("fasta_softmask_intervals.py") == False:
        raise RuntimeError("ERROR: fasta_softmask_intervals.py" + \
               " not in PATH.\n Ensure that the latest lastz is installed" + \
               " and lastz/tools is in path\n\n")
    
    # make sure fragment size is even so they can overlap by exactly one half. 
    if options.fragment % 2:
        options.fragment += 1
    
    # make temporary working directory in same path as output
    currentDir = os.path.dirname(args[1])
    tempDir = tempfile.mkdtemp(dir=currentDir)
    maskInfoFile = os.path.join(tempDir, "maskFile.dat")
    targetFile = os.path.join(tempDir, "target.fa")
    
    #Make temporary target file, if more than one file
    if len(targetFiles) > 0:
        system("cat %s > %s" % (" ".join(targetFiles)), targetFile)
    else:
        system("ln %s %s" % targetFiles[0], targetFile)

    try:
        # chop up input fasta file into into fragments of specified size.  fragments overlap by 
        # half their length. 
        fragCmdLine = 'cat ' + queryFile + ' | cactus_fasta_fragments.py ' + '--fragment=' + \
                        str(options.fragment) + ' --step=' + str(options.fragment / 2)
        
        # lastz each fragment against the entire input sequence.  Each time a fragment aligns to a base
        # in the sequence, that base's match count is incremented.  
        # Note, we multiply period by two because we specify that every base is covered by two
        # overlapping fragments
        
        ##Not currently true
        # Also note: repeats already masked in the input sequence are ignored (as by default in lastz).  
        # This behaviour can be changed by something like replacing [multiple] with [multiple,unmask]'
        lastzCmdLine = options.lastzCmd + ' ' + targetFile + \
        '[multiple,unmask][nameparse=darkspace] /dev/stdin[unmask][nameparse=darkspace] ' + \
        options.lastzOptions + \
        (' --querydepth=keep,nowarn:%i --format=general:name1,zstart1,end1,name2,zstart2+,end2+ --markend ' % \
         (options.period * 2 + 1))

        #This runs Bob's covered intervals program, which combins the lastz alignment info into intervals of the query.
        coveredIntervalsCmdLine = "cactus_covered_intervals  --queryoffsets --origin=one > %s" % maskInfoFile

        system(fragCmdLine + ' | ' + lastzCmdLine + ' | ' + coveredIntervalsCmdLine)

        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
        softMaskCmdLine = 'cat ' + queryFile + ' | cactus_fasta_softmask_intervals.py --origin=1 ' + maskInfoFile + \
                            ' > ' + outputFile
    
        system(softMaskCmdLine)
    
    except Exception, e:
        # delete the temporary files
        shutil.rmtree(tempDir)
        raise e

    # delete the temporary files
    shutil.rmtree(tempDir)
        
if __name__ == '__main__':
    exit(main())
