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

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <input FASTA sequence file> <output FASTA sequence file>"
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
                     default="12")
    
    parser.add_option("--lastzOpts", dest="lastzOptions", type="string",
                      help="lastz options for repeat identification",
                      default="")
    
    parser.add_option("--lastzCmd", dest="lastzCmd", type="string",
                      help="lastz executable",
                      default="lastz")
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        return 1
    
    inputFile = args[0]
    outputFile = args[1]
    
    assert os.path.isfile(inputFile)
    assert options.fragment > 1
    
    # make sure fragment size is even so they can overlap by exactly one half. 
    if options.fragment % 2:
        options.fragment += 1
    
    # make temporary working directory in same path as output
    currentDir = os.path.dirname(args[1])
    tempDir = tempfile.mkdtemp(dir=currentDir)
    fragFile = os.path.join(tempDir, "fragFile.fa")
    maskInfoFile = os.path.join(tempDir, "maskFile.dat")

    try:
        # chop up input fasta file into into fragments of specified size.  fragments overlap by 
        # half their length. 
        fragCmdLine = 'cat ' + inputFile + ' | fasta_fragments.py ' + '--fragment=' + \
                        str(options.fragment) + ' --step=' + str(options.fragment / 2) + ' > ' + fragFile
        
        assert os.system(fragCmdLine) == 0
        
        # lastz each fragment against the entire input sequence.  Each time a fragment aligns to a base
        # in the sequence, that base's match count is incremented.  base's whose match count exceeds the 
        # input period are softmasked.  
        # Note, we multiply period by two because we specify that every base is covered by two
        # overlapping fragments
        # Also note: repeats already masked in the input sequence are ignored (as by default in lastz).  
        # This behaviour can be changed by something like replacing [multiple] with [multiple,unmask]'
        lastzCmdLine = options.lastzCmd + ' ' + inputFile + '[multiple] ' + fragFile + ' ' + options.lastzOptions + \
                        ' --masking=' + str(options.period * 2) + ' --outputmasking+:soft=' + maskInfoFile + \
                        ' --format=none' 

        assert os.system(lastzCmdLine) == 0
        
        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
        softMaskCmdLine = 'cat ' + inputFile + ' | fasta_softmask_intervals.py --origin=1 ' + maskInfoFile + \
                            ' > ' + outputFile
    
        assert os.system(softMaskCmdLine) == 0
    
    except Exception, e:
        # delete the temporary files
        shutil.rmtree(tempDir)
        raise e

    # delete the temporary files
    shutil.rmtree(tempDir)
        
if __name__ == '__main__':
    main()
