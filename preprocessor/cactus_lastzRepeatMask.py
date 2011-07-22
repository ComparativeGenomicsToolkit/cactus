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
                     default="12")
    
    parser.add_option("--lastzOpts", dest="lastzOptions", type="string",
                      help="lastz options for repeat identification",
                      default="")
    
    parser.add_option("--lastzCmd", dest="lastzCmd", type="string",
                      help="lastz executable",
                      default="lastz")
    
    options, args = parser.parse_args()
    
    if len(args) != 3:
        parser.print_help()
        return 1
    
    queryFile = args[0]
    targetFile = args[1]
    outputFile = args[2]
    
    assert os.path.isfile(queryFile)
    assert os.path.isfile(targetFile)
    assert options.fragment > 1
    
    if testExec("fasta_fragments.py") == False or testExec("fasta_softmask_intervals.py") == False:
        print "ERROR: fasta_fragments.py or fasta_softmask_intervals.py" + \
               " not in PATH.\n Ensure that the latest lastz is installed" + \
               " and lastz/tools is in path\n\n"
        raise
    
    # make sure fragment size is even so they can overlap by exactly one half. 
    if options.fragment % 2:
        options.fragment += 1
    
    # make temporary working directory in same path as output
    currentDir = os.path.dirname(args[1])
    tempDir = tempfile.mkdtemp(dir=currentDir)
    fragFile = os.path.join(tempDir, "fragFile.fa")
    maskInfoFile = os.path.join(tempDir, "maskFile.dat")
    cleanTargetFile = os.path.join(tempDir, "cleanTargetFile.fa")
    
    
    # strip out the |1| strings from header or they will be eaten by lastz!!
    pipeCode = "__#%x__" % random.randint(0, 16777215)
    assert os.system("sed -e \"s/|1|/%s/g\" %s > %s" % (pipeCode, targetFile, cleanTargetFile)) == 0

    try:
        # chop up input fasta file into into fragments of specified size.  fragments overlap by 
        # half their length. 
        fragCmdLine = 'cat ' + queryFile + ' | fasta_fragments.py ' + '--fragment=' + \
                        str(options.fragment) + ' --step=' + str(options.fragment / 2) + ' > ' + fragFile
        
        assert os.system(fragCmdLine) == 0
        
        # lastz each fragment against the entire input sequence.  Each time a fragment aligns to a base
        # in the sequence, that base's match count is incremented.  base's whose match count exceeds the 
        # input period are softmasked.  
        # Note, we multiply period by two because we specify that every base is covered by two
        # overlapping fragments
        # Also note: repeats already masked in the input sequence are ignored (as by default in lastz).  
        # This behaviour can be changed by something like replacing [multiple] with [multiple,unmask]'
        lastzCmdLine = options.lastzCmd + ' ' + cleanTargetFile + '[multiple] ' + fragFile + ' ' + options.lastzOptions + \
                        ' --masking=' + str(options.period * 2) + ' --outputmasking+:soft=' + maskInfoFile + \
                        ' --format=none' 

        assert os.system(lastzCmdLine) == 0
        
        # reinsert the |1|'s
        assert os.system("sed -i -e \"s/%s/|1|/g\" %s" % (pipeCode, maskInfoFile)) == 0
        
        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
        softMaskCmdLine = 'cat ' + targetFile + ' | fasta_softmask_intervals.py --origin=1 ' + maskInfoFile + \
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
