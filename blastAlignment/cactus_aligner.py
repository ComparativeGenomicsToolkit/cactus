#!/usr/bin/env python

#Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script for computing alignments for a reconstruction problem.
"""

import os
from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from cactus.blastAlignment.cactus_batch import makeBlastFromOptions
from cactus.blastAlignment.cactus_batch import makeStandardBlastOptions

from cactus.blastAlignment.cactus_alignerTestAligner import MakeBlastsLoader as MakeBlastsTest
    
class MakeSequences(Target):
    """Take a reconstruction problem and generate the sequences to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, flowerName, resultsFile, blastOptions):
        Target.__init__(self, time=0.0099)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.resultsFile = resultsFile
        self.blastOptions = blastOptions
        
    def run(self):
        ##########################################
        #Setup the temp files
        ##########################################
        
        tempSeqFile = getTempFile(rootDir=self.getGlobalTempDir())
        tempResultsFile = getTempFile(rootDir=self.getGlobalTempDir())
        
        logger.info("Built temporary files")
        
        ##########################################
        #Construct the sequences file for doing all against all blast.
        ##########################################
        
        system("cactus_aligner '%s' %s %s" % (self.cactusDisk, self.flowerName, tempSeqFile))
        
        logger.info("Got the sequence files to align")
        
        ##########################################
        #Make blast target
        ##########################################
        
        self.addChildTarget(self.blastOptions.makeBlastOptions([ tempSeqFile ], tempResultsFile))
        logger.info("Added child target okay")
        
        ##########################################
        #Setup follow on coordinates
        ##########################################
        
        self.setFollowOnTarget(ModifyBlasts(tempSeqFile, tempResultsFile, self.resultsFile))
        logger.info("Created modify blasts target")
    
class ModifyBlasts(Target):
    """Modifies the alignments file so that the sequences have the correct coordinates.
    """
    
    def __init__(self, tempSeqFile, tempResultsFile, resultsFile):
        Target.__init__(self, time=0.062069)
        self.tempSeqFile = tempSeqFile
        self.tempResultsFile = tempResultsFile
        self.resultsFile = resultsFile    
    
    def run(self):
        os.remove(self.tempSeqFile)
        logger.info("Removed the temporary fasta file for the blast step")
        
        ##########################################
        #Translate the coordinates
        ##########################################
        
        tempFile = os.path.join(self.getLocalTempDir(), "temp.txt")
        fileHandle = open(tempFile, 'w')
        fileHandle.write(self.tempResultsFile)
        fileHandle.close()
        system("cactus_batch_convertCoordinates %s %s" % (tempFile, self.resultsFile))
        logger.info("Translated the coordinates of the alignments to the final file: %s", self.resultsFile)

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    
    parser.add_option("--cactusDisk", dest="cactusDisk", 
                      help="The path to the flower-disk")
    
    parser.add_option("--flowerName", dest="flowerName", 
                      help="The name of the flower in which to get the sequences to align")
    
    parser.add_option("--useDummy", dest="useDummy", action="store_true",
                      help="Use a dummy blast aligner target (for testing)",
                      default=False)
    
    parser.add_option("--resultsFile", dest="resultsFile", type="string",
                      help="The file to put the alignments in")
        
    parsedOptions, args = parser.parse_args()
        
    assert len(args) == 0
     
    blastOptions = makeBlastFromOptions(makeStandardBlastOptions())
    if parsedOptions.useDummy:
        blastOptions = MakeBlastsTest()
    
    firstTarget = MakeSequences(parsedOptions.cactusDisk, 
                                parsedOptions.flowerName, 
                                parsedOptions.resultsFile, 
                                blastOptions)
    
    Stack(firstTarget).startJobTree(parsedOptions)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blastAlignment.cactus_aligner import *
    _test()
    main()
