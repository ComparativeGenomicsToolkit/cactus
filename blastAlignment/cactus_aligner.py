#!/usr/bin/env python

"""Script for computing alignments for a reconstruction problem.
"""

import os

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system
from workflow.jobTree.scriptTree.target import Target

from cactus.blastAlignment.cactus_batch import makeBlastFromOptions
from cactus.blastAlignment.cactus_batch import makeStandardBlastOptions

from cactus.blastAlignment.cactus_alignerTestAligner import MakeBlastsLoader as MakeBlastsTest
    
class MakeSequences(Target):
    """Take a reconstruction problem and generate the sequences to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, netDisk, netName, resultsFile, blastOptions):
        Target.__init__(self)
        self.netDisk = netDisk
        self.netName = netName
        self.resultsFile = resultsFile
        self.blastOptions = blastOptions
        
    def run(self, localTempDir, globalTempDir):
        ##########################################
        #Setup the temp files
        ##########################################
        
        tempSeqFile = getTempFile(rootDir=globalTempDir)
        tempResultsFile = getTempFile(rootDir=globalTempDir)
        
        logger.info("Built temporary files")
        
        ##########################################
        #Construct the sequences file for doing all against all blast.
        ##########################################
        
        system("cactus_aligner %s %s %s" % (self.netDisk, self.netName, tempSeqFile))
        
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
        Target.__init__(self)
        self.tempSeqFile = tempSeqFile
        self.tempResultsFile = tempResultsFile
        self.resultsFile = resultsFile    
        
    def cleanup(self, localTempDir, globalTempDir):
        os.remove(self.tempSeqFile)
        logger.info("Removed the temporary fasta file for the blast step")
    
    def run(self, localTempDir, globalTempDir):
        ##########################################
        #Translate the coordinates
        ##########################################
        
        system("cactus_batch_convertCoordinates %s %s" % (self.tempResultsFile, self.resultsFile))
        logger.info("Translated the coordinates of the alignments to the final file: %s", self.resultsFile)

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = getBasicOptionParser("usage: %prog [options]", "%prog 0.1")
    
    parser.add_option("--job", dest="jobFile", 
                      help="Job file containing command to run")
    
    parser.add_option("--netDisk", dest="netDisk", 
                      help="The path to the net-disk")
    
    parser.add_option("--netName", dest="netName", 
                      help="The name of the net in which to get the sequences to align")
    
    parser.add_option("--useDummy", dest="useDummy", action="store_true",
                      help="Use a dummy blast aligner target (for testing)",
                      default=False)
    
    parser.add_option("--resultsFile", dest="resultsFile", type="string",
                      help="The file to put the alignments in")
        
    parsedOptions, args = parseBasicOptions(parser)
        
    assert len(args) == 0
    logger.info("Parsed arguments")
     
    blastOptions = makeBlastFromOptions(makeStandardBlastOptions())
    if parsedOptions.useDummy:
        blastOptions = MakeBlastsTest()
    
    firstTarget = MakeSequences(parsedOptions.netDisk, 
                                parsedOptions.netName, 
                                parsedOptions.resultsFile, 
                                blastOptions)
    firstTarget.execute(parsedOptions.jobFile)
    
    logger.info("Ran the first target okay")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
