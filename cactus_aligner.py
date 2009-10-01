#!/usr/bin/env python

"""Script for computing alignments for a reconstruction problem.
"""

import xml.etree.ElementTree as ET
import os

from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import fastaEncodeHeader
from sonLib.bioio import fastaDecodeHeader
from sonLib.bioio import cigarRead
from sonLib.bioio import cigarWrite
from sonLib.bioio import system
from workflow.jobTree.scriptTree.target import Target

from pecan2.pecan2_batch import pecan2BatchWrapperTopLevel

def cactusAlignerTestAligner(job, sequenceFiles, resultsFile):
    """Function to create a cactus_alignerTestAligner.TestAligner target.
    """
    from cactus.cactus_alignerTestAligner import MakeBlasts
    target = MakeBlasts(job, sequenceFiles, resultsFile)
    logger.info("Constructed the cactus_alignerTestAligner target")
    return target

class MakeSequencesOptions:
    """A class to contain options which can be passed to the MakeSequences target.
    """
    def __init__(self, makeBlastTarget=pecan2BatchWrapperTopLevel):
        self.makeBlastTarget = makeBlastTarget
    
class MakeSequences(Target):
    """Take a reconstruction problem and generate the sequences to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, job, netDisk, 
                 netName, resultsFile, options):
        self.netDisk = netDisk
        self.netName = netName
        self.resultsFile = resultsFile
        self.options = options
        Target.__init__(self, job, None)
        
    def run(self, job):
        ##########################################
        #Setup the temp files
        ##########################################
        
        tempSeqFile = getTempFile(rootDir=job.attrib["global_temp_dir"])
        tempResultsFile = getTempFile(rootDir=job.attrib["global_temp_dir"])
        
        logger.info("Built temporary files")
        
        ##########################################
        #Construct the sequences file for doing all against all blast.
        ##########################################
        
        system("cactus_aligner %s %s %s" % (self.netDisk, self.netName, tempSeqFile))
        
        logger.info("Got the sequence files to align")
        
        ##########################################
        #Make blast target
        ##########################################
        
        self.addChildTarget(self.options.makeBlastTarget(job, [ tempSeqFile ], tempResultsFile))
        logger.info("Added child target okay")
        
        ##########################################
        #Setup follow on coordinates
        ##########################################
        
        ModifyBlasts(job, self, tempSeqFile, tempResultsFile, self.resultsFile)
        logger.info("Created modify blasts target")
    
class ModifyBlasts(Target):
    """Modifies the alignments file so that the sequences have the correct coordinates.
    """
    
    def __init__(self, job, previousTarget, tempSeqFile, tempResultsFile, resultsFile):
        self.tempSeqFile = tempSeqFile
        self.tempResultsFile = tempResultsFile
        self.resultsFile = resultsFile
        Target.__init__(self, job, previousTarget)      
        
    def cleanup(self, job):
        os.remove(self.tempSeqFile)
        logger.info("Removed the temporary fasta file for the blast step")
    
    def run(self, job):
        ##########################################
        #Translate the coordinates
        ##########################################
        
        fileHandle = open(self.tempResultsFile, 'r')
        fileHandle2 = open(self.resultsFile, 'w')
        
        for pairwiseAlignment in cigarRead(fileHandle):
            #Adjust contig1 coordinates
            attributes = fastaDecodeHeader(pairwiseAlignment.contig1)
            start = int(attributes[-1])
            pairwiseAlignment.contig1 = fastaEncodeHeader(attributes[:-1])
            pairwiseAlignment.start1 += start
            pairwiseAlignment.end1 += start
            #Adjust contig2 coordinates
            attributes = fastaDecodeHeader(pairwiseAlignment.contig2)
            start = int(attributes[-1])
            pairwiseAlignment.contig2 = fastaEncodeHeader(attributes[:-1])
            pairwiseAlignment.start2 += start
            pairwiseAlignment.end2 += start
            
            cigarWrite(fileHandle2, pairwiseAlignment)
        
        fileHandle.close()
        fileHandle2.close()
        
        logger.info("Processed alignments okay")
        
        ##########################################
        #Make final cleanup target
        ##########################################
        
        ##
        #We allow jobTree to clean this up, as it doesn't seem worth creating a target just for this small amount of data.
        ##
        #CleanupTarget(job, self, [ self.tempResultsFile ])

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    options = MakeSequencesOptions()
    
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
    
    if parsedOptions.useDummy:
        options.makeBlastTarget = cactusAlignerTestAligner
    
    job = ET.parse(parsedOptions.jobFile).getroot()
    
    firstTarget = MakeSequences(job, parsedOptions.netDisk, 
                                parsedOptions.netName, 
                                parsedOptions.resultsFile, options)
    firstTarget.execute(parsedOptions.jobFile)
    
    logger.info("Ran the first target okay")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()