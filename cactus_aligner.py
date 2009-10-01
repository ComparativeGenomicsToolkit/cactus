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
from workflow.jobTree.scriptTree.target import Target
from workflow.jobTree.scriptTree.target import CleanupTarget

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
    def __init__(self, job, absolutePathPrefix, 
                 reconstructionProblem, resultsFile, options):
        self.absolutePathPrefix = absolutePathPrefix
        self.reconstructionProblem = reconstructionProblem
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
        #Read the top level reconstruction tags and the file containing the reconstruction problem node
        ##########################################
        
        reconstructionProblemTag = ET.parse(os.path.join(self.absolutePathPrefix, self.reconstructionProblem)).getroot()
        
        logger.info("Parsed the input reconstruction problem")
        
        ##########################################
        #Construct the sequences file for doing all against all blast.
        ##########################################
        
        strings = reconstructionProblemTag.find("strings")
        sequences = reconstructionProblemTag.find("sequences")
        sequencesMap = {}
        for sequence in sequences.findall("sequence"):
            sequencesMap[sequence.attrib["contig"]] = sequence
        
        fileHandle = open(tempSeqFile, 'w')
        
        for string in strings.findall("string"):
            contig = string.attrib["contig"]
            assert contig != None
            start = int(string.attrib["start"])
            length = int(string.attrib["length"])
            sequenceTag = sequencesMap[contig]
            
            fileHandle.write(">%s\n" % fastaEncodeHeader([ contig, start ])) #Do not use the checks of the fastaWrite function
            fileHandle2 = open(os.path.join(self.absolutePathPrefix, sequenceTag.attrib["sequence_file"]), 'r')
            fileHandle2.seek(start)
            seq = fileHandle2.read(length)
            fileHandle2.close()
            assert len(seq) == length
            fileHandle.write(seq)
            seq = ""
            fileHandle.write("\n")
        
        fileHandle.close()
        
        logger.info("Written")
        
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
    
    parser.add_option("--absolutePathPrefix", dest="absolutePathPrefix", 
                      help="The path to the root of the reconstruction tree problem")
    
    parser.add_option("--reconstructionProblem", dest="reconstructionProblem", 
                      help="The file containing the reconstruction problem")
    
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
    
    firstTarget = MakeSequences(job, parsedOptions.absolutePathPrefix, 
                                parsedOptions.reconstructionProblem, 
                                parsedOptions.resultsFile, options)
    firstTarget.execute(parsedOptions.jobFile)
    
    logger.info("Ran the first target okay")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()