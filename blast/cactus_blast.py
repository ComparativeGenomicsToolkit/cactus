#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
import sys
from optparse import OptionParser
from sonLib.bioio import TempFileTree
from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch
from sonLib.bioio import getLogLevelString
from sonLib.bioio import makeSubDir
from sonLib.bioio import catFiles
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

class BlastOptions:
    def __init__(self, chunkSize=10000000, overlapSize=10000, 
                 lastzArguments="", compressFiles=True, realign=False, realignArguments="",
                 minimumSequenceLength=1, memory=sys.maxint):
        """Class defining options for blast
        """
        self.chunkSize = chunkSize
        self.overlapSize = overlapSize
        
        if realign:
            self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] | cactus_realign %s SEQ_FILE_1 SEQ_FILE_2 > CIGARS_FILE"  % (lastzArguments, realignArguments) 
        else:
            self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE"  % lastzArguments 
        if realign:
            self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial | cactus_realign %s SEQ_FILE > CIGARS_FILE" % (lastzArguments, realignArguments)
        else:
            self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial > CIGARS_FILE" % lastzArguments
        self.compressFiles = compressFiles
        self.minimumSequenceLength = 10
        self.memory = memory
    
class BlastFlower(Target):
    """Take a reconstruction problem and generate the sequences in chunks to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, flowerName, finalResultsFile, blastOptions):
        Target.__init__(self)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 2
        
    def run(self):
        chunksDir = makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks"))
        chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkFlowerSequences %s '%s' %s %i %i %i %s" % \
                                                          (getLogLevelString(), self.cactusDisk, self.flowerName, 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          self.blastOptions.minimumSequenceLength,
                                                          chunksDir)).split("\n") if chunk != "" ]
        logger.info("Broken up the flowers into individual 'chunk' files")
        self.addChildTarget(MakeBlastsAllAgainstAll(self.blastOptions, chunks, self.finalResultsFile))
        
class BlastSequencesAllAgainstAll(Target):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFiles1, finalResultsFile, blastOptions):
        Target.__init__(self)
        self.sequenceFiles1 = sequenceFiles1
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 1
    
    def getChunks(self, sequenceFiles, chunksDir):
        return [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                                          (getLogLevelString(), 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          chunksDir,
                                                          " ".join(sequenceFiles))).split("\n") if chunk != "" ]
        
    def run(self):
        chunks = self.getChunks(self.sequenceFiles1, makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks")))
        logger.info("Broken up the sequence files into individual 'chunk' files")
        self.addChildTarget(MakeBlastsAllAgainstAll(self.blastOptions, chunks, self.finalResultsFile))
        
class MakeBlastsAllAgainstAll(Target):
    """Breaks up the inputs into bits and builds a bunch of alignment jobs.
    """
    def __init__(self, blastOptions, chunks, finalResultsFile):
        Target.__init__(self)
        self.blastOptions = blastOptions
        self.chunks = chunks
        self.finalResultsFile = finalResultsFile
        
    def run(self):
        #Avoid compression if just one chunk
        self.blastOptions.compressFiles = self.blastOptions.compressFiles and len(self.chunks) > 2
        selfResultsDir = makeSubDir(os.path.join(self.getGlobalTempDir(), "selfResults"))
        resultsFiles = []
        for i in xrange(len(self.chunks)):
            resultsFile = os.path.join(selfResultsDir, str(i))
            resultsFiles.append(resultsFile)
            self.addChildTarget(RunSelfBlast(self.blastOptions, self.chunks[i], resultsFile))
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        self.setFollowOnTarget(MakeBlastsAllAgainstAll2(self.blastOptions, self.chunks, resultsFiles, self.finalResultsFile))
    
class MakeBlastsAllAgainstAll2(MakeBlastsAllAgainstAll):
        def __init__(self, blastOptions, chunks, resultsFiles, finalResultsFile):
            MakeBlastsAllAgainstAll.__init__(self, blastOptions, chunks, finalResultsFile)
            self.resultsFiles = resultsFiles
           
        def run(self):
            tempFileTree = TempFileTree(os.path.join(self.getGlobalTempDir(), "allAgainstAllResults"))
            #Make the list of blast jobs.
            for i in xrange(0, len(self.chunks)):
                for j in xrange(i+1, len(self.chunks)):
                    resultsFile = tempFileTree.getTempFile()
                    self.resultsFiles.append(resultsFile)
                    self.addChildTarget(RunBlast(self.blastOptions, self.chunks[i], self.chunks[j], resultsFile))
            logger.info("Made the list of all-against-all blasts")
            #Set up the job to collate all the results
            self.setFollowOnTarget(CollateBlasts(self.finalResultsFile, self.resultsFiles))
            
class BlastSequencesAgainstEachOther(BlastSequencesAllAgainstAll):
    """Take two sets of sequences, chunks them up and blasts one set against the other.
    """
    def __init__(self, sequenceFiles1, sequenceFiles2, finalResultsFile, blastOptions):
        BlastSequencesAllAgainstAll.__init__(self, sequenceFiles1, finalResultsFile, blastOptions)
        self.sequenceFiles2 = sequenceFiles2
        
    def run(self):
        chunks1 = self.getChunks(self.sequenceFiles1, makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks1")))
        chunks2 = self.getChunks(self.sequenceFiles2, makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks2")))
        tempFileTree = TempFileTree(os.path.join(self.getGlobalTempDir(), "allAgainstAllResults"))
        resultsFiles = []
        #Make the list of blast jobs.
        for chunk1 in chunks1:
            for chunk2 in chunks2:
                resultsFile = tempFileTree.getTempFile()
                resultsFiles.append(resultsFile)
                #TODO: Make the compression work
                self.blastOptions.compressFiles = False
                self.addChildTarget(RunBlast(self.blastOptions, chunk1, chunk2, resultsFile))
        logger.info("Made the list of blasts")
        #Set up the job to collate all the results
        self.setFollowOnTarget(CollateBlasts(self.finalResultsFile, resultsFiles))
        
def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    system("bzip2 --keep --fast %s" % fileName)
        
class RunSelfBlast(Target):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFile, resultsFile):
        Target.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFile = seqFile
        self.resultsFile = resultsFile
    
    def run(self):   
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
        command = self.blastOptions.selfBlastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE", self.seqFile)
        system(command)
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, self.resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        if self.blastOptions.compressFiles:
            compressFastaFile(self.seqFile)
        logger.info("Ran the self blast okay")

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    system("bunzip2 --stdout %s > %s" % (fileName, tempFileName))
    return tempFileName
    
class RunBlast(Target):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFile1, seqFile2, resultsFile):
        Target.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFile1 = seqFile1
        self.seqFile2 = seqFile2
        self.resultsFile = resultsFile
    
    def run(self):
        if self.blastOptions.compressFiles:
            self.seqFile1 = decompressFastaFile(self.seqFile1 + ".bz2", os.path.join(self.getLocalTempDir(), "1.fa"))
            self.seqFile2 = decompressFastaFile(self.seqFile2 + ".bz2", os.path.join(self.getLocalTempDir(), "2.fa"))
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
        command = self.blastOptions.blastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE_1", self.seqFile1).replace("SEQ_FILE_2", self.seqFile2)
        system(command)
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, self.resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        logger.info("Ran the blast okay")
    
class CollateBlasts(Target):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, finalResultsFile, resultsFiles):
        Target.__init__(self)
        self.finalResultsFile = finalResultsFile
        self.resultsFiles = resultsFiles
    
    def run(self):
        catFiles(self.resultsFiles, self.finalResultsFile)
        logger.info("Collated the alignments to the file: %s",  self.finalResultsFile)
        
class SortCigarAlignmentsInPlace(Target):
    """Sorts an alignment file in place.
    """
    def __init__(self, cigarFile):
        Target.__init__(self)
        self.cigarFile = cigarFile
    
    def run(self):
        #Need to add sorting by start coordinate for coverage calculation
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
        system("cactus_blast_sortAlignments %s %s %i" % (getLogLevelString(), self.cigarFile, tempResultsFile))
        logger.info("Sorted the alignments okay")
        system("mv %s %s" % (tempResultsFile, self.cigarFile))

class CalculateAlignmentCoverage(Target):
    """Calculates the coverage of an genome according to an alignment.
    """
    def __init__(self, sortedCigarFile, sequenceFiles, outputCoverageFile):
        Target.__init__(self)
        self.sortedCigarFile = sortedCigarFile
        self.sequenceFiles = sequenceFiles
        self.outputCoverageFile = outputCoverageFile
    
    def run(self):
        #TODO - I'd suggest the coverage file be in a standard format, i.e. BED?
        #It should report for each interval of the reference the number of distinct places the interval is aligned.
        #I'd probably do this by looking at gapless alignments - i.e. the "matched" segments of the alignment.
        #The input should first be sorted by sequence position, so this can be done by iterating over the alignments in the order of the sequence.
        
        pass
    
class TrimGenome(Target):
    """Trim genome according to coverage file.
    """
    def __init__(self, sequenceFiles, coverageFile, outputSequenceFile, parameters):
        Target.__init__(self)
        self.sequenceFiles = sequenceFiles
        self.coverageFile = coverageFile
        self.outputSequenceFile = outputSequenceFile
        self.parameters = parameters
    
    def run(self):
        #TODO - I'd suggest the parameters be:
        #Coverage threshold - the amount of coverage to decide to include or exclude a subsequence (integer)
        #Exclude/include parameter- flag which decides if we include region with a given or higher coverage threshold or exclude them.
        #Window length - the length of the window to integrate the coverage over, this will act as both a smoothing parameter and ensure we don't report tiny regions).
        #Avg/median parameter - allow use of median or avg. coverage.
        pass

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    blastOptions = BlastOptions()
    
    #output stuff
    parser.add_option("--cigars", dest="cigarFile", 
                      help="File to write cigars in",
                      default="cigarFile.txt")
    
    parser.add_option("--chunkSize", dest="chunkSize", type="int", 
                     help="The size of chunks passed to lastz (must be at least twice as big as overlap)",
                     default=blastOptions.chunkSize)
    
    parser.add_option("--overlapSize", dest="overlapSize", type="int",
                     help="The size of the overlap between the chunks passed to lastz (min size 2)",
                     default=blastOptions.overlapSize)
    
    parser.add_option("--blastString", dest="blastString", type="string",
                     help="The default string used to call the blast program. \
Must contain three strings: SEQ_FILE_1, SEQ_FILE_2 and CIGARS_FILE which will be \
replaced with the two sequence files and the results file, respectively",
                     default=blastOptions.blastString)
    
    parser.add_option("--selfBlastString", dest="selfBlastString", type="string",
                     help="The default string used to call the blast program for self alignment. \
Must contain three strings: SEQ_FILE and CIGARS_FILE which will be \
replaced with the the sequence file and the results file, respectively",
                     default=blastOptions.selfBlastString)
   
    parser.add_option("--compressFiles", dest="compressFiles", action="store_false",
                      help="Turn of bz2 based file compression of sequences for I/O transfer", 
                      default=blastOptions.compressFiles)
    
    parser.add_option("--lastzMemory", dest="memory", type="int",
                      help="Lastz memory (in bytes)", 
                      default=blastOptions.memory)
    
    parser.add_option("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")
    
    parser.add_option("--targetSequenceFiles", dest="targetSequenceFiles", type="string",
                     help="Sequences to align against the input sequences against. If these are not provided then the input sequences are aligned against each other.",
                     default=None)
    
    options, args = parser.parse_args()
    if options.test:
        _test()
    
    if options.targetSequenceFiles == None:
        firstTarget = BlastSequencesAllAgainstAll(args, options.cigarFile, options)
    else:
        firstTarget = BlastSequencesAgainstEachOther(args, options.targetSequenceFiles.split(), options.cigarFile, options)
    Stack(firstTarget).startJobTree(options)

def _test():
    import doctest
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blast.cactus_blast import *
    main()
