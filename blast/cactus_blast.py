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
from bz2 import BZ2File
from sonLib.bioio import TempFileTree
from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch
from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getLogLevelString
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

class BlastOptions:
    def __init__(self, chunkSize=10000000, overlapSize=10000, 
                 lastzArguments="", compressFiles=True, 
                 minimumSequenceLength=1, memory=sys.maxint):
        """Class defining options for blast
        """
        self.chunkSize = chunkSize
        self.overlapSize = overlapSize
        self.blastString = "lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE"  % lastzArguments 
        self.selfBlastString = "lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[multiple][nameparse=darkspace] --notrivial > CIGARS_FILE" % lastzArguments
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
        
    def run(self):
        chunksDir = os.path.join(self.getGlobalTempDir(), "chunks")
        if not os.path.exists(chunksDir):
            os.mkdir(chunksDir)
        chunks = [ line.split()[0] for line in popenCatch("cactus_blast_chunkFlowerSequences %s '%s' %s %i %i %s" % \
                                                          (getLogLevelString(), self.cactusDisk, self.flowerName, 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          self.blastOptions.minimumSequenceLength,
                                                          chunksDir)) ]
        logger.info("Broken up the flowers into individual 'chunk' files")
        self.addChildTarget(MakeBlasts(self.blastOptions, chunks, self.finalResultsFile))
        
class BlastSequences(Target):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFiles, finalResultsFile, blastOptions):
        Target.__init__(self)
        self.sequenceFiles = sequenceFiles
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
        
    def run(self):
        chunksDir = os.path.join(self.getGlobalTempDir(), "chunks")
        if not os.path.exists(chunksDir):
            os.mkdir(chunksDir)
        chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                                          (getLogLevelString(), 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          chunksDir,
                                                          " ".join(self.sequenceFiles))).split("\n") if chunk != "" ]
        logger.info("Broken up the sequence files into individual 'chunk' files")
        self.addChildTarget(MakeBlasts(self.blastOptions, chunks, self.finalResultsFile))
        
class MakeBlasts(Target):
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
        selfResultsDir = os.path.join(self.getGlobalTempDir(), "selfResults")
        if not os.path.exists(selfResultsDir):
            os.mkdir(selfResultsDir)
        resultsFiles = []
        for i in xrange(len(self.chunks)):
            resultsFile = os.path.join(selfResultsDir, str(i))
            resultsFiles.append(resultsFile)
            self.addChildTarget(RunSelfBlast(self.blastOptions, self.chunks[i], resultsFile))
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        self.setFollowOnTarget(MakeBlasts2(self.blastOptions, self.chunks, resultFiles, self.finalResultFile))
    
class MakeBlasts2(MakeBlasts):
        def __init__(self, blastOptions, chunks, resultsFiles, finalResultsFile):
            MakeBlasts.__init__(self, blastOptions, chunks, finalResultsFile)
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
        
def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    fileHandleOut = BZ2File(fileName + ".bz2", 'w')
    fileHandleIn = open(fileName, 'r')
    for fastaHeader, seq in fastaRead(fileHandleIn):
        fastaWrite(fileHandleOut, fastaHeader, seq)
    fileHandleIn.close()
    fileHandleOut.close()
        
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
        command = selfBlastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE", self.seqFile)
        system(command)
        system("cactus_blast_convertCoordinates %s %s" % (tempResultsFile, self.resultsFile))
        if self.blastOptions.compressFiles:
            compressFastaFile(self.seqFile)
        logger.info("Ran the self blast okay")

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    fileHandleOut = open(tempFileName, 'w')
    fileHandleIn = BZ2File(fileName, 'r')
    for fastaHeader, seq in fastaRead(fileHandleIn):
        fastaWrite(fileHandleOut, fastaHeader, seq)
    fileHandleIn.close()
    fileHandleOut.close()
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
        command = blastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE_1", self.seqFile1).replace("SEQ_FILE_2", self.seqFile2)
        system(command)
        system("cactus_blast_convertCoordinates %s %s" % (tempResultsFile, self.resultsFile))
        logger.info("Ran the blast okay")

def catFiles(filesToCat, catFile):
    """Cats a bunch of files into one file. Ensures a no more than maxCat files
    are concatenated at each step.
    """
    maxCat = 25
    system("cat %s > %s" % (" ".join(filesToCat[:maxCat]), catFile))
    filesToCat = filesToCat[maxCat:]
    while len(filesToCat) > 0:
        system("cat %s >> %s" % (" ".join(filesToCat[:maxCat]), catFile))
        filesToCat = filesToCat[maxCat:]
    
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
                      default=sys.maxint)
    
    options, args = parser.parse_args()
    
    firstTarget = BlastSequences(args, options.cigarFile, options)
    Stack(firstTarget).startJobTree(options)

def _test():
    import doctest 
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blast.cactus_blast import *
    _test()
    main()
