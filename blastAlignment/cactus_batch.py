#!/usr/bin/env python

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
from bz2 import BZ2File

from sonLib.bioio import TempFileTree
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getBasicOptionParser
from sonLib.bioio import parseBasicOptions
from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite

from workflow.jobTree.scriptTree.target import Target
from workflow.jobTree.scriptTree.stack import Stack

class MakeBlastOptions:
    def __init__(self, chunkSize, overlapSize, 
                 blastString, selfBlastString, chunksPerJob, compressFiles):
        """Method makes options which can be passed to the to the make blasts target.
        """
        self.chunkSize = chunkSize
        self.overlapSize = overlapSize
        self.blastString = blastString
        self.selfBlastString = selfBlastString
        self.chunksPerJob = 1
        self.compressFiles = compressFiles

def makeStandardBlastOptions():
    """Function to create options for a pecan2_batch.MakeBlasts target for middle level 
    alignments (20 MB range)
    """
    chunkSize = 10000000
    overlapSize = 10000
    chunksPerJob = 1
    compressFiles = True
    blastString="lastz --format=cigar SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE"
    selfBlastString="lastz --format=cigar SEQ_FILE[nameparse=darkspace] --self > CIGARS_FILE"
    return MakeBlastOptions(chunkSize=chunkSize, overlapSize=overlapSize,
                                blastString=blastString, selfBlastString=selfBlastString,
                                chunksPerJob=chunksPerJob, compressFiles=compressFiles)

class makeBlastFromOptions:
    def __init__(self, blastOptions):
        self.blastOptions = blastOptions
    
    def makeBlastOptions(self, sequenceFiles, resultsFile):
        return MakeBlasts(self.blastOptions, sequenceFiles, resultsFile)

class MakeBlasts(Target):
    """Breaks up the inputs into bits a builds a bunch of alignment jobs.
    """
    def __init__(self, options, sequences, finalResultsFile):
        Target.__init__(self)
        assert options.chunkSize > options.overlapSize
        assert options.overlapSize >= 2
        assert options.chunksPerJob >= 1
        self.options = options
        self.sequences = sequences
        self.finalResultsFile = finalResultsFile
        
    def run(self, localTempDir, globalTempDir):
        ##########################################
        #Setup the temp file structure.
        ##########################################
        
        tempFileTreeDir = getTempDirectory(globalTempDir)
        tempFileTree = TempFileTree(tempFileTreeDir)
        
        logger.info("Setup the temporary files")
        
        ##########################################
        #Break up the fasta sequences into overlapping chunks.
        ##########################################
        
        #We break up a sequence into a series of chunks
        #Bridging the breaks between chunks (only between adjacency chunks) we create smaller overlapping, call these bridges.
        
        def processSequences(sequenceFiles, tempFilesDir):
            chunkFiles = getTempFile(suffix=".txt", rootDir=localTempDir)
            bridgeFiles = getTempFile(suffix=".txt", rootDir=localTempDir)
            system("cactus_batch_chunkSequences %s %s %i %i %s %i %s" % \
                   (chunkFiles, bridgeFiles, 
                    self.options.chunkSize, self.options.overlapSize,
                    tempFilesDir, self.options.compressFiles, " ".join(sequenceFiles)))
            
            def readSequenceData(sequenceDataFile):
                l = []
                fileHandle = open(sequenceDataFile, 'r')
                line = fileHandle.readline()
                while line != '':
                    seqFile, length = line.split()
                    length = int(length)
                    l.append((length, seqFile))
                    line = fileHandle.readline()
                fileHandle.close()
                return l
            tempSeqFiles = readSequenceData(chunkFiles)
            tempOverlapSeqFiles = readSequenceData(bridgeFiles)
            os.remove(chunkFiles)
            os.remove(bridgeFiles)
            return tempSeqFiles, tempOverlapSeqFiles
        
        tempSeqFilesDir = getTempDirectory(globalTempDir)
        tempSeqFiles, tempOverlapSeqFiles = processSequences(self.sequences, tempSeqFilesDir)
        logger.info("Broken up the sequence files into individual contigs files")
    
        ##########################################
        #Make all against all blast jobs lists for non overlapping files.
        ##########################################
        
        #We align each chunk + bridge against every chunk + bridge.
        #We align each chunk + bridge against itself.
        
        ###Quick Math to work out how to partition jobs for minimum I/O.
        #N is number of jobs.
        #G is size of total input sequences.
        #x and y are the sizes of the query and target sequences.
        #N = G/x * G/y
        #The amount of data to be moved to to the nodes (alpha) is proportional to
        #alpha = N * (x + y)
        #As
        #y = G^2/N*x
        #alpha = Nx + G^2/x
        
        #For a given N, this number is minimal when x and y are equal.
        ##Numbers for 10 genomes....
        
        #G = 10 x 3000 MB (million bases) = 30 GB (billion bases)
        #N = 10000
        
        #x = y = 300MB
        #alpha = 6,000,000 megabases = 6 terabases of data to move
        
        #x = 600MB, y = 150MB
        #alpha = 7,500,000 MB
        
        #x = 1000MB, y = 90 MB
        #alpha = 10,900,000 MB
        
        #N = 5000
        #x = y = 424.3
        #alpha = 4,243,000 MB
        
        #N = 2000
        #x = y = 670.9 MB
        #alpha = 2,683,600 MB
        
        #N = 1000
        #x = y = 948.5 MB
        #alpha = 1,897,366 MB
        
        ##Bzip2 compresses fasta to approx 2bits per base.. so reduce above alpha number by factor of 4 to get bits.
        def getChunks(tempSeqFiles):
            """Gets all seq files until their combined length exceeds the chunk size, or the list is exhausted.
            """
            l = []
            while len(tempSeqFiles) > 0 and sum([ i[0] for i in l ]) < self.options.chunkSize:
                l.append(tempSeqFiles.pop())
            assert sum([ i[0] for i in l ]) <= self.options.chunkSize + self.options.overlapSize
            return [ i[1] for i in l ]
        
        #Get the list of chunks
        l = tempSeqFiles + tempOverlapSeqFiles
        l.reverse() #This line ensures the list of sequences is added biggest to smallest. 
        chunks = []
        while len(l) > 0:
            chunks.append(getChunks(l))
               
        resultsFiles = []
    
        def makeBlastJobs(chunks1, chunks2):
            """Makes blast jobs.
            """
            while len(chunks2) > 0:
                resultsFile = tempFileTree.getTempFile()
                resultsFiles.append(resultsFile)
                l = chunks2[:self.options.chunksPerJob]
                self.addChildTarget(RunBlast(self.options, chunks1[:], l, resultsFile))
                chunks2 = chunks2[self.options.chunksPerJob:]
                
        def makeSelfBlastJobs(seqFiles):
            """Makes self blast job.
            """
            resultsFile = tempFileTree.getTempFile()
            resultsFiles.append(resultsFile)
            self.addChildTarget(RunSelfBlast(self.options, seqFiles, resultsFile))
        
        #Make the list of self blast jobs.
        for chunk in chunks:
            makeSelfBlastJobs(chunk)
        
        #Make the list of blast jobs.
        while len(chunks) > 0:
            l = chunks[:self.options.chunksPerJob]
            chunks = chunks[self.options.chunksPerJob:]
            makeBlastJobs(l, chunks)
            
        logger.info("Made the list of child targets")
            
        ##########################################
        #Make follow on job to collate results
        ##########################################

        self.setFollowOnTarget(CollateBlasts(self.options, self.finalResultsFile, resultsFiles, tempFileTreeDir, tempSeqFilesDir))
        
        logger.info("Made the follow on job")
        
def executeBlast(seqFile1, seqFile2, resultsFile, blastString):
    """Run actual blast command.
    """
    open(resultsFile, 'w').close() #For safety
    command = blastString.replace("CIGARS_FILE", resultsFile).replace("SEQ_FILE_1", seqFile1).replace("SEQ_FILE_2", seqFile2)
    system(command)
    logger.info("Ran a blast command okay")
    
def executeSelfBlast(seqFile, resultsFile, selfBlastString):
    """Run the actual self blast command.
    """
    open(resultsFile, 'w').close() #For safety
    command = selfBlastString.replace("CIGARS_FILE", resultsFile).replace("SEQ_FILE", seqFile)
    system(command)
    logger.info("Ran the self blast command okay")
    
def readFastaTemporaryFiles(fileNames, tempDir, compressFiles):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    tempFileName = getTempFile(suffix=".fa", rootDir=tempDir)
    fileHandle = open(tempFileName, 'w')
    for fileName in fileNames:
        if compressFiles:
            fileHandle2 = BZ2File(fileName, 'r')
        else:
            fileHandle2 = open(fileName, 'r')
        for fastaHeader, seq in fastaRead(fileHandle2):
            fastaWrite(fileHandle, fastaHeader, seq)
        fileHandle2.close()
    fileHandle.close()
    return tempFileName

def catFiles(filesToCat, catFile):
    """Cats a bunch of files into one file. Ensures a no more than MAX_CAT files
    are concatenated at each step.
    """
    MAX_CAT = 25
    system("cat %s > %s" % (" ".join(filesToCat[:MAX_CAT]), catFile))
    filesToCat = filesToCat[MAX_CAT:]
    while len(filesToCat) > 0:
        system("cat %s >> %s" % (" ".join(filesToCat[:MAX_CAT]), catFile))
        filesToCat = filesToCat[MAX_CAT:]

class RunBlast(Target):
    """Runs blast as a job.
    """
    def __init__(self, options, seqFilesChunks1, seqFilesChunks2, resultsFile):
        Target.__init__(self)
        self.options = options
        self.seqFilesChunks1 = seqFilesChunks1
        self.seqFilesChunks2 = seqFilesChunks2
        self.resultsFile = resultsFile
    
    def run(self, localTempDir, globalTempDir):
        self.seqFilesChunks1 = [ readFastaTemporaryFiles(seqFiles, localTempDir, self.options.compressFiles) for seqFiles in self.seqFilesChunks1 ]
        self.seqFilesChunks2 = [ readFastaTemporaryFiles(seqFiles, localTempDir, self.options.compressFiles) for seqFiles in self.seqFilesChunks2 ]
        
        logger.info("Created the temporary sequence files and copied input files to the local temporary directory")
        
        tempResultsFiles = []
        for tempSeqFile1 in self.seqFilesChunks1:
            for tempSeqFile2 in self.seqFilesChunks2:
                tempResultsFile = getTempFile(suffix=".cigars", rootDir=localTempDir)
                tempResultsFiles.append(tempResultsFile)
                logger.info("Got a temporary results file")
                executeBlast(tempSeqFile1, tempSeqFile2, tempResultsFile, self.options.blastString)
                
        #Write stuff back to the central dirs.
        catFiles(tempResultsFiles, self.resultsFile)
        
        logger.info("Copied back the results files")
        
        for tempSeqFile in self.seqFilesChunks1 + self.seqFilesChunks2 + tempResultsFiles:
            os.remove(tempSeqFile)
                
        logger.info("Removed the temporary sequence files")
    
class RunSelfBlast(Target):
    """Runs blast as a job.
    """
    def __init__(self, options, seqFiles, resultsFile):
        Target.__init__(self)
        self.options = options
        self.seqFiles = seqFiles
        self.resultsFile = resultsFile
    
    def run(self, localTempDir, globalTempDir):
        #Get temp file to put sequences and results in locally.
        for i in xrange(len(self.seqFiles)):
            self.seqFiles[i] = readFastaTemporaryFiles(self.seqFiles[i:i+1], localTempDir, self.options.compressFiles)
        
        tempResultsFile1 = getTempFile(suffix=".cigars", rootDir=localTempDir)
        tempResultsFile2 = getTempFile(suffix=".cigars", rootDir=localTempDir)
        tempSeqFile1 = getTempFile(suffix=".fa", rootDir=localTempDir)
        tempSeqFile2 = getTempFile(suffix=".fa", rootDir=localTempDir)
        
        logger.info("Created the temporary files")
        
        def fn(seqFiles):
            """If there are 2 or more seq files, the function divides the input into two and aligns one against the other,
            then recurses on the two subdivisions.
            If there is only on seq file it does a self alignment of the seq file against itself.
            """
            if len(seqFiles) >= 2:
                seqFiles1 = seqFiles[:len(seqFiles)/2]
                seqFiles2 = seqFiles[len(seqFiles)/2:]
                #Do blast job
                catFiles(seqFiles1, tempSeqFile1)
                catFiles(seqFiles2, tempSeqFile2)
                
                executeBlast(tempSeqFile1, tempSeqFile2, tempResultsFile1, self.options.blastString)
                system("cat %s >> %s" % (tempResultsFile1, tempResultsFile2))
                logger.info("Copied input file to the local temporary directory")
                #Now recurse
                fn(seqFiles1)
                fn(seqFiles2)
                logger.info("Recursed on split sequences successfully")
            else:
                assert len(seqFiles) == 1
                #Do self blast
                executeSelfBlast(seqFiles[0], tempResultsFile1, self.options.selfBlastString)
                system("cat %s >> %s" % (tempResultsFile1, tempResultsFile2))
                logger.info("Copied input file to the local temporary directory")
        open(tempResultsFile2, 'w').close() #Defensive, to ensue tempResultsFile2 is empty
        fn(self.seqFiles)
        
        #Write stuff back to the central dirs.
        system("cp %s %s" % (tempResultsFile2, self.resultsFile))
        logger.info("Copied the results back")
            
        #Cleanup
        for tempFile in self.seqFiles + [ tempSeqFile1, tempSeqFile2, tempResultsFile1, tempResultsFile2 ]:
            os.remove(tempFile)
        
        logger.info("Cleaned up")
    
class CollateBlasts(Target):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, options, finalResultsFile, resultsFiles, tempFileTreeDir, tempSeqFilesDir):
        Target.__init__(self)
        self.options = options
        self.finalResultsFile = finalResultsFile
        self.resultsFiles = resultsFiles
        self.tempFileTreeDir = tempFileTreeDir
        self.tempSeqFilesDir = tempSeqFilesDir
    
    def run(self, localTempDir, globalTempDir):
        ##########################################
        #Collate the results.
        ##########################################
        
        system("cactus_batch_convertCoordinates %s %s" % (" ".join(self.resultsFiles), self.finalResultsFile))
        logger.info("Collated the alignments to the file: %s",  self.finalResultsFile)

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    parser = getBasicOptionParser("usage: %prog [options] contigFilexN", "%prog 0.1")
    Stack.addJobTreeOptions(parser)
    options = makeStandardBlastOptions()
    
    #output stuff
    parser.add_option("--cigars", dest="cigarFile", 
                      help="File to write cigars in",
                      default="cigarFile.txt")
    
    parser.add_option("--chunkSize", dest="chunkSize", type="int",
                     help="The size of chunks passed to lastz (must be at least twice as big as overlap)",
                     default=options.chunkSize)
    
    parser.add_option("--overlapSize", dest="overlapSize", type="int",
                     help="The size of the overlap between the chunks passed to lastz (min size 2)",
                     default=options.overlapSize)
    
    parser.add_option("--blastString", dest="blastString", type="string",
                     help="The default string used to call the blast program. \
Must contain three strings: SEQ_FILE_1, SEQ_FILE_2 and CIGARS_FILE which will be \
replaced with the two sequence files and the results file, respectively",
                     default=options.blastString)
    
    parser.add_option("--selfBlastString", dest="selfBlastString", type="string",
                     help="The default string used to call the blast program for self alignment. \
Must contain three strings: SEQ_FILE and CIGARS_FILE which will be \
replaced with the the sequence file and the results file, respectively",
                     default=options.selfBlastString)
    
    parser.add_option("--chunksPerJob", dest="chunksPerJob", type="int",
                      help="The number of blast chunks to align per job. Every chunk is aligned against every other chunk, \
this allows each job to more than one chunk comparison per job, which will save on I/O.", default=options.chunksPerJob)
    
    parser.add_option("--compressFiles", dest="compressFiles", action="store_false",
                      help="Turn of bz2 based file compression of sequences for I/O transfer", 
                      default=options.compressFiles)
    
    options, args = parseBasicOptions(parser)
    logger.info("Parsed arguments")
    
    firstTarget = MakeBlasts(options, args, options.cigarFile)
    Stack(firstTarget).startJobTree(options)
    
    logger.info("Ran the first target okay")

def _test():
    import doctest 
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blastAlignment.cactus_batch import *
    _test()
    main()
