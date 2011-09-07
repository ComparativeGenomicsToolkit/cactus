#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
import errno
from optparse import OptionParser
from bz2 import BZ2File

from sonLib.bioio import TempFileTree
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getLogLevelString

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from cactus.blastAlignment.cactus_batch import decompressFastaFile

def compressFastaFile(fileName, tempDir, compressFiles):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    if compressFiles:
        tempFileName = getTempFile(suffix=".fa.bz2", rootDir=tempDir)
        fileHandle = BZ2File(tempFileName, 'w')
        fileHandle2 = open(fileName, 'r')
        for fastaHeader, seq in fastaRead(fileHandle2):
            fastaWrite(fileHandle, fastaHeader, seq)
        fileHandle2.close()
        fileHandle.close()
        return tempFileName
    return fileName

def fileList(path):
    """ return a list of files in the directory
    """
    if os.path.isdir(path):
        contents = os.listdir(path)
        files = []
        for i in contents:
            if i[0] != '.':
                fpath = os.path.join(path, i)
                if os.path.isfile(fpath):
                    files.append(fpath)
        return files
    else:
        return [path]

class PreprocessorOptions:
    def __init__(self, chunkSize, chunksPerJob, overlapSize, compressFiles, cmdLine):
        self.chunkSize = chunkSize
        self.chunksPerJob = chunksPerJob
        self.overlapSize = overlapSize
        self.compressFiles = compressFiles
        self.cmdLine = cmdLine

class PreprocessChunks(Target):
    """ locally preprocess some fasta chunks, output then copied back to input
    """
    def __init__(self, prepOptions, seqPath, chunkList):
        Target.__init__(self)
        self.prepOptions = prepOptions 
        self.seqPath = seqPath
        self.chunkList = chunkList
    
    def run(self):
        #Full sequence file can be quite big: only decompress it into the local
        #path if it's actually used by the preprocessor
        localSequencePath = ""
        if self.prepOptions.cmdLine.find("QUERY_FILE") >= 0:            
            localSequencePath = decompressFastaFile(self.seqPath,
                                                    self.getLocalTempDir(), 
                                                    self.prepOptions.compressFiles)
        for chunk in self.chunkList:
            localChunkPath = decompressFastaFile(chunk, self.getLocalTempDir(),
                                                 self.prepOptions.compressFiles)
            prepChunkPath = getTempFile(rootDir=self.getLocalTempDir())
            tempPath = getTempFile(rootDir=self.getLocalTempDir())
            
            cmdline = self.prepOptions.cmdLine.replace("QUERY_FILE", "\"" + localSequencePath + "\"")
            cmdline = cmdline.replace("TARGET_FILE", "\"" + localChunkPath + "\"")
            cmdline = cmdline.replace("OUT_FILE", "\"" + prepChunkPath + "\"")
            cmdline = cmdline.replace("TEMP_FILE", "\"" + tempPath + "\"")
            
            logger.info("Preprocessor exec " + cmdline)
            assert os.system(cmdline) == 0
           
            compressedChunk = compressFastaFile(prepChunkPath, self.getLocalTempDir(),
                                                self.prepOptions.compressFiles)
            
            assert os.system("mv %s %s" % (compressedChunk, chunk)) == 0

class MergeChunks(Target):
    """ merge a list of chunks into a fasta file
    """
    def __init__(self, prepOptions, chunkListPath, outSequencePath):
        Target.__init__(self)
        self.prepOptions = prepOptions 
        self.chunkListPath = chunkListPath
        self.outSequencePath = outSequencePath
    
    def run(self):
        baseDir = os.path.dirname(self.outSequencePath)
        
        #somewhat threadsafe 
        try:
            os.makedirs(baseDir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise e

        sysRet = system("cactus_batch_mergeChunks %s %s %i" % \
                        (self.chunkListPath, self.outSequencePath, self.prepOptions.compressFiles))
        
        assert sysRet == 0
 
class PreprocessSequence(Target):
    """ cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequencePath, outSequencePath):
        Target.__init__(self)
        self.prepOptions = prepOptions 
        self.inSequencePath = inSequencePath
        self.outSequencePath = outSequencePath
    
    # Chunk the input path (inSequencePath).  return a path containing a list of
    # input sequence chunks (chunkListPath)
    def makeChunkList(self):
        chunkListPath = getTempFile(suffix=".txt", rootDir=self.getGlobalTempDir())
        inSeqListPath = getTempFile(suffix=".txt", rootDir=self.getLocalTempDir())
        inSeqListHandle = open(inSeqListPath, "w")
        inSeqListHandle.write(self.inSequencePath + "\n")
        inSeqListHandle.close()
        chunkDirectory = getTempDirectory(self.getGlobalTempDir())
        
        sysRet = system("cactus_batch_chunkSequences %s %i %i %s %i %s" % \
                        (chunkListPath, self.prepOptions.chunkSize, self.prepOptions.overlapSize,
                         chunkDirectory, self.prepOptions.compressFiles, inSeqListPath))   
        assert sysRet == 0
        
        return chunkListPath
    
    def run(self):
        #make compressed version of the sequence
        #Full sequence file can be quite big: only compress it into the local
        #path if it's actually used by the preprocessor
        seqPath = ""
        if self.prepOptions.cmdLine.find("QUERY_FILE") >= 0:        
            seqPath = compressFastaFile(self.inSequencePath, self.getGlobalTempDir(), 
                                    self.prepOptions.compressFiles)
        
        logger.info("Preparing sequence for preprocessing")
        # chunk it up
        chunkListPath = self.makeChunkList()
        # read each line of chunk file, stripping the stupid \n characters
        chunkList = map(lambda x: x.split()[0], open(chunkListPath, "r").readlines())
     
        # for every chunksPerJob chunks in list
        for i in range(0, len(chunkList), self.prepOptions.chunksPerJob):
            chunkSubList = chunkList[i : i + self.prepOptions.chunksPerJob]
            self.addChildTarget(PreprocessChunks(self.prepOptions, seqPath, chunkSubList))

        # follow on to merge chunks
        self.setFollowOnTarget(MergeChunks(self.prepOptions, chunkListPath, self.outSequencePath))

class BatchPreprocessor(Target):
    def __init__(self, options, globalInSequences, inSequences, outDirBase, iteration = 0):
        Target.__init__(self, time=0.0002)
        self.options = options 
        self.globalInSequences = globalInSequences
        self.inSequences = inSequences
        self.outDirBase = outDirBase
        self.prepOptions = None
        self.iteration = iteration
       
    def run(self):
        # Parse the "preprocessor" config xml element
        prepNodes = self.options.config.findall("preprocessor")
        
        assert self.iteration < len(prepNodes)
        
        prepNode = prepNodes[self.iteration]
    
        self.prepOptions = PreprocessorOptions(int(prepNode.get("chunkSize", default="2147483647")),
                                          int(prepNode.get("chunksPerJob", default="1")),
                                          int(prepNode.get("overlapSize", default="10")),
                                          prepNode.get("compressFiles", default="True").lower() == "true",
                                          prepNode.attrib["preprocessorString"])
        
        #iterate over each input fasta file
        inSeqFiles = []
        map(inSeqFiles.extend, map(fileList, self.inSequences))
        assert len(inSeqFiles) >= len(self.inSequences)
        
        globalInSeqFiles = []
        map(globalInSeqFiles.extend, map(fileList, self.globalInSequences))
        assert len(globalInSeqFiles) == len(inSeqFiles)
        
        #output to temporary directory unless we are on the last iteration
        lastIteration = self.iteration == len(prepNodes) - 1
        outDir = self.outDirBase
        if lastIteration == False:
            outDir = getTempDirectory(self.getGlobalTempDir())
            
        outSeqFiles = map(lambda x: outDir  + "/" + x, globalInSeqFiles)
        outSeq = iter(outSeqFiles)
        assert len(outSeqFiles) == len(inSeqFiles)
        
        for inSeq in inSeqFiles:
            self.addChildTarget(PreprocessSequence(self.prepOptions, inSeq, outSeq.next()))
        
        if lastIteration == False:
            self.setFollowOnTarget(BatchPreprocessor(self.options, self.globalInSequences, outSeqFiles, 
                                                     self.outDirBase, self.iteration + 1))
            