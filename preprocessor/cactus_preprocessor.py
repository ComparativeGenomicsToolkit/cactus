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
import copy

from sonLib.bioio import TempFileTree
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getLogLevelString
from sonLib.bioio import newickTreeParser


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
    """ return a list of files in the directory, or just the
        directory name if its empty
    """
    if os.path.isdir(path):
        contents = os.listdir(path)
        files = []
        for i in contents:
            if i[0] != '.':
                fpath = os.path.join(path, i)
                if os.path.isfile(fpath):
                    files.append(fpath)
        if len(files) > 0:
            return files
    return [path]

class PreprocessorHelper:
    def __init__(self, cactusWorkflowArguments, sequences):     
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.sequences = sequences
        self.fileEventMap = self.__computeFileEventMap()
    
    def getFilteredXmlElems(self, sequence):
        prepNodes = self.cactusWorkflowArguments.configNode.findall("preprocessor")
        filteredNodes = []
        event = self.fileEventMap[sequence]
        leafEvents = getattr(self.cactusWorkflowArguments, 'globalLeafEventSet', set([event]))         
        for node in prepNodes:
            scope = node.get("scope", default="leaves").lower()
            if event in leafEvents:
                if scope != 'internal':
                    filteredNodes.append(node)
            elif scope != 'leaves':
                filteredNodes.append(node)
        return filteredNodes
    
    # link each fasta file to an event name and store 
    # relies on sequences aways being read in same order
    def __computeFileEventMap(self):
        seqIterator = iter(self.sequences)
        eventMap = dict()
        tree = newickTreeParser(self.cactusWorkflowArguments.speciesTree)
        dfStack = [tree]
        while dfStack:
            node = dfStack.pop(-1)
            if node is not None and not node.internal:
                eventMap[seqIterator.next()] = node.iD 
            else:
                dfStack.append(node.right)
                dfStack.append(node.left) 
        assert len(eventMap) == len(self.sequences)
        fileEventMap = dict()
        for seq in self.sequences:
            event = eventMap[seq]            
            for file in fileList(seq):
                assert file not in fileEventMap
                fileEventMap[file] = event
            fileEventMap[seq] = event
        return fileEventMap
            

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
    def __init__(self, prepOptions, seqPath, chunkList, event):
        Target.__init__(self)
        self.prepOptions = prepOptions 
        self.seqPath = seqPath
        self.chunkList = chunkList
        self.event = event
    
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
            cmdline = cmdline.replace("EVENT_STRING", self.event)
            
            logger.info("Preprocessor exec " + cmdline)
            system(cmdline)
           
            compressedChunk = compressFastaFile(prepChunkPath, self.getLocalTempDir(),
                                                self.prepOptions.compressFiles)
            
            system("mv %s %s" % (compressedChunk, chunk))

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

        system("cactus_batch_mergeChunks %s %s %i" % \
               (self.chunkListPath, self.outSequencePath, self.prepOptions.compressFiles))
        
 
class PreprocessSequence(Target):
    """ cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequencePath, outSequencePath, event):
        Target.__init__(self)
        self.prepOptions = prepOptions 
        self.inSequencePath = inSequencePath
        self.outSequencePath = outSequencePath
        self.event = event
    
    # Chunk the input path (inSequencePath).  return a path containing a list of
    # input sequence chunks (chunkListPath)
    def makeChunkList(self):
        chunkListPath = getTempFile(suffix=".txt", rootDir=self.getGlobalTempDir())
        inSeqListPath = getTempFile(suffix=".txt", rootDir=self.getLocalTempDir())
        inSeqListHandle = open(inSeqListPath, "w")
        inSeqListHandle.write(self.inSequencePath + "\n")
        inSeqListHandle.close()
        chunkDirectory = getTempDirectory(self.getGlobalTempDir())
        
        system("cactus_batch_chunkSequences %s %i %i %s %i %s" % \
               (chunkListPath, self.prepOptions.chunkSize, self.prepOptions.overlapSize,
                chunkDirectory, self.prepOptions.compressFiles, inSeqListPath))   
        
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
            self.addChildTarget(PreprocessChunks(self.prepOptions, seqPath, chunkSubList, self.event))

        # follow on to merge chunks
        self.setFollowOnTarget(MergeChunks(self.prepOptions, chunkListPath, self.outSequencePath))

class BatchPreprocessor(Target):
    def __init__(self, cactusWorkflowArguments, event, prepXmlElems, inSequence, globalOutSequence, iteration = 0):
        Target.__init__(self, time=0.0002)
        self.cactusWorkflowArguments = cactusWorkflowArguments 
        self.event = event
        self.prepXmlElems = prepXmlElems
        self.inSequence = inSequence
        self.globalOutSequence = globalOutSequence
        self.iteration = iteration
              
    def run(self):
        # Parse the "preprocessor" config xml element     
        assert self.iteration < len(self.prepXmlElems)
        
        prepNode = self.prepXmlElems[self.iteration]
    
        prepOptions = PreprocessorOptions(int(prepNode.get("chunkSize", default="2147483647")),
                                          int(prepNode.get("chunksPerJob", default="1")),
                                          int(prepNode.get("overlapSize", default="10")),
                                          prepNode.get("compressFiles", default="True").lower() == "true",
                                          prepNode.attrib["preprocessorString"])
        
        #output to temporary directory unless we are on the last iteration
        lastIteration = self.iteration == len(self.prepXmlElems) - 1        
        if lastIteration == False:
            outDir = getTempDirectory(self.getGlobalTempDir())
            outSeq = os.path.join(outDir, os.path.split(self.inSequence)[1])
        else:
            outSeq = self.globalOutSequence      
            outDir = self.globalOutSequence                    
            if not os.path.isdir(self.inSequence):
                outDir = os.path.split(outDir)[0]
        
        #iterate over each input fasta file
        inSeqFiles = fileList(self.inSequence)
        outSeqFiles = []
        for inSeqFile in inSeqFiles:
            fileName = os.path.split(inSeqFile)[1]
            outSeqFiles.append(os.path.join(outDir, fileName))      
        outSeqFile = iter(outSeqFiles)
        assert len(outSeqFiles) == len(inSeqFiles)
        
        # either process a sequence, or propagate an empty directory
        for inSeqFile in inSeqFiles:
            if not os.path.isdir(inSeqFile):
                self.addChildTarget(PreprocessSequence(prepOptions, inSeqFile, outSeqFile.next(), self.event))
            else:
                os.makedirs(outSeqFile.next())

        if lastIteration == False:
            self.setFollowOnTarget(BatchPreprocessor(self.cactusWorkflowArguments, self.event, self.prepXmlElems, outSeq,
                                                     self.globalOutSequence, self.iteration + 1))
            
