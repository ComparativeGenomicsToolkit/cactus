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

from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch
from sonLib.bioio import getLogLevelString
from sonLib.bioio import newickTreeParser
from sonLib.bioio import makeSubDir
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from cactus.blast.cactus_blast import catFiles

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
            if node.get("preprocessorString", default=None) is not None:
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
    def __init__(self, chunkSize, cmdLine, memory, cpu):
        self.chunkSize = chunkSize
        self.cmdLine = cmdLine
        self.memory = memory
        self.cpu = cpu

class PreprocessChunk(Target):
    """ locally preprocess a fasta chunk, output then copied back to input
    """
    def __init__(self, prepOptions, seqPath, chunk, event):
        Target.__init__(self, memory=prepOptions.memory, cpu=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.seqPath = seqPath
        self.chunk = chunk
        self.event = event
    
    def run(self):
        prepChunkPath = os.path.join(self.getLocalTempDir(), "prepChunk")
        tempPath = os.path.join(self.getLocalTempDir(), "tempPath")
        
        cmdline = self.prepOptions.cmdLine.replace("QUERY_FILE", "\"" + self.seqPath + "\"")
        cmdline = cmdline.replace("TARGET_FILE", "\"" + self.chunk + "\"")
        cmdline = cmdline.replace("OUT_FILE", "\"" + prepChunkPath + "\"")
        cmdline = cmdline.replace("TEMP_FILE", "\"" + tempPath + "\"")
        cmdline = cmdline.replace("EVENT_STRING", self.event)
        
        logger.info("Preprocessor exec " + cmdline)
        system(cmdline)

class MergeChunks(Target):
    """ merge a list of chunks into a fasta file
    """
    def __init__(self, prepOptions, chunkList, outSequencePath):
        Target.__init__(self, cpu=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.chunkList = chunkList
        self.outSequencePath = outSequencePath
    
    def run(self):
        system("cactus_batch_mergeChunks %s %s" % \
               (self.outSequencePath, " ".join(self.chunkList)))
 
class PreprocessSequence(Target):
    """ cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequencePath, outSequencePath, event):
        Target.__init__(self, cpu=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.inSequencePath = inSequencePath
        self.outSequencePath = outSequencePath
        self.event = event
    
    def run(self):        
        logger.info("Preparing sequence for preprocessing")
        # chunk it up
        chunkDirectory = makeSubDir(os.path.join(self.getGlobalTempDir(), "preprocessChunks"))
        chunkList = [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i 0 %s %s" % \
               (getLogLevelString(), self.prepOptions.chunkSize,
                chunkDirectory, self.inSequencePath)).split("\n") if chunk != "" ]   
        for chunk in chunkList:
            self.addChildTarget(PreprocessChunk(self.prepOptions, self.inSequencePath, chunk, self.event))
        # follow on to merge chunks
        self.setFollowOnTarget(MergeChunks(self.prepOptions, chunkList, self.outSequencePath))

class BatchPreprocessor(Target):
    def __init__(self, cactusWorkflowArguments, event, prepXmlElems, inSequence, 
                 globalOutSequence, memory, cpu, iteration = 0):
        Target.__init__(self, time=0.0002)
        self.cactusWorkflowArguments = cactusWorkflowArguments 
        self.event = event
        self.prepXmlElems = prepXmlElems
        self.inSequence = inSequence
        self.globalOutSequence = globalOutSequence
        self.memory = memory
        self.cpu = cpu
        self.iteration = iteration
              
    def run(self):
        # Parse the "preprocessor" config xml element     
        assert self.iteration < len(self.prepXmlElems)
        
        prepNode = self.prepXmlElems[self.iteration]
        prepOptions = PreprocessorOptions(int(prepNode.get("chunkSize", default="2147483647")),
                                          prepNode.attrib["preprocessorString"],
                                          int(self.memory),
                                          int(self.cpu))
        
        if os.path.isdir(self.inSequence):
            tempFile = os.path.join(self.getGlobalTempDir(), "catSeq.fa")
            catFiles([ os.path.join(self.inSequence, f) for f in os.listdir(self.inSequence)], tempFile)
            self.inSequence = tempFile
        
        #output to temporary directory unless we are on the last iteration
        lastIteration = self.iteration == len(self.prepXmlElems) - 1
        if lastIteration == False:
            outSeq = os.path.join(self.getGlobalTempDir(), str(self.iteration))
        else:
            outSeq = self.globalOutSequence
        
        self.addChildTarget(PreprocessSequence(prepOptions, self.inSequence, outSeq, self.event)) 
        
        if lastIteration == False:
            self.setFollowOnTarget(BatchPreprocessor(self.cactusWorkflowArguments, self.event, 
                                                     self.prepXmlElems, outSeq,
                                                     self.globalOutSequence, self.memory, 
                                                     self.cpu, self.iteration + 1))   
