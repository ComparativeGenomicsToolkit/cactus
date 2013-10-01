#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
import sys
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
from cactus.shared.common import getOptionalAttrib

class PreprocessorOptions:
    def __init__(self, chunkSize, cmdLine, memory, cpu, check):
        self.chunkSize = chunkSize
        self.cmdLine = cmdLine
        self.memory = memory
        self.cpu = cpu
        self.check = check

class PreprocessChunk(Target):
    """ locally preprocess a fasta chunk, output then copied back to input
    """
    def __init__(self, prepOptions, seqPath, inChunk, outChunk):
        Target.__init__(self, memory=prepOptions.memory, cpu=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.seqPath = seqPath
        self.inChunk = inChunk
        self.outChunk = outChunk
    
    def run(self):
        tempPath = os.path.join(self.getLocalTempDir(), "tempPath")
        cmdline = self.prepOptions.cmdLine.replace("GENOME_FILE", "\"" + self.seqPath + "\"")
        cmdline = cmdline.replace("IN_FILE", "\"" + self.inChunk + "\"")
        cmdline = cmdline.replace("OUT_FILE", "\"" + self.outChunk + "\"")
        cmdline = cmdline.replace("TEMP_FILE", "\"" + tempPath + "\"")
        logger.info("Preprocessor exec " + cmdline)
        system(cmdline)
        if self.prepOptions.check:
            system("cp %s %s" % (self.inChunk, self.outChunk))

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
    """Cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequencePath, outSequencePath):
        Target.__init__(self, cpu=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.inSequencePath = inSequencePath
        self.outSequencePath = outSequencePath
    
    def run(self):        
        logger.info("Preparing sequence for preprocessing")
        # chunk it up
        inChunkDirectory = makeSubDir(os.path.join(self.getGlobalTempDir(), "preprocessChunksIn"))
        inChunkList = [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i 0 %s %s" % \
               (getLogLevelString(), self.prepOptions.chunkSize,
                inChunkDirectory, self.inSequencePath)).split("\n") if chunk != "" ]   
        outChunkDirectory = makeSubDir(os.path.join(self.getGlobalTempDir(), "preprocessChunksOut"))
        outChunkList = [] 
        #For each input chunk we create an output chunk, it is the output chunks that get concatenated together.
        for i in xrange(len(inChunkList)):
            outChunkList.append(os.path.join(outChunkDirectory, "chunk_%i" % i))
            self.addChildTarget(PreprocessChunk(self.prepOptions, self.inSequencePath, inChunkList[i], outChunkList[i]))
        # follow on to merge chunks
        self.setFollowOnTarget(MergeChunks(self.prepOptions, outChunkList, self.outSequencePath))

class BatchPreprocessor(Target):
    def __init__(self, prepXmlElems, inSequence, 
                 globalOutSequence, iteration = 0):
        Target.__init__(self, time=0.0002) 
        self.prepXmlElems = prepXmlElems
        self.inSequence = inSequence
        self.globalOutSequence = globalOutSequence
        prepNode = self.prepXmlElems[iteration]
        self.memory = getOptionalAttrib(prepNode, "memory", typeFn=int, default=sys.maxint)
        self.cpu = getOptionalAttrib(prepNode, "cpu", typeFn=int, default=sys.maxint)
        self.iteration = iteration
              
    def run(self):
        # Parse the "preprocessor" config xml element     
        assert self.iteration < len(self.prepXmlElems)
        
        prepNode = self.prepXmlElems[self.iteration]
        prepOptions = PreprocessorOptions(int(prepNode.get("chunkSize", default="-1")),
                                          prepNode.attrib["preprocessorString"],
                                          int(self.memory),
                                          int(self.cpu),
                                          bool(int(prepNode.get("check", default="0"))),)
        
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
        
        if prepOptions.chunkSize <= 0: #In this first case we don't need to break up the sequence
            self.addChildTarget(PreprocessChunk(prepOptions, self.inSequence, self.inSequence, outSeq))
        else:
            self.addChildTarget(PreprocessSequence(prepOptions, self.inSequence, outSeq)) 
        
        if lastIteration == False:
            self.setFollowOnTarget(BatchPreprocessor(self.prepXmlElems, outSeq,
                                                     self.globalOutSequence, self.iteration + 1))   
            

############################################################
############################################################
############################################################
##The preprocessor phase, which modifies the input sequences
############################################################
############################################################
############################################################

def getOutputSequenceFile(outputSequenceDir, inputSequenceFile):
    return os.path.join(outputSequenceDir, os.path.split(inputSequenceFile)[1] + ".preprocessed")        

class CactusPreprocessor(Target):
    """Modifies the input genomes, doing things like masking/checking, etc.
    """
    def __init__(self, inputSequences, outputSequenceDir, configNode):
        Target.__init__(self)
        self.inputSequences = inputSequences
        self.outputSequenceDir = outputSequenceDir
        self.configNode = configNode  
    def run(self):
        if not os.path.isdir(self.outputSequenceDir):
            os.mkdir(self.outputSequenceDir)
        for sequence in self.inputSequences:
            outputSequenceFile = getOutputSequenceFile(self.outputSequenceDir, sequence)
            assert sequence != outputSequenceFile
            if not os.path.isfile(outputSequenceFile): #Only create the output sequence if it doesn't already exist. This prevents reprocessing if the sequence is used in multiple places between runs.
                prepXmlElems = self.configNode.findall("preprocessor")
                if len(prepXmlElems) == 0: #Just ln the file to the output dir
                    system("ln %s %s" % (sequence, outputSequenceFile))
                else:
                    logger.info("Adding child batch_preprocessor target")
                    self.addChildTarget(BatchPreprocessor(prepXmlElems, 
                                                          sequence, outputSequenceFile, 
                                                          0))            

