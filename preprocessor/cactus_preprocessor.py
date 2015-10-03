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
import math
import errno
from argparse import ArgumentParser
from bz2 import BZ2File
import copy
import xml.etree.ElementTree as ET

from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch, popenPush
from sonLib.bioio import getLogLevelString
from sonLib.bioio import newickTreeParser
from sonLib.bioio import makeSubDir
from sonLib.bioio import catFiles, getTempFile, getTempDirectory
from toil.job import Job

from cactus.shared.common import getOptionalAttrib, runCactusAnalyseAssembly
from cactus.shared.common import WritePermanentFile
from toil.lib.bioio import setLoggingFromOptions
from cactus.shared.configWrapper import ConfigWrapper

class PreprocessorOptions:
    def __init__(self, chunkSize, cmdLine, memory, cpu, check, proportionToSample):
        self.chunkSize = chunkSize
        self.cmdLine = cmdLine
        self.memory = memory
        self.cpu = cpu
        self.check = check
        self.proportionToSample=proportionToSample

class PreprocessChunk(Job):
    """ locally preprocess a fasta chunk, output then copied back to input
    """
    def __init__(self, prepOptions, seqIDs, proportionSampled, inChunkID):
        Job.__init__(self, memory=prepOptions.memory, cores=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.seqIDs = seqIDs
        self.inChunkID = inChunkID
        self.proportionSampled = proportionSampled
    
    def run(self, fileStore):
        inChunk = fileStore.readGlobalFile(self.inChunkID)
        outChunk = fileStore.getLocalTempFile()
        seqPaths = [fileStore.readGlobalFile(fileID) for fileID in self.seqIDs]
        cmdline = self.prepOptions.cmdLine.replace("IN_FILE", "\"" + inChunk + "\"")
        cmdline = cmdline.replace("OUT_FILE", "\"" + outChunk + "\"")
        cmdline = cmdline.replace("TEMP_DIR", "\"" + fileStore.getLocalTempDir() + "\"")
        cmdline = cmdline.replace("PROPORTION_SAMPLED", str(self.proportionSampled))
        logger.info("Preprocessor exec " + cmdline)
        #print "command", cmdline
        #sys.exit(1)
        popenPush(cmdline, " ".join(seqPaths))
        if self.prepOptions.check:
            system("cp %s %s" % (inChunk, outChunk))
        return fileStore.writeGlobalFile(outChunk)

class MergeChunks(Job):
    """ merge a list of chunks into a fasta file
    """
    def __init__(self, prepOptions, chunkIDList):
        Job.__init__(self, cores=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.chunkIDList = chunkIDList
    
    def run(self, fileStore):
        chunkList = [fileStore.readGlobalFile(fileID) for fileID in self.chunkIDList]
        outSequencePath = fileStore.getLocalTempFile()
        popenPush("cactus_batch_mergeChunks > %s" % outSequencePath, " ".join(chunkList))
        return fileStore.writeGlobalFile(outSequencePath)
 
class PreprocessSequence(Job):
    """Cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequenceID):
        Job.__init__(self, cores=prepOptions.cpu)
        self.prepOptions = prepOptions 
        self.inSequenceID = inSequenceID
    
    def run(self, fileStore):        
        logger.info("Preparing sequence for preprocessing")
        # chunk it up
        inSequence = fileStore.readGlobalFile(self.inSequenceID)
        inChunkDirectory = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        inChunkList = [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i 0 %s %s" % \
               (getLogLevelString(), self.prepOptions.chunkSize,
                inChunkDirectory, inSequence)).split("\n") if chunk != "" ]
        inChunkIDList = [fileStore.writeGlobalFile(chunk) for chunk in inChunkList]
        outChunkIDList = [] 
        #For each input chunk we create an output chunk, it is the output chunks that get concatenated together.
        for i in xrange(len(inChunkList)):
            #Calculate the number of chunks to use
            inChunkNumber = int(max(1, math.ceil(len(inChunkList) * self.prepOptions.proportionToSample)))
            assert inChunkNumber <= len(inChunkList) and inChunkNumber > 0
            #Now get the list of chunks flanking and including the current chunk
            j = max(0, i - inChunkNumber/2)
            inChunkIDs = inChunkIDList[j:j+inChunkNumber]
            if len(inChunkIDs) < inChunkNumber: #This logic is like making the list circular
                inChunkIDs += inChunkIDList[:inChunkNumber-len(inChunkIDs)]
            assert len(inChunkIDs) == inChunkNumber
            outChunkIDList.append(self.addChild(PreprocessChunk(self.prepOptions, inChunkIDs, float(inChunkNumber)/len(inChunkIDList), inChunkIDList[i])).rv())
        # follow on to merge chunks
        return self.addFollowOn(MergeChunks(self.prepOptions, outChunkIDList)).rv()

class BatchPreprocessor(Job):
    def __init__(self, prepXmlElems, inSequenceID, iteration = 0):
        Job.__init__(self) 
        self.prepXmlElems = prepXmlElems
        self.inSequenceID = inSequenceID
        prepNode = self.prepXmlElems[iteration]
        self.memory = getOptionalAttrib(prepNode, "memory", typeFn=int, default=None)
        self.cores = getOptionalAttrib(prepNode, "cpu", typeFn=int, default=None)
        self.iteration = iteration
              
    def run(self, fileStore):
        # Parse the "preprocessor" config xml element     
        assert self.iteration < len(self.prepXmlElems)
        
        prepNode = self.prepXmlElems[self.iteration]
        prepOptions = PreprocessorOptions(int(prepNode.get("chunkSize", default="-1")),
                                          prepNode.attrib["preprocessorString"],
                                          self.memory,
                                          self.cores,
                                          bool(int(prepNode.get("check", default="0"))),
                                          getOptionalAttrib(prepNode, "proportionToSample", typeFn=float, default=1.0))
        
        #output to temporary directory unless we are on the last iteration
        lastIteration = self.iteration == len(self.prepXmlElems) - 1
        
        if prepOptions.chunkSize <= 0: #In this first case we don't need to break up the sequence
            outSeqID = self.addChild(PreprocessChunk(prepOptions, [ self.inSequenceID ], 1.0, self.inSequenceID)).rv()
        else:
            outSeqID = self.addChild(PreprocessSequence(prepOptions, self.inSequenceID)).rv()
        
        if lastIteration == False:
            return self.addFollowOn(BatchPreprocessor(self.prepXmlElems, outSeqID, self.iteration + 1)).rv()
        else:
            return self.addFollowOn(BatchPreprocessorEnd(outSeqID)).rv()

class BatchPreprocessorEnd(Job):
    def __init__(self,  globalOutSequenceID):
        Job.__init__(self) 
        self.globalOutSequenceID = globalOutSequenceID
        
    def run(self, fileStore):
        globalOutSequence = fileStore.readGlobalFile(self.globalOutSequenceID)
        analysisString = runCactusAnalyseAssembly(globalOutSequence)
        fileStore.logToMaster("After preprocessing assembly we got the following stats: %s" % analysisString)
        return self.globalOutSequenceID

############################################################
############################################################
############################################################
##The preprocessor phase, which modifies the input sequences
############################################################
############################################################
############################################################

class CactusPreprocessor(Job):
    """Modifies the input genomes, doing things like masking/checking, etc.
    """
    def __init__(self, inputSequences, outputSequences, configNode):
        Job.__init__(self)
        self.inputSequences = inputSequences
        self.outputSequences = outputSequences
        assert len(self.inputSequences) == len(self.outputSequences) #If these are not the same length then we have a problem
        self.configNode = configNode  
    
    def run(self, fileStore):
        for inputSequenceFileOrDirectory, outputSequenceFile in zip(self.inputSequences, self.outputSequences):
            if not os.path.isfile(outputSequenceFile): #Only create the output sequence if it doesn't already exist. This prevents reprocessing if the sequence is used in multiple places between runs.
                self.addChild(CactusPreprocessor2(inputSequenceFileOrDirectory, outputSequenceFile, self.configNode))
  
    @staticmethod
    def getOutputSequenceFiles(inputSequences, outputSequenceDir):
        """Function to get unambiguous file names for each input sequence in the output sequence dir. 
        """
        if not os.path.isdir(outputSequenceDir):
            os.mkdir(outputSequenceDir)
        return [ os.path.join(outputSequenceDir, inputSequences[i].split("/")[-1] + "_%i" % i) for i in xrange(len(inputSequences)) ]
        #return [ os.path.join(outputSequenceDir, "_".join(inputSequence.split("/"))) for inputSequence in inputSequences ]
  
class CactusPreprocessor2(Job):
    def __init__(self, inputSequenceFileOrDirectory, outputSequenceFile, configNode):
        Job.__init__(self)
        self.inputSequenceFileOrDirectory = inputSequenceFileOrDirectory
        self.outputSequenceFile = outputSequenceFile
        self.configNode = configNode
        
    def run(self, fileStore):
        #If the files are in a sub-dir then rip them out.
        if os.path.isdir(self.inputSequenceFileOrDirectory):
            tempFile = fileStore.getLocalTempFile()
            catFiles([ os.path.join(self.inputSequenceFileOrDirectory, f) for f in os.listdir(self.inputSequenceFileOrDirectory)], tempFile)
            inputSequenceFile = tempFile
        else:
            inputSequenceFile = self.inputSequenceFileOrDirectory
            
        assert inputSequenceFile != self.outputSequenceFile
        
        prepXmlElems = self.configNode.findall("preprocessor")
        
        analysisString = runCactusAnalyseAssembly(inputSequenceFile)
        fileStore.logToMaster("Before running any preprocessing on the assembly: %s got following stats (assembly may be listed as temp file if input sequences from a directory): %s" % \
                         (self.inputSequenceFileOrDirectory, analysisString))
        
        if len(prepXmlElems) == 0: #Just cp the file to the output file
            system("cp %s %s" % (inputSequenceFile, self.outputSequenceFile))
        else:
            logger.info("Adding child batch_preprocessor target")
            inputSequenceID = fileStore.writeGlobalFile(inputSequenceFile)
            outputSequenceID = self.addChild(BatchPreprocessor(prepXmlElems, inputSequenceID, 0)).rv()
            self.addFollowOn(WritePermanentFile(outputSequenceID, self.outputSequenceFile))
                    
def main():
    usage = "usage: %prog outputSequenceDir configXMLFile inputSequenceFastaFilesxN [options]"
    parser = ArgumentParser(usage=usage)
    Job.Runner.addToilOptions(parser)
    parser.add_argument("--outputSequenceDir", dest="outputSequenceDir", type=str)
    parser.add_argument("--configFile", dest="configFile", type=str)
    parser.add_argument("--inputSequences", dest="inputSequences", type=str)
    
    options = parser.parse_args()
    setLoggingFromOptions(options)
    
    if not (options.outputSequenceDir and options.configFile and options.inputSequences):
        raise RuntimeError("Too few input arguments")
    
    
    #Replace any constants
    configNode = ET.parse(options.configFile).getroot()
    inputSequences = options.inputSequences.split()
    outputSequences = CactusPreprocessor.getOutputSequenceFiles(inputSequences, options.outputSequenceDir)
    if configNode.find("constants") != None:
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals()
    
    Job.Runner.startToil(CactusPreprocessor(inputSequences, outputSequences, configNode), options)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.preprocessor.cactus_preprocessor import *
    main()
