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

from toil.lib.bioio import logger
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString

from sonLib.bioio import catFiles, getTempFile, getTempDirectory, popenCatch, popenPush, newickTreeParser
from toil.job import Job
from toil.common import Toil
from cactus.shared.common import cactus_call

from cactus.shared.common import getOptionalAttrib, runCactusAnalyseAssembly
from toil.lib.bioio import setLoggingFromOptions
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import nameValue
from cactus.shared.common import runGetChunks
from cactus.shared.common import makeURL

from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import lastzRepeatMaskJob

class PreprocessorOptions:
    def __init__(self, chunkSize, memory, cpu, check, proportionToSample, unmask,
                 preprocessJob, checkAssemblyHub=None, lastzOptions=None, minPeriod=None):
        self.chunkSize = chunkSize
        self.memory = memory
        self.cpu = cpu
        self.check = check
        self.proportionToSample=proportionToSample
        self.unmask = unmask
        self.preprocessJob = preprocessJob
        self.checkAssemblyHub = checkAssemblyHub
        self.lastzOptions = lastzOptions
        self.minPeriod = minPeriod

class PreprocessChunk(Job):
    """ locally preprocess a fasta chunk, output then copied back to input
    """
    def __init__(self, prepOptions, seqIDs, proportionSampled, inChunkID, disk=None):
        Job.__init__(self, memory=prepOptions.memory, cores=prepOptions.cpu, disk=disk)
        self.prepOptions = prepOptions 
        self.seqIDs = seqIDs
        self.inChunkID = inChunkID
    
    def run(self, fileStore):
        outChunkID = None
        if self.prepOptions.preprocessJob == "checkUniqueHeaders":
            inChunk = fileStore.readGlobalFile(self.inChunkID)
            outChunk = fileStore.getLocalTempFile()
            seqPaths = [fileStore.readGlobalFile(fileID) for fileID in self.seqIDs]
            seqString = " ".join(seqPaths)
            cactus_call(tool="cactus", stdin_string=seqString,
                        parameters=["cactus_checkUniqueHeaders.py",
                                    nameValue("checkAssemblyHub", self.prepOptions.checkAssemblyHub, bool),
                                    outChunk])
            outChunkID = fileStore.writeGlobalFile(outChunk)
        elif self.prepOptions.preprocessJob == "lastzRepeatMask":
            outChunkID = self.addChildJobFn(lastzRepeatMaskJob, queryID=self.inChunkID, targetIDs=self.seqIDs,
                                      proportionSampled=self.prepOptions.proportionToSample,
                                      minPeriod=self.prepOptions.minPeriod).rv()
        
        
        if self.prepOptions.check:
            return self.inChunkID
        else:
            return outChunkID

class MergeChunks(Job):
    def __init__(self, prepOptions, chunkIDList):
        Job.__init__(self)
        self.prepOptions = prepOptions
        self.chunkIDList = chunkIDList

    def run(self, fileStore):
        chunkIDListSize = sum([chunkID.size for chunkID in self.chunkIDList])
        return self.addFollowOn(MergeChunks2(self.prepOptions, self.chunkIDList, disk=3*chunkIDListSize)).rv()


class MergeChunks2(Job):
    """ merge a list of chunks into a fasta file
    """
    def __init__(self, prepOptions, chunkIDList, disk=None):
        Job.__init__(self, cores=prepOptions.cpu, memory=prepOptions.memory, disk=disk)
        self.prepOptions = prepOptions 
        self.chunkIDList = chunkIDList
    
    def run(self, fileStore):
        chunkList = [fileStore.readGlobalFile(fileID) for fileID in self.chunkIDList]

        #Docker expects paths relative to the work dir
        chunkList = [os.path.basename(chunk) for chunk in chunkList]
        outSequencePath = fileStore.getLocalTempFile()
        cactus_call(tool="cactus", outfile=outSequencePath, stdin_string=" ".join(chunkList),
                    parameters=["cactus_batch_mergeChunks"])
        return fileStore.writeGlobalFile(outSequencePath)
 
class PreprocessSequence(Job):
    """Cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequenceID, disk=None):
        Job.__init__(self, cores=prepOptions.cpu, memory=prepOptions.memory, disk=disk)
        self.prepOptions = prepOptions 
        self.inSequenceID = inSequenceID
    
    def run(self, fileStore):        
        logger.info("Preparing sequence for preprocessing")
        # chunk it up
        inSequence = fileStore.readGlobalFile(self.inSequenceID)
        inChunkDirectory = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        inChunkList = runGetChunks(sequenceFiles=[inSequence], chunksDir=inChunkDirectory,
                                   chunkSize=self.prepOptions.chunkSize,
                                   overlapSize=0)
        inChunkList = [os.path.abspath(path) for path in inChunkList]
        logger.info("Chunks = %s" % inChunkList)
        logger.info("Chunks dir = %s" % os.listdir(inChunkDirectory))
                
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
            diskForPreprocessChunk = 2*(sum([inChunkID.size for inChunkID in inChunkIDs]) + 
                    inChunkIDList[i].size)
            outChunkIDList.append(self.addChild(PreprocessChunk(self.prepOptions, inChunkIDs, float(inChunkNumber)/len(inChunkIDList), inChunkIDList[i], disk=diskForPreprocessChunk)).rv())
        # follow on to merge chunks
        return self.addFollowOn(MergeChunks(self.prepOptions, outChunkIDList)).rv()

def unmaskFasta(inFasta, outFasta):
    """Uppercase a fasta file (removing the soft-masking)."""
    with open(outFasta, 'w') as out:
        for line in open(inFasta):
            if len(line) > 0 and line[0] != ">":
                out.write(line.upper())
            else:
                out.write(line)

class BatchPreprocessor(Job):
    def __init__(self, prepXmlElems, inSequenceID, iteration = 0):
        Job.__init__(self) 
        self.prepXmlElems = prepXmlElems
        self.inSequenceID = inSequenceID
        prepNode = self.prepXmlElems[iteration]
        self._memory = getOptionalAttrib(prepNode, "memory", typeFn=int, default=None)
        self._cores = getOptionalAttrib(prepNode, "cpu", typeFn=int, default=None)
        self.iteration = iteration
              
    def run(self, fileStore):
        # Parse the "preprocessor" config xml element     
        assert self.iteration < len(self.prepXmlElems)
        
        prepNode = self.prepXmlElems[self.iteration]
        prepOptions = PreprocessorOptions(chunkSize = int(prepNode.get("chunkSize", default="-1")),
                                          preprocessJob=prepNode.attrib["preprocessJob"],
                                          memory = int(prepNode.get("memory", default=0)),
                                          cpu = int(prepNode.get("cpu", default=1)),
                                          check = bool(int(prepNode.get("check", default="0"))),
                                          proportionToSample = getOptionalAttrib(prepNode, "proportionToSample", typeFn=float, default=1.0),
                                          unmask = getOptionalAttrib(prepNode, "unmask", typeFn=bool, default=False),
                                          lastzOptions = getOptionalAttrib(prepNode, "lastzOpts", default=""),
                                          minPeriod = getOptionalAttrib(prepNode, "minPeriod", typeFn=int, default="0"),
                                          checkAssemblyHub = getOptionalAttrib(prepNode, "checkAssemblyHub", typeFn=bool, default=False))
        
        #output to temporary directory unless we are on the last iteration
        lastIteration = self.iteration == len(self.prepXmlElems) - 1


        if prepOptions.unmask:
            inSequence = fileStore.readGlobalFile(self.inSequenceID)
            unmaskedInputFile = fileStore.getLocalTempFile()
            unmaskFasta(inSequence, unmaskedInputFile)
            self.inSequenceID = fileStore.writeGlobalFile(inSequence)
            
        inSequenceSize = 4*1024*1024*1024
        if prepOptions.chunkSize <= 0: #In this first case we don't need to break up the sequence
            #Estimate the size of a genome since we can't get the sizes of files imported
            #through Toil.importFile yet.
            outSeqID = self.addChild(PreprocessChunk(prepOptions, [ self.inSequenceID ], 1.0, self.inSequenceID, disk=3*inSequenceSize)).rv()
        else:
            outSeqID = self.addChild(PreprocessSequence(prepOptions, self.inSequenceID, disk=2*inSequenceSize)).rv()
        
        if lastIteration == False:
            return self.addFollowOn(BatchPreprocessor(self.prepXmlElems, outSeqID, self.iteration + 1)).rv()
        else:
            return self.addFollowOn(BatchPreprocessorEnd(outSeqID, disk=2*inSequenceSize)).rv()



class BatchPreprocessorEnd(Job):
    def __init__(self,  globalOutSequenceID, disk=None):
        Job.__init__(self, disk=disk) 
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
    def __init__(self, inputSequenceIDs, configNode):
        Job.__init__(self)
        self.inputSequenceIDs = inputSequenceIDs
        self.configNode = configNode  
    
    def run(self, fileStore):
        outputSequenceIDs = []
        for inputSequenceID in self.inputSequenceIDs:
            outputSequenceIDs.append(self.addChild(CactusPreprocessor2(inputSequenceID, self.configNode)).rv())
        return outputSequenceIDs
  
    @staticmethod
    def getOutputSequenceFiles(inputSequences, outputSequenceDir):
        """Function to get unambiguous file names for each input sequence in the output sequence dir. 
        """
        if not os.path.isdir(outputSequenceDir):
            os.mkdir(outputSequenceDir)
        return [ os.path.join(outputSequenceDir, inputSequences[i].split("/")[-1] + "_%i" % i) for i in xrange(len(inputSequences)) ]
  
class CactusPreprocessor2(Job):
    def __init__(self, inputSequenceID, configNode):
        Job.__init__(self)
        self.inputSequenceID = inputSequenceID
        self.configNode = configNode
        
    def run(self, fileStore):
        inputSequenceFile = fileStore.readGlobalFile(self.inputSequenceID)
        prepXmlElems = self.configNode.findall("preprocessor")
        
        #analysisString = runCactusAnalyseAssembly(inputSequenceFile)
        #fileStore.logToMaster("Before running any preprocessing on the assembly: %s got following stats (assembly may be listed as temp file if input sequences from a directory): %s" % \
        #                 (self.inputSequenceID, analysisString))
        
        if len(prepXmlElems) == 0: #Just cp the file to the output file
            return self.inputSequenceID
        else:
            logger.info("Adding child batch_preprocessor target")
            return self.addChild(BatchPreprocessor(prepXmlElems, self.inputSequenceID, 0)).rv()

def stageWorkflow(outputSequenceDir, configFile, inputSequences, toil):
    #Replace any constants
    configNode = ET.parse(configFile).getroot()
    outputSequences = CactusPreprocessor.getOutputSequenceFiles(inputSequences, outputSequenceDir)
    if configNode.find("constants") != None:
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals()
    inputSequenceIDs = [toil.importFile(makeURL(seq)) for seq in inputSequences]
    outputSequenceIDs = toil.start(CactusPreprocessor(inputSequenceIDs, configNode))
    for seqID, path in zip(outputSequenceIDs, outputSequences):
        toil.exportFile(seqID, makeURL(path))

def runCactusPreprocessor(outputSequenceDir, configFile, inputSequences, toilDir):
    toilOptions = Job.Runner.getDefaultOptions(toilDir)
    with Toil(toilOptions) as toil:
        stageWorkflow(outputSequenceDir, configFile, inputSequences, toil)
                    
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
    with Toil(options) as toil:
        stageWorkflow(outputSquenceDir=options.outputSequenceDir, configFile=options.configFile, inputSequences=options.inputSequences.split(), toil=toil)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.preprocessor.cactus_preprocessor import *
    main()
