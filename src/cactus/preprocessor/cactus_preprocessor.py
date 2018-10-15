#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
import math
from argparse import ArgumentParser
import xml.etree.ElementTree as ET

from toil.lib.bioio import logger

from sonLib.bioio import getTempDirectory
from toil.common import Toil
from toil.job import Job
from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import runGetChunks
from cactus.shared.common import makeURL
from cactus.shared.common import readGlobalFileWithoutCache
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper

from toil.lib.bioio import setLoggingFromOptions

from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import LastzRepeatMaskJob
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import RepeatMaskOptions

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

class PreprocessChunk(RoundedJob):
    """locally preprocess a fasta chunk, output then copied back to input"""
    def __init__(self, prepOptions, seqIDs, proportionSampled, inChunkID):
        disk = sum([seqID.size for seqID in seqIDs]) + 3*inChunkID.size
        RoundedJob.__init__(self, memory=prepOptions.memory, cores=prepOptions.cpu, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions 
        self.seqIDs = seqIDs
        self.inChunkID = inChunkID

    def run(self, fileStore):
        outChunkID = None
        if self.prepOptions.preprocessJob == "checkUniqueHeaders":
            inChunk = fileStore.readGlobalFile(self.inChunkID)
            seqPaths = [fileStore.readGlobalFile(fileID) for fileID in self.seqIDs]
            seqString = " ".join(seqPaths)
            args = [inChunk]
            if self.prepOptions.checkAssemblyHub:
                args += ["--checkAssemblyHub"]
            cactus_call(stdin_string=seqString,
                        parameters=["cactus_checkUniqueHeaders.py"] + args)
            outChunkID = self.inChunkID
        elif self.prepOptions.preprocessJob == "lastzRepeatMask":
            repeatMaskOptions = RepeatMaskOptions(proportionSampled=self.prepOptions.proportionToSample,
                    minPeriod=self.prepOptions.minPeriod)
            outChunkID = self.addChild(LastzRepeatMaskJob(repeatMaskOptions=repeatMaskOptions, 
                    queryID=self.inChunkID, targetIDs=self.seqIDs)).rv()
        elif self.prepOptions.preprocessJob == "none":
            outChunkID = self.inChunkID

        return outChunkID

class MergeChunks(RoundedJob):
    def __init__(self, prepOptions, chunkIDList):
        RoundedJob.__init__(self, preemptable=True)
        self.prepOptions = prepOptions
        self.chunkIDList = chunkIDList

    def run(self, fileStore):
        return self.addFollowOn(MergeChunks2(self.prepOptions, self.chunkIDList)).rv()

class MergeChunks2(RoundedJob):
    """merge a list of chunks into a fasta file"""
    def __init__(self, prepOptions, chunkIDList):
        disk = 2*sum([chunkID.size for chunkID in chunkIDList])
        RoundedJob.__init__(self, cores=prepOptions.cpu, memory=prepOptions.memory, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions 
        self.chunkIDList = chunkIDList

    def run(self, fileStore):
        chunkList = [readGlobalFileWithoutCache(fileStore, fileID) for fileID in self.chunkIDList]

        #Docker expects paths relative to the work dir
        chunkList = [os.path.basename(chunk) for chunk in chunkList]
        outSequencePath = fileStore.getLocalTempFile()
        cactus_call(outfile=outSequencePath, stdin_string=" ".join(chunkList),
                    parameters=["cactus_batch_mergeChunks"])
        return fileStore.writeGlobalFile(outSequencePath)

class PreprocessSequence(RoundedJob):
    """Cut a sequence into chunks, process, then merge
    """
    def __init__(self, prepOptions, inSequenceID, chunksToCompute=None):
        disk = 3*inSequenceID.size if hasattr(inSequenceID, "size") else None
        RoundedJob.__init__(self, cores=prepOptions.cpu, memory=prepOptions.memory, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions 
        self.inSequenceID = inSequenceID
        self.chunksToCompute = chunksToCompute

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
        if not self.chunksToCompute:
            self.chunksToCompute = range(len(inChunkList))
        for i in self.chunksToCompute:
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

def unmaskFasta(inFasta, outFasta):
    """Uppercase a fasta file (removing the soft-masking)."""
    with open(outFasta, 'w') as out:
        for line in open(inFasta):
            if len(line) > 0 and line[0] != ">":
                out.write(line.upper())
            else:
                out.write(line)

class BatchPreprocessor(RoundedJob):
    def __init__(self, prepXmlElems, inSequenceID, iteration = 0):
        self.prepXmlElems = prepXmlElems
        self.inSequenceID = inSequenceID
        self.iteration = iteration
        RoundedJob.__init__(self, preemptable=True)
              
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
        
        lastIteration = self.iteration == len(self.prepXmlElems) - 1

        if prepOptions.unmask:
            inSequence = fileStore.readGlobalFile(self.inSequenceID)
            unmaskedInputFile = fileStore.getLocalTempFile()
            unmaskFasta(inSequence, unmaskedInputFile)
            self.inSequenceID = fileStore.writeGlobalFile(inSequence)
            
        if prepOptions.chunkSize <= 0: #In this first case we don't need to break up the sequence
            outSeqID = self.addChild(PreprocessChunk(prepOptions, [ self.inSequenceID ], 1.0, self.inSequenceID)).rv()
        else:
            outSeqID = self.addChild(PreprocessSequence(prepOptions, self.inSequenceID)).rv()
        
        if lastIteration == False:
            return self.addFollowOn(BatchPreprocessor(self.prepXmlElems, outSeqID, self.iteration + 1)).rv()
        else:
            return outSeqID

############################################################
############################################################
############################################################
##The preprocessor phase, which modifies the input sequences
############################################################
############################################################
############################################################

class CactusPreprocessor(RoundedJob):
    """Modifies the input genomes, doing things like masking/checking, etc.
    """
    def __init__(self, inputSequenceIDs, configNode):
        RoundedJob.__init__(self, disk=sum([id.size for id in inputSequenceIDs]), preemptable=True)
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
        if not outputSequenceDir.startswith("s3://") and not os.path.isdir(outputSequenceDir):
            os.mkdir(outputSequenceDir)
        return [ os.path.join(outputSequenceDir, inputSequences[i].split("/")[-1] + "_%i" % i) for i in xrange(len(inputSequences)) ]

class CactusPreprocessor2(RoundedJob):
    def __init__(self, inputSequenceID, configNode):
        RoundedJob.__init__(self, preemptable=True)
        self.inputSequenceID = inputSequenceID
        self.configNode = configNode
        
    def run(self, fileStore):
        prepXmlElems = self.configNode.findall("preprocessor")

        if len(prepXmlElems) == 0: #Just cp the file to the output file
            return self.inputSequenceID
        else:
            logger.info("Adding child batch_preprocessor target")
            return self.addChild(BatchPreprocessor(prepXmlElems, self.inputSequenceID, 0)).rv()

def stageWorkflow(outputSequenceDir, configFile, inputSequences, toil, restart=False):
    #Replace any constants
    configNode = ET.parse(configFile).getroot()
    outputSequences = CactusPreprocessor.getOutputSequenceFiles(inputSequences, outputSequenceDir)
    if configNode.find("constants") != None:
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals()
    if not restart:
        inputSequenceIDs = [toil.importFile(makeURL(seq)) for seq in inputSequences]
        outputSequenceIDs = toil.start(CactusPreprocessor(inputSequenceIDs, configNode))
    else:
        outputSequenceIDs = toil.restart()
    for seqID, path in zip(outputSequenceIDs, outputSequences):
        toil.exportFile(seqID, makeURL(path))

def runCactusPreprocessor(outputSequenceDir, configFile, inputSequences, toilDir):
    toilOptions = Job.Runner.getDefaultOptions(toilDir)
    toilOptions.logLevel = "INFO"
    toilOptions.disableCaching = True
    with Toil(toilOptions) as toil:
        stageWorkflow(outputSequenceDir, configFile, inputSequences, toil)

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument("outputSequenceDir", help='Directory where the processed sequences will be placed')
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("inputSequences", nargs='+', help='input FASTA file(s)')

    options = parser.parse_args()
    setLoggingFromOptions(options)

    with Toil(options) as toil:
        stageWorkflow(outputSequenceDir=options.outputSequenceDir, configFile=options.configFile, inputSequences=options.inputSequences, toil=toil, restart=options.restart)

if __name__ == '__main__':
    main()
