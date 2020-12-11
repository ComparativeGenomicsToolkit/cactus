#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""

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
from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import enableDumpStack
from toil.lib.bioio import setLoggingFromOptions
from toil.realtimeLogger import RealtimeLogger

from cactus.shared.common import cactus_override_toil_options
from cactus.preprocessor.checkUniqueHeaders import checkUniqueHeaders
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import LastzRepeatMaskJob
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import RepeatMaskOptions
from cactus.preprocessor.dnabrnnMasking import DnabrnnMaskJob

class PreprocessorOptions:
    def __init__(self, chunkSize, memory, cpu, check, proportionToSample, unmask,
                 preprocessJob, checkAssemblyHub=None, lastzOptions=None, minPeriod=None,
                 gpuLastz=False, dnabrnnOpts=None, dnabrnnLength=None, hardmask=None):
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
        self.gpuLastz = gpuLastz
        self.gpuLastzInterval = self.chunkSize
        if self.gpuLastz:
            self.chunkSize = 0
        self.dnabrnnOpts = dnabrnnOpts
        self.dnabrnnLength = dnabrnnLength
        self.hardmask = hardmask

class CheckUniqueHeaders(RoundedJob):
    """
    Check that the headers of the input file meet certain naming requirements.
    """
    def __init__(self, prepOptions, inChunkID):
        disk = inChunkID.size
        RoundedJob.__init__(self, memory=prepOptions.memory, cores=prepOptions.cpu, disk=disk,
                     preemptable=True)
        self.prepOptions = prepOptions
        self.inChunkID = inChunkID

    def run(self, fileStore):
        inChunk = fileStore.readGlobalFile(self.inChunkID)
        with open(inChunk) as inFile:
            checkUniqueHeaders(inFile, checkAssemblyHub=self.prepOptions.checkAssemblyHub)
        # We re-write the file here so that the output's lifecycle
        # matches the other chunked jobs, which usually write a new
        # chunk.
        return fileStore.writeGlobalFile(inChunk)

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

    def getChunkedJobForCurrentStage(self, seqIDs, proportionSampled, inChunkID):
        """
        Give the chunked work to the appropriate job.
        """
        if self.prepOptions.preprocessJob == "checkUniqueHeaders":
            return CheckUniqueHeaders(self.prepOptions, inChunkID)
        elif self.prepOptions.preprocessJob == "lastzRepeatMask":
            repeatMaskOptions = RepeatMaskOptions(proportionSampled=proportionSampled,
                                                  minPeriod=self.prepOptions.minPeriod,
                                                  lastzOpts=self.prepOptions.lastzOptions,
                                                  gpuLastz=self.prepOptions.gpuLastz,
                                                  gpuLastzInterval=self.prepOptions.gpuLastzInterval)
            return LastzRepeatMaskJob(repeatMaskOptions=repeatMaskOptions,
                                      queryID=inChunkID,
                                      targetIDs=seqIDs)
        elif self.prepOptions.preprocessJob == "dna-brnn":
            return DnabrnnMaskJob(inChunkID,
                                  minLength=self.prepOptions.dnabrnnLength,
                                  dnabrnnOpts=self.prepOptions.dnabrnnOpts,
                                  hardmask=self.prepOptions.hardmask)
        else:
            raise RuntimeError("Unknown preprocess job %s" % self.prepOptions.preprocessJob)

    def run(self, fileStore):
        logger.info("Preparing sequence for preprocessing")

        inSequence = fileStore.readGlobalFile(self.inSequenceID)

        if self.prepOptions.chunkSize <= 0:
            # In this first case we don't need to break up the sequence
            chunked = False
            inChunkList = [inSequence]
        else:
            # chunk it up
            chunked = True
            inChunkDirectory = getTempDirectory(rootDir=fileStore.getLocalTempDir())
            inChunkList = runGetChunks(sequenceFiles=[inSequence], chunksDir=inChunkDirectory,
                                       chunkSize=self.prepOptions.chunkSize,
                                       overlapSize=0)
            inChunkList = [os.path.abspath(path) for path in inChunkList]
        logger.info("Chunks = %s" % inChunkList)

        inChunkIDList = [fileStore.writeGlobalFile(chunk, cleanup=True) for chunk in inChunkList]
        outChunkIDList = []
        #For each input chunk we create an output chunk, it is the output chunks that get concatenated together.
        if not self.chunksToCompute:
            self.chunksToCompute = list(range(len(inChunkList)))
        for i in self.chunksToCompute:
            #Calculate the number of chunks to use
            inChunkNumber = int(max(1, math.ceil(len(inChunkList) * self.prepOptions.proportionToSample)))
            assert inChunkNumber <= len(inChunkList) and inChunkNumber > 0
            #Now get the list of chunks flanking and including the current chunk
            j = max(0, i - inChunkNumber//2)
            inChunkIDs = inChunkIDList[j:j+inChunkNumber]
            if len(inChunkIDs) < inChunkNumber: #This logic is like making the list circular
                inChunkIDs += inChunkIDList[:inChunkNumber-len(inChunkIDs)]
            assert len(inChunkIDs) == inChunkNumber
            if self.prepOptions.gpuLastz:
                # when using lastz, we pass through the proportion directly to segalign
                proportionSampled = self.prepOptions.proportionToSample
            else:
                # otherwise, it's taken from the ratio of chunks
                proportionSampled = float(inChunkNumber)/len(inChunkIDList)
            outChunkIDList.append(self.addChild(self.getChunkedJobForCurrentStage(inChunkIDs, proportionSampled, inChunkIDList[i])).rv())

        if chunked:
            # Merge results of the chunking process back into a genome-wide file
            return self.addFollowOn(MergeChunks(self.prepOptions, outChunkIDList)).rv()
        else:
            # Didn't chunk--we have a genome-wide fasta file
            return outChunkIDList[0]

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

        lastIteration = self.iteration == len(self.prepXmlElems) - 1

        prepNode = self.prepXmlElems[self.iteration]
        if getOptionalAttrib(prepNode, "active", typeFn = bool, default=True):
            prepOptions = PreprocessorOptions(chunkSize = int(prepNode.get("chunkSize", default="-1")),
                                              preprocessJob=prepNode.attrib["preprocessJob"],
                                              memory = int(prepNode.get("memory", default=0)),
                                              cpu = int(prepNode.get("cpu", default=1)),
                                              check = bool(int(prepNode.get("check", default="0"))),
                                              proportionToSample = getOptionalAttrib(prepNode, "proportionToSample", typeFn=float, default=1.0),
                                              unmask = getOptionalAttrib(prepNode, "unmask", typeFn=bool, default=False),
                                              lastzOptions = getOptionalAttrib(prepNode, "lastzOpts", default=""),
                                              minPeriod = getOptionalAttrib(prepNode, "minPeriod", typeFn=int, default=0),
                                              checkAssemblyHub = getOptionalAttrib(prepNode, "checkAssemblyHub", typeFn=bool, default=False),
                                              gpuLastz = getOptionalAttrib(prepNode, "gpuLastz", typeFn=bool, default=False),
                                              dnabrnnOpts = getOptionalAttrib(prepNode, "dna-brnnOpts", default=""),
                                              dnabrnnLength = getOptionalAttrib(prepNode, "minLength", typeFn=int, default=1),
                                              hardmask = getOptionalAttrib(prepNode, "hardmask", typeFn=bool, default=False))

            if prepOptions.unmask:
                inSequence = fileStore.readGlobalFile(self.inSequenceID)
                unmaskedInputFile = fileStore.getLocalTempFile()
                unmaskFasta(inSequence, unmaskedInputFile)
                self.inSequenceID = fileStore.writeGlobalFile(unmaskedInputFile)

            outSeqID = self.addChild(PreprocessSequence(prepOptions, self.inSequenceID)).rv()
        else:
            logger.info("Skipping inactive preprocessor {}".format(prepNode.attrib["preprocessJob"]))
            outSeqID = self.inSequenceID

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
        RoundedJob.__init__(self, disk=sum([id.size for id in inputSequenceIDs if hasattr(id, 'size')]), preemptable=True)
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
        return [ os.path.join(outputSequenceDir, inputSequences[i].split("/")[-1] + "_%i" % i) for i in range(len(inputSequences)) ]

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

def stageWorkflow(outputSequenceDir, configFile, inputSequences, toil, restart=False, outputSequences = [], maskAlpha=False):
    #Replace any constants
    configNode = ET.parse(configFile).getroot()
    if not outputSequences:
        outputSequences = CactusPreprocessor.getOutputSequenceFiles(inputSequences, outputSequenceDir)
    else:
        assert len(outputSequences) == len(inputSequences)
    if configNode.find("constants") != None:
        ConfigWrapper(configNode).substituteAllPredefinedConstantsWithLiterals()
    if maskAlpha:
        ConfigWrapper(configNode).setPreprocessorActive("lastzRepeatMask", False)
        ConfigWrapper(configNode).setPreprocessorActive("dna-brnn", True)
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
    parser.add_argument("inSeqFile", type=str, nargs='?', default=None, help = "Input Seq file")
    parser.add_argument("outSeqFile", type=str, nargs='?', default=None, help = "Output Seq file (ex generated with cactus-prepare)")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--inputNames", nargs='*', help='input genome names (not paths) to preprocess (all leaves from Input Seq file if none specified)')
    parser.add_argument("--inPaths", nargs='*', help='Space-separated list of input fasta paths (to be used in place of --inSeqFile')
    parser.add_argument("--outPaths", nargs='*', help='Space-separated list of output fasta paths (one for each inPath, used in place of --outSeqFile)')
    parser.add_argument("--maskAlpha", action='store_true', help='Use dna-brnn instead of lastz for repeatmasking')
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()
    setupBinaries(options)
    setLoggingFromOptions(options)
    enableDumpStack()

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    # we have two modes: operate directly on paths or rely on the seqfiles.  they cannot be mixed
    if options.inSeqFile or options.outSeqFile:
        if not options.inSeqFile or not options.outSeqFile or options.inPaths or options.outPaths:
            raise RuntimeError('--inSeqFile must be used in conjunction with --outSeqFile and not with --inPaths nor --outPaths')
    elif options.inPaths or options.outPaths:
        if not options.inPaths or not options.outPaths or options.inSeqFile or options.outSeqFile or options.inputNames:
            raise RuntimeError('--inPaths must be used in conjunction with --outPaths and not with --inSeqFile, --outSeqFile nor --inputNames')
        if len(options.inPaths) != len(options.outPaths):
            raise RuntimeError('--inPaths and --outPaths must have the same number of arguments')
    else:
        raise RuntimeError('--inSeqFile/--outSeqFile/--inputNames or --inPaths/--outPaths required to specify input')


    inSeqPaths = []
    outSeqPaths = []

    # mine the paths out of the seqfiles
    if options.inSeqFile:
        inSeqFile = SeqFile(options.inSeqFile)
        outSeqFile = SeqFile(options.outSeqFile)

        inNames = options.inputNames
        if not inNames:
            inNames = [inSeqFile.tree.getName(node) for node in inSeqFile.tree.getLeaves()]


        for inName in inNames:
            if inName not in inSeqFile.pathMap or inName not in outSeqFile.pathMap:
                raise RuntimeError('{} not present in input and output Seq files'.format(inNmae))
            inPath = inSeqFile.pathMap[inName]
            outPath = outSeqFile.pathMap[inName]
            if os.path.isdir(inPath):
                try:
                    os.makedirs(outPath)
                except:
                    pass
                assert os.path.isdir(inPath) == os.path.isdir(outPath)
                inSeqPaths += [os.path.join(inPath, seqPath) for seqPath in os.listdir(inPath)]
                outSeqPaths += [os.path.join(outPath, seqPath) for seqPath in os.listdir(inPath)]
            else:
                inSeqPaths += [inPath]
                outSeqPaths += [outPath]

    # we got path names directly from the command line
    else:
        inSeqPaths = options.inPaths
        outSeqPaths = options.outPaths

    with Toil(options) as toil:
        stageWorkflow(outputSequenceDir=None, configFile=options.configFile, inputSequences=inSeqPaths, toil=toil, restart=options.restart, outputSequences=outSeqPaths,
                      maskAlpha=options.maskAlpha)

if __name__ == '__main__':
    main()
