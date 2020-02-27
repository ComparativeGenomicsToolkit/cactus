#!/usr/bin/env python3
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the toil framework to parallelise the blasts.
"""
import os
import shutil
from toil.lib.bioio import logger
from toil.lib.bioio import system

from sonLib.bioio import catFiles, nameValue, popenCatch, getTempDirectory

from cactus.shared.common import RoundedJob
from cactus.shared.common import cactus_call
from cactus.shared.common import runLastz, runSelfLastz
from cactus.shared.common import runCactusRealign, runCactusSelfRealign
from cactus.shared.common import runGetChunks
from cactus.shared.common import readGlobalFileWithoutCache
from cactus.shared.common import ChildTreeJob
from cactus.blast.upconvertCoordinates import upconvertCoords
from cactus.blast.trimSequences import trimSequences

class BlastOptions(object):
    def __init__(self, chunkSize=10000000, overlapSize=10000,
                 lastzArguments="", compressFiles=True, realign=False, realignArguments="",
                 minimumSequenceLength=1, memory=None,
                 smallDisk = None,
                 largeDisk = None,
                 # Trim options for trimming ingroup seqs:
                 trimFlanking=10, trimMinSize=20,
                 trimWindowSize=10, trimThreshold=1,
                 trimOutgroupDepth=1,
                 # Trim options for trimming outgroup seqs (options
                 # other than flanking sequence will need a check to
                 # remove alignments that don't qualify)
                 # HACK: outgroup flanking is only set so high by
                 # default because it's needed for the tests (which
                 # don't use realign.)
                 trimOutgroupFlanking=2000,
                 keepParalogs=False):
        """Class defining options for blast
        """
        self.chunkSize = chunkSize
        self.overlapSize = overlapSize

        self.realignArguments = realignArguments
        self.lastzArguments = lastzArguments
        self.realign = realign

        self.compressFiles = compressFiles
        self.minimumSequenceLength = 10
        self.memory = memory
        self.smallDisk = smallDisk
        self.largeDisk = largeDisk
        self.trimFlanking = trimFlanking
        self.trimMinSize = trimMinSize
        self.trimThreshold = trimThreshold
        self.trimWindowSize = trimWindowSize
        self.trimOutgroupDepth = trimOutgroupDepth
        self.trimOutgroupFlanking = trimOutgroupFlanking
        self.keepParalogs = keepParalogs

class BlastSequencesAllAgainstAll(RoundedJob):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFileIDs1, blastOptions):
        disk = 4*sum([seqFileID.size for seqFileID in sequenceFileIDs1])
        cores = 1
        memory = blastOptions.memory

        super(BlastSequencesAllAgainstAll, self).__init__(disk=disk, cores=cores, memory=memory, preemptable=True)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.blastOptions = blastOptions
        self.blastOptions.compressFiles = False
        self.blastOptions.roundsOfCoordinateConversion = 1

    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        chunks = runGetChunks(sequenceFiles=sequenceFiles1,
                              chunksDir=getTempDirectory(rootDir=fileStore.getLocalTempDir()),
                              chunkSize=self.blastOptions.chunkSize, overlapSize=self.blastOptions.overlapSize)
        if len(chunks) == 0:
            raise Exception("no chunks produced for files: {} ".format(sequenceFiles1))
        logger.info("Broken up the sequence files into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk, cleanup=True) for chunk in chunks]

        diagonalResultsID = self.addChild(MakeSelfBlasts(self.blastOptions, chunkIDs)).rv()
        offDiagonalResultsID = self.addChild(MakeOffDiagonalBlasts(self.blastOptions, chunkIDs)).rv()
        logger.debug("Collating the blasts after blasting all-against-all")
        return self.addFollowOn(CollateBlasts(self.blastOptions, [diagonalResultsID, offDiagonalResultsID])).rv()

class MakeSelfBlasts(ChildTreeJob):
    """Breaks up the inputs into bits and builds a bunch of alignment jobs.
    """
    def __init__(self, blastOptions, chunkIDs):
        super(MakeSelfBlasts, self).__init__(preemptable=True)
        self.blastOptions = blastOptions
        self.chunkIDs = chunkIDs

    def run(self, fileStore):
        logger.info("Chunk IDs: %s" % self.chunkIDs)
        #Avoid compression if just one chunk
        self.blastOptions.compressFiles = self.blastOptions.compressFiles and len(self.chunkIDs) > 2
        resultsIDs = []
        for i in range(len(self.chunkIDs)):
            resultsIDs.append(self.addChild(RunSelfBlast(self.blastOptions, self.chunkIDs[i])).rv())
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        logger.debug("Collating self blasts.")
        logger.info("Blast file IDs: %s" % resultsIDs)
        return self.addFollowOn(CollateBlasts(self.blastOptions, resultsIDs)).rv()

class MakeOffDiagonalBlasts(ChildTreeJob):
        def __init__(self, blastOptions, chunkIDs):
            super(MakeOffDiagonalBlasts, self).__init__(preemptable=True)
            self.chunkIDs = chunkIDs
            self.blastOptions = blastOptions
            self.blastOptions.compressFiles = False

        def run(self, fileStore):
            resultsIDs = []
            #Make the list of blast jobs.
            for i in range(0, len(self.chunkIDs)):
                for j in range(i+1, len(self.chunkIDs)):
                    resultsIDs.append(self.addChild(RunBlast(blastOptions=self.blastOptions, seqFileID1=self.chunkIDs[i], seqFileID2=self.chunkIDs[j])).rv())

            return self.addFollowOn(CollateBlasts(self.blastOptions, resultsIDs)).rv()

class BlastSequencesAgainstEachOther(ChildTreeJob):
    """Take two sets of sequences, chunks them up and blasts one set against the other.
    """
    def __init__(self, sequenceFileIDs1, sequenceFileIDs2, blastOptions):
        disk = 3*(sum([seqID.size for seqID in sequenceFileIDs1]) + sum([seqID.size for seqID in sequenceFileIDs2]))
        cores = 1
        memory = blastOptions.memory

        super(BlastSequencesAgainstEachOther, self).__init__(disk=disk, cores=cores, memory=memory, preemptable=True)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.sequenceFileIDs2 = sequenceFileIDs2
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1

    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        sequenceFiles2 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs2]
        chunks1 = runGetChunks(sequenceFiles=sequenceFiles1, chunksDir=getTempDirectory(rootDir=fileStore.getLocalTempDir()), chunkSize=self.blastOptions.chunkSize, overlapSize=self.blastOptions.overlapSize)
        chunks2 = runGetChunks(sequenceFiles=sequenceFiles2, chunksDir=getTempDirectory(rootDir=fileStore.getLocalTempDir()), chunkSize=self.blastOptions.chunkSize, overlapSize=self.blastOptions.overlapSize)
        chunkIDs1 = [fileStore.writeGlobalFile(chunk, cleanup=True) for chunk in chunks1]
        chunkIDs2 = [fileStore.writeGlobalFile(chunk, cleanup=True) for chunk in chunks2]
        resultsIDs = []
        #Make the list of blast jobs.
        for chunkID1 in chunkIDs1:
            for chunkID2 in chunkIDs2:
                #TODO: Make the compression work
                self.blastOptions.compressFiles = False
                resultsIDs.append(self.addChild(RunBlast(self.blastOptions, chunkID1, chunkID2)).rv())
        logger.info("Made the list of blasts")
        #Set up the job to collate all the results
        return self.addFollowOn(CollateBlasts(self.blastOptions, resultsIDs)).rv()

class BlastIngroupsAndOutgroups(RoundedJob):
    """Blast ingroup sequences against each other, and against the given
    outgroup sequences in succession. The next outgroup is only
    aligned against the regions that are not found in the previous
    outgroup.
    """
    def __init__(self, blastOptions, ingroupNames, ingroupSequenceIDs,
                 outgroupNames, outgroupSequenceIDs):
        super(BlastIngroupsAndOutgroups, self).__init__(memory=blastOptions.memory, preemptable=True)
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        self.ingroupNames = ingroupNames
        self.outgroupNames = outgroupNames
        self.ingroupSequenceIDs = ingroupSequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs

    def run(self, fileStore):
        fileStore.logToMaster("Blasting ingroups vs outgroups. "
                              "Ingroup genomes: %s, outgroup genomes: %s" \
                              % (", ".join(self.ingroupNames), ", ".join(self.outgroupNames)))

        ingroupAlignmentsID = self.addChild(BlastSequencesAllAgainstAll(self.ingroupSequenceIDs,
                                                        blastOptions=self.blastOptions)).rv()
        if len(self.outgroupSequenceIDs) > 0:
            blastFirstOutgroupJob = self.addChild(BlastFirstOutgroup(
                ingroupNames=self.ingroupNames,
                untrimmedSequenceIDs=self.ingroupSequenceIDs,
                sequenceIDs=self.ingroupSequenceIDs,
                outgroupNames=self.outgroupNames,
                outgroupSequenceIDs=self.outgroupSequenceIDs,
                outgroupFragmentIDs=[],
                outgroupResultsID=None,
                blastOptions=self.blastOptions,
                outgroupNumber=1,
                ingroupCoverageIDs=[]))
            outgroupAlignmentsID = blastFirstOutgroupJob.rv(0)
            outgroupFragmentIDs = blastFirstOutgroupJob.rv(1)
            ingroupCoverageIDs = blastFirstOutgroupJob.rv(2)
            alignmentsID = self.addFollowOn(CollateBlasts(blastOptions=self.blastOptions, resultsFileIDs=[ingroupAlignmentsID, outgroupAlignmentsID])).rv()
        else:
            alignmentsID = ingroupAlignmentsID
            outgroupFragmentIDs = []
            ingroupCoverageIDs = []

        return (alignmentsID, outgroupFragmentIDs, ingroupCoverageIDs)

class BlastFirstOutgroup(RoundedJob):
    """Blast the given sequence(s) against the first of a succession of
    outgroups, only aligning fragments that haven't aligned to the
    previous outgroups. Then recurse on the other outgroups.
    """
    def __init__(self, ingroupNames, untrimmedSequenceIDs, sequenceIDs,
                 outgroupNames, outgroupSequenceIDs, outgroupFragmentIDs,
                 outgroupResultsID, blastOptions, outgroupNumber,
                 ingroupCoverageIDs):
        super(BlastFirstOutgroup, self).__init__(memory=blastOptions.memory, preemptable=True)
        self.ingroupNames = ingroupNames
        self.untrimmedSequenceIDs = untrimmedSequenceIDs
        self.sequenceIDs = sequenceIDs
        self.outgroupNames = outgroupNames
        self.outgroupSequenceIDs = outgroupSequenceIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.outgroupResultsID = outgroupResultsID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber
        self.ingroupCoverageIDs = ingroupCoverageIDs

    def run(self, fileStore):
        logger.info("Blasting ingroup sequences to outgroup %s",
                    self.outgroupNames[self.outgroupNumber - 1])
        alignmentsID = self.addChild(BlastSequencesAgainstEachOther(
            self.sequenceIDs,
            [self.outgroupSequenceIDs[0]],
            self.blastOptions)).rv()
        trimRecurseJob = self.addFollowOn(TrimAndRecurseOnOutgroups(
            ingroupNames=self.ingroupNames,
            untrimmedSequenceIDs=self.untrimmedSequenceIDs,
            sequenceIDs=self.sequenceIDs,
            outgroupNames=self.outgroupNames,
            outgroupSequenceIDs=self.outgroupSequenceIDs,
            outgroupFragmentIDs=self.outgroupFragmentIDs,
            mostRecentResultsID=alignmentsID,
            outgroupResultsID=self.outgroupResultsID,
            blastOptions=self.blastOptions,
            outgroupNumber=self.outgroupNumber,
            ingroupCoverageIDs=self.ingroupCoverageIDs))
        outgroupAlignmentsID = trimRecurseJob.rv(0)
        outgroupFragmentIDs = trimRecurseJob.rv(1)
        ingroupCoverageIDs = trimRecurseJob.rv(2)
        return (outgroupAlignmentsID, outgroupFragmentIDs, ingroupCoverageIDs)

class TrimAndRecurseOnOutgroups(RoundedJob):
    def __init__(self, ingroupNames, untrimmedSequenceIDs, sequenceIDs,
                 outgroupNames, outgroupSequenceIDs, outgroupFragmentIDs,
                 mostRecentResultsID, outgroupResultsID,
                 blastOptions, outgroupNumber, ingroupCoverageIDs):
        super(TrimAndRecurseOnOutgroups, self).__init__(preemptable=True)
        self.ingroupNames = ingroupNames
        self.untrimmedSequenceIDs = untrimmedSequenceIDs
        self.sequenceIDs = sequenceIDs
        self.outgroupNames = outgroupNames
        self.outgroupSequenceIDs = outgroupSequenceIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.mostRecentResultsID = mostRecentResultsID
        self.outgroupResultsID = outgroupResultsID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber
        self.ingroupCoverageIDs = ingroupCoverageIDs

    def run(self, fileStore):
        # Trim outgroup, convert outgroup coordinates, and add to
        # outgroup fragments dir

        outgroupSequenceFiles = [fileStore.readGlobalFile(fileID) for fileID in self.outgroupSequenceIDs]
        mostRecentResultsFile = fileStore.readGlobalFile(self.mostRecentResultsID)
        trimmedOutgroup = fileStore.getLocalTempFile()
        outgroupCoverage = fileStore.getLocalTempFile()
        calculateCoverage(outgroupSequenceFiles[0],
                          mostRecentResultsFile, outgroupCoverage)
        # The windowSize and threshold are fixed at 1: anything more
        # and we will run into problems with alignments that aren't
        # covered in a matching trimmed sequence.
        trimSequences(outgroupSequenceFiles[0], outgroupCoverage,
                      trimmedOutgroup, flanking=self.blastOptions.trimOutgroupFlanking,
                      windowSize=1, threshold=1)
        outgroupConvertedResultsFile = fileStore.getLocalTempFile()
        with open(outgroupConvertedResultsFile, 'w') as f:
            upconvertCoords(cigarPath=mostRecentResultsFile,
                            fastaPath=trimmedOutgroup,
                            contigNum=1,
                            outputFile=f)

        self.outgroupFragmentIDs.append(fileStore.writeGlobalFile(trimmedOutgroup))
        sequenceFiles = [fileStore.readGlobalFile(path) for path in self.sequenceIDs]
        untrimmedSequenceFiles = [fileStore.readGlobalFile(path) for path in self.untrimmedSequenceIDs]

        # Report coverage of the latest outgroup on the trimmed ingroups.
        for trimmedIngroupSequence, ingroupSequence, ingroupName in zip(sequenceFiles, untrimmedSequenceFiles, self.ingroupNames):
            tmpIngroupCoverage = fileStore.getLocalTempFile()
            calculateCoverage(trimmedIngroupSequence, mostRecentResultsFile,
                              tmpIngroupCoverage)
            fileStore.logToMaster("Coverage on %s from outgroup #%d, %s: %s%% (current ingroup length %d, untrimmed length %d). Outgroup trimmed to %d bp from %d" % (ingroupName, self.outgroupNumber, self.outgroupNames[self.outgroupNumber - 1], percentCoverage(trimmedIngroupSequence, tmpIngroupCoverage), sequenceLength(trimmedIngroupSequence), sequenceLength(ingroupSequence), sequenceLength(trimmedOutgroup), sequenceLength(outgroupSequenceFiles[0])))

        # Convert the alignments' ingroup coordinates.
        ingroupConvertedResultsFile = fileStore.getLocalTempFile()
        if self.sequenceIDs == self.untrimmedSequenceIDs:
            # No need to convert ingroup coordinates on first run.
            shutil.copy(outgroupConvertedResultsFile,
                        ingroupConvertedResultsFile)
        else:
            cactus_call(parameters=["cactus_blast_convertCoordinates",
                                    "--onlyContig1",
                                    outgroupConvertedResultsFile,
                                    ingroupConvertedResultsFile,
                                    "1"])
        # Append the latest results to the accumulated outgroup coverage file
        if self.outgroupResultsID:
            outgroupResultsFile = fileStore.readGlobalFile(self.outgroupResultsID, mutable=True)
        else:
            outgroupResultsFile = fileStore.getLocalTempFile()
        with open(ingroupConvertedResultsFile) as results:
            with open(outgroupResultsFile, 'a') as output:
                output.write(results.read())

        self.outgroupResultsID = fileStore.writeGlobalFile(outgroupResultsFile)

        # Report coverage of the all outgroup alignments so far on the ingroups.
        ingroupCoverageFiles = []
        self.ingroupCoverageIDs = []
        for ingroupSequence, ingroupName in zip(untrimmedSequenceFiles, self.ingroupNames):
            ingroupCoverageFile = fileStore.getLocalTempFile()
            calculateCoverage(sequenceFile=ingroupSequence, cigarFile=outgroupResultsFile,
                              outputFile=ingroupCoverageFile, depthById=self.blastOptions.trimOutgroupDepth > 1)
            ingroupCoverageFiles.append(ingroupCoverageFile)
            self.ingroupCoverageIDs.append(fileStore.writeGlobalFile(ingroupCoverageFile))
            fileStore.logToMaster("Cumulative coverage of %d outgroups on ingroup %s: %s" % (self.outgroupNumber, ingroupName, percentCoverage(ingroupSequence, ingroupCoverageFile)))

        if len(self.outgroupSequenceIDs) > 1:
            # Trim ingroup seqs and recurse on the next outgroup.
            trimmedSeqs = []
            # Use the accumulated results so far to trim away the
            # aligned parts of the ingroups.
            for i, sequenceFile in enumerate(untrimmedSequenceFiles):
                outgroupCoverageFile = ingroupCoverageFiles[i]
                selfCoverageFile = fileStore.getLocalTempFile()
                coverageFile = fileStore.getLocalTempFile()
                if self.blastOptions.keepParalogs:
                    subtractBed(outgroupCoverageFile, selfCoverageFile, coverageFile)
                else:
                    coverageFile = outgroupCoverageFile

                trimmed = fileStore.getLocalTempFile()
                trimSequences(sequenceFile, coverageFile, trimmed,
                              complement=True, flanking=self.blastOptions.trimFlanking,
                              minSize=self.blastOptions.trimMinSize,
                              threshold=self.blastOptions.trimThreshold,
                              windowSize=self.blastOptions.trimWindowSize,
                              depth=self.blastOptions.trimOutgroupDepth)
                trimmedSeqs.append(trimmed)
            trimmedSeqIDs = [fileStore.writeGlobalFile(path, cleanup=True) for path in trimmedSeqs]
            return self.addChild(BlastFirstOutgroup(
                ingroupNames=self.ingroupNames,
                untrimmedSequenceIDs=self.untrimmedSequenceIDs,
                sequenceIDs=trimmedSeqIDs,
                outgroupNames=self.outgroupNames,
                outgroupSequenceIDs=self.outgroupSequenceIDs[1:],
                outgroupFragmentIDs=self.outgroupFragmentIDs,
                outgroupResultsID=self.outgroupResultsID,
                blastOptions=self.blastOptions,
                outgroupNumber=self.outgroupNumber + 1,
                ingroupCoverageIDs=self.ingroupCoverageIDs)).rv()
        else:
            # Finally, put the ingroups and outgroups results together
            return (self.outgroupResultsID, self.outgroupFragmentIDs, self.ingroupCoverageIDs)

def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    system("bzip2 --keep --fast %s" % fileName)
    return fileName + ".bz2"

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    system("bunzip2 --stdout %s > %s" % (fileName, tempFileName))
    return tempFileName

class RunSelfBlast(RoundedJob):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID):
        disk = 3*seqFileID.size
        memory = 3*seqFileID.size

        super(RunSelfBlast, self).__init__(memory=memory, disk=disk, preemptable=True)
        self.blastOptions = blastOptions
        self.seqFileID = seqFileID

    def run(self, fileStore):
        blastResultsFile = fileStore.getLocalTempFile()
        seqFile = fileStore.readGlobalFile(self.seqFileID)
        runSelfLastz(seqFile, blastResultsFile, lastzArguments=self.blastOptions.lastzArguments)
        if self.blastOptions.realign:
            realignResultsFile = fileStore.getLocalTempFile()
            runCactusSelfRealign(seqFile, inputAlignmentsFile=blastResultsFile,
                                 outputAlignmentsFile=realignResultsFile,
                                 realignArguments=self.blastOptions.realignArguments)
            blastResultsFile = realignResultsFile
        resultsFile = fileStore.getLocalTempFile()
        cactus_call(parameters=["cactus_blast_convertCoordinates",
                                blastResultsFile,
                                resultsFile,
                                str(self.blastOptions.roundsOfCoordinateConversion)])
        if self.blastOptions.compressFiles:
            #TODO: This throws away the compressed file
            seqFile = compressFastaFile(seqFile)
        logger.info("Ran the self blast okay")
        return fileStore.writeGlobalFile(resultsFile)

class RunBlast(RoundedJob):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID1, seqFileID2):
        if hasattr(seqFileID1, "size") and hasattr(seqFileID2, "size"):
            disk = 2*(seqFileID1.size + seqFileID2.size)
            memory = 2*(seqFileID1.size + seqFileID2.size)
        else:
            disk = None
            memory = None
        super(RunBlast, self).__init__(memory=memory, disk=disk, preemptable=True)
        self.blastOptions = blastOptions
        self.seqFileID1 = seqFileID1
        self.seqFileID2 = seqFileID2

    def run(self, fileStore):
        seqFile1 = fileStore.readGlobalFile(self.seqFileID1)
        seqFile2 = fileStore.readGlobalFile(self.seqFileID2)
        if self.blastOptions.compressFiles:
            seqFile1 = decompressFastaFile(seqFile1, fileStore.getLocalTempFile())
            seqFile2 = decompressFastaFile(seqFile2, fileStore.getLocalTempFile())
        blastResultsFile = fileStore.getLocalTempFile()

        runLastz(seqFile1, seqFile2, blastResultsFile, lastzArguments = self.blastOptions.lastzArguments)
        if self.blastOptions.realign:
            realignResultsFile = fileStore.getLocalTempFile()
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile=blastResultsFile,
                             outputAlignmentsFile=realignResultsFile,
                             realignArguments=self.blastOptions.realignArguments)
            blastResultsFile = realignResultsFile

        resultsFile = fileStore.getLocalTempFile()
        cactus_call(parameters=["cactus_blast_convertCoordinates",
                                blastResultsFile,
                                resultsFile,
                                str(self.blastOptions.roundsOfCoordinateConversion)])
        logger.info("Ran the blast okay")
        return fileStore.writeGlobalFile(resultsFile)

class CollateBlasts(RoundedJob):
    def __init__(self, blastOptions, resultsFileIDs):
        super(CollateBlasts, self).__init__(preemptable=True)
        self.blastOptions = blastOptions
        self.resultsFileIDs = resultsFileIDs

    def run(self, fileStore):
        return self.addFollowOn(CollateBlasts2(self.blastOptions, self.resultsFileIDs)).rv()

class CollateBlasts2(ChildTreeJob):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, blastOptions, resultsFileIDs):
        disk = 8*sum([alignmentID.size for alignmentID in resultsFileIDs])
        memory = blastOptions.memory
        super(CollateBlasts2, self).__init__(memory=memory, disk=disk, preemptable=True)
        self.resultsFileIDs = resultsFileIDs
        # it's slow to run fileStore.deleteGlobalFile, so we do it in parallel batches
        self.delete_batch_size = 50

    def run(self, fileStore):
        logger.info("Results IDs: %s" % self.resultsFileIDs)
        resultsFiles = [readGlobalFileWithoutCache(fileStore, fileID) for fileID in self.resultsFileIDs]
        collatedResultsFile = fileStore.getLocalTempFile()
        catFiles(resultsFiles, collatedResultsFile)
        logger.info("Collated the alignments to the file: %s",  collatedResultsFile)
        collatedResultsID = fileStore.writeGlobalFile(collatedResultsFile)
        for i in range(0, len(self.resultsFileIDs), self.delete_batch_size):
            self.addChild(DeleteFileIDs(self.resultsFileIDs[i:i+self.delete_batch_size]))        
        return collatedResultsID

class DeleteFileIDs(RoundedJob):
    """Deletes some files from the file store
    """
    def __init__(self, fileIDs):
        super(DeleteFileIDs, self).__init__(preemptable=True)
        self.fileIDs = fileIDs
        
    def run(self, fileStore):
        for fileID in self.fileIDs:
            fileStore.deleteGlobalFile(fileID)        

def sequenceLength(sequenceFile):
    """Get the total # of bp from a fasta file."""
    seqLength = 0
    for line in open(sequenceFile):
        line = line.strip()
        if line == '' or line[0] == '>':
            continue
        seqLength += len(line)
    return seqLength

def percentCoverage(sequenceFile, coverageFile):
    """Get the % coverage of a sequence from a coverage file."""
    sequenceLen = sequenceLength(sequenceFile)
    if sequenceLen == 0:
        return 0
    coverage = popenCatch("awk '{ total += $3 - $2 } END { print total }' %s" % coverageFile)
    if coverage.strip() == '': # No coverage lines
        return 0
    return 100*float(coverage)/sequenceLen

def calculateCoverage(sequenceFile, cigarFile, outputFile, fromGenome=None, depthById=False, work_dir=None):
    logger.info("Calculating coverage of cigar file %s on %s, writing to %s" % (
        cigarFile, sequenceFile, outputFile))
    args = [sequenceFile, cigarFile]
    if fromGenome is not None:
        args += ["--from", fromGenome]
    if depthById:
        args += ["--depthById"]
    cactus_call(outfile=outputFile, work_dir=work_dir,
                parameters=["cactus_coverage"] + args)

def subtractBed(bed1, bed2, destBed):
    """Subtract two non-bed12 beds"""
    # tmp. don't really want to use bedtools
    if os.path.getsize(bed1) == 0 or os.path.getsize(bed2) == 0:
        # bedtools will complain on zero-size beds
        os.rename(bed1, destBed)
    else:
        cactus_call(outfile=destBed,
                    parameters=["subtract",
                                "-a", bed1,
                                "-b", bed2])
