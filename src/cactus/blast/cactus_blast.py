#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
import sys
from argparse import ArgumentParser
from toil.lib.bioio import logger
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString
from toil.lib.bioio import getTempFile

from sonLib.bioio import catFiles, nameValue, popenCatch, getTempDirectory

from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.common import runLastz, runSelfLastz
from cactus.shared.common import runCactusRealign, runCactusSelfRealign
from cactus.shared.common import runGetChunks

class BlastOptions:
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
        
class BlastFlower(Job):
    """Take a reconstruction problem and generate the sequences in chunks to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, cactusSequencesID, flowerName, blastOptions):
        disk = 2*cactusSequencesID.size
        cores = 1
        memory = blastOptions.memory
        
        Job.__init__(self, memory=memory, cores=cores, disk=disk)
        self.cactusDisk = cactusDisk
        self.cactusSequencesID = cactusSequencesID
        self.flowerName = flowerName
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 2
        
    def run(self, fileStore):
        chunksDir = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        cactusSequencesPath = fileStore.readGlobalFile(self.cactusSequencesID)
        chunks = [ chunk for chunk in cactus_call(tool="cactus", check_output=True,
                                                  parameters=["cactus_blast_chunkFlowerSequences"],
                                                  option_string="%s '%s' '%s' %s %i %i %i %s" % \
                                                          (getLogLevelString(), self.cactusDisk, cactusSequencesPath, self.flowerName, 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          self.blastOptions.minimumSequenceLength,
                                                           os.path.basename(chunksDir))).split("\n") if chunk != "" ]
        logger.info("Broken up the flowers into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk, cleanup=False) for chunk in chunks]
        selfResultsID = self.addChild(MakeSelfBlasts(self.blastOptions, chunkIDs)).rv()
        offDiagonalResultsID = self.addChild(MakeOffDiagonalBlasts(blastOptions=self.blastOptions, chunkIDs=chunkIDs)).rv()
        return self.addFollowOn(CollateBlasts(self.blastOptions, [selfResultsID, offDiagonalResultsID])).rv()
    
class BlastSequencesAllAgainstAll(Job):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFileIDs1, blastOptions):
        disk = 4*sum([seqFileID.size for seqFileID in sequenceFileIDs1])
        cores = 1
        memory = blastOptions.memory
        
        Job.__init__(self, disk=disk)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.blastOptions = blastOptions
        self.blastOptions.compressFiles = False
        self.blastOptions.roundsOfCoordinateConversion = 1
        
    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        chunks = runGetChunks(sequenceFiles=sequenceFiles1, chunksDir=getTempDirectory(rootDir=fileStore.getLocalTempDir()), chunkSize = self.blastOptions.chunkSize, overlapSize=self.blastOptions.overlapSize)
        assert len(chunks) > 0
        logger.info("Broken up the sequence files into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk, cleanup=False) for chunk in chunks]

        diagonalResultsID = self.addChild(MakeSelfBlasts(self.blastOptions, chunkIDs)).rv()
        offDiagonalResultsID = self.addChild(MakeOffDiagonalBlasts(self.blastOptions, chunkIDs)).rv()
        logger.debug("Collating the blasts after blasting all-against-all")
        return self.addFollowOn(CollateBlasts(self.blastOptions, [diagonalResultsID, offDiagonalResultsID])).rv()
        
class MakeSelfBlasts(Job):
    """Breaks up the inputs into bits and builds a bunch of alignment jobs.
    """
    def __init__(self, blastOptions, chunkIDs):
        Job.__init__(self)
        self.blastOptions = blastOptions
        self.chunkIDs = chunkIDs
        
    def run(self, fileStore):
        logger.info("Chunk IDs: %s" % self.chunkIDs)
        #Avoid compression if just one chunk
        self.blastOptions.compressFiles = self.blastOptions.compressFiles and len(self.chunkIDs) > 2
        resultsIDs = []
        for i in xrange(len(self.chunkIDs)):
            resultsIDs.append(self.addChild(RunSelfBlast(self.blastOptions, self.chunkIDs[i])).rv())
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        logger.debug("Collating self blasts.")
        logger.info("Blast file IDs: %s" % resultsIDs)
        return self.addFollowOn(CollateBlasts(self.blastOptions, resultsIDs)).rv()
    
class MakeOffDiagonalBlasts(Job):
        def __init__(self, blastOptions, chunkIDs):
            Job.__init__(self)
            self.chunkIDs = chunkIDs
            self.blastOptions = blastOptions
            self.blastOptions.compressFiles = False

        def run(self, fileStore):
            resultsIDs = []
            #Make the list of blast jobs.
            for i in xrange(0, len(self.chunkIDs)):
                for j in xrange(i+1, len(self.chunkIDs)):
                    self.resultsIDs.append(self.addChild(RunBlast(blastOptions=self.blastOptions, seqFileID1=self.chunkIDs[i], seqFileID2=self.chunkIDs[j])).rv())


            return self.addFollowOn(CollateBlasts(self.blastOptions, resultsIDs)).rv()

            
class BlastSequencesAgainstEachOther(Job):
    """Take two sets of sequences, chunks them up and blasts one set against the other.
    """
    def __init__(self, sequenceFileIDs1, sequenceFileIDs2, blastOptions):
        disk = 3*(sum([seqID.size for seqID in sequenceFileIDs1]) + sum([seqID.size for seqID in sequenceFileIDs2]))
        cores = 1
        memory = blastOptions.memory
        
        Job.__init__(self, disk=disk)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.sequenceFileIDs2 = sequenceFileIDs2
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        
    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        sequenceFiles2 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs2]
        chunks1 = runGetChunks(sequenceFiles=sequenceFiles1, chunksDir=getTempDirectory(rootDir=fileStore.getLocalTempDir()), chunkSize=self.blastOptions.chunkSize, overlapSize=self.blastOptions.overlapSize)
        chunks2 = runGetChunks(sequenceFiles=sequenceFiles2, chunksDir=getTempDirectory(rootDir=fileStore.getLocalTempDir()), chunkSize=self.blastOptions.chunkSize, overlapSize=self.blastOptions.overlapSize)
        logger.info("Chunks1 = %s" % chunks1)
        chunkIDs1 = [fileStore.writeGlobalFile(chunk, cleanup=False) for chunk in chunks1]
        chunkIDs2 = [fileStore.writeGlobalFile(chunk, cleanup=False) for chunk in chunks2]
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
        
class BlastIngroupsAndOutgroups(Job):
    """Blast ingroup sequences against each other, and against the given
    outgroup sequences in succession. The next outgroup is only
    aligned against the regions that are not found in the previous
    outgroup.
    """
    def __init__(self, blastOptions, ingroupSequenceIDs,
                 outgroupSequenceIDs):
        Job.__init__(self, memory = blastOptions.memory)
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        self.ingroupSequenceIDs = ingroupSequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs

    def run(self, fileStore):
        fileStore.logToMaster("Blasting ingroups vs outgroups")
        
        ingroupAlignmentsID = self.addChild(BlastSequencesAllAgainstAll(self.ingroupSequenceIDs,
                                                        blastOptions=self.blastOptions)).rv()

        blastFirstOutgroupJob = self.addFollowOn(BlastFirstOutgroup(untrimmedSequenceIDs=self.ingroupSequenceIDs,
                                        sequenceIDs=self.ingroupSequenceIDs,
                                        outgroupSequenceIDs=self.outgroupSequenceIDs,
                                        outgroupFragmentIDs=[],
                                        ingroupResultsID=ingroupAlignmentsID,
                                        outgroupResultsID=None,
                                        blastOptions=self.blastOptions,
                                        outgroupNumber=1,
                                        ingroupCoverageIDs=[]))
        alignmentsID = blastFirstOutgroupJob.rv(0)
        outgroupFragmentIDs = blastFirstOutgroupJob.rv(1)
        ingroupCoverageIDs = blastFirstOutgroupJob.rv(2)
        

        return (alignmentsID, outgroupFragmentIDs, ingroupCoverageIDs)

class BlastFirstOutgroup(Job):
    """Blast the given sequence(s) against the first of a succession of
    outgroups, only aligning fragments that haven't aligned to the
    previous outgroups. Then recurse on the other outgroups.
    """
    def __init__(self, untrimmedSequenceIDs, sequenceIDs,
                 outgroupSequenceIDs, outgroupFragmentIDs,
                 ingroupResultsID, outgroupResultsID,
                 blastOptions, outgroupNumber, ingroupCoverageIDs):
        Job.__init__(self, memory=blastOptions.memory)
        self.untrimmedSequenceIDs = untrimmedSequenceIDs
        self.sequenceIDs = sequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.ingroupResultsID = ingroupResultsID
        self.outgroupResultsID = outgroupResultsID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber
        self.ingroupCoverageIDs = ingroupCoverageIDs

    def run(self, fileStore):
        logger.info("Blasting ingroup sequences to outgroup")
        ingroupAlignmentsID = self.addChild(BlastSequencesAgainstEachOther(self.sequenceIDs,
                                                           [self.outgroupSequenceIDs[0]],
                                                           self.blastOptions)).rv()
        trimRecurseJob = self.addFollowOn(TrimAndRecurseOnOutgroups(untrimmedSequenceIDs=self.untrimmedSequenceIDs,
                                                         sequenceIDs=self.sequenceIDs,
                                                         outgroupSequenceIDs=self.outgroupSequenceIDs,
                                                         outgroupFragmentIDs=self.outgroupFragmentIDs,
                                                         mostRecentResultsID=ingroupAlignmentsID,
                                                         ingroupResultsID=self.ingroupResultsID,
                                                         outgroupResultsID=self.outgroupResultsID,
                                                         blastOptions=self.blastOptions,
                                                         outgroupNumber=self.outgroupNumber,
                                                         ingroupCoverageIDs=self.ingroupCoverageIDs))
        alignmentsID = trimRecurseJob.rv(0)
        outgroupFragmentIDs = trimRecurseJob.rv(1)
        ingroupCoverageIDs = trimRecurseJob.rv(2)
        return (alignmentsID, outgroupFragmentIDs, ingroupCoverageIDs)
        
        
        
        

class TrimAndRecurseOnOutgroups(Job):
    def __init__(self, untrimmedSequenceIDs, sequenceIDs,
                 outgroupSequenceIDs, outgroupFragmentIDs,
                 mostRecentResultsID, ingroupResultsID, outgroupResultsID,
                 blastOptions,
                 outgroupNumber, ingroupCoverageIDs):
        Job.__init__(self)
        self.untrimmedSequenceIDs = untrimmedSequenceIDs
        self.sequenceIDs = sequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.mostRecentResultsID = mostRecentResultsID
        self.ingroupResultsID = ingroupResultsID
        self.outgroupResultsID = outgroupResultsID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber
        self.ingroupCoverageIDs = ingroupCoverageIDs

    def run(self, fileStore):
        # TODO: split up this function, it's getting unwieldy

        # Trim outgroup, convert outgroup coordinates, and add to
        # outgroup fragments dir

        outgroupSequenceFiles = [fileStore.readGlobalFile(fileID) for fileID in self.outgroupSequenceIDs]
        fileStore.logToMaster("Outgroup sequence files: %s" % outgroupSequenceFiles)
        mostRecentResultsFile = fileStore.readGlobalFile(self.mostRecentResultsID)
        trimmedOutgroup = fileStore.getLocalTempFile()
        outgroupCoverage = fileStore.getLocalTempFile()
        calculateCoverage(outgroupSequenceFiles[0],
                          mostRecentResultsFile, outgroupCoverage)
        # The windowSize and threshold are fixed at 1: anything more
        # and we will run into problems with alignments that aren't
        # covered in a matching trimmed sequence.
        trimGenome(outgroupSequenceFiles[0], outgroupCoverage,
                   trimmedOutgroup, flanking=self.blastOptions.trimOutgroupFlanking,
                   windowSize=1, threshold=1)
        outgroupConvertedResultsFile = fileStore.getLocalTempFile()
        cactus_call(tool="cactus", outfile=outgroupConvertedResultsFile,
                    parameters=["cactus_upconvertCoordinates.py",
                                trimmedOutgroup,
                                mostRecentResultsFile,
                                1])

        self.outgroupFragmentIDs.append(fileStore.writeGlobalFile(trimmedOutgroup))
        sequenceFiles = [fileStore.readGlobalFile(path) for path in self.sequenceIDs]
        untrimmedSequenceFiles = [fileStore.readGlobalFile(path) for path in self.untrimmedSequenceIDs]

        # Report coverage of the latest outgroup on the trimmed ingroups.
        for trimmedIngroupSequence, ingroupSequence in zip(sequenceFiles, untrimmedSequenceFiles):
            tmpIngroupCoverage = fileStore.getLocalTempFile()
            calculateCoverage(trimmedIngroupSequence, mostRecentResultsFile,
                              tmpIngroupCoverage)
            fileStore.logToMaster("Coverage on %s from outgroup #%d, %s: %s%% (current ingroup length %d, untrimmed length %d). Outgroup trimmed to %d bp from %d" % (os.path.basename(ingroupSequence), self.outgroupNumber, os.path.basename(outgroupSequenceFiles[0]), percentCoverage(trimmedIngroupSequence, tmpIngroupCoverage), sequenceLength(trimmedIngroupSequence), sequenceLength(ingroupSequence), sequenceLength(trimmedOutgroup), sequenceLength(outgroupSequenceFiles[0])))


        # Convert the alignments' ingroup coordinates.
        ingroupConvertedResultsFile = fileStore.getLocalTempFile()
        if self.sequenceIDs == self.untrimmedSequenceIDs:
            # No need to convert ingroup coordinates on first run.
            fileStore.logToMaster("Copying outgroup results file to ingroup: %s %s" % (outgroupConvertedResultsFile, ingroupConvertedResultsFile))
            system("cp %s %s" % (outgroupConvertedResultsFile,
                                 ingroupConvertedResultsFile))
        else:
            cactus_call(tool="cactus",
                        parameters=["cactus_blast_convertCoordinates",
                                    "--onlyContig1",
                                    outgroupConvertedResultsFile,
                                    ingroupConvertedResultsFile,
                                    1])
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
        for ingroupSequence in untrimmedSequenceFiles:
            ingroupCoverageFile = fileStore.getLocalTempFile()
            calculateCoverage(sequenceFile=ingroupSequence, cigarFile=outgroupResultsFile,
                              outputFile=ingroupCoverageFile, depthById=self.blastOptions.trimOutgroupDepth > 1)
            ingroupCoverageFiles.append(ingroupCoverageFile)
            self.ingroupCoverageIDs.append(fileStore.writeGlobalFile(ingroupCoverageFile))
            fileStore.logToMaster("Cumulative coverage of %d outgroups on ingroup %s: %s" % (self.outgroupNumber, os.path.basename(ingroupSequence), percentCoverage(ingroupSequence, ingroupCoverageFile)))

        # Trim ingroup seqs and recurse on the next outgroup.

        # TODO: Optionally look at coverage on ingroup vs. outgroup,
        # and if coverage is >1 among ingroups but 1 in outgroups,
        # look for it in the next outgroup as well. Would require
        # doing self blast first and sending the alignments here.
        # (Probably needs an extra option in cactus coverage to only
        # count self-alignments, since we need to cut the ingroup
        # sequences in question and not something aligning to both of
        # them.)
        # Could also just ignore the coverage on the outgroup to
        # start, since the fraction of duplicated sequence will be
        # relatively small.
        if self.ingroupResultsID:
            ingroupResultsFile = fileStore.readGlobalFile(self.ingroupResultsID, mutable=True)
        else:
            ingroupResultsFile = fileStore.getLocalTempFile()
        if len(self.outgroupSequenceIDs) > 1:
            trimmedSeqs = []
            # Use the accumulated results so far to trim away the
            # aligned parts of the ingroups.
            for i, sequenceFile in enumerate(untrimmedSequenceFiles):
                outgroupCoverageFile = ingroupCoverageFiles[i]
                selfCoverageFile = fileStore.getLocalTempFile()
                calculateCoverage(sequenceFile, ingroupResultsFile,
                                  selfCoverageFile, fromGenome=sequenceFile)
                self.ingroupResultsID = fileStore.writeGlobalFile(ingroupResultsFile)
                fileStore.logToMaster("Self-coverage on sequence %s: %s%%" % (os.path.basename(sequenceFile), percentCoverage(sequenceFile, selfCoverageFile)))
                coverageFile = fileStore.getLocalTempFile()
                if self.blastOptions.keepParalogs:
                    subtractBed(outgroupCoverageFile, selfCoverageFile, coverageFile)
                else:
                    coverageFile = outgroupCoverageFile

                trimmed = fileStore.getLocalTempFile()
                trimGenome(sequenceFile, coverageFile, trimmed,
                           complement=True, flanking=self.blastOptions.trimFlanking,
                           minSize=self.blastOptions.trimMinSize,
                           threshold=self.blastOptions.trimThreshold,
                           windowSize=self.blastOptions.trimWindowSize,
                           depth=self.blastOptions.trimOutgroupDepth)
                trimmedSeqs.append(trimmed)
            trimmedSeqIDs = [fileStore.writeGlobalFile(path) for path in trimmedSeqs]
            return self.addChild(BlastFirstOutgroup(untrimmedSequenceIDs=self.untrimmedSequenceIDs,
                                                   sequenceIDs=trimmedSeqIDs,
                                                   outgroupSequenceIDs=self.outgroupSequenceIDs[1:],
                                                   outgroupFragmentIDs=self.outgroupFragmentIDs,
                                                   ingroupResultsID=self.ingroupResultsID,
                                                   outgroupResultsID=self.outgroupResultsID,
                                                   blastOptions=self.blastOptions,
                                                   outgroupNumber=self.outgroupNumber + 1,
                                                   ingroupCoverageIDs=self.ingroupCoverageIDs)).rv()
            #self.outgroupResultsID = blastFirstOutgroupJob.rv(0)
            #self.outgroupFragmentIDs = blastFirstOutgroupJob.rv(1)
            #self.ingroupCoverageIDs = blastFirstOutgroupJob.rv(2)
            
        else:
            # Finally, put the ingroups and outgroups results together
            alignmentsID = self.addFollowOn(CollateBlasts(blastOptions=self.blastOptions, resultsFileIDs=[self.ingroupResultsID,
                                                  self.outgroupResultsID])).rv()
            return (alignmentsID, self.outgroupFragmentIDs, self.ingroupCoverageIDs)
def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    system("bzip2 --keep --fast %s" % fileName)
    return fileName + ".bz2"
        
class RunSelfBlast(Job):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID):
        disk = 3*seqFileID.size
        memory = 3*seqFileID.size
        
        Job.__init__(self, memory=memory, disk=disk)
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
        cactus_call(tool="cactus",
                    parameters=["cactus_blast_convertCoordinates",
                                blastResultsFile,
                                resultsFile,
                                self.blastOptions.roundsOfCoordinateConversion])
        if self.blastOptions.compressFiles:
            #TODO: This throws away the compressed file
            seqFile = compressFastaFile(seqFile)
        logger.info("Ran the self blast okay")
        return fileStore.writeGlobalFile(resultsFile, cleanup=False)

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    system("bunzip2 --stdout %s > %s" % (fileName, tempFileName))
    return tempFileName
    
class RunBlast(Job):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID1, seqFileID2):
        if hasattr(seqFileID1, "size") and hasattr(seqFileID2, "size"):
            disk = 2*(seqFileID1.size + seqFileID2.size)
            memory = 2*(seqFileID1.size + seqFileID2.size)
        else:
            disk = None
            memory = None
        Job.__init__(self, memory=memory, disk=disk)
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
        cactus_call(tool="cactus",
                    parameters=["cactus_blast_convertCoordinates",
                                blastResultsFile,
                                resultsFile,
                                self.blastOptions.roundsOfCoordinateConversion])
        logger.info("Ran the blast okay")
        return fileStore.writeGlobalFile(resultsFile, cleanup=False)

class CollateBlasts(Job):
    def __init__(self, blastOptions, resultsFileIDs):
        Job.__init__(self)
        self.blastOptions = blastOptions
        self.resultsFileIDs = resultsFileIDs

    def run(self, fileStore):
        return self.addFollowOn(CollateBlasts2(self.blastOptions, self.resultsFileIDs)).rv()
class CollateBlasts2(Job):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, blastOptions, resultsFileIDs):
        disk = 8*sum([alignmentID.size for alignmentID in resultsFileIDs])
        cores = 1
        memory = blastOptions.memory
        Job.__init__(self, memory = memory, disk=disk)
        self.resultsFileIDs = resultsFileIDs
    
    def run(self, fileStore):
        logger.info("Results IDs: %s" % self.resultsFileIDs)
        resultsFiles = [fileStore.readGlobalFile(fileID) for fileID in self.resultsFileIDs]
        collatedResultsFile = fileStore.getLocalTempFile()
        catFiles(resultsFiles, collatedResultsFile)
        logger.info("Collated the alignments to the file: %s",  collatedResultsFile)
        collatedResultsID = fileStore.writeGlobalFile(collatedResultsFile, cleanup=False)
        return collatedResultsID
        
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
    fromGenome = nameValue("from", fromGenome).split()
    cactus_call(tool="cactus", outfile=outputFile, work_dir=work_dir,
                parameters=["cactus_coverage",
                            sequenceFile,
                            cigarFile] +
                            fromGenome +
                            [nameValue("depthById", depthById, bool)])

def trimGenome(sequenceFile, coverageFile, outputFile, complement=False,
               flanking=0, minSize=1, windowSize=10, threshold=1, depth=None):
    cactus_call(tool="cactus", outfile=outputFile,
                parameters=["cactus_trimSequences.py",
                            nameValue("complement", complement, valueType=bool),
                            nameValue("flanking", flanking), nameValue("minSize", minSize),
                            nameValue("windowSize", windowSize), nameValue("threshold", threshold),
                            nameValue("depth", depth), sequenceFile, coverageFile])

def subtractBed(bed1, bed2, destBed):
    """Subtract two non-bed12 beds"""
    # tmp. don't really want to use bedtools
    if os.path.getsize(bed1) == 0 or os.path.getsize(bed2) == 0:
        # bedtools will complain on zero-size beds
        os.rename(bed1, destBed)
    else:
        cactus_call(tool="bedtools", outfile=destBed,
                    parameters=["subtract",
                                "-a", bed1,
                                "-b", bed2])

def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    blastOptions = BlastOptions()
    
    #output stuff
    parser.add_option("--cigars", dest="cigarFile", 
                      help="File to write cigars in",
                      default="cigarFile.txt")
    
    parser.add_option("--chunkSize", dest="chunkSize", type="int", 
                     help="The size of chunks passed to lastz (must be at least twice as big as overlap)",
                     default=blastOptions.chunkSize)
    
    parser.add_option("--overlapSize", dest="overlapSize", type="int",
                     help="The size of the overlap between the chunks passed to lastz (min size 2)",
                     default=blastOptions.overlapSize)
    
    parser.add_option("--blastString", dest="blastString", type="string",
                     help="The default string used to call the blast program. \
Must contain three strings: SEQ_FILE_1, SEQ_FILE_2 and CIGARS_FILE which will be \
replaced with the two sequence files and the results file, respectively",
                     default=blastOptions.blastString)
    
    parser.add_option("--selfBlastString", dest="selfBlastString", type="string",
                     help="The default string used to call the blast program for self alignment. \
Must contain three strings: SEQ_FILE and CIGARS_FILE which will be \
replaced with the the sequence file and the results file, respectively",
                     default=blastOptions.selfBlastString)
   
    parser.add_option("--compressFiles", dest="compressFiles", action="store_false",
                      help="Turn of bz2 based file compression of sequences for I/O transfer", 
                      default=blastOptions.compressFiles)
    
    parser.add_option("--lastzMemory", dest="memory", type="int",
                      help="Lastz memory (in bytes)", 
                      default=blastOptions.memory)
    
    parser.add_option("--trimFlanking", type=int, help="Amount of flanking sequence to leave on trimmed ingroup sequences", default=blastOptions.trimFlanking)
    parser.add_option("--trimMinSize", type=int, help="Minimum size, before adding flanking sequence, of ingroup sequence to align against the next outgroup", default=blastOptions.trimMinSize)
    parser.add_option("--trimThreshold", type=int, help="Coverage threshold for an ingroup region to not be aligned against the next outgroup", default=blastOptions.trimThreshold)
    parser.add_option("--trimWindowSize", type=int, help="Windowing size to integrate ingroup coverage over", default=blastOptions.trimWindowSize)
    parser.add_option("--trimOutgroupFlanking", type=int, help="Amount of flanking sequence to leave on trimmed outgroup sequences", default=blastOptions.trimOutgroupFlanking)
    parser.add_option("--trimOutgroupDepth", type=int, help="Trim regions away "
                      "after this many separate outgroups have been aligned to "
                      "the region", default=1)
    
    parser.add_option("--keepParalogs", action="store_true", help="Never trim away any ingroup sequence that aligns somewhere else in the ingroup", default=blastOptions.keepParalogs)

    parser.add_option("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")
    
    parser.add_option("--targetSequenceFiles", dest="targetSequenceFiles", type="string",
                     help="Sequences to align against the input sequences against. If these are not provided then the input sequences are aligned against each other.",
                     default=None)

    parser.add_option("--ingroups", type=str, default=None,
                      help="Ingroups to align (comma-separated) (--outgroups "
                      "must be provided as well")

    parser.add_option("--outgroups", type=str, default=None,
                      help="Outgroups to align (comma-separated) (--ingroups "
                      "must be provided as well")

    parser.add_option("--outgroupFragmentsDir", type=str,
                      default="outgroupFragments/", help= "Directory to "
                      "store outgroup fragments in")

    parser.add_option("--ingroupCoverageDir", type=str,
                      help="Directory to store ingroup coverage beds in "
                      "(only works in ingroups vs. outgroups mode")

    options, args = parser.parse_args()
    if options.test:
        _test()

    if (options.ingroups is not None) ^ (options.outgroups is not None):
        raise RuntimeError("--ingroups and --outgroups must be provided "
                           "together")
    if options.ingroups:
        firstTarget = BlastIngroupsAndOutgroups(options,
                                                options.ingroups.split(','),
                                                options.outgroups.split(','),
                                                options.cigarFile,
                                                options.outgroupFragmentsDir,
                                                options.ingroupCoverageDir)
    elif options.targetSequenceFiles == None:
        firstTarget = BlastSequencesAllAgainstAll(args, options.cigarFile, options)
    else:
        firstTarget = BlastSequencesAgainstEachOther(args, options.targetSequenceFiles.split(), options.cigarFile, options)
    Stack(firstTarget).startJobTree(options)

def _test():
    import doctest
    return doctest.testmod()
