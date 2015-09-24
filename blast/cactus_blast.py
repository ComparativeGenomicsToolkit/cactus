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
from sonLib.bioio import TempFileTree
from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch
from sonLib.bioio import getLogLevelString
from sonLib.bioio import makeSubDir
from sonLib.bioio import catFiles
from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import nameValue
from toil.job import Job
from cactus.shared.common import WritePermanentFile


class BlastOptions:
    def __init__(self, chunkSize=10000000, overlapSize=10000, 
                 lastzArguments="", compressFiles=True, realign=False, realignArguments="",
                 minimumSequenceLength=1, memory=None,
                 # Trim options for trimming ingroup seqs:
                 trimFlanking=10, trimMinSize=20,
                 trimWindowSize=10, trimThreshold=1,
                 # Trim options for trimming outgroup seqs (options
                 # other than flanking sequence will need a check to
                 # remove alignments that don't qualify)
                 # HACK: outgroup flanking is only set so high by
                 # default because it's needed for the tests (which
                 # don't use realign.)
                 trimOutgroupFlanking=2000):
        """Class defining options for blast
        """
        self.chunkSize = chunkSize
        self.overlapSize = overlapSize
        
        if realign:
            self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] | cactus_realign %s SEQ_FILE_1 SEQ_FILE_2 > CIGARS_FILE"  % (lastzArguments, realignArguments) 
        else:
            self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE"  % lastzArguments 
        if realign:
            self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial | cactus_realign %s SEQ_FILE > CIGARS_FILE" % (lastzArguments, realignArguments)
        else:
            self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial > CIGARS_FILE" % lastzArguments
        self.compressFiles = compressFiles
        self.minimumSequenceLength = 10
        self.memory = memory
        self.trimFlanking = trimFlanking
        self.trimMinSize = trimMinSize
        self.trimThreshold = trimThreshold
        self.trimWindowSize = trimWindowSize
        self.trimOutgroupFlanking = trimOutgroupFlanking
        
class BlastFlowerWrapper(Job):
    def __init__(self, cactusDisk, flowerName, finalResultsFile, blastOptions):
        Job.__init__(self)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
    def run(self, fileStore):
        resultsID = self.addChild(BlastFlower(self.cactusDisk, self.flowerName, self.blastOptions)).rv()
        self.addFollowOn(WritePermanentFile(resultsID, self.finalResultsFile))

class BlastFlower(Job):
    """Take a reconstruction problem and generate the sequences in chunks to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, flowerName, blastOptions):
        Job.__init__(self)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 2
        
    def run(self, fileStore):
        chunksDir = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkFlowerSequences %s '%s' %s %i %i %i %s" % \
                                                          (getLogLevelString(), self.cactusDisk, self.flowerName, 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          self.blastOptions.minimumSequenceLength,
                                                          chunksDir)).split("\n") if chunk != "" ]
        logger.info("Broken up the flowers into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk) for chunk in chunks]
        selfResultsID = self.addChild(MakeSelfBlasts(self.blastOptions, chunkIDs)).rv()
        offDiagonalResultsID = self.addChild(MakeOffDiagonalBlasts(self.blastOptions, chunkIDs)).rv()
        return self.addFollowOn(CollateBlasts([selfResultsID, offDiagonalResultsID])).rv()
    
class BlastSequencesAllAgainstAllWrapper(Job):
    def __init__(self, sequenceFiles1, outputFile, blastOptions):
        Job.__init__(self)
        self.sequenceFiles1 = sequenceFiles1
        self.outputFile = outputFile
        self.blastOptions = blastOptions
    def run(self, fileStore):
        sequenceIDs1 = [fileStore.writeGlobalFile(seq) for seq in self.sequenceFiles1]
        outputID = self.addChild(BlastSequencesAllAgainstAll(sequenceIDs1, self.blastOptions)).rv()
        self.addFollowOn(WritePermanentFile(outputID, self.outputFile))

class BlastSequencesAllAgainstAll(Job):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFileIDs1, blastOptions):
        Job.__init__(self)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.blastOptions = blastOptions
        self.blastOptions.compressFiles = False
        blastOptions.roundsOfCoordinateConversion = 1
        
    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        chunks = getChunks(sequenceFiles1, getTempDirectory(rootDir=fileStore.getLocalTempDir()), self.blastOptions)
        logger.info("Broken up the sequence files into individual 'chunk' files")
        chunkIDs = [fileStore.writeGlobalFile(chunk) for chunk in chunks]
        diagonalResultsID = self.addChild(MakeSelfBlasts(self.blastOptions, chunkIDs)).rv()
        offDiagonalResultsID = self.addChild(MakeOffDiagonalBlasts(self.blastOptions, chunkIDs)).rv()
        logger.debug("Collating the blasts after blasting all-against-all")
        return self.addFollowOn(CollateBlasts([diagonalResultsID, offDiagonalResultsID])).rv()
        
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
        selfResultsDir = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        resultsIDs = []
        for i in xrange(len(self.chunkIDs)):
            resultsIDs.append(self.addChild(RunSelfBlast(self.blastOptions, self.chunkIDs[i])).rv())
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        logger.debug("Collating self blasts.")
        logger.info("Blast file IDs: %s" % resultsIDs)
        return self.addFollowOn(CollateBlasts(resultsIDs)).rv()
    
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
                    resultsIDs.append(self.addChild(RunBlast(self.blastOptions, self.chunkIDs[i], self.chunkIDs[j])).rv())
            logger.info("Made the list of all-against-all blasts")
            #Set up the job to collate all the results
            logger.debug("Collating off-diagonal blasts")
            return self.addFollowOn(CollateBlasts(resultsIDs)).rv()
class BlastSequencesAgainstEachOtherWrapper(Job):
    def __init__(self, sequenceFiles1, sequenceFiles2, cigarFile, blastOptions):
        Job.__init__(self)
        self.sequenceFiles1 = sequenceFiles1
        self.sequenceFiles2 = sequenceFiles2
        self.cigarFile = cigarFile
        self.blastOptions = blastOptions
    def run(self, fileStore):
        seqIDs1 = [fileStore.writeGlobalFile(seq) for seq in self.sequenceFiles1]
        seqIDs2 = [fileStore.writeGlobalFile(seq) for seq in self.sequenceFiles2]
        cigarID = self.addChild(BlastSequencesAgainstEachOther(seqIDs1, seqIDs2, self.blastOptions)).rv()
        self.addFollowOn(WritePermanentFile(cigarID, self.cigarFile))
            
class BlastSequencesAgainstEachOther(Job):
    """Take two sets of sequences, chunks them up and blasts one set against the other.
    """
    def __init__(self, sequenceFileIDs1, sequenceFileIDs2, blastOptions):
        Job.__init__(self)
        self.sequenceFileIDs1 = sequenceFileIDs1
        self.sequenceFileIDs2 = sequenceFileIDs2
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        
    def run(self, fileStore):
        sequenceFiles1 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs1]
        sequenceFiles2 = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceFileIDs2]
        chunks1 = getChunks(sequenceFiles1, getTempDirectory(rootDir=fileStore.getLocalTempDir()), self.blastOptions)
        chunks2 = getChunks(sequenceFiles2, getTempDirectory(rootDir=fileStore.getLocalTempDir()), self.blastOptions)
        chunkIDs1 = [fileStore.writeGlobalFile(chunk) for chunk in chunks1]
        chunkIDs2 = [fileStore.writeGlobalFile(chunk) for chunk in chunks2]
        resultsIDs = []
        #Make the list of blast jobs.
        for chunkID1 in chunkIDs1:
            for chunkID2 in chunkIDs2:
                #TODO: Make the compression work
                self.blastOptions.compressFiles = False
                resultsIDs.append(self.addChild(RunBlast(self.blastOptions, chunkID1, chunkID2)).rv())
        logger.info("Made the list of blasts")
        #Set up the job to collate all the results
        return self.addFollowOn(CollateBlasts(resultsIDs)).rv()
class BlastIngroupsAndOutgroupsWrapper(Job):
    def __init__(self, blastOptions, ingroups, outgroups, cigarFile, outgroupFragmentsDir):
        Job.__init__(self)
        self.blastOptions = blastOptions
        self.outgroups = outgroups
        self.ingroups = ingroups
        self.cigarFile = cigarFile
        self.outgroupFragmentsDir = outgroupFragmentsDir
    def run(self, fileStore):
        fileStore.logToMaster("Blasting ingroups vs outgroups to file %s" % (self.cigarFile))
        ingroupIDs = [fileStore.writeGlobalFile(seq) for seq in self.ingroups]
        outgroupIDs = [fileStore.writeGlobalFile(seq) for seq in self.outgroups]
        try:
            os.makedirs(self.outgroupFragmentsDir)
        except os.error:
            # Directory already exists
            pass
        results = self.addChild(BlastIngroupsAndOutgroups(self.blastOptions, ingroupIDs, outgroupIDs)).rv()
        
        self.addFollowOn(WriteIngroupAndOutgroupResults(results, self.outgroups, self.cigarFile, self.outgroupFragmentsDir))
class WriteIngroupAndOutgroupResults(Job):
    def __init__(self, results, outgroups, cigarFile, outgroupFragmentsDir):
        Job.__init__(self)
        self.results = results
        self.outgroups = outgroups
        self.cigarFile = cigarFile
        self.outgroupFragmentsDir = outgroupFragmentsDir
        
    def run(self, fileStore):
        ingroupResultsID = self.results["ingroupResultsID"]
        outgroupFragmentsAndOutgroupResults = self.results["outgroupFragmentsAndOutgroupResults"]
        outgroupFragmentIDs = outgroupFragmentsAndOutgroupResults["outgroupFragmentIDs"]
        outgroupResultsID = outgroupFragmentsAndOutgroupResults["outputID"]
        logger.info("Outgroup fragment IDs: %s" % outgroupFragmentIDs)
        logger.info("outgroupResultsID: %s" % outgroupResultsID)
        logger.info("ingroupResultsID: %s" % ingroupResultsID)

        outgroupFragments = [fileStore.readGlobalFile(seq) for seq in outgroupFragmentIDs]
        ingroupResults = fileStore.readGlobalFile(ingroupResultsID)
        outgroupResults = fileStore.readGlobalFile(outgroupResultsID)

        catFiles([ingroupResults, outgroupResults], self.cigarFile)
        logger.info("Writing BlastIngroupsAndOutgroups cigar file to: %s" % self.cigarFile)
        assert len(outgroupFragments) == len(self.outgroups)
        for (outgroupFragment, outgroupName) in zip(outgroupFragments, self.outgroups):
            system("cp %s %s" % (outgroupFragment, os.path.join(self.outgroupFragmentsDir, os.path.basename(outgroupName))))
            
        
class BlastIngroupsAndOutgroups(Job):
    """Blast ingroup sequences against each other, and against the given
    outgroup sequences in succession. The next outgroup is only
    aligned against the regions that are not found in the previous
    outgroup.
    """
    def __init__(self, blastOptions, ingroupSequenceIDs, outgroupSequenceIDs):
        Job.__init__(self, memory = blastOptions.memory)
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        self.ingroupSequenceIDs = ingroupSequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs

    def run(self, fileStore):
        ingroupResultsID = self.addChild(BlastSequencesAllAgainstAll(self.ingroupSequenceIDs, self.blastOptions)).rv()

        outgroupFragmentsAndOutgroupResults = self.addChild(BlastFirstOutgroup(untrimmedSequenceIDs=self.ingroupSequenceIDs,
                                               sequenceIDs=self.ingroupSequenceIDs,
                                               outgroupSequenceIDs=self.outgroupSequenceIDs,
                                               outgroupFragmentIDs=[],
                                               outputID=None,
                                               blastOptions=self.blastOptions,
                                               outgroupNumber=1)).rv()
        return {"ingroupResultsID":ingroupResultsID, "outgroupFragmentsAndOutgroupResults":outgroupFragmentsAndOutgroupResults}
        

class BlastFirstOutgroup(Job):
    """Blast the given sequence(s) against the first of a succession of
    outgroups, only aligning fragments that haven't aligned to the
    previous outgroups. Then recurse on the other outgroups. Returns the 
    blast cigar file and a a list of the trimmed outgroups.
    """
    def __init__(self, untrimmedSequenceIDs, sequenceIDs, outgroupSequenceIDs, outgroupFragmentIDs, outputID,
                 blastOptions, outgroupNumber):
        Job.__init__(self, memory=blastOptions.memory)
        self.untrimmedSequenceIDs = untrimmedSequenceIDs
        self.sequenceIDs = sequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.outputID = outputID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber

    def run(self, fileStore):
        logger.info("Blasting ingroup sequences %s to outgroup %s" % (self.sequenceIDs, self.outgroupSequenceIDs[0]))
        blastResultsID = self.addChild(BlastSequencesAgainstEachOther(self.sequenceIDs,
                                                           [self.outgroupSequenceIDs[0]],
                                                           self.blastOptions)).rv()
        return self.addFollowOn(TrimAndRecurseOnOutgroups(untrimmedSequenceIDs=self.untrimmedSequenceIDs,
                                                          sequenceIDs=self.sequenceIDs,
                                                          outgroupSequenceIDs=self.outgroupSequenceIDs,
                                                          mostRecentResultsID=blastResultsID,
                                                          outgroupFragmentIDs=self.outgroupFragmentIDs,
                                                          outputID=self.outputID,
                                                          blastOptions=self.blastOptions,
                                                          outgroupNumber=self.outgroupNumber)).rv()

class TrimAndRecurseOnOutgroups(Job):
    def __init__(self, untrimmedSequenceIDs, sequenceIDs,
                 outgroupSequenceIDs, mostRecentResultsID, outgroupFragmentIDs, outputID, blastOptions, outgroupNumber):
        Job.__init__(self)
        self.untrimmedSequenceIDs = untrimmedSequenceIDs
        self.sequenceIDs = sequenceIDs
        self.outgroupSequenceIDs = outgroupSequenceIDs
        self.mostRecentResultsID = mostRecentResultsID
        self.outgroupFragmentIDs = outgroupFragmentIDs
        self.outputID = outputID
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber

    def run(self, fileStore):
        # Trim outgroup, convert outgroup coordinates, and add to
        # outgroup fragments dir
        untrimmedSequences = [fileStore.readGlobalFile(fileID) for fileID in self.untrimmedSequenceIDs]
        outgroupSequences = [fileStore.readGlobalFile(fileID) for fileID in self.outgroupSequenceIDs]
        sequences = [fileStore.readGlobalFile(fileID) for fileID in self.sequenceIDs]
        mostRecentResults = fileStore.readGlobalFile(self.mostRecentResultsID)
        if self.outputID:
            outputFile = fileStore.readGlobalFile(self.outputID)
        else:
            outputFile = fileStore.getLocalTempFile()
        trimmedOutgroup = fileStore.getLocalTempFile()
        outgroupCoverage = fileStore.getLocalTempFile()
        calculateCoverage(outgroupSequences[0],
                          mostRecentResults, outgroupCoverage)
        # The windowSize and threshold are fixed at 1: anything more
        # and we will run into problems with alignments that aren't
        # covered in a matching trimmed sequence.
        trimGenome(outgroupSequences[0], outgroupCoverage,
                   trimmedOutgroup, flanking=self.blastOptions.trimOutgroupFlanking,
                   windowSize=1, threshold=1)
        outgroupConvertedResultsFile = fileStore.getLocalTempFile()
        system("cactus_upconvertCoordinates.py %s %s 1 > %s" %\
               (trimmedOutgroup, mostRecentResults,
                outgroupConvertedResultsFile))
        self.outgroupFragmentIDs.append(fileStore.writeGlobalFile(trimmedOutgroup))


        # Report coverage of the latest outgroup on the trimmed ingroups.
        for trimmedIngroupSequence, ingroupSequence in zip(sequences, untrimmedSequences):
            tmpIngroupCoverage = fileStore.getLocalTempFile()
            calculateCoverage(trimmedIngroupSequence, mostRecentResults,
                              tmpIngroupCoverage)
            fileStore.logToMaster("Coverage on %s from outgroup #%d, %s: %s%% (current ingroup length %d, untrimmed length %d). Outgroup trimmed to %d bp from %d" % (os.path.basename(ingroupSequence), self.outgroupNumber, os.path.basename(outgroupSequences[0]), percentCoverage(trimmedIngroupSequence, tmpIngroupCoverage), sequenceLength(trimmedIngroupSequence), sequenceLength(ingroupSequence), sequenceLength(trimmedOutgroup), sequenceLength(outgroupSequences[0])))


        # Convert the alignments' ingroup coordinates.
        ingroupConvertedResultsFile = fileStore.getLocalTempFile()
        if self.sequenceIDs == self.untrimmedSequenceIDs:
            # No need to convert ingroup coordinates on first run.
            system("cp %s %s" % (outgroupConvertedResultsFile,
                                 ingroupConvertedResultsFile))
        else:
            system("cactus_blast_convertCoordinates --onlyContig1 %s %s 1" % (
                outgroupConvertedResultsFile, ingroupConvertedResultsFile))
        
        # Append the latest results to the accumulated outgroup coverage file
        with open(ingroupConvertedResultsFile) as results:
            with open(outputFile, 'a') as output:
                output.write(results.read())
        self.outputID = fileStore.writeGlobalFile(outputFile)

        # Report coverage of the all outgroup alignments so far on the ingroups.
        ingroupCoverageFiles = []
        for ingroupSequence in untrimmedSequences:
            tmpIngroupCoverage = fileStore.getLocalTempFile()
            calculateCoverage(ingroupSequence, outputFile,
                              tmpIngroupCoverage)
            fileStore.logToMaster("Cumulative coverage of %d outgroups on ingroup %s: %s" % (self.outgroupNumber, os.path.basename(ingroupSequence), percentCoverage(ingroupSequence, tmpIngroupCoverage)))
            ingroupCoverageFiles.append(tmpIngroupCoverage)

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
        if len(outgroupSequences) > 1:
            trimmedSeqs = []
            # Use the accumulated results so far to trim away the
            # aligned parts of the ingroups.
            for i, sequenceFile in enumerate(untrimmedSequences):
                coverageFile = ingroupCoverageFiles[i]

                trimmed = fileStore.getLocalTempFile()
                trimGenome(sequenceFile, coverageFile, trimmed,
                           complement=True, flanking=self.blastOptions.trimFlanking,
                           minSize=self.blastOptions.trimMinSize,
                           threshold=self.blastOptions.trimThreshold,
                           windowSize=self.blastOptions.trimWindowSize)
                trimmedSeqs.append(trimmed)
            trimmedSeqIDs = [fileStore.writeGlobalFile(seq) for seq in trimmedSeqs]
            return self.addChild(BlastFirstOutgroup(untrimmedSequenceIDs=self.untrimmedSequenceIDs,
                                                   sequenceIDs=trimmedSeqIDs,
                                                   outgroupSequenceIDs=self.outgroupSequenceIDs[1:],
                                                   outgroupFragmentIDs=self.outgroupFragmentIDs,
                                                   outputID=self.outputID,
                                                   blastOptions=self.blastOptions,
                                                   outgroupNumber=self.outgroupNumber + 1)).rv()
        else:
            return {"outgroupFragmentIDs":self.outgroupFragmentIDs, "outputID":self.outputID}
                

def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    system("bzip2 --keep --fast %s" % fileName)
    return fileName + ".bz2"
        
class RunSelfBlast(Job):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID):
        Job.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFileID = seqFileID
    
    def run(self, fileStore):   
        tempResultsFile = fileStore.getLocalTempFile()
        seqFile = fileStore.readGlobalFile(self.seqFileID)
        command = self.blastOptions.selfBlastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE", seqFile)
        system(command)
        resultsFile = fileStore.getLocalTempFile()
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        if self.blastOptions.compressFiles:
            #TODO: This throws away the compressed file
            seqFile = compressFastaFile(seqFile)
        logger.info("Ran the self blast okay")
        return fileStore.writeGlobalFile(resultsFile)

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    system("bunzip2 --stdout %s > %s" % (fileName, tempFileName))
    return tempFileName
    
class RunBlast(Job):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFileID1, seqFileID2):
        Job.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFileID1 = seqFileID1
        self.seqFileID2 = seqFileID2
    
    def run(self, fileStore):
        seqFile1 = fileStore.readGlobalFile(self.seqFileID1)
        seqFile2 = fileStore.readGlobalFile(self.seqFileID2)
        if self.blastOptions.compressFiles:
            seqFile1 = decompressFastaFile(seqFile1, fileStore.getLocalTempFile())
            seqFile2 = decompressFastaFile(seqFile2, fileStore.getLocalTempFile())
        tempResultsFile = fileStore.getLocalTempFile()
        command = self.blastOptions.blastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE_1", seqFile1).replace("SEQ_FILE_2", seqFile2)
        system(command)
        resultsFile = fileStore.getLocalTempFile()
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        logger.info("Ran the blast okay")
        return fileStore.writeGlobalFile(resultsFile)

class CollateBlasts(Job):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, resultsFileIDs):
        Job.__init__(self)
        self.resultsFileIDs = resultsFileIDs
    
    def run(self, fileStore):
        logger.info("Results IDs: %s" % self.resultsFileIDs)
        resultsFiles = [fileStore.readGlobalFile(fileID) for fileID in self.resultsFileIDs]
        collatedResultsFile = fileStore.getLocalTempFile()
        catFiles(resultsFiles, collatedResultsFile)
        logger.info("Collated the alignments to the file: %s",  collatedResultsFile)
        return fileStore.writeGlobalFile(collatedResultsFile)
        
        
class SortCigarAlignmentsInPlace(Job):
    """Sorts an alignment file in place.
    """
    def __init__(self, cigarFile):
        Job.__init__(self)
        self.cigarFile = cigarFile
    
    def run(self, fileStore):
        tempResultsFile = os.path.join(fileStore.getLocalTempDir(), "tempResults.cig")
        system("cactus_blast_sortAlignments %s %s %i" % (getLogLevelString(), self.cigarFile, tempResultsFile))
        logger.info("Sorted the alignments okay")
        system("mv %s %s" % (tempResultsFile, self.cigarFile))

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

def calculateCoverage(sequenceFile, cigarFile, outputFile):
    logger.info("Calculating coverage of cigar file %s on %s, writing to %s" % (
        cigarFile, sequenceFile, outputFile))
    system("cactus_coverage %s %s > %s" % (sequenceFile,
                                           cigarFile,
                                           outputFile))

def trimGenome(sequenceFile, coverageFile, outputFile, complement=False,
               flanking=0, minSize=1, windowSize=10, threshold=1):
    system("cactus_trimSequences.py %s %s %s %s %s %s %s > %s" % (
        nameValue("complement", complement, valueType=bool),
        nameValue("flanking", flanking), nameValue("minSize", minSize),
        nameValue("windowSize", windowSize), nameValue("threshold", threshold),
        sequenceFile, coverageFile, outputFile))
def getChunks(sequenceFiles, chunksDir, blastOptions):
    return [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                           (getLogLevelString(), 
                                            blastOptions.chunkSize, 
                                            blastOptions.overlapSize,
                                            chunksDir,
                                            " ".join(sequenceFiles))).split("\n") if chunk != "" ]
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    blastOptions = BlastOptions()
    
    #output stuff
    parser.add_argument("--seqFiles", dest="seqFiles",
                        help="Input sequence files.",
                        default=None)
    parser.add_argument("--cigars", dest="cigarFile", 
                      help="File to write cigars in",
                      default="cigarFile.txt")
    
    parser.add_argument("--chunkSize", dest="chunkSize", type=int, 
                     help="The size of chunks passed to lastz (must be at least twice as big as overlap)",
                     default=blastOptions.chunkSize)
    
    parser.add_argument("--overlapSize", dest="overlapSize", type=int,
                     help="The size of the overlap between the chunks passed to lastz (min size 2)",
                     default=blastOptions.overlapSize)
    
    parser.add_argument("--blastString", dest="blastString", type=str,
                     help="The default string used to call the blast program. \
Must contain three strings: SEQ_FILE_1, SEQ_FILE_2 and CIGARS_FILE which will be \
replaced with the two sequence files and the results file, respectively",
                     default=blastOptions.blastString)
    
    parser.add_argument("--selfBlastString", dest="selfBlastString", type=str,
                     help="The default string used to call the blast program for self alignment. \
Must contain three strings: SEQ_FILE and CIGARS_FILE which will be \
replaced with the the sequence file and the results file, respectively",
                     default=blastOptions.selfBlastString)
   
    parser.add_argument("--compressFiles", dest="compressFiles", action="store_false",
                      help="Turn of bz2 based file compression of sequences for I/O transfer", 
                      default=blastOptions.compressFiles)
    
    parser.add_argument("--lastzMemory", dest="memory", type=int,
                      help="Lastz memory (in bytes)", 
                      default=blastOptions.memory)
    
    parser.add_argument("--trimFlanking", type=int, help="Amount of flanking sequence to leave on trimmed ingroup sequences", default=blastOptions.trimFlanking)
    parser.add_argument("--trimMinSize", type=int, help="Minimum size, before adding flanking sequence, of ingroup sequence to align against the next outgroup", default=blastOptions.trimMinSize)
    parser.add_argument("--trimThreshold", type=int, help="Coverage threshold for an ingroup region to not be aligned against the next outgroup", default=blastOptions.trimThreshold)
    parser.add_argument("--trimWindowSize", type=int, help="Windowing size to integrate ingroup coverage over", default=blastOptions.trimWindowSize)
    parser.add_argument("--trimOutgroupFlanking", type=int, help="Amount of flanking sequence to leave on trimmed outgroup sequences", default=blastOptions.trimOutgroupFlanking)
    

    parser.add_argument("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")
    
    parser.add_argument("--targetSequenceFiles", dest="targetSequenceFiles", type=str,
                     help="Sequences to align against the input sequences against. If these are not provided then the input sequences are aligned against each other.",
                     default=None)

    parser.add_argument("--ingroups", type=str, default=None,
                      help="Ingroups to align (comma-separated) (--outgroups "
                      "must be provided as well")

    parser.add_argument("--outgroups", type=str, default=None,
                      help="Outgroups to align (comma-separated) (--ingroups "
                      "must be provided as well")

    parser.add_argument("--outgroupFragmentsDir", type=str,
                      default="outgroupFragments/", help= "Directory to "
                      "store outgroup fragments in")

    options = parser.parse_args()
    if options.test:
        _test()

    if (options.ingroups is not None) ^ (options.outgroups is not None):
        raise RuntimeError("--ingroups and --outgroups must be provided "
                           "together")
    if options.ingroups:
        firstJob = BlastIngroupsAndOutgroupsWrapper(options,
                                                options.ingroups.split(','),
                                                options.outgroups.split(','),
                                                options.cigarFile,
                                                options.outgroupFragmentsDir)
    elif options.targetSequenceFiles == None:
        firstJob = BlastSequencesAllAgainstAllWrapper(options.seqFiles.split(), options.cigarFile, options)
    else:
        firstJob = BlastSequencesAgainstEachOtherWrapper(options.seqFiles.split(), options.targetSequenceFiles.split(), options.cigarFile, options)
    Job.Runner.startToil(firstJob, options)

def _test():
    import doctest
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blast.cactus_blast import *
    main()
