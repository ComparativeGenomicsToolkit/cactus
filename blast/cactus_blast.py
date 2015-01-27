#!/usr/bin/env python
#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script for running an all against all (including self) set of alignments on a set of input
sequences. Uses the jobTree framework to parallelise the blasts.
"""
import os
import sys
from optparse import OptionParser
from sonLib.bioio import TempFileTree
from sonLib.bioio import logger
from sonLib.bioio import system, popenCatch
from sonLib.bioio import getLogLevelString
from sonLib.bioio import makeSubDir
from sonLib.bioio import catFiles
from sonLib.bioio import getTempFile
from sonLib.bioio import nameValue
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

class BlastOptions:
    def __init__(self, chunkSize=10000000, overlapSize=10000, 
                 lastzArguments="", lastArguments = "", compressFiles=True, realign=False, realignArguments="",
                 minimumSequenceLength=1, memory=sys.maxint,
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
        
        #switch from lastZ to LAST. LAST only outputs MAF format, so convert to CIGAR
        if realign:
        	self.blastString = "cactus_lastdb %s temp SEQ_FILE_1 && cactus_lastal temp SEQ_FILE_2 | maf2cigar | cactus_realign %s SEQ_FILE_1 SEQ_FILE_2 > CIGARS_FILE" % (lastArguments, realignArguments) 
        else:
        	self.blastString = "cactus_lastdb %s temp SEQ_FILE_1 && cactus_lastal temp SEQ_FILE_2 | maf2cigar > CIGARS_FILE" % lastArguments
        if realign:
        	self.selfBlastString = "cactus_lastdb %s temp SEQ_FILE && cactus_lastal temp SEQ_FILE | maf2cigar | cactus_realign %s SEQ_FILE > CIGARS_FILE" % (lastArguments, realignArguments)
        else:
			self.selfBlastString = "cactus_lastdb %s temp SEQ_FILE && cactus_lastal temp SEQ_FILE | maf2cigar > CIGARS_FILE" % lastArguments

        #if realign:
        #    self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] | cactus_realign %s SEQ_FILE_1 SEQ_FILE_2 > CIGARS_FILE"  % (lastzArguments, realignArguments) 
        #else:
        #    self.blastString = "cactus_lastz --format=cigar %s SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE"  % lastzArguments 
        #if realign:
        #    self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial | cactus_realign %s SEQ_FILE > CIGARS_FILE" % (lastzArguments, realignArguments)
        #else:
        #    self.selfBlastString = "cactus_lastz --format=cigar %s SEQ_FILE[multiple][nameparse=darkspace] SEQ_FILE[nameparse=darkspace] --notrivial > CIGARS_FILE" % lastzArguments
        
        self.compressFiles = compressFiles
        self.minimumSequenceLength = 10
        self.memory = memory
        self.trimFlanking = trimFlanking
        self.trimMinSize = trimMinSize
        self.trimThreshold = trimThreshold
        self.trimWindowSize = trimWindowSize
        self.trimOutgroupFlanking = trimOutgroupFlanking

class BlastFlower(Target):
    """Take a reconstruction problem and generate the sequences in chunks to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, flowerName, finalResultsFile, blastOptions):
        Target.__init__(self)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 2
        
    def run(self):
        chunksDir = makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks"))
        chunks = [ chunk for chunk in popenCatch("cactus_blast_chunkFlowerSequences %s '%s' %s %i %i %i %s" % \
                                                          (getLogLevelString(), self.cactusDisk, self.flowerName, 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          self.blastOptions.minimumSequenceLength,
                                                          chunksDir)).split("\n") if chunk != "" ]
        logger.info("Broken up the flowers into individual 'chunk' files")
        self.addChildTarget(MakeBlastsAllAgainstAll(self.blastOptions, chunks, self.finalResultsFile))
        
class BlastSequencesAllAgainstAll(Target):
    """Take a set of sequences, chunks them up and blasts them.
    """
    def __init__(self, sequenceFiles1, finalResultsFile, blastOptions):
        Target.__init__(self)
        self.sequenceFiles1 = sequenceFiles1
        self.finalResultsFile = finalResultsFile
        self.blastOptions = blastOptions
        blastOptions.roundsOfCoordinateConversion = 1
    
    def getChunks(self, sequenceFiles, chunksDir):
        return [ chunk for chunk in popenCatch("cactus_blast_chunkSequences %s %i %i %s %s" % \
                                                          (getLogLevelString(), 
                                                          self.blastOptions.chunkSize, 
                                                          self.blastOptions.overlapSize,
                                                          chunksDir,
                                                          " ".join(sequenceFiles))).split("\n") if chunk != "" ]
        
    def run(self):
        chunks = self.getChunks(self.sequenceFiles1, makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks")))
        logger.info("Broken up the sequence files into individual 'chunk' files")
        self.addChildTarget(MakeBlastsAllAgainstAll(self.blastOptions, chunks, self.finalResultsFile))
        
class MakeBlastsAllAgainstAll(Target):
    """Breaks up the inputs into bits and builds a bunch of alignment jobs.
    """
    def __init__(self, blastOptions, chunks, finalResultsFile):
        Target.__init__(self)
        self.blastOptions = blastOptions
        self.chunks = chunks
        self.finalResultsFile = finalResultsFile
        
    def run(self):
        #Avoid compression if just one chunk
        self.blastOptions.compressFiles = self.blastOptions.compressFiles and len(self.chunks) > 2
        selfResultsDir = makeSubDir(os.path.join(self.getGlobalTempDir(), "selfResults"))
        resultsFiles = []
        for i in xrange(len(self.chunks)):
            resultsFile = os.path.join(selfResultsDir, str(i))
            resultsFiles.append(resultsFile)
            self.addChildTarget(RunSelfBlast(self.blastOptions, self.chunks[i], resultsFile))
        logger.info("Made the list of self blasts")
        #Setup job to make all-against-all blasts
        self.setFollowOnTarget(MakeBlastsAllAgainstAll2(self.blastOptions, self.chunks, resultsFiles, self.finalResultsFile))
    
class MakeBlastsAllAgainstAll2(MakeBlastsAllAgainstAll):
        def __init__(self, blastOptions, chunks, resultsFiles, finalResultsFile):
            MakeBlastsAllAgainstAll.__init__(self, blastOptions, chunks, finalResultsFile)
            self.resultsFiles = resultsFiles
           
        def run(self):
            tempFileTree = TempFileTree(os.path.join(self.getGlobalTempDir(), "allAgainstAllResults"))
            #Make the list of blast jobs.
            for i in xrange(0, len(self.chunks)):
                for j in xrange(i+1, len(self.chunks)):
                    resultsFile = tempFileTree.getTempFile()
                    self.resultsFiles.append(resultsFile)
                    self.addChildTarget(RunBlast(self.blastOptions, self.chunks[i], self.chunks[j], resultsFile))
            logger.info("Made the list of all-against-all blasts")
            #Set up the job to collate all the results
            self.setFollowOnTarget(CollateBlasts(self.finalResultsFile, self.resultsFiles))
            
class BlastSequencesAgainstEachOther(BlastSequencesAllAgainstAll):
    """Take two sets of sequences, chunks them up and blasts one set against the other.
    """
    def __init__(self, sequenceFiles1, sequenceFiles2, finalResultsFile, blastOptions):
        BlastSequencesAllAgainstAll.__init__(self, sequenceFiles1, finalResultsFile, blastOptions)
        self.sequenceFiles2 = sequenceFiles2
        
    def run(self):
        chunks1 = self.getChunks(self.sequenceFiles1, makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks1")))
        chunks2 = self.getChunks(self.sequenceFiles2, makeSubDir(os.path.join(self.getGlobalTempDir(), "chunks2")))
        tempFileTree = TempFileTree(os.path.join(self.getGlobalTempDir(), "allAgainstAllResults"))
        resultsFiles = []
        #Make the list of blast jobs.
        for chunk1 in chunks1:
            for chunk2 in chunks2:
                resultsFile = tempFileTree.getTempFile()
                resultsFiles.append(resultsFile)
                #TODO: Make the compression work
                self.blastOptions.compressFiles = False
                self.addChildTarget(RunBlast(self.blastOptions, chunk1, chunk2, resultsFile))
        logger.info("Made the list of blasts")
        #Set up the job to collate all the results
        self.setFollowOnTarget(CollateBlasts(self.finalResultsFile, resultsFiles))

class BlastIngroupsAndOutgroups(Target):
    """Blast ingroup sequences against each other, and against the given
    outgroup sequences in succession. The next outgroup is only
    aligned against the regions that are not found in the previous
    outgroup.
    """
    def __init__(self, blastOptions, ingroupSequenceFiles,
                 outgroupSequenceFiles, finalResultsFile,
                 outgroupFragmentsDir):
        Target.__init__(self, memory = blastOptions.memory)
        self.blastOptions = blastOptions
        self.blastOptions.roundsOfCoordinateConversion = 1
        self.ingroupSequenceFiles = ingroupSequenceFiles
        self.outgroupSequenceFiles = outgroupSequenceFiles
        self.finalResultsFile = finalResultsFile
        self.outgroupFragmentsDir = outgroupFragmentsDir

    def run(self):
        self.logToMaster("Blasting ingroups vs outgroups to file %s" % (self.finalResultsFile))
        try:
            os.makedirs(self.outgroupFragmentsDir)
        except os.error:
            # Directory already exists
            pass
        ingroupResults = getTempFile(rootDir=self.getGlobalTempDir())
        self.addChildTarget(BlastSequencesAllAgainstAll(self.ingroupSequenceFiles,
                                                        ingroupResults,
                                                        self.blastOptions))
        outgroupResults = getTempFile(rootDir=self.getGlobalTempDir())
        self.addChildTarget(BlastFirstOutgroup(self.ingroupSequenceFiles,
                                               self.ingroupSequenceFiles,
                                               self.outgroupSequenceFiles,
                                               self.outgroupFragmentsDir,
                                               outgroupResults,
                                               self.blastOptions, 1))
        self.setFollowOnTarget(CollateBlasts(self.finalResultsFile,
                                             [ingroupResults, outgroupResults]))

class BlastFirstOutgroup(Target):
    """Blast the given sequence(s) against the first of a succession of
    outgroups, only aligning fragments that haven't aligned to the
    previous outgroups. Then recurse on the other outgroups.
    """
    def __init__(self, untrimmedSequenceFiles, sequenceFiles,
                 outgroupSequenceFiles, outgroupFragmentsDir, outputFile,
                 blastOptions, outgroupNumber):
        Target.__init__(self, memory=blastOptions.memory)
        self.untrimmedSequenceFiles = untrimmedSequenceFiles
        self.sequenceFiles = sequenceFiles
        self.outgroupSequenceFiles = outgroupSequenceFiles
        self.outgroupFragmentsDir = outgroupFragmentsDir
        self.outputFile = outputFile
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber

    def run(self):
        logger.info("Blasting ingroup sequences %s to outgroup %s" % (self.sequenceFiles, self.outgroupSequenceFiles[0]))
        blastResults = getTempFile(rootDir=self.getGlobalTempDir())
        self.addChildTarget(BlastSequencesAgainstEachOther(self.sequenceFiles,
                                                           [self.outgroupSequenceFiles[0]],
                                                           blastResults,
                                                           self.blastOptions))
        self.setFollowOnTarget(TrimAndRecurseOnOutgroups(self.untrimmedSequenceFiles,
                                                         self.sequenceFiles,
                                                         self.outgroupSequenceFiles,
                                                         self.outgroupFragmentsDir,
                                                         blastResults,
                                                         self.outputFile,
                                                         self.blastOptions,
                                                         self.outgroupNumber))

class TrimAndRecurseOnOutgroups(Target):
    def __init__(self, untrimmedSequenceFiles, sequenceFiles,
                 outgroupSequenceFiles, outgroupFragmentsDir,
                 mostRecentResultsFile, outputFile, blastOptions,
                 outgroupNumber):
        Target.__init__(self)
        self.untrimmedSequenceFiles = untrimmedSequenceFiles
        self.sequenceFiles = sequenceFiles
        self.outgroupSequenceFiles = outgroupSequenceFiles
        self.outgroupFragmentsDir = outgroupFragmentsDir
        self.mostRecentResultsFile = mostRecentResultsFile
        self.outputFile = outputFile
        self.blastOptions = blastOptions
        self.outgroupNumber = outgroupNumber

    def run(self):
        # Trim outgroup, convert outgroup coordinates, and add to
        # outgroup fragments dir
        trimmedOutgroup = getTempFile(rootDir=self.getGlobalTempDir())
        outgroupCoverage = getTempFile(rootDir=self.getGlobalTempDir())
        calculateCoverage(self.outgroupSequenceFiles[0],
                          self.mostRecentResultsFile, outgroupCoverage)
        # The windowSize and threshold are fixed at 1: anything more
        # and we will run into problems with alignments that aren't
        # covered in a matching trimmed sequence.
        trimGenome(self.outgroupSequenceFiles[0], outgroupCoverage,
                   trimmedOutgroup, flanking=self.blastOptions.trimOutgroupFlanking,
                   windowSize=1, threshold=1)
        outgroupConvertedResultsFile = getTempFile(rootDir=self.getGlobalTempDir())
        system("cactus_upconvertCoordinates.py %s %s 1 > %s" %\
               (trimmedOutgroup, self.mostRecentResultsFile,
                outgroupConvertedResultsFile))
        system("cat %s > %s" % (trimmedOutgroup, os.path.join(self.outgroupFragmentsDir, os.path.basename(self.outgroupSequenceFiles[0]))))

        # Report coverage of the latest outgroup on the trimmed ingroups.
        for trimmedIngroupSequence, ingroupSequence in zip(self.sequenceFiles, self.untrimmedSequenceFiles):
            tmpIngroupCoverage = getTempFile(rootDir=self.getGlobalTempDir())
            calculateCoverage(trimmedIngroupSequence, self.mostRecentResultsFile,
                              tmpIngroupCoverage)
            self.logToMaster("Coverage on %s from outgroup #%d, %s: %s%% (current ingroup length %d, untrimmed length %d). Outgroup trimmed to %d bp from %d" % (os.path.basename(ingroupSequence), self.outgroupNumber, os.path.basename(self.outgroupSequenceFiles[0]), percentCoverage(trimmedIngroupSequence, tmpIngroupCoverage), sequenceLength(trimmedIngroupSequence), sequenceLength(ingroupSequence), sequenceLength(trimmedOutgroup), sequenceLength(self.outgroupSequenceFiles[0])))


        # Convert the alignments' ingroup coordinates.
        ingroupConvertedResultsFile = getTempFile(rootDir=self.getGlobalTempDir())
        if self.sequenceFiles == self.untrimmedSequenceFiles:
            # No need to convert ingroup coordinates on first run.
            system("cp %s %s" % (outgroupConvertedResultsFile,
                                 ingroupConvertedResultsFile))
        else:
            system("cactus_blast_convertCoordinates --onlyContig1 %s %s 1" % (
                outgroupConvertedResultsFile, ingroupConvertedResultsFile))
        
        # Append the latest results to the accumulated outgroup coverage file
        with open(ingroupConvertedResultsFile) as results:
            with open(self.outputFile, 'a') as output:
                output.write(results.read())
        os.remove(outgroupConvertedResultsFile)
        os.remove(ingroupConvertedResultsFile)
        os.remove(outgroupCoverage)
        os.remove(trimmedOutgroup)

        # Report coverage of the all outgroup alignments so far on the ingroups.
        ingroupCoverageFiles = []
        for ingroupSequence in self.untrimmedSequenceFiles:
            tmpIngroupCoverage = getTempFile(rootDir=self.getGlobalTempDir())
            calculateCoverage(ingroupSequence, self.outputFile,
                              tmpIngroupCoverage)
            self.logToMaster("Cumulative coverage of %d outgroups on ingroup %s: %s" % (self.outgroupNumber, os.path.basename(ingroupSequence), percentCoverage(ingroupSequence, tmpIngroupCoverage)))
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
        if len(self.outgroupSequenceFiles) > 1:
            trimmedSeqs = []
            # Use the accumulated results so far to trim away the
            # aligned parts of the ingroups.
            for i, sequenceFile in enumerate(self.untrimmedSequenceFiles):
                coverageFile = ingroupCoverageFiles[i]

                trimmed = getTempFile(rootDir=self.getGlobalTempDir())
                trimGenome(sequenceFile, coverageFile, trimmed,
                           complement=True, flanking=self.blastOptions.trimFlanking,
                           minSize=self.blastOptions.trimMinSize,
                           threshold=self.blastOptions.trimThreshold,
                           windowSize=self.blastOptions.trimWindowSize)
                trimmedSeqs.append(trimmed)
            self.addChildTarget(BlastFirstOutgroup(self.untrimmedSequenceFiles,
                                                   trimmedSeqs,
                                                   self.outgroupSequenceFiles[1:],
                                                   self.outgroupFragmentsDir,
                                                   self.outputFile,
                                                   self.blastOptions,
                                                   self.outgroupNumber + 1))

def compressFastaFile(fileName):
    """Compress a fasta file.
    """
    system("bzip2 --keep --fast %s" % fileName)
        
class RunSelfBlast(Target):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFile, resultsFile):
        Target.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFile = seqFile
        self.resultsFile = resultsFile
    
    def run(self):   
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
        command = self.blastOptions.selfBlastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE", self.seqFile)
        system(command)
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, self.resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        if self.blastOptions.compressFiles:
            compressFastaFile(self.seqFile)
        logger.info("Ran the self blast okay")

def decompressFastaFile(fileName, tempFileName):
    """Copies the file from the central dir to a temporary file, returning the temp file name.
    """
    system("bunzip2 --stdout %s > %s" % (fileName, tempFileName))
    return tempFileName
    
class RunBlast(Target):
    """Runs blast as a job.
    """
    def __init__(self, blastOptions, seqFile1, seqFile2, resultsFile):
        Target.__init__(self, memory=blastOptions.memory)
        self.blastOptions = blastOptions
        self.seqFile1 = seqFile1
        self.seqFile2 = seqFile2
        self.resultsFile = resultsFile
    
    def run(self):
        if self.blastOptions.compressFiles:
            self.seqFile1 = decompressFastaFile(self.seqFile1 + ".bz2", os.path.join(self.getLocalTempDir(), "1.fa"))
            self.seqFile2 = decompressFastaFile(self.seqFile2 + ".bz2", os.path.join(self.getLocalTempDir(), "2.fa"))
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
        command = self.blastOptions.blastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE_1", self.seqFile1).replace("SEQ_FILE_2", self.seqFile2)
		#print >> sys.stderr, "Running blast"
        system(command)
        system("cactus_blast_convertCoordinates %s %s %i" % (tempResultsFile, self.resultsFile, self.blastOptions.roundsOfCoordinateConversion))
        logger.info("Ran the blast okay")

class CollateBlasts(Target):
    """Collates all the blasts into a single alignments file.
    """
    def __init__(self, finalResultsFile, resultsFiles):
        Target.__init__(self)
        self.finalResultsFile = finalResultsFile
        self.resultsFiles = resultsFiles
    
    def run(self):
        catFiles(self.resultsFiles, self.finalResultsFile)
        logger.info("Collated the alignments to the file: %s",  self.finalResultsFile)
        
class SortCigarAlignmentsInPlace(Target):
    """Sorts an alignment file in place.
    """
    def __init__(self, cigarFile):
        Target.__init__(self)
        self.cigarFile = cigarFile
    
    def run(self):
        tempResultsFile = os.path.join(self.getLocalTempDir(), "tempResults.cig")
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
                                                options.outgroupFragmentsDir)
    elif options.targetSequenceFiles == None:
        firstTarget = BlastSequencesAllAgainstAll(args, options.cigarFile, options)
    else:
        firstTarget = BlastSequencesAgainstEachOther(args, options.targetSequenceFiles.split(), options.cigarFile, options)
    Stack(firstTarget).startJobTree(options)

def _test():
    import doctest
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blast.cactus_blast import *
    main()
