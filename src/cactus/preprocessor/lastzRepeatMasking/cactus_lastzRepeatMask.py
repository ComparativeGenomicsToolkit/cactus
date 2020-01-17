#!/usr/bin/env python3

## USE LASTZ TO SOFTMASK REPEATS OF A GIVEN FASTA SEQUENCE FILE.

import os

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob

class RepeatMaskOptions:
    def __init__(self,
            fragment=200,
            minPeriod=10,
            lastzOpts="",
            unmaskInput=False,
            unmaskOutput=False,
            proportionSampled=1.0):
        self.fragment = fragment
        self.minPeriod = minPeriod
        self.lastzOpts = lastzOpts
        self.unmaskInput = unmaskInput
        self.unmaskOutput = unmaskOutput
        self.proportionSampled = proportionSampled

        self.period = max(1, round(self.proportionSampled * self.minPeriod))

        # make sure fragment size is even so they can overlap by exactly one half.
        if self.fragment % 2:
            self.fragment += 1


class LastzRepeatMaskJob(RoundedJob):
    def __init__(self, repeatMaskOptions, queryID, targetIDs):
        targetsSize = sum(targetID.size for targetID in targetIDs)
        memory = 4*1024*1024*1024
        disk = 2*(queryID.size + targetsSize)
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.repeatMaskOptions = repeatMaskOptions
        self.queryID = queryID
        self.targetIDs = targetIDs

    def getFragments(self, fileStore, queryFile):
        """
        Chop up the query fasta into fragments of a certain size, overlapping by half their length.
        """
        fragments = fileStore.getLocalTempFile()
        cactus_call(infile=queryFile, outfile=fragments,
                    parameters=["cactus_fasta_fragments.py",
                                "--fragment=%s" % str(self.repeatMaskOptions.fragment),
                                "--step=%s" % (str(self.repeatMaskOptions.fragment // 2)),
                                "--origin=zero"])
        return fragments

    def alignFastaFragments(self, fileStore, targetFiles, fragments):
        """
        Align each query fragment against all the target chunks, stopping
        early to avoid exponential blowup if too many alignments are found.
        """
        target = fileStore.getLocalTempFile()
        catFiles(targetFiles, target)
        lastZSequenceHandling  = ['%s[multiple][nameparse=darkspace]' % os.path.basename(target), '%s[nameparse=darkspace]' % os.path.basename(fragments)]
        if self.repeatMaskOptions.unmaskInput:
            lastZSequenceHandling  = ['%s[multiple,unmask][nameparse=darkspace]' % os.path.basename(target), '%s[unmask][nameparse=darkspace]' % os.path.basename(fragments)]
        alignment = fileStore.getLocalTempFile()
        # Each time a fragment aligns to a base in the sequence, that
        # base's match count is incremented.  the plus three for the
        # period parameter is a fudge to ensure sufficient alignments
        # are found
        cactus_call(outfile=alignment,
                    parameters=["cPecanLastz"] + lastZSequenceHandling +
                                self.repeatMaskOptions.lastzOpts.split() +
                                ["--querydepth=keep,nowarn:%i" % (self.repeatMaskOptions.period+3),
                                 "--format=general:name1,zstart1,end1,name2,zstart2+,end2+",
                                 "--markend"])
        return alignment

    def maskCoveredIntervals(self, fileStore, queryFile, alignment):
        """
        Mask the query fasta using the alignments to the target. Anything with more alignments than the period gets masked.
        """
        #This runs Bob's covered intervals program, which combines the lastz alignment info into intervals of the query.
        maskInfo = fileStore.getLocalTempFile()
        cactus_call(infile=alignment, outfile=maskInfo,
                    parameters=["cactus_covered_intervals",
                                "--queryoffsets",
                                "--origin=one",
                                # * 2 takes into account the effect of the overlap
                                "M=%s" % (int(self.repeatMaskOptions.period*2))])

        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output.
        args = ["--origin=one"]
        if self.repeatMaskOptions.unmaskOutput:
            args.append("--unmask")
        args.append(maskInfo)
        maskedQuery = fileStore.getLocalTempFile()
        cactus_call(infile=queryFile, outfile=maskedQuery,
                    parameters=["cactus_fasta_softmask_intervals.py"] + args)
        return maskedQuery

    def run(self, fileStore):
        """
        Using sampled target fragments, mask repetitive regions of the query.
        """
        assert len(self.targetIDs) >= 1
        assert self.repeatMaskOptions.fragment > 1
        queryFile = fileStore.readGlobalFile(self.queryID)
        targetFiles = [fileStore.readGlobalFile(fileID) for fileID in self.targetIDs]

        fragments = self.getFragments(fileStore, queryFile)
        alignment = self.alignFastaFragments(fileStore, targetFiles, fragments)
        maskedQuery = self.maskCoveredIntervals(fileStore, queryFile, alignment)
        return fileStore.writeGlobalFile(maskedQuery)
