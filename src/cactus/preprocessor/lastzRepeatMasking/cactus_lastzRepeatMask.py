#!/usr/bin/env python

## USE LASTZ TO SOFTMASK REPEATS OF A GIVEN FASTA SEQUENCE FILE.  

############################################################
##  NOTE:   NEW VERSION OF LASTZ REQUIRED (MINIMUM 1.02.40)
##          lastz/tools MUST BE IN SYSTEM PATH
############################################################

import os
import re
import sys
import tempfile
import shutil
import random

from argparse import ArgumentParser
from sonLib.bioio import system, catFiles

from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import readGlobalFileWithoutCache

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


class AlignFastaFragments(RoundedJob):
    def __init__(self, repeatMaskOptions, fragmentsID, targetIDs):
        if hasattr(fragmentsID, "size"):
            targetsSize = sum(targetID.size for targetID in targetIDs)
            memory = 3500000000
            disk = 2*(fragmentsID.size + targetsSize)
        else:
            memory = None
            disk = None
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.repeatMaskOptions = repeatMaskOptions
        self.fragmentsID = fragmentsID
        self.targetIDs = targetIDs

    def run(self, fileStore):
        # Align each fragment against a chunk of the input sequence.  Each time a fragment aligns to a base
        # in the sequence, that base's match count is incremented.
        # the plus three for the period parameter is a fudge to ensure sufficient alignments are found
        fragments = fileStore.readGlobalFile(self.fragmentsID)
        targetFiles = [fileStore.readGlobalFile(fileID) for fileID in self.targetIDs]
        target = fileStore.getLocalTempFile()
        catFiles(targetFiles, target)
        lastZSequenceHandling  = ['%s[multiple][nameparse=darkspace]' % os.path.basename(target), '%s[nameparse=darkspace]' % os.path.basename(fragments)]
        if self.repeatMaskOptions.unmaskInput:
            lastZSequenceHandling  = ['%s[multiple,unmask][nameparse=darkspace]' % os.path.basename(target), '%s[unmask][nameparse=darkspace]' % os.path.basename(fragments)]
        alignment = fileStore.getLocalTempFile()
        cactus_call(outfile=alignment,
                    parameters=["cPecanLastz"] + lastZSequenceHandling +
                                self.repeatMaskOptions.lastzOpts.split() +
                                ["--querydepth=keep,nowarn:%i" % (self.repeatMaskOptions.period+3),
                                 "--format=general:name1,zstart1,end1,name2,zstart2+,end2+",
                                 "--markend"])
        return fileStore.writeGlobalFile(alignment)

class MaskCoveredIntervals(RoundedJob):
    def __init__(self, repeatMaskOptions, alignmentsID, queryID):
        RoundedJob.__init__(self, preemptable=True)
        self.repeatMaskOptions = repeatMaskOptions
        self.alignmentsID = alignmentsID
        self.queryID = queryID

    def run(self, fileStore):
        #This runs Bob's covered intervals program, which combins the lastz alignment info into intervals of the query.
        alignments = fileStore.readGlobalFile(self.alignmentsID)
        query = fileStore.readGlobalFile(self.queryID)
        maskInfo = fileStore.getLocalTempFile()
        cactus_call(infile=alignments, outfile=maskInfo,
                    parameters=["cactus_covered_intervals",
                                "--queryoffsets",
                                "--origin=one",
                                "M=%s" % (int(self.repeatMaskOptions.period*2))])

        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
        args = ["--origin=one"]
        if self.repeatMaskOptions.unmaskOutput:
            args.append("--unmask")
        args.append(maskInfo)
        maskedQuery = fileStore.getLocalTempFile()
        cactus_call(infile=query, outfile=maskedQuery,
                    parameters=["cactus_fasta_softmask_intervals.py"] + args)
        tmp = fileStore.writeGlobalFile(maskedQuery)
        return tmp

class LastzRepeatMaskJob(RoundedJob):
    def __init__(self, repeatMaskOptions, queryID, targetIDs):
        RoundedJob.__init__(self, preemptable=True)
        self.repeatMaskOptions = repeatMaskOptions
        self.queryID = queryID
        self.targetIDs = targetIDs

    def run(self, fileStore):
        assert len(self.targetIDs) >= 1
        assert self.repeatMaskOptions.fragment > 1
        queryFile = fileStore.readGlobalFile(self.queryID)

        # chop up input fasta file into into fragments of specified size.  fragments overlap by 
        # half their length.
        fragOutput = fileStore.getLocalTempFile()
        cactus_call(infile=queryFile, outfile=fragOutput,
                    parameters=["cactus_fasta_fragments.py",
                                "--fragment=%s" % str(self.repeatMaskOptions.fragment),
                                "--step=%s" % (str(self.repeatMaskOptions.fragment /2)),
                                "--origin=zero"])
        fragmentsID = fileStore.writeGlobalFile(fragOutput)

        alignmentJob = self.addChild(AlignFastaFragments(repeatMaskOptions=self.repeatMaskOptions, 
                    fragmentsID=fragmentsID, targetIDs=self.targetIDs))

        maskCoveredIntervalsJob = self.addChild(MaskCoveredIntervals(repeatMaskOptions=self.repeatMaskOptions, alignmentsID=alignmentJob.rv(), queryID=self.queryID))
        alignmentJob.addFollowOn(maskCoveredIntervalsJob)

        return maskCoveredIntervalsJob.rv()
