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
from sonLib.bioio import system

from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call

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


class AlignFastaFragments(Job):
    def __init__(self, repeatMaskOptions, fragmentsID, targetID):
        if hasattr(targetID, "size"):
            memory = 2*(fragmentsID.size + targetID.size)
            disk = 2*(fragmentsID.size + targetID.size)
        else:
            memory = None
            disk = None
        Job.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.repeatMaskOptions = repeatMaskOptions
        self.fragmentsID = fragmentsID
        self.targetID = targetID
    def run(self, fileStore):
        # Align each fragment against a chunk of the input sequence.  Each time a fragment aligns to a base
        # in the sequence, that base's match count is incremented.
        # the plus three for the period parameter is a fudge to ensure sufficient alignments are found

        fragments = fileStore.readGlobalFile(self.fragmentsID)
        target = fileStore.readGlobalFile(self.targetID)
        lastZSequenceHandling  = '%s[multiple][nameparse=darkspace] %s[nameparse=darkspace] ' % (os.path.basename(target), os.path.basename(fragments))
        if self.repeatMaskOptions.unmaskInput:
            lastZSequenceHandling  = '%s[multiple,unmask][nameparse=darkspace] %s[unmask][nameparse=darkspace] ' % (os.path.basename(target), os.path.basename(fragments))
        alignment = fileStore.getLocalTempFile()
        cactus_call(outfile=alignment,
                    parameters=["cPecanLastz", lastZSequenceHandling,
                                self.repeatMaskOptions.lastzOpts,
                                "--querydepth=keep,nowarn:%i --format=general:name1,zstart1,end1,name2,zstart2+,end2+ --markend" % (self.repeatMaskOptions.period+3)])
        return fileStore.writeGlobalFile(alignment)

class CollateAlignments(Job):
    def __init__(self, alignmentIDs):
        if hasattr(alignmentIDs[0], "size"):
            disk = 2*sum([alignmentID.size for alignmentID in alignmentIDs])
        else:
            disk = None
        Job.__init__(self, disk=disk, preemptable=True)
        self.alignmentIDs = alignmentIDs
    def run(self, fileStore):
        alignments = [fileStore.readGlobalFile(alignmentID) for alignmentID in self.alignmentIDs]

        #Sort the alignments by the start position of the alignment in the chunk
        #being repeat-masked. These will be out of order due to the parallelization
        sortedAlignments = fileStore.getLocalTempFile()
        system("cat %s | awk '@include \"join\";{split($4,a,\"_\"); $5 += a[length(a)]; $6 += a[length(a)]; $4 = join(a, 1, length(a) - 1, \"_\"); print $0}' | sort -k4,4 -k5,5n > %s" % (" ".join(alignments), sortedAlignments))
        return fileStore.writeGlobalFile(sortedAlignments)

class MaskCoveredIntervals(Job):
    def __init__(self, repeatMaskOptions, alignmentsID, queryID):
        Job.__init__(self, preemptable=True)
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
                                "--origin=one",
                                "M=%s" % (int(self.repeatMaskOptions.period*2))])

        #open(maskInfoFile, "w").close()
        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
        unmaskString = ""
        if self.repeatMaskOptions.unmaskOutput:
            unmaskString = "--unmask"
        maskedQuery = fileStore.getLocalTempFile()
        cactus_call(infile=query, outfile=maskedQuery,
                    parameters=["cactus_fasta_softmask_intervals.py",
                                "--origin=one",
                                unmaskString,
                                maskInfo])
        return fileStore.writeGlobalFile(maskedQuery)

class LastzRepeatMaskJob(Job):
    def __init__(self, repeatMaskOptions, queryID, targetIDs):
        Job.__init__(self, preemptable=True)
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

        alignmentJobs = []
        alignmentIDs = []
        for targetID in self.targetIDs:
            alignmentJob = self.addChild(AlignFastaFragments(repeatMaskOptions=self.repeatMaskOptions, 
                    fragmentsID=fragmentsID, targetID=targetID))
            alignmentIDs.append(alignmentJob.rv())
            alignmentJobs.append(alignmentJob)

        mergeAlignmentsJob = self.addChild(CollateAlignments(alignmentIDs=alignmentIDs))
        for alignmentJob in alignmentJobs:
            alignmentJob.addFollowOn(mergeAlignmentsJob)

        maskCoveredIntervalsJob = self.addChild(MaskCoveredIntervals(repeatMaskOptions=self.repeatMaskOptions, alignmentsID = mergeAlignmentsJob.rv(), queryID=self.queryID))
        mergeAlignmentsJob.addFollowOn(maskCoveredIntervalsJob)

        return maskCoveredIntervalsJob.rv()
