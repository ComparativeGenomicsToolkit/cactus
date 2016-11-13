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
from cactus.shared.bioio import catFiles

from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call


def lastzRepeatMaskJob(job, queryID, targetIDs,
                       fragment=200,
                       minPeriod=10,
                       lastzOpts="",
                       unmaskInput=False,
                       unmaskOutput=False,
                       proportionSampled=1.0):
                       
                       
    
    assert len(targetIDs) >= 1
    assert fragment > 1
    queryFile = job.fileStore.readGlobalFile(queryID)
    targetFiles = [job.fileStore.readGlobalFile(fileID) for fileID in targetIDs]
    
    #Adjust the period parameter using the amount of genome sampled
    period = max(1, round(proportionSampled * minPeriod))
    
    # make sure fragment size is even so they can overlap by exactly one half. 
    if fragment % 2:
        fragment += 1
    
    maskInfoFile = job.fileStore.getLocalTempFile()
    targetFile = job.fileStore.getLocalTempFile()

    #Make temporary target file, if more than one file
    catFiles(targetFiles, targetFile)

    # chop up input fasta file into into fragments of specified size.  fragments overlap by 
    # half their length.
    fragOutput = job.fileStore.getLocalTempFile()
    cactus_call(tool="quay.io/adderan/cactus", infile=queryFile, outfile=fragOutput,
                parameters=["cactus_fasta_fragments.py",
                            "--fragment=%s" % str(fragment),
                            "--step=%s" % (str(fragment /2)),
                            "--origin=zero"])

    # lastz each fragment against the entire input sequence.  Each time a fragment aligns to a base
    # in the sequence, that base's match count is incremented.
    # the plus three for the period parameter is a fudge to ensure sufficient alignments are found
    lastZSequenceHandling  = '[multiple][nameparse=darkspace] /dev/stdin[nameparse=darkspace] '
    if unmaskInput:
        lastZSequenceHandling  = '[multiple,unmask][nameparse=darkspace] /dev/stdin[unmask][nameparse=darkspace] '
    lastzOutput = job.fileStore.getLocalTempFile()
    cactus_call(tool="quay.io/adderan/cpecan-lastz", infile=fragOutput, outfile=lastzOutput,
                parameters=["%s%s" % (os.path.basename(targetFile), lastZSequenceHandling),
                            lastzOpts,
                            "--querydepth=keep,nowarn:%i --format=general:name1,zstart1,end1,name2,zstart2+,end2+ --markend" % (period+3)])


    #This runs Bob's covered intervals program, which combins the lastz alignment info into intervals of the query.
    cactus_call(tool="quay.io/adderan/cactus", infile=lastzOutput, outfile=maskInfoFile,
                parameters=["cactus_covered_intervals",
                            "--queryoffsets",
                            "--origin=one",
                            "M=%s" % (int(period*2))])

    #open(maskInfoFile, "w").close()
    # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
    # we finish by applying these intervals to the input file, to produce the final, softmasked output. 
    unmaskString = ""
    if unmaskOutput:
        unmaskString = "--unmask"
    outputFile = job.fileStore.getLocalTempFile()
    cactus_call(tool="quay.io/adderan/cactus", infile=queryFile, outfile=outputFile,
                parameters=["cactus_fasta_softmask_intervals.py",
                            "--origin=one",
                            unmaskString,
                            maskInfoFile])
    return job.fileStore.writeGlobalFile(outputFile)
