#!/usr/bin/env python
#Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""Script for modifying the scores of pairwise alignments to reflect mapping qualities and then (optionally)
filtering the alignments to remove lower probability alignments.

- Input: set of sequences S and set of pairwise alignments T
- Output: modified set of alignments T' in which scores are replaced with mapping qualities, 
  and optionally filtered to keep only the single most probable alignment per position (the primary alignments).
  This involves chopping up alignments in T to avoid partial overlaps.

- Overview of procedure (top level in Python in this script):
        - Add mirror alignments to T  and ensure alignments are reported with repsect to positive strand of first sequence 
        (this ensures that each alignment is considered on both sequences 
        to which it aligns): C subscript: cactus_mirrorAndOrientAlignments.c
        - Sort alignments in T by coordinates on S: Unix sort
        - Split alignments in T so that they don't partially overlap on S: C subscript: cactus_splitAlignmentOverlaps
            - Each alignment defines an interval on a sequence in S
            - Split alignments into sub-alignments so for any two alignments in the set 
            if they overlap they have the same interval. 
        - Calculate mapping qualities for each alignments and optionally filter alignments, 
        for example to only keep the primary alignment: C subscript: cactus_calculateMappingQualities

"""
import os
import shutil
from toil.lib.bioio import logger
from toil.lib.bioio import system

from cactus.shared.common import cactus_call

def countLines(inputFile):
    with open(inputFile, 'r') as f:
        return sum(1 for line in f)

def mappingQualityRescoring(job, inputAlignmentFileID, 
                            minimumMapQValue, maxAlignmentsPerSite, logLevel):
    """
    Function to rescore and filter alignments by calculating the mapping quality of sub-alignments
    """
    # Get temporary files
    tempAlignmentFile = job.fileStore.getLocalTempFile()
    tempAlignmentFile2 = job.fileStore.getLocalTempFile()
    
    inputAlignmentFile = job.fileStore.readGlobalFile(inputAlignmentFileID)
    
    job.fileStore.logToMaster("Input cigar file has %s lines" % countLines(inputAlignmentFile))
    
    # Mirror and orient alignments
    cactus_call(parameters=["cactus_mirrorAndOrientAlignments",
                             logLevel,
                             inputAlignmentFile,
                             tempAlignmentFile ])
    
    #job.fileStore.logToMaster("Mirrored and oriented cigar file has %s lines" % countLines(tempAlignmentFile))
    
    # Sort
    cactus_call(parameters=[ "cactus_blast_sortAlignmentsByQuery",
                             logLevel,
                             tempAlignmentFile,
                             tempAlignmentFile2])
    
    #job.fileStore.logToMaster("Sorted cigar file has %s lines" % countLines(tempAlignmentFile2))
    
    # Split overlaps 
    cactus_call(parameters=["cactus_splitAlignmentOverlaps",
                             logLevel,
                             tempAlignmentFile2,
                             tempAlignmentFile ])
    
    job.fileStore.logToMaster("Split cigar file has %s lines" % countLines(tempAlignmentFile))
    
    # Calculate mapping qualities
    cactus_call(parameters=["cactus_calculateMappingQualities",
                             logLevel,
                             tempAlignmentFile,
                             tempAlignmentFile2,
                             str(maxAlignmentsPerSite), str(minimumMapQValue) ])
    
    job.fileStore.logToMaster("Filtered, non-overlapping cigar file has %s lines" % countLines(tempAlignmentFile2))
    
    # Now write back alignments results file and return
    return job.fileStore.writeGlobalFile(tempAlignmentFile2)