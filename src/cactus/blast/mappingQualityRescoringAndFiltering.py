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

def mappingQualityRescoring(job, inputAlignmentFile, outputAlignmentFile, cactusWorkflowArguments):
    """
    Function to rescore and filter alignments
    """
    # Get temporary files
    tempAlignmentFile, tempAlignmentFile2 = job.fileStore.getLocalTemporaryFile(), 
    job.fileStore.getLocalTemporaryFile()
    
    # Mirror alignments
    cactus_call(parameters=["cactus_mirrorAndOrientAlignments",
                             getLogLevelString(),
                             inputAlignmentFile,
                             tempAlignmentFile ])
    
    # Sort
    cactus_call(parameters=[ "cactus_blast_sortAlignmentsByQuery",
                             getLogLevelString(),
                             tempAlignmentFile,
                             tempAlignmentFile2])
    
    # Split overlaps 
    cactus_call(parameters=["cactus_splitAlignmentOverlaps",
                             getLogLevelString(),
                             tempAlignmentFile2,
                             tempAlignmentFile ])
    
    # Calculate mapping qualities
    cactus_call(parameters=["cactus_calculateMappingQualities",
                             getLogLevelString(),
                             tempAlignmentFile,
                             outputAlignmentFile ])
    