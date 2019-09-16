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
from cactus.shared.common import cactus_call

def countLines(inputFile):
    with open(inputFile, 'r') as f:
        return sum(1 for line in f)

def mappingQualityRescoring(job, inputAlignmentFileID, 
                            minimumMapQValue, maxAlignmentsPerSite, alpha, logLevel):
    """
    Function to rescore and filter alignments by calculating the mapping quality of sub-alignments
    
    Returns primary alignments and secondary alignments in two separate files.
    """
    inputAlignmentFile = job.fileStore.readGlobalFile(inputAlignmentFileID)
    
    job.fileStore.logToMaster("Input cigar file has %s lines" % countLines(inputAlignmentFile))
    
    # Get temporary file
    assert maxAlignmentsPerSite >= 1
    tempAlignmentFiles = [job.fileStore.getLocalTempFile() for i in xrange(maxAlignmentsPerSite)]
    
    # Mirror and orient alignments, sort, split overlaps and calculate mapping qualities
    cactus_call(parameters=[["cat", inputAlignmentFile],
                            ["cactus_mirrorAndOrientAlignments", logLevel],
                            ["sort", "-T{}".format(job.fileStore.getLocalTempDir()), "-k6,6", "-k7,7n", "-k8,8n"], # This sorts by coordinate
                            ["uniq"], # This eliminates any annoying duplicates if lastz reports the alignment in both orientations
                            ["cactus_splitAlignmentOverlaps", logLevel],
                            ["cactus_calculateMappingQualities", logLevel, str(maxAlignmentsPerSite),
                             str(minimumMapQValue), str(alpha)] + tempAlignmentFiles])

    # Merge together the output files in order
    secondaryTempAlignmentFile = job.fileStore.getLocalTempFile()
    if len(tempAlignmentFiles) > 1:
        cactus_call(parameters=[["cat" ] + tempAlignmentFiles[1:]], outfile=secondaryTempAlignmentFile)

    job.fileStore.logToMaster("Filtered, non-overlapping primary cigar file has %s lines" % countLines(tempAlignmentFiles[0]))
    job.fileStore.logToMaster("Filtered, non-overlapping secondary cigar file has %s lines" % countLines(secondaryTempAlignmentFile))

    # Now write back alignments results file and return
    return job.fileStore.writeGlobalFile(tempAlignmentFiles[0]), job.fileStore.writeGlobalFile(secondaryTempAlignmentFile)
