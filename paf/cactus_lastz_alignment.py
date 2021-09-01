#!/usr/bin/env python3

"""
Generate the alignment file the needs to be input to cactus_consolidated

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.common import Toil
from toil.job import Job
from toil.lib.bioio import system
from cactus.shared.common import makeURL

def run_lastz(job, genome_A, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    alignment_file = job.fileStore.getLocalTempFile()

    # Get the params
    ##

    # Generate the alignment
    system("lastz {} {} --format=paf:minimap2 > {}".format(job.fileStore.readGlobalFile(genome_A),
                                                           job.fileStore.readGlobalFile(genome_B),
                                                           alignment_file))

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)

def chain_alignments(job, alignments_A_B, alignments_A_C, alignments_B_C):
    # Create a local temporary file to put the alignments in.
    output_alignments_file = job.fileStore.getLocalTempFile()

    # Run the chaining
    system("cactus_chain {} {} {} {}".format(job.fileStore.readGlobalFile(alignments_A_B),
                                              job.fileStore.readGlobalFile(alignments_A_C),
                                              job.fileStore.readGlobalFile(alignments_B_C),
                                              output_alignments_file))

    return job.fileStore.writeGlobalFile(output_alignments_file)

def run_lastz_alignments(job, genome_A, genome_B, genome_C, params,
                         distance_A_B=1.0, distance_A_C=1.0, distance_B_C=1.0):
    # Generate the alignments
    alignments_A_B = job.addChildJobFn(run_lastz, genome_A, genome_B, distance_A_B, params).rv()  # Generate alignments A to B
    alignments_A_C = job.addChildJobFn(run_lastz, genome_A, genome_C, distance_A_C, params).rv()  # Generate alignments A to C
    alignments_B_C = job.addChildJobFn(run_lastz, genome_B, genome_C, distance_B_C, params).rv()  # Generate alignments B to C
    # Now do the chaining
    return job.addFollowOnJobFn(chain_alignments, alignments_A_B, alignments_A_C, alignments_B_C).rv()

if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("genome_A", help="Ingroup genomes A file [fasta]")
    parser.add_argument("genome_B", help="Ingroup genomes B file [fasta]")
    parser.add_argument("genome_C", help="Outgroup genomes C file [fasta]")
    parser.add_argument("output_file", help="Output alignments file [paf]")
    parser.add_argument("--params", help="Cactus parameters file")

    options = parser.parse_args()
    options.clean = "always"
    with Toil(options) as toil:
        output_file_ID = toil.start(Job.wrapJobFn(run_lastz_alignments,
                                    toil.importFile(makeURL(options.genome_A)),
                                    toil.importFile(makeURL(options.genome_B)),
                                    toil.importFile(makeURL(options.genome_C)),
                                    options))
        toil.exportFile(output_file_ID, makeURL(options.output_file))
