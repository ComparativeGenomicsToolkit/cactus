#!/usr/bin/env python3

"""
Generate the alignment file that needs to be input to cactus_consolidated

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.job import Job
from toil.statsAndLogging import logger
from toil.lib.bioio import getLogLevelString
from sonLib.bioio import newickTreeParser
import os
from cactus.paf.paf import get_leaf_event_pairs, get_leaves, get_node, get_distances
from cactus.shared.common import cactus_call, getOptionalAttrib


def run_lastz(job, genome_A, genome_B, distance, params):    
    # Create a local temporary file to put the alignments in.
    work_dir = job.fileStore.getLocalTempDir()
    alignment_file = os.path.join(work_dir, 'aln.paf')

    # Get the input files
    genome_a_file = os.path.join(work_dir, 'a.fa')
    genome_b_file = os.path.join(work_dir, 'b.fa')
    job.fileStore.readGlobalFile(genome_A, genome_a_file)
    job.fileStore.readGlobalFile(genome_B, genome_b_file)

    # Get the params to do the alignment
    lastz_params_node = params.find("blast")
    lastz_divergence_node = lastz_params_node.find("divergence")
    divergences = params.find("constants").find("divergences")
    for i in "one", "two", "three", "four", "five":
        if distance <= float(divergences.attrib[i]):
            lastz_params = lastz_divergence_node.attrib[i]
            break
    else:
        lastz_params = lastz_divergence_node.attrib["default"]
    logger.info("For distance {} for genomes {}, {} using {} lastz parameters".format(distance, genome_A,
                                                                                      genome_B, lastz_params))
    
    if getOptionalAttrib(lastz_params_node, 'gpuLastz', typeFn=bool, default=False):
        lastz_bin = 'run_segalign'
        suffix_a, suffix_b = '', ''
    else:
        lastz_bin = 'lastz'
        suffix_a = '[multiple][nameparse=darkspace]'
        suffix_b = '[nameparse=darkspace]'

    # Generate the alignment
    lastz_cmd = [lastz_bin,
                 '{}{}'.format(os.path.basename(genome_a_file), suffix_a),
                 '{}{}'.format(os.path.basename(genome_b_file), suffix_b),
                 '--format=paf:minimap2'] + lastz_params.split(' ')
    # note: it's very important to set the work_dir here, because cactus_call is not able to
    # sort out the mount directory by itself, presumably due to square brackets...
    cactus_call(parameters=lastz_cmd, outfile=alignment_file, work_dir=work_dir)

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)


def run_minimap2(job, genome_A, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    alignment_file = job.fileStore.getLocalTempFile()

    minimap2_params = params.find("blast").attrib["minimap2_params"]
    
    # Generate the alignment
    minimap2_cmd = ['minimap2',
                    '-c',
                    job.fileStore.readGlobalFile(genome_A),
                    job.fileStore.readGlobalFile(genome_B)] + minimap2_params.split(' ')
    cactus_call(parameters=minimap2_cmd, outfile=alignment_file)

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)


def combine_chunks(job, chunked_alignment_files):
    # Make combined alignments file
    alignment_file = job.fileStore.getLocalTempFile()
    for chunk in chunked_alignment_files: # Append each of the chunked alignment files into one file
        cactus_call(parameters=['paf_dechunk', '-i', job.fileStore.readGlobalFile(chunk)],
                    outfile=alignment_file, outappend=True)
        job.fileStore.deleteGlobalFile(chunk) # Cleanup the old files
    return job.fileStore.writeGlobalFile(alignment_file)  # Return the final alignments file copied to the jobstore


def make_chunked_alignments(job, genome_a, genome_b, distance, params):
    lastz_params_node = params.find("blast")
    if getOptionalAttrib(lastz_params_node, 'gpuLastz', typeFn=bool, default=False):
        # wga-gpu has a 6G limit, so we always override
        lastz_params_node.attrib['chunkSize'] = '6000000000'
    lastz_cores = getOptionalAttrib(lastz_params_node, 'cpu', typeFn=int, default=None)
    
    def make_chunks(genome):
        output_chunks_dir = job.fileStore.getLocalTempDir()
        fasta_chunk_cmd = ['fasta_chunk',
                           '-c', params.find("blast").attrib["chunkSize"],
                           '-o', params.find("blast").attrib["overlapSize"],
                           '--dir', output_chunks_dir,
                           job.fileStore.readGlobalFile(genome)]
        cactus_call(parameters=fasta_chunk_cmd)
        return [job.fileStore.writeGlobalFile(os.path.join(output_chunks_dir, chunk), cleanup=True) for chunk in os.listdir(output_chunks_dir)]
    # Chunk each input genome
    chunks_a = make_chunks(genome_a)
    chunks_b = make_chunks(genome_b)

    # Align all chunks from genome_A against all chunks from genome_B
    chunked_alignment_files = []
    for chunk_a in chunks_a:
        for chunk_b in chunks_b:
            mappers = {"lastz":run_lastz, "minimap2":run_minimap2}
            mappingFn = mappers[params.find("blast").attrib["mapper"]]
            chunked_alignment_files.append(job.addChildJobFn(mappingFn, chunk_a, chunk_b, distance, params,
                                                             cores=lastz_cores, disk=4*(chunk_a.size+chunk_b.size)).rv())

    return job.addFollowOnJobFn(combine_chunks, chunked_alignment_files).rv()  # Combine the chunked alignment files


def make_ingroup_to_outgroup_alignments(job, ingroup, outgroup_events, event_tree, event_names_to_sequences, params):
    #  a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    #  align ingroup to first outgroup to produce paf alignments
    outgroup = outgroup_events[0] # The first outgroup
    logger.info("Building alignment between ingroup event: {} and outgroup event: {}".format(ingroup.iD, outgroup.iD))
    alignment = root_job.addChildJobFn(make_chunked_alignments, event_names_to_sequences[outgroup.iD],
                                       event_names_to_sequences[ingroup.iD], 1.0, params,
                                       disk=2*(event_names_to_sequences[ingroup.iD].size+event_names_to_sequences[outgroup.iD].size)).rv()

    #  post process the alignments and recursively generate alignments to remaining outgroups
    return root_job.addFollowOnJobFn(make_ingroup_to_outgroup_alignments2, alignment, ingroup, outgroup_events[1:], event_tree,
                                     event_names_to_sequences, params).rv() if len(outgroup_events) > 1 else alignment


def make_ingroup_to_outgroup_alignments2(job, alignments, ingroup, outgroup_events, event_tree,
                                         event_names_to_sequences, params):
    # a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    # identify all ingroup sub-sequences that remain unaligned longer than a threshold as follows:

    # run paf_to_bed to create a bed of aligned coverage
    bed_file = job.fileStore.getLocalTempFile()  # Get a temporary file to store the bed file in
    messages = cactus_call(parameters=['paf_to_bed', "--excludeAligned", "--binary",
                                       "--minSize", params.find("blast").attrib["trimMinSize"],
                                       "-i", job.fileStore.readGlobalFile(alignments),
                                       "--logLevel", getLogLevelString()],
                                       outfile=bed_file, returnStdErr=True)
    job.fileStore.logToMaster("paf_to_bed event:{}\n{}".format(ingroup.iD, messages[:-1]))  # Log paf_to_bed

    # run fasta_extract to extract unaligned sequences longer than a threshold creating a reduced subset of A
    seq_file = job.fileStore.getLocalTempFile()  # Get a temporary file to store the subsequences in
    ingroup_seq_file = job.fileStore.readGlobalFile(event_names_to_sequences[ingroup.iD])  # The ingroup sequences
    messages = cactus_call(parameters=['fasta_extract', "-i", bed_file, ingroup_seq_file,
                                       "--flank", params.find("blast").attrib["trimFlanking"],
                                       "--logLevel", getLogLevelString()],
                           outfile=seq_file, returnStdErr=True)
    job.fileStore.logToMaster("fasta_extract event:{}\n{}".format(ingroup.iD, messages[:-1]))  # Log fasta_extract

    # replace the ingroup sequences with remaining sequences
    event_names_to_sequences[ingroup.iD] = job.fileStore.writeGlobalFile(seq_file)

    # recursively make alignments with the remaining outgroups
    alignments2 = root_job.addChildJobFn(make_ingroup_to_outgroup_alignments, ingroup, outgroup_events, event_tree,
                                         event_names_to_sequences, params).rv()

    return root_job.addFollowOnJobFn(make_ingroup_to_outgroup_alignments3, ingroup, event_names_to_sequences[ingroup.iD],
                                     alignments, alignments2).rv()


def make_ingroup_to_outgroup_alignments3(job, ingroup, ingroup_seq_file, alignments, alignments2):
    # Merge the together two alignment files

    alignments = job.fileStore.readGlobalFile(alignments)  # Copy the global alignment files locally
    alignments2 = job.fileStore.readGlobalFile(alignments2)

    # use fasta_dechunk to correct the subsequence coordinates of alignments2
    alignments2_corrected = job.fileStore.getLocalTempFile()
    messages = cactus_call(parameters=['paf_dechunk', "-i", alignments2, "--query", "--logLevel", getLogLevelString()],
                           outfile=alignments2_corrected, returnStdErr=True)
    job.fileStore.logToMaster("paf_dechunk event:{}\n{}".format(ingroup.iD, messages[:-1]))  # Log paf_dechunk

    # merge the two alignment files together
    merged_alignments = job.fileStore.getLocalTempFile()
    cactus_call(parameters=['cat', alignments, alignments2_corrected], outfile=merged_alignments)  # Don't bother to log

    # Delete the file containing the subset of ingroup sequences as we don't need it any longer and it takes up space
    job.fileStore.deleteGlobalFile(ingroup_seq_file)

    return job.fileStore.writeGlobalFile(merged_alignments)


def chain_alignments(job, alignment_files, reference_event_name, params):
    # Create a local temporary file to put the alignments in.
    output_alignments_file = job.fileStore.getLocalTempFile()

    # Get temporary file to store concatenated, chained alignments
    chained_alignment_file = job.fileStore.getLocalTempFile()

    # Run the chaining
    for i in alignment_files:
        i = job.fileStore.readGlobalFile(i)  # Copy the global alignment file locally
        j = job.fileStore.getLocalTempFile()  # Get a temporary file to store local alignments and their inverse in
        cactus_call(parameters=['cat', i], outfile=j, outappend=True)
        cactus_call(parameters=['paf_invert', "-i", i], outfile=j, outappend=True)  # Not bothering to log this one
        messages = cactus_call(parameters=['paf_chain', "-i", j,
                                           "--maxGapLength", params.find("blast").attrib["chainMaxGapLength"],
                                           "--chainGapOpen", params.find("blast").attrib["chainGapOpen"],
                                           "--chainGapExtend", params.find("blast").attrib["chainGapExtend"],
                                           "--trimFraction", params.find("blast").attrib["chainTrimFraction"],
                                           "--logLevel", getLogLevelString()],
                               outfile=chained_alignment_file, outappend=True, returnStdErr=True)
        job.fileStore.logToMaster("paf_chain {}\n{}".format(reference_event_name, messages[:-1]))  # Log paf_chain

    # Now tile
    messages = cactus_call(parameters=['paf_tile', "-i", chained_alignment_file, "--logLevel", getLogLevelString()],
                           outfile=output_alignments_file, returnStdErr=True)
    job.fileStore.logToMaster("paf_tile event:{}\n{}".format(reference_event_name, messages[:-1]))  # Log paf_tile

    # Cleanup the old alignment files
    for i in alignment_files:
        job.fileStore.deleteGlobalFile(i)

    return job.fileStore.writeGlobalFile(output_alignments_file)  # Copy back


def make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    # a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    job.fileStore.logToMaster("Parsing species tree: {}".format(event_tree_string))
    event_tree = newickTreeParser(event_tree_string)
    distances = get_distances(event_tree)  # Distances between all pairs of nodes

    ancestor_event = get_node(event_tree, ancestor_event_string)
    ingroup_events = get_leaves(ancestor_event) # Get the set of ingroup events
    job.fileStore.logToMaster("Got ingroup events: {} for ancestor event: {}".format(" ".join([i.iD for i in ingroup_events]),
                                                                       ancestor_event_string))

    outgroup_events = [event for event in get_leaves(event_tree) if event not in ingroup_events]  # Set of outgroups
    outgroup_events.sort(key=lambda outgroup: distances[ancestor_event, outgroup])  # Sort from closest to furthest
    job.fileStore.logToMaster("Got outgroup events: {} for ancestor event: {}".format(" ".join([i.iD for i in outgroup_events]),
                                                                        ancestor_event_string))

    # Calculate the total sequence size
    total_sequence_size = sum(event_names_to_sequences[event.iD].size for event in get_leaves(event_tree))

    # for each pair of ingroups make alignments
    ingroup_alignments = []
    for ingroup, ingroup2, distance_a_b in get_leaf_event_pairs(ancestor_event):
        job.fileStore.logToMaster("Building alignment between event: {} (ingroup) and event: {} (ingroup)".format(ingroup.iD, ingroup2.iD))
        ingroup_alignments.append(root_job.addChildJobFn(make_chunked_alignments, event_names_to_sequences[ingroup.iD],
                                                         event_names_to_sequences[ingroup2.iD], distance_a_b, params,
                                                         disk=2*total_sequence_size).rv())

    # for each ingroup make alignments to the outgroups
    outgroup_alignments = [root_job.addChildJobFn(make_ingroup_to_outgroup_alignments, ingroup, outgroup_events,
                                                  event_tree, event_names_to_sequences, params).rv()
                           for ingroup in ingroup_events] if len(outgroup_events) > 0 else []

    # Now do the chaining
    return root_job.addFollowOnJobFn(chain_alignments, ingroup_alignments + outgroup_alignments,
                                     ancestor_event_string, params, disk=2*total_sequence_size).rv()


# Todo: Sort out right level of logging

# Todo: Write a function to trim the outgroup sequences by coverage
# For each outgroup: Calculate coverage on the outgroup using paf_to_bed
# Todo: ensure no overlap between bed_intervals and allow for flanks

# Extract the aligned subsequences using fasta_extract
# Todo: add option to ignore any region not in input sequences

# Convert alignments to refer to subsequences

# Write back file


