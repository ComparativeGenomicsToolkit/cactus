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
from cactus.paf.paf import get_leaf_event_pairs, get_subtree_nodes, get_leaves, get_node
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
            mappers = { "lastz":run_lastz, "minimap2":run_minimap2}
            mappingFn = mappers[params.find("blast").attrib["mapper"]]
            chunked_alignment_files.append(job.addChildJobFn(mappingFn, chunk_a, chunk_b, distance, params,
                                                             cores=lastz_cores, disk=4*(chunk_a.size+chunk_b.size)).rv())

    return job.addFollowOnJobFn(combine_chunks, chunked_alignment_files).rv()  # Combine the chunked alignment files


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
                                "--logLevel", getLogLevelString()], outfile=chained_alignment_file, outappend=True, returnStdErr=True)
        job.fileStore.logToMaster("paf_chain {}\n{}".format(reference_event_name, messages[:-1]))  # Log paf_chain messages

    # Now tile
    messages = cactus_call(parameters=['paf_tile', "-i", chained_alignment_file, "--logLevel", getLogLevelString()],
                           outfile=output_alignments_file, returnStdErr=True)
    job.fileStore.logToMaster("paf_tile event:{}\n{}".format(reference_event_name, messages[:-1]))  # Log paf_tile messages

    # Cleanup the old alignment files
    for i in alignment_files:
        job.fileStore.deleteGlobalFile(i)

    return job.fileStore.writeGlobalFile(output_alignments_file)  # Copy back


def make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    # a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)
    logger.info("Parsing species tree: {}".format(event_tree_string))
    event_tree = newickTreeParser(event_tree_string)

    # Make a map of event names to nodes in the event tree
    event_names_to_events = { node.iD : node for node in get_subtree_nodes(event_tree) }

    # Check we have a sequence file for all events
    total_sequence_size = 0
    for event in get_leaves(event_tree):
        if event.iD not in event_names_to_sequences:
            raise RuntimeError("No sequence found for event (aka node) {}".format(event.iD))
        total_sequence_size += event_names_to_sequences[event.iD].size

    ancestor_event = get_node(event_tree, ancestor_event_string)
    ingroup_events = get_leaves(ancestor_event) # Get the set of ingroup events
    logger.info("Got ingroup events: {} for ancestor event: {}".format(" ".join([ i.iD for i in ingroup_events ]), ancestor_event_string))

    # Get pairs of sequences in the tree and their MRCA
    alignments = []
    for event_a, event_b, distance_a_b in get_leaf_event_pairs(event_tree):
        if event_a in ingroup_events or event_b in ingroup_events: # If either is an ingroup we align them
            logger.info("Building alignment between event: "
                        "{} (ingroup:{}) and event: {} (ingroup:{})".format(event_a.iD, event_a in ingroup_events,
                                                                            event_b.iD, event_b in ingroup_events))
            alignment = root_job.addChildJobFn(make_chunked_alignments, event_names_to_sequences[event_a.iD],
                                               event_names_to_sequences[event_b.iD], distance_a_b, params,
                                               disk=2*total_sequence_size).rv()
            alignments.append(alignment)

    # Now do the chaining
    return root_job.addFollowOnJobFn(chain_alignments, alignments, ancestor_event_string, params,
                                     disk=2*total_sequence_size).rv()
