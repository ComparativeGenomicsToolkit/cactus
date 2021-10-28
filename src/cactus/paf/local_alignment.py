#!/usr/bin/env python3

"""
Generate the alignment file that needs to be input to cactus_consolidated

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.lib.bioio import system
from toil.statsAndLogging import logger
from sonLib.bioio import newickTreeParser
import os
from cactus.paf.paf import get_leaf_event_pairs, get_subtree_nodes, get_leaves, get_node

def run_lastz(job, genome_A, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    alignment_file = job.fileStore.getLocalTempFile()

    # Get the params to do the alignment
    lastz_params_node = params.find("blast").find("divergence")
    divergences = params.find("constants").find("divergences")
    for i in "one", "two", "three", "four", "five":
        if distance <= float(divergences.attrib[i]):
            lastz_params = lastz_params_node.attrib[i]
            break
    else:
        lastz_params = lastz_params_node.attrib["default"]
    logger.info("For distance {} for genomes {}, {} using {} lastz parameters".format(distance, genome_A,
                                                                                      genome_B, lastz_params))

    # Generate the alignment
    system("lastz {}[multiple][nameparse=darkspace] {}[nameparse=darkspace] {} --format=paf:minimap2 > {}".format(job.fileStore.readGlobalFile(genome_A),
                                                           job.fileStore.readGlobalFile(genome_B),
                                                           lastz_params, alignment_file))

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)

def run_minimap2(job, genome_A, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    alignment_file = job.fileStore.getLocalTempFile()

    # Generate the alignment
    system("minimap2 -c {} {} {} > {}".format(job.fileStore.readGlobalFile(genome_A),
                                                                           job.fileStore.readGlobalFile(genome_B),
                                                                           params.find("blast").attrib["minimap2_params"],
                                                                           alignment_file))  # Minimap2 must be installed

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)


def combine_chunks(job, chunked_alignment_files):
    # Make combined alignments file
    alignment_file = job.fileStore.getLocalTempFile()
    for chunk in chunked_alignment_files: # Append each of the chunked alignment files into one file
        system("paf_dechunk -i {} >> {}".format(job.fileStore.readGlobalFile(chunk), alignment_file))
        job.fileStore.deleteGlobalFile(chunk) # Cleanup the old files
    return job.fileStore.writeGlobalFile(alignment_file)  # Return the final alignments file copied to the jobstore


def make_chunked_alignments(job, genome_a, genome_b, distance, params):
    def make_chunks(genome):
        output_chunks_dir = job.fileStore.getLocalTempDir()
        system("fasta_chunk -c {} -o {} --dir {} {}".format(int(params.find("blast").attrib["chunkSize"]),
                                                            int(params.find("blast").attrib["overlapSize"]),
                                                            output_chunks_dir, job.fileStore.readGlobalFile(genome)))
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
            chunked_alignment_files.append(job.addChildJobFn(mappingFn, chunk_a, chunk_b, distance, params).rv())

    return job.addFollowOnJobFn(combine_chunks, chunked_alignment_files).rv()  # Combine the chunked alignment files


def chain_alignments(job, alignment_files):
    # Create a local temporary file to put the alignments in.
    output_alignments_file = job.fileStore.getLocalTempFile()

    # Copy the alignment files locally
    local_alignment_files = [job.fileStore.readGlobalFile(i) for i in alignment_files]

    # Run the chaining
    system("cactus_chain {} > {}".format(" ".join(local_alignment_files), output_alignments_file))

    # Cleanup the old alignment files
    for i in alignment_files:
        job.fileStore.deleteGlobalFile(i)

    return job.fileStore.writeGlobalFile(output_alignments_file)  # Copy back


def make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    logger.info("Parsing species tree: {}".format(event_tree_string))
    event_tree = newickTreeParser(event_tree_string)

    # Make a map of event names to nodes in the event tree
    event_names_to_events = { node.iD : node for node in get_subtree_nodes(event_tree) }

    # Check we have a sequence file for all events
    for event in get_leaves(event_tree):
        if event.iD not in event_names_to_sequences:
            raise RuntimeError("No sequence found for event (aka node) {}".format(event.iD))

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
            alignment = job.addChildJobFn(make_chunked_alignments, event_names_to_sequences[event_a.iD],
                                          event_names_to_sequences[event_b.iD], distance_a_b, params).rv()
            alignments.append(alignment)

    # Now do the chaining
    return job.addFollowOnJobFn(chain_alignments, alignments).rv()
