#!/usr/bin/env python3

"""
Generate the alignment file that needs to be input to cactus_consolidated

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.common import Toil
from toil.job import Job
from toil.lib.bioio import system
from toil.statsAndLogging import logger
from cactus.shared.common import makeURL
from sonLib.bioio import newickTreeParser
import xml.etree.ElementTree as ET
from impl.paf import get_leaf_event_pairs, get_subtree_nodes, get_leaves, get_node

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
    system("lastz {} {} {} --format=paf:minimap2 > {}".format(job.fileStore.readGlobalFile(genome_A),
                                                           job.fileStore.readGlobalFile(genome_B),
                                                           lastz_params, alignment_file))

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file)

def chain_alignments(job, alignment_files):
    # Create a local temporary file to put the alignments in.
    output_alignments_file = job.fileStore.getLocalTempFile()

    # Copy the alignment files locally
    local_alignment_files = [ job.fileStore.readGlobalFile(i) for i in alignment_files ]

    # Run the chaining
    system("cactus_chain {} > {}".format(" ".join(local_alignment_files), output_alignments_file))

    return job.fileStore.writeGlobalFile(output_alignments_file)  # Copy back

def make_alignments(job, event_tree, events_to_sequences, ancestor_event, params):
    ingroup_events = get_leaves(ancestor_event) # Get the set of ingroup events
    logger.info("Got ingroup events: {}".format(" ".join([ i.iD for i in ingroup_events ])))

    # Get pairs of sequences in the tree and their MRCA
    alignments = []
    for event_a, event_b, distance_a_b in get_leaf_event_pairs(event_tree):
        if event_a in ingroup_events or event_b in ingroup_events: # If either is an ingroup we align them
            logger.info("Building alignment between event: "
                        "{} (ingroup:{}) and event: {} (ingroup:{})".format(event_a.iD, event_a in ingroup_events,
                                                                            event_b.iD, event_b in ingroup_events))
            alignment = job.addChildJobFn(run_lastz, events_to_sequences[event_a],
                                          events_to_sequences[event_b], distance_a_b, params).rv()
            alignments.append(alignment)

    # Now do the chaining
    return job.addFollowOnJobFn(chain_alignments, alignments).rv()

if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("sequences", help="[eventName fastaFile]xN: [Required] The sequence file for each event(node) "
                                          "in the species tree")
    parser.add_argument("--speciesTree", help="[Required] The species tree", required=True)
    parser.add_argument("--referenceEvent", help="[Required] The internal node (aka event) in the species tree that "
                        " we are attempting to align back to", required=True)
    parser.add_argument("--outputFile", help="Output alignments file [paf]", default="output.paf")
    parser.add_argument("--params", help="[Required] Cactus parameters file", required=True)

    options = parser.parse_args()
    options.clean = "always"
    with Toil(options) as toil:
        # Parse the config file
        logger.info("Parsing parameters file: {}".format(options.params))
        params = ET.parse(options.params).getroot()

        # Parse the event tree
        logger.info("Parsing species tree: {}".format(options.speciesTree))
        event_tree = newickTreeParser(options.speciesTree)

        # Make a map of event names to nodes in the event tree
        event_names_to_events = { node.iD : node for node in get_subtree_nodes(event_tree) }

        # Build a map of events to sequence files
        options.sequences = options.sequences.split() # Tokenize into words
        logger.info("Parsing sequence files for nodes in the input tree: {}".format(options.sequences))
        if len(options.sequences) % 2 != 0: # Not even
            raise RuntimeError("Mismatched number of events to sequences")
        events_to_sequences = {}
        for i in range(0, len(options.sequences), 2):
            event_name, sequence_file = options.sequences[i], options.sequences[i+1]
            logger.info("Parsing sequence file: {} for species node: {}\n".format(sequence_file, event_name))
            event = event_names_to_events[event_name]
            if event is None:
                raise RuntimeError("Event (aka node) {} not found in the "
                                   "given species tree: {}".format(event_name, options.speciesTree))
            events_to_sequences[event] = toil.importFile(makeURL(sequence_file))

        # Run the actual toil workflow to create the alignments
        output_file_ID = toil.start(Job.wrapJobFn(make_alignments, event_tree, events_to_sequences,
                                                  get_node(event_tree, options.referenceEvent), params))

        # Now export the final file
        toil.exportFile(output_file_ID, makeURL(options.outputFile))
