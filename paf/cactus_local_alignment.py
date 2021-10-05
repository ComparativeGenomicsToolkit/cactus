#!/usr/bin/env python3

"""
Generate the alignment file that needs to be input to cactus_consolidated

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.common import Toil
from toil.job import Job
from toil.statsAndLogging import logger
from cactus.shared.common import makeURL
import xml.etree.ElementTree as ET
from cactus.paf.local_alignment import make_paf_alignments

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

        # Build a map of events to sequence files
        options.sequences = options.sequences.split() # Tokenize into words
        logger.info("Parsing sequence files for nodes in the input tree: {}".format(options.sequences))
        if len(options.sequences) % 2 != 0: # Not even
            raise RuntimeError("Mismatched number of events to sequences")
        event_names_to_sequences = {}
        for i in range(0, len(options.sequences), 2):
            event_name, sequence_file = options.sequences[i], options.sequences[i+1]
            logger.info("Parsing sequence file: {} for species node: {}\n".format(sequence_file, event_name))
            event_names_to_sequences[event_name] = toil.importFile(makeURL(sequence_file))

        # Run the actual toil workflow to create the alignments
        output_file_ID = toil.start(Job.wrapJobFn(make_paf_alignments, options.speciesTree, event_names_to_sequences,
                                                  options.referenceEvent, params))

        # Now export the final file
        toil.exportFile(output_file_ID, makeURL(options.outputFile))
