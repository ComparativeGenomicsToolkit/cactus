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
from cactus.paf.paf import get_event_pairs, get_leaves, get_node, get_distances
from cactus.shared.common import cactus_call, getOptionalAttrib
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers

def run_lastz(job, name_A, genome_A, name_B, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    work_dir = job.fileStore.getLocalTempDir()
    alignment_file = os.path.join(work_dir, '{}_{}.paf'.format(name_A, name_B))

    # Get the input files
    genome_a_file = os.path.join(work_dir, '{}.fa'.format(name_A))
    genome_b_file = os.path.join(work_dir, '{}.fa'.format(name_B))
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

    gpu = getOptionalAttrib(lastz_params_node, 'gpu', typeFn=int, default=0)
    if gpu:
        lastz_bin = 'run_segalign'
        suffix_a, suffix_b = '', ''
        assert gpu > 0
        lastz_params += ' --num_gpu {}'.format(gpu)
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
    segalign_messages = cactus_call(parameters=lastz_cmd, outfile=alignment_file, work_dir=work_dir, returnStdErr=gpu>0, gpus=gpu)

    if gpu:
        # run_segalign can crash and still exit 0, so it's worth taking a moment to check the log for errors
        segalign_messages = segalign_messages.lower()
        for line in segalign_messages.split("\n"):
            if not line.startswith("signals delivered"):
                for keyword in ['terminate', 'error', 'fail', 'assert', 'signal', 'abort', 'segmentation', 'sigsegv', 'kill']:
                    if keyword in line and 'signals' not in line:
                        job.fileStore.logToMaster("Segalign offending line: " + line)  # Log the messages
                        raise RuntimeError('{} exited 0 but keyword "{}" found in stderr'.format(lastz_cmd, keyword))

    # Add in mismatches so we can trim alignments later
    alignment_file_with_mismatches = os.path.join(work_dir, '{}_{}_mismatches.paf'.format(name_A, name_B))
    cactus_call(parameters=['paffy', 'add_mismatches', "-i", alignment_file, genome_a_file, genome_b_file ],
                outfile=alignment_file_with_mismatches)

    # Return the alignment file
    return job.fileStore.writeGlobalFile(alignment_file_with_mismatches)


def run_minimap2(job, name_A, genome_A, name_B, genome_B, distance, params):
    # Create a local temporary file to put the alignments in.
    work_dir = job.fileStore.getLocalTempDir()
    alignment_file = os.path.join(work_dir, '{}_{}.paf'.format(name_A, name_B))

    # Get the input files
    genome_a_file = os.path.join(work_dir, '{}.fa'.format(name_A))
    genome_b_file = os.path.join(work_dir, '{}.fa'.format(name_B))
    job.fileStore.readGlobalFile(genome_A, genome_a_file)
    job.fileStore.readGlobalFile(genome_B, genome_b_file)

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
        cactus_call(parameters=['paffy', 'dechunk', '-i', job.fileStore.readGlobalFile(chunk)],
                    outfile=alignment_file, outappend=True)
        job.fileStore.deleteGlobalFile(chunk) # Cleanup the old files
    return job.fileStore.writeGlobalFile(alignment_file)  # Return the final alignments file copied to the jobstore


def make_chunked_alignments(job, event_a, genome_a, event_b, genome_b, distance, params):
    lastz_params_node = params.find("blast")
    gpu = getOptionalAttrib(lastz_params_node, 'gpu', typeFn=int, default=0)
    if gpu:
        # wga-gpu has a 6G limit, so we always override
        lastz_params_node.attrib['chunkSize'] = '6000000000'
    lastz_cores = getOptionalAttrib(lastz_params_node, 'cpu', typeFn=int, default=None)
    
    def make_chunks(genome):
        output_chunks_dir = job.fileStore.getLocalTempDir()
        fasta_chunk_cmd = ['faffy', 'chunk',
                           '-c', params.find("blast").attrib["chunkSize"],
                           '-o', params.find("blast").attrib["overlapSize"],
                           '--dir', output_chunks_dir,
                           job.fileStore.readGlobalFile(genome)]
        cactus_call(parameters=fasta_chunk_cmd)
        return [job.fileStore.writeGlobalFile(os.path.join(output_chunks_dir, chunk), cleanup=True)
                for chunk in os.listdir(output_chunks_dir)]
    # Chunk each input genome
    chunks_a = make_chunks(genome_a)
    chunks_b = make_chunks(genome_b)

    # Align all chunks from genome_A against all chunks from genome_B
    accelerators = 'cuda:{}'.format(gpu) if gpu else None
    chunked_alignment_files = []
    for i, chunk_a in enumerate(chunks_a):
        for j, chunk_b in enumerate(chunks_b):
            mappers = { "lastz":run_lastz, "minimap2":run_minimap2}
            mappingFn = mappers[params.find("blast").attrib["mapper"]]
            chunked_alignment_files.append(job.addChildJobFn(mappingFn, '{}_{}'.format(event_a, i), chunk_a,
                                                             '{}_{}'.format(event_b, j), chunk_b, distance, params,
                                                             cores=lastz_cores, disk=4*(chunk_a.size+chunk_b.size),
                                                             memory=4*(chunk_a.size+chunk_b.size),
                                                             accelerators=accelerators).rv())

    return job.addFollowOnJobFn(combine_chunks, chunked_alignment_files).rv()  # Combine the chunked alignment files


def make_ingroup_to_outgroup_alignments_1(job, ingroup_event, outgroup_events, event_names_to_sequences, distances, params):
    #  a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    #  align ingroup to first outgroup to produce paf alignments
    outgroup = outgroup_events[0] # The first outgroup
    logger.info("Building alignment between ingroup event: {} and outgroup event: {}".format(ingroup_event.iD, outgroup.iD))
    alignment = root_job.addChildJobFn(make_chunked_alignments,
                                       outgroup.iD, event_names_to_sequences[outgroup.iD],
                                       ingroup_event.iD, event_names_to_sequences[ingroup_event.iD], distances[ingroup_event, outgroup], params,
                                       disk=2*(event_names_to_sequences[ingroup_event.iD].size+event_names_to_sequences[outgroup.iD].size)).rv()

    #  post process the alignments and recursively generate alignments to remaining outgroups
    return root_job.addFollowOnJobFn(make_ingroup_to_outgroup_alignments_2, alignment, ingroup_event, outgroup_events[1:],
                                     event_names_to_sequences, distances, params).rv() if len(outgroup_events) > 1 else alignment


def make_ingroup_to_outgroup_alignments_2(job, alignments, ingroup_event, outgroup_events,
                                          event_names_to_sequences, distances, params):
    # a job should never set its own follow-on, so we hang everything off root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    # identify all ingroup sub-sequences that remain unaligned longer than a threshold as follows:

    # run paffy to_bed to create a bed of aligned coverage
    bed_file = job.fileStore.getLocalTempFile()  # Get a temporary file to store the bed file in
    messages = cactus_call(parameters=['paffy', 'to_bed', "--excludeAligned", "--binary",
                                       "--minSize", params.find("blast").attrib["trimMinSize"],
                                       "-i", job.fileStore.readGlobalFile(alignments),
                                       "--logLevel", getLogLevelString()],
                                       outfile=bed_file, returnStdErr=True)
    logger.info("paffy to_bed event:{}\n{}".format(ingroup_event.iD, messages[:-1]))  # Log paffy to_bed

    # run faffy extract to extract unaligned sequences longer than a threshold creating a reduced subset of A
    seq_file = job.fileStore.getLocalTempFile()  # Get a temporary file to store the subsequences in
    ingroup_seq_file = job.fileStore.readGlobalFile(event_names_to_sequences[ingroup_event.iD])  # The ingroup sequences
    messages = cactus_call(parameters=['faffy', 'extract', "-i", bed_file, ingroup_seq_file,
                                       "--flank", params.find("blast").attrib["trimFlanking"],
                                       "--logLevel", getLogLevelString()],
                           outfile=seq_file, returnStdErr=True)
    logger.info("faffy extract event:{}\n{}".format(ingroup_event.iD, messages[:-1]))  # Log faffy extract

    # replace the ingroup sequences with remaining sequences
    event_names_to_sequences[ingroup_event.iD] = job.fileStore.writeGlobalFile(seq_file)

    # recursively make alignments with the remaining outgroups
    alignments2 = root_job.addChildJobFn(make_ingroup_to_outgroup_alignments_1, ingroup_event, outgroup_events,
                                         event_names_to_sequences, distances, params).rv()

    return root_job.addFollowOnJobFn(make_ingroup_to_outgroup_alignments_3, ingroup_event, event_names_to_sequences[ingroup_event.iD],
                                     alignments, alignments2).rv()


def make_ingroup_to_outgroup_alignments_3(job, ingroup_event, ingroup_seq_file, alignments, alignments2):
    # Merge the together two alignment files

    alignments = job.fileStore.readGlobalFile(alignments)  # Copy the global alignment files locally
    alignments2 = job.fileStore.readGlobalFile(alignments2)

    # use paffy dechunk to correct the subsequence coordinates of alignments2
    alignments2_corrected = job.fileStore.getLocalTempFile()
    messages = cactus_call(parameters=['paffy', 'dechunk', "-i", alignments2, "--query", "--logLevel", getLogLevelString()],
                           outfile=alignments2_corrected, returnStdErr=True)
    logger.info("paffy dechunk event:{}\n{}".format(ingroup_event.iD, messages[:-1]))  # Log paffy dechunk

    # merge the two alignment files together
    merged_alignments = job.fileStore.getLocalTempFile()
    cactus_call(parameters=['cat', alignments, alignments2_corrected], outfile=merged_alignments)  # Don't bother to log

    # Delete the file containing the subset of ingroup sequences as we don't need it any longer and it takes up space
    job.fileStore.deleteGlobalFile(ingroup_seq_file)

    return job.fileStore.writeGlobalFile(merged_alignments)


def chain_alignments(job, alignment_files, reference_event_name, params):
    # The following recapitulates the pipeline showing in paffy/tests/paf_pipeline_test.sh

    # Create a local temporary file to put the alignments in.
    output_alignments_file = job.fileStore.getLocalTempFile()

    # Get temporary files to store intermediate alignments in
    j = job.fileStore.getLocalTempFile()
    k = job.fileStore.getLocalTempFile()

    # Run the chaining
    for i in alignment_files:
        # Copy the alignments from the job store
        i = job.fileStore.readGlobalFile(i)  # Copy the global alignment file locally

        # Get the forward and reverse versions of each alignment for symmetry with chaining
        cactus_call(parameters=['cat', i], outfile=j)
        cactus_call(parameters=['paffy', 'invert', "-i", i], outfile=j, outappend=True)  # Not bothering to log this one

        # Now chain the alignments
        messages = cactus_call(parameters=['paffy', 'chain', "-i", j,
                                           "--maxGapLength", params.find("blast").attrib["chainMaxGapLength"],
                                           "--chainGapOpen", params.find("blast").attrib["chainGapOpen"],
                                           "--chainGapExtend", params.find("blast").attrib["chainGapExtend"],
                                           "--trimFraction", params.find("blast").attrib["chainTrimFraction"],
                                           "--logLevel", getLogLevelString()],
                               outfile=k, outappend=True, returnStdErr=True)
        logger.info("paffy chain {}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy chain

    # Now tile to select the primary alignments
    messages = cactus_call(parameters=['paffy', 'tile', "-i", k, "--logLevel", getLogLevelString()],
                           outfile=j, returnStdErr=True)
    logger.info("paffy tile event:{}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy tile

    # Trim the poorly aligned tails off the ends of the alignments after tiling - doing so before tiling
    # can create gaps in alignment chains which allows in spurious chains
    messages = cactus_call(parameters=['paffy', 'trim', "-i", j,
                                       "--trimIdentity", params.find("blast").attrib["pafTrimIdentity"]],
                           outfile=k, returnStdErr=True)
    logger.info("paffy trim {}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy trim

    # Filter to primary alignments
    messages = cactus_call(parameters=['paffy', 'filter', "-i", k, "--maxTileLevel", "1"],
                           outfile=j, returnStdErr=True)
    logger.info("paffy filter to get primaries {}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy filter

    # Do we want to include the secondary alignments
    use_secondary_alignments = int(params.find("blast").attrib["outputSecondaryAlignments"])  # We should really switch to
    # the empty string being false instead of 0

    if use_secondary_alignments:
        # Filter to secondary alignments and put in the final output file
        messages = cactus_call(parameters=['paffy', 'filter', "-i", k, "--maxTileLevel", "1", '-x'],
                               outfile=output_alignments_file, returnStdErr=True)
        logger.info("paffy filter to get secondaries {}\n{}".format(reference_event_name, messages[:-1]))

    # Rechain the "primary" alignments so we can see how good the chains of the primary alignments are
    messages = cactus_call(parameters=['paffy', 'chain', "-i", j,
                                       "--maxGapLength", params.find("blast").attrib["chainMaxGapLength"],
                                       "--chainGapOpen", params.find("blast").attrib["chainGapOpen"],
                                       "--chainGapExtend", params.find("blast").attrib["chainGapExtend"],
                                       "--trimFraction", params.find("blast").attrib["chainTrimFraction"],
                                       "--logLevel", getLogLevelString()],
                           outfile=k, returnStdErr=True)
    logger.info("paffy chain of primary alignments {}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy chain

    # Filter primary alignments not in good chains
    messages = cactus_call(parameters=['paffy', 'filter', "-i", k,
                                       "--minChainScore", params.find("blast").attrib["minPrimaryChainScore"]],
                           outfile=output_alignments_file, outappend=True, returnStdErr=True)
    logger.info("paffy filter {}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy filter

    if use_secondary_alignments:
        # Get the primaries we've filtered and switch them to secondaries in the final output
        messages = cactus_call(parameters=['paffy', 'filter', "-i", k, "-x",
                                           "--minChainScore", params.find("blast").attrib["minPrimaryChainScore"]],
                               outfile=j, returnStdErr=True)
        logger.info("paffy filter {}\n{}".format(reference_event_name, messages[:-1]))  # Log paffy filter
        messages = cactus_call(parameters=['sed', 's/tp:A:P/tp:A:S/', k], outfile=j, returnStdErr=True)
        logger.info("switching low score primaries to secondaries {}\n{}".format(reference_event_name, messages[:-1]))
        messages = cactus_call(parameters=['sed', 's/tl:i:1/tl:i:2/', j], outfile=output_alignments_file,
                               outappend=True, returnStdErr=True)
        logger.info("switching changing tile level of low score primaries {}\n{}".format(reference_event_name, messages[:-1]))

    # Cleanup the old alignment files
    for i in alignment_files:
        job.fileStore.deleteGlobalFile(i)

    return job.fileStore.writeGlobalFile(output_alignments_file)  # Copy back

def sanitize_then_make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    sanitize_job = job.addChildJobFn(sanitize_fasta_headers, event_names_to_sequences)
    paf_job = sanitize_job.addFollowOnJobFn(make_paf_alignments, event_tree_string, sanitize_job.rv(), ancestor_event_string, params)
    return paf_job.rv()

def make_paf_alignments(job, event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    # a job should never set its own follow-on, so we hang everything off the root_job here to encapsulate
    root_job = Job()
    job.addChild(root_job)

    logger.info("Parsing species tree: {}".format(event_tree_string))
    event_tree = newickTreeParser(event_tree_string)

    ancestor_event = get_node(event_tree, ancestor_event_string)
    ingroup_events = get_leaves(ancestor_event) # Get the set of ingroup events
    # to be consistent with pre-refactor (and work with updating tests), we include the root when its id is input
    # and just treat it as an ingroup
    if ancestor_event.iD in event_names_to_sequences and event_names_to_sequences[ancestor_event.iD]:
        ingroup_events.append(ancestor_event)
    outgroup_events = [event for event in get_leaves(event_tree) if event not in ingroup_events]  # Set of outgroups
    logger.info("Got ingroup events: {} for ancestor event: {}".format(" ".join([i.iD for i in ingroup_events]),
                                                                       ancestor_event_string))

    # Calculate the total sequence size
    total_sequence_size = sum(event_names_to_sequences[event.iD].size for event in get_leaves(event_tree))

    # for each pair of ingroups make alignments
    ingroup_alignments = []
    for ingroup, ingroup2, distance_a_b in get_event_pairs(ancestor_event, ingroup_events):
        logger.info("Building alignment between event: {} (ingroup) and event: {} (ingroup)".format(ingroup.iD, ingroup2.iD))
        ingroup_alignments.append(root_job.addChildJobFn(make_chunked_alignments,
                                                         ingroup.iD, event_names_to_sequences[ingroup.iD],
                                                         ingroup2.iD, event_names_to_sequences[ingroup2.iD], distance_a_b, params,
                                                         disk=2*total_sequence_size).rv())

    # Get the outgroup events
    distances = get_distances(event_tree)  # Distances between all pairs of nodes
    outgroup_events.sort(key=lambda outgroup: distances[ancestor_event, outgroup])  # Sort from closest to furthest
    logger.info("Got outgroup events: {} for ancestor event: {}".format(" ".join([i.iD for i in outgroup_events]),
                                                                        ancestor_event.iD))

    # for each ingroup make alignments to the outgroups
    if int(params.find("blast").attrib["trimIngroups"]):  # Trim the ingroup sequences
        outgroup_alignments = [root_job.addChildJobFn(make_ingroup_to_outgroup_alignments_1, ingroup, outgroup_events,
                                                      dict(event_names_to_sequences), distances, params).rv()
                                for ingroup in ingroup_events] if len(outgroup_events) > 0 else []
    else:
        outgroup_alignments = [root_job.addChildJobFn(make_chunked_alignments,
                                                      ingroup.iD, event_names_to_sequences[ingroup.iD],
                                                      outgroup.iD, event_names_to_sequences[outgroup.iD],
                                                      distances[ingroup, outgroup], params,
                                                      disk=2*total_sequence_size).rv()
                               for ingroup in ingroup_events for outgroup in outgroup_events]

    # Now do the chaining
    return root_job.addFollowOnJobFn(chain_alignments, ingroup_alignments + outgroup_alignments,
                                     ancestor_event_string, params, disk=2*total_sequence_size).rv()


def trim_unaligned_sequences(job, sequences, alignments, params):
    alignments = job.fileStore.readGlobalFile(alignments)  # Download the alignments

    # Make a bed of aligned sequence
    # run paffy to_bed to create a bed of aligned coverage
    bed_file = job.fileStore.getLocalTempFile()  # Get a temporary file to store the bed file in
    messages = cactus_call(parameters=['paffy', 'to_bed', "--binary", "--excludeUnaligned", "--includeInverted",
                                       '-i', alignments, "--logLevel", getLogLevelString()],
                           outfile=bed_file, returnStdErr=True)
    # Log paffy to_bed
    logger.info("paffy to_bed:\n{}".format(messages[:-1]))

    trimmed_sequence_files = []
    for sequence in sequences:
        # run faffy extract to extract unaligned sequences longer than a threshold creating a reduced subset of A
        seq_file = job.fileStore.readGlobalFile(sequence)  # Load original sequence file
        trimmed_seq_file = job.fileStore.getLocalTempFile()  # Get a temporary file to store the extracted subsequences
        messages = cactus_call(parameters=['faffy', 'extract', "-i", bed_file, seq_file, "--skipMissing", "--minSize", "1",
                                           "--flank", params.find("blast").attrib["trimOutgroupFlanking"],
                                           "--logLevel", getLogLevelString()],
                               outfile=trimmed_seq_file, returnStdErr=True)
        logger.info("faffy extract \n{}".format(messages[:-1]))  # Log faffy extract
        trimmed_sequence_files.append(trimmed_seq_file)

    # Now convert the alignments to refer to the reduced sequences
    trimmed_alignments = job.fileStore.getLocalTempFile()  # Get a temporary file to store the "trimmed" alignments
    messages = cactus_call(parameters=['paffy', 'upconvert', "-i", alignments, "--logLevel", getLogLevelString()] +
                           trimmed_sequence_files, outfile=trimmed_alignments, returnStdErr=True)
    logger.info("paffy upconvert\n{}".format(messages[:-1]))  # Log

    return [job.fileStore.writeGlobalFile(i) for i in trimmed_sequence_files], \
        job.fileStore.writeGlobalFile(trimmed_alignments)  # Return the trimmed sequence files and trimmed alignments

