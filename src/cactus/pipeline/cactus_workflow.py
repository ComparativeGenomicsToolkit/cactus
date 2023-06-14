#!/usr/bin/env python3

#Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Script strings together all the components to make the basic pipeline for reconstruction.
"""

import os
import sys
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString
from toil.realtimeLogger import RealtimeLogger
from toil.lib.conversions import bytes2human
from sonLib.nxnewick import NXNewick
from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import findRequiredNode, getOptionalAttrib

############################################################
############################################################
############################################################
##The consolidate phase, which runs the setup, caf,
## bar, reference and cactus to hal algorithms in one job
## on a multi-node machine
############################################################
############################################################
############################################################

def cactus_cons_with_resources(job, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                               cons_cores = None, cons_memory = None, intermediate_results_url = None, chrom_name = None):
    ''' run cactus_consolidated as a child job, requesting resources based on input sizes '''

    # compute resource requirements
    total_sequence_size = sum([seq_id.size for seq_id in seq_id_map.values()])
    disk = 5 * total_sequence_size + 2 * paf_id.size

    # constant factor
    mem = 1e8
    # abPOA needs a bunch of memory for its table, even for tiny alignments
    if getOptionalAttrib(findRequiredNode(config_node, 'bar'), 'partialOrderAlignment', typeFn=bool, default=True):
        mem = 4e9
    # add function of paf size
    mem += 3 * paf_id.size
    # and quadratic in sequence size (if doable in kb)
    if total_sequence_size > 1024:
        total_sequence_size_kb = total_sequence_size / 1024.
        if len(seq_id_map) < 6:
            mem += 1024. * (total_sequence_size_kb ** 1.75)
        else:
            # probably in pangenome mode: go a bit easier
            mem += 1024. * (total_sequence_size_kb ** 1.30)
    else:
        mem += 4 * total_sequence_size        
    # but we cap things at 512 Gigs
    mem=int(min(mem, 512e09))

    RealtimeLogger.info('Estimating cactus_consolidated({}) memory as {} bytes from {} sequences with total-sequence-size {} and paf-size {}'.format(ancestor_event, mem, len(seq_id_map), total_sequence_size, paf_id.size))

    if cons_memory is not None and cons_memory < mem:
        RealtimeLogger.info('Overriding cactus_conslidated memory estimate of {} with {} from --consMemory'.format(
            bytes2human(mem), bytes2human(cons_memory)))
        mem = cons_memory

    cons_job = job.addChildJobFn(cactus_cons, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                                 intermediate_results_url=intermediate_results_url, chrom_name=chrom_name, cores = cons_cores,
                                 memory=mem, disk=disk)
    return cons_job.rv()

def cactus_cons(job, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                intermediate_results_url = None, chrom_name = None):
    ''' run cactus_consolidated '''

    # Build up a genome -> fasta map.
    work_dir = job.fileStore.getLocalTempDir()
    seq_path_map = {}
    for event, seq_id in seq_id_map.items():
        seq_path = os.path.join(work_dir, '{}.fa'.format(event))
        job.fileStore.readGlobalFile(seq_id, seq_path)
        seq_path_map[event] = seq_path

    outgroups = og_map[ancestor_event] if ancestor_event in og_map else []

    # Get the alignments file
    paf_path = os.path.join(work_dir, f'{ancestor_event}.paf')
    job.fileStore.readGlobalFile(paf_id, paf_path)

    # Split the alignments file into primary and secondary
    primary_alignment_file = os.path.join(work_dir, f'{ancestor_event}_primary.paf')
    system(f"grep -v 'tp:A:S' {paf_path} > {primary_alignment_file} || true")  # Alignments that are not-secondaries

    # Optionally parse our secondary alignments
    use_secondary_alignments = int(config_node.find("blast").attrib["outputSecondaryAlignments"])  # We should really switch to
    # the empty string being false instead of 0
    assert use_secondary_alignments == 0 or use_secondary_alignments == 1
    if use_secondary_alignments:
        secondary_alignment_file = os.path.join(work_dir, f'{ancestor_event}_secondary.paf')
        system(f"grep 'tp:A:S' {paf_path} > {secondary_alignment_file} || true")  # Alignments that are secondaries

    # Optionally copy the alignments to a specified location for debug purposes
    if config_node.find("caf").attrib["writeInputAlignmentsTo"]:
        system(f'cp {paf_path} {config_node.find("caf").attrib["writeInputAlignmentsTo"]}/{ancestor_event}.paf')

    # Temporary place to store the output c2h file
    tmpHal = os.path.join(work_dir, f'{ancestor_event}.c2h')
    tmpFasta = os.path.join(work_dir, f'{ancestor_event}.c2h.fa')
    tmpRef = os.path.join(work_dir, f'{ancestor_event}.ref')

    tmpConfig = os.path.join(work_dir, f'{ancestor_event}.config.xml')
    ConfigWrapper(config_node).writeXML(tmpConfig)

    # We pass in the genome->sequence map as a series of paired arguments: [genome, faPath]*N.
    pairs = [[genome, faPath] for genome, faPath in list(seq_path_map.items())]
    args = ["--sequences", " ".join([item for sublist in pairs for item in sublist]),
            "--speciesTree", NXNewick().writeString(tree), "--logLevel", getLogLevelString(),
            "--alignments", primary_alignment_file, "--params", tmpConfig, "--outputFile", tmpHal,
            "--outputHalFastaFile", tmpFasta, "--outputReferenceFile", tmpRef, "--outgroupEvents", " ".join(outgroups),
            "--referenceEvent", ancestor_event, "--threads", str(job.cores)]
    if use_secondary_alignments:  # Optionally add the secondary alignments
        args += ["--secondaryAlignments", secondary_alignment_file]

    messages = cactus_call(check_output=True, returnStdErr=True,
                           realtimeStderrPrefix=f'cactus_consolidated({chrom_name if chrom_name else ancestor_event})',
                           parameters=["cactus_consolidated"] + args,
                           job_memory=job.memory)[1]  # Get just the standard error output

    # if cactus was run with --realTimeLogging, cactus_call will print out conslidated's stderr messages as they happen
    # (and not return anything)
    # otherwise, if run without --realTimeLogging, cactus_call will return (but not print) stderr messages
    if messages:
        job.fileStore.logToMaster(f"cactus_consolidated event:{ancestor_event}\n{messages}")  # Log the messages
    else:
        job.fileStore.logToMaster("Ran cactus consolidated okay")

    # Write the temporary output files to the final output
    # At top level--have the final .c2h file
    halID = job.fileStore.writeGlobalFile(tmpHal)
    fastaID = job.fileStore.writeGlobalFile(tmpFasta)
    referenceID = job.fileStore.writeGlobalFile(tmpRef)

    if intermediate_results_url is not None:
        # The user requested to keep the c2h files in a separate place. Export it there.
        url = intermediate_results_url + ".c2h"
        job.fileStore.exportFile(halID, makeURL(url))

        # The user requested to keep the hal fasta files in a separate place. Export it there.
        url = intermediate_results_url + ".hal.fa"
        job.fileStore.exportFile(fastaID, makeURL(url))

        # The user requested to keep the reference fasta files in a separate place. Export it there.
        url = intermediate_results_url + ".reference.fa"
        job.fileStore.exportFile(referenceID, makeURL(url))

    return (ancestor_event, halID, fastaID, referenceID)


if __name__ == '__main__':
    runCactusWorkflow(sys.argv)
