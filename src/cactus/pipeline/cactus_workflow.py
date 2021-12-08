#!/usr/bin/env python3

#Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Script strings together all the components to make the basic pipeline for reconstruction.
"""

import os
import sys
import xml.etree.ElementTree as ET
import math
from argparse import ArgumentParser
from toil.lib.bioio import system
from toil.lib.bioio import getTempFile
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.lib.bioio import getLogLevelString
from sonLib.bioio import catFiles
from sonLib.nxnewick import NXNewick

from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import runCactusConsolidated
from cactus.shared.common import findRequiredNode

from cactus.shared.configWrapper import ConfigWrapper


############################################################
############################################################
############################################################
##The consolidate phase, which runs the setup, caf,
## bar, reference and cactus to hal algorithms in one job
## on a multi-node machine
############################################################
############################################################
############################################################

def cactus_cons_with_resources(job, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id, cons_cores = None, intermediate_results_url = None):
    ''' run cactus_consolidated as a child job, requesting resources based on input sizes '''

    # compute resource requirements
    total_sequence_size = 0
    for seq_id in seq_id_map.values():
        seq_path = job.fileStore.readGlobalFile(seq_id, mutable=True, cache=False)
        total_sequence_size += os.stat(seq_path).st_size
        os.remove(seq_path)

    disk = int(3 * total_sequence_size)
 
    # this is the old caf job's memory function
    memoryPoly = [1.80395944e+01, 7.96042247e+07]
    memoryCap = 120e09
    mem = 0
    for degree, coefficient in enumerate(reversed(memoryPoly)):
        mem += coefficient * (total_sequence_size**degree)
    mem = int(min(mem, memoryCap))

    cons_job = job.addChildJobFn(cactus_cons, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                                 intermediate_results_url=intermediate_results_url, cores = cons_cores,
                                 memory=mem, disk=disk)
    return cons_job.rv()
   
def cactus_cons(job, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id, intermediate_results_url = None):
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
    paf_path = os.path.join(work_dir, '{}.paf'.format(ancestor_event))
    job.fileStore.readGlobalFile(paf_id, paf_path)
    
    # Split the alignments file into primary and secondary
    primary_alignment_file = os.path.join(work_dir, '{}_primary.paf'.format(ancestor_event))
    system("grep -v 'tp:A:S' {} > {} || true".format(paf_path, primary_alignment_file))
    secondary_alignment_file = os.path.join(work_dir, '{}_secondary.paf'.format(ancestor_event))
    system("grep 'tp:A:S' {} > {} || true".format(paf_path, secondary_alignment_file))

    # Temporary place to store the output c2h file
    tmpHal = os.path.join(work_dir, '{}.c2h'.format(ancestor_event))
    tmpFasta = os.path.join(work_dir, '{}.c2h.fa'.format(ancestor_event))
    tmpRef = os.path.join(work_dir, '{}.ref'.format(ancestor_event))

    tmpConfig = os.path.join(work_dir, '{}.config.xml'.format(ancestor_event))
    ConfigWrapper(config_node).writeXML(tmpConfig)

    messages = runCactusConsolidated(seqMap=seq_path_map,
                                     newickTreeString=NXNewick().writeString(tree),
                                     cactusParams=tmpConfig,
                                     alignmentsFile=primary_alignment_file,
                                     outputFile=tmpHal,
                                     outputHalFastaFile=tmpFasta,
                                     outputReferenceFile=tmpRef,
                                     secondaryAlignmentsFile=secondary_alignment_file,
                                     constraintAlignmentsFile=None,
                                     logLevel=getLogLevelString(),
                                     outgroupEvents=outgroups,
                                     referenceEvent=ancestor_event,
                                     cores=job.cores)

    # Log back any messages
    for message in messages:
        fileStore.logToMaster(message)

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
