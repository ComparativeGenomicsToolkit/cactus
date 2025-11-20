#!/usr/bin/env python3
"""Checks headers are all unique.
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from toil.job import Job
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import cactus_call, getOptionalAttrib
from toil.realtimeLogger import RealtimeLogger

def unmask_contigs_all(job, event_names_to_sequences, ingroup_events, params):
    """ unmask the given events and return and updated fasta ids for each ingroup as list """
    root_job = Job()
    job.addChild(root_job)
    output_ids = []
    for event in ingroup_events:
        unmask_job = root_job.addChildJobFn(unmask_contigs_one, event, event_names_to_sequences[event], params,
                                            memory=event_names_to_sequences[event].size * 2,
                                            disk=event_names_to_sequences[event].size * 5)
        output_ids.append(unmask_job.rv())
       
    return output_ids

def unmask_contigs_one(job, event, fasta_id, params):
    """ unmask a single fasta """
    lastz_params_node = params.find("blast")
    unmask_params_node = lastz_params_node.find("unmask")
    unmask_threshold = getOptionalAttrib(unmask_params_node, 'threshold', typeFn=float, default=1.0)
    unmask_threshold_bp = getOptionalAttrib(unmask_params_node, 'threshold_bp', typeFn=int, default=0)
    unmask_action = getOptionalAttrib(lastz_params_node.find("unmask"), 'action', typeFn=str, default='none')
    #assert unmask_action in ['unmask', 'remask']
    # todo: remask not yet implemented
    assert unmask_action in ['unmask']

    work_dir = job.fileStore.getLocalTempDir()
    fasta_path = os.path.join(work_dir, '{}.fa'.format(event))
    job.fileStore.readGlobalFile(fasta_id, fasta_path)
    bed_path = os.path.join(work_dir, '{}.unmasked.bed'.format(event))

    unmasked_bp = 0
    with open(bed_path, 'w') as bed_file:
        for seq_record in SeqIO.parse(fasta_path, 'fasta'):
            header = seq_record.description
            seq = seq_record.seq
            seq_str = str(seq)  # Convert to string once for faster iteration
            threshold = max(unmask_threshold_bp, int(unmask_threshold * len(seq)))
            unmasked_bases = 0
            too_masked = True
            for c in seq_str:
                if c.isupper():
                    unmasked_bases += 1
                    if unmasked_bases >= threshold:
                        too_masked = False
                        break
            if too_masked:
                bed_file.write('{}\t0\t{}\n'.format(header, len(seq)))
                unmasked_bp += len(seq)

    if unmasked_bp > 0:
        RealtimeLogger.info('Unmasked {} bp of contigs for {} using thresholds {}, {}'.format(
            unmasked_bp, event, unmask_threshold, unmask_threshold_bp))
        unmasked_fasta_path = os.path.join(work_dir, '{}.unmask.fa'.format(event))
        cactus_call(parameters=['cactus_fasta_softmask_intervals.py', '--origin=zero', '--unmask', bed_path],
                    infile=fasta_path, outfile=unmasked_fasta_path)
        fasta_id = job.fileStore.writeGlobalFile(unmasked_fasta_path)

    return fasta_id


    
    
    
    
