#!/usr/bin/env python3
"""Unmasks or remasks contigs that are too heavily masked.
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from toil.job import Job
from cactus.shared.common import cactus_clamp_memory, catFiles, clean_jobstore_files
from cactus.shared.common import cactus_call, getOptionalAttrib
from toil.realtimeLogger import RealtimeLogger
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

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
    assert unmask_action in ['unmask', 'remask']

    work_dir = job.fileStore.getLocalTempDir()
    fasta_path = os.path.join(work_dir, '{}.fa'.format(event))
    job.fileStore.readGlobalFile(fasta_id, fasta_path)
    if unmask_action == 'unmask':
        bed_path = os.path.join(work_dir, '{}.unmasked.bed'.format(event))
        bed_file = open(bed_path, 'w')
    else:
        remask_input_fasta_path = os.path.join(work_dir, '{}.to-remask.fasta'.format(event))
        remainder_fasta_path = os.path.join(work_dir, '{}.remainder.fasta'.format(event))
        remask_input_fasta_file = open(remask_input_fasta_path, 'w')
        remainder_fasta_file = open(remainder_fasta_path, 'w')

    unmasked_bp = 0
    for seq_record in SeqIO.parse(fasta_path, 'fasta'):
        header = seq_record.description
        seq = seq_record.seq
        seq_str = str(seq)  # Convert to string once for faster iteration
        threshold = max(unmask_threshold_bp, int(unmask_threshold * len(seq_str)))
        threshold = min(threshold, len(seq_str))
        unmasked_bases = 0
        too_masked = True
        for c in seq_str:
            if c.isupper():
                unmasked_bases += 1
                if unmasked_bases >= threshold:
                    too_masked = False
                    break
        
        if too_masked:
            unmasked_bp += len(seq_str)
            if unmask_action == 'unmask':
                bed_file.write('{}\t0\t{}\n'.format(seq_record.id, len(seq_str)))
            else:
                SeqIO.write(seq_record, remask_input_fasta_file, 'fasta')
        elif unmask_action != 'unmask':
            SeqIO.write(seq_record, remainder_fasta_file, 'fasta')
    if unmask_action == 'unmask':
        bed_file.close()
    else:
        remask_input_fasta_file.close()
        remainder_fasta_file.close()                

    if unmasked_bp > 0:
        if unmask_action == 'unmask':
            RealtimeLogger.info('Unmasked {} bp of contigs for {} using thresholds {}, {}'.format(
                unmasked_bp, event, unmask_threshold, unmask_threshold_bp))
            unmasked_fasta_path = os.path.join(work_dir, '{}.unmask.fa'.format(event))
            cactus_call(parameters=['cactus_fasta_softmask_intervals.py', '--origin=zero', '--unmask', bed_path],
                        infile=fasta_path, outfile=unmasked_fasta_path)
            fasta_id = job.fileStore.writeGlobalFile(unmasked_fasta_path)
        else:
            RealtimeLogger.info('Remasking {} bp of contigs for {} using thresholds {}, {}'.format(
                unmasked_bp, event, unmask_threshold, unmask_threshold_bp))
            for node in params.findall("preprocessor"):
                if getOptionalAttrib(node, "preprocessJob") in ['cutHeaders', 'checkUniqueHeaders']:
                    node.attrib['active'] = '0'
                else:
                    node.attrib['unmask'] = '1'                    
            remask_fasta_id = job.fileStore.writeGlobalFile(remask_input_fasta_path)
            remainder_fasta_id = job.fileStore.writeGlobalFile(remainder_fasta_path)
            # run the preprocessor
            pp_job = job.addChild(CactusPreprocessor([remask_fasta_id], params, eventNames=[event]))
            pp_id = pp_job.rv(0)
            fa_merge_job = pp_job.addFollowOnJobFn(merge_fa, remainder_fasta_id, pp_id)
            fa_merge_job.addFollowOnJobFn(clean_jobstore_files, file_ids=[remask_fasta_id])
            fasta_id = fa_merge_job.rv()

    return fasta_id

def merge_fa(job, fasta1_id, fasta2_id):
    """ cat some fastas together """
    work_dir = job.fileStore.getLocalTempDir()

    if fasta1_id.size == 0:
        return fasta2_id
    elif fasta2_id.size == 0:
        return fasta1_id

    fa_1 = os.path.join(work_dir, '1.fa')
    fa_2 = os.path.join(work_dir, '2.fa')
    fa_merge = os.path.join(work_dir, 'merged.fa')
    job.fileStore.readGlobalFile(fasta1_id, fa_1)
    job.fileStore.readGlobalFile(fasta2_id, fa_2)
    catFiles([fa_1, fa_2], fa_merge)
    job.fileStore.deleteGlobalFile(fasta1_id)
    job.fileStore.deleteGlobalFile(fasta2_id)

    return job.fileStore.writeGlobalFile(fa_merge)


    
    
    
    
