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
    min_length = getOptionalAttrib(lastz_params_node.find("unmask"), 'min_length', typeFn=int, default=0)
    assert unmask_action in ['unmask', 'remask']

    work_dir = job.fileStore.getLocalTempDir()
    fasta_path = os.path.join(work_dir, '{}.fa'.format(event))
    job.fileStore.readGlobalFile(fasta_id, fasta_path)

    unmasked_contigs = {} # map name to length
    for seq_record in SeqIO.parse(fasta_path, 'fasta'):
        seq = seq_record.seq
        seq_str = str(seq)  # Convert to string once for faster iteration
        if len(seq_str) < min_length:
            continue
        
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
            unmasked_contigs[seq_record.id] = len(seq_str)

    if unmasked_contigs:
        bed_path = os.path.join(work_dir, '{}.unmasked.bed'.format(event))
        with open(bed_path, 'w') as bed_file:
            for contig, length in unmasked_contigs.items():
                bed_file.write('{}\t0\t{}\n'.format(contig, length))
        
        if unmask_action == 'unmask':
            unmasked_fasta_path = os.path.join(work_dir, '{}.unmask.fa'.format(event))
            cactus_call(parameters=['cactus_fasta_softmask_intervals.py', '--origin=zero', '--unmask', bed_path],
                    infile=fasta_path, outfile=unmasked_fasta_path)
            job.fileStore.deleteGlobalFile(fasta_id)
            fasta_id = job.fileStore.writeGlobalFile(unmasked_fasta_path)
            RealtimeLogger.info('Unmasked {} bp of contigs for {} using thresholds {}, {}'.format(
                sum(x for x in unmasked_contigs.values()), event, unmask_threshold, unmask_threshold_bp))
        else:
            active = {}
            for node in params.findall("preprocessor"):
                # deactivate everything but masking
                if getOptionalAttrib(node, "preprocessJob") not in  ['lastzRepeatMask', 'red', 'fastan']:
                    node.attrib['active'] = '0'
                # figure out which masking needs to unmask first
                elif getOptionalAttrib(node, 'active', typeFn=bool, default=False):
                    active[getOptionalAttrib(node, "preprocessJob")] = node
                    node.attrib['unmask'] = '0'
            if active:
                # note: order here must be same as in preprocessor, since we only want to call unmask once
                for prep in ['lastzRepeatMask', 'red', 'fastan']:
                    if prep in active:
                        active[prep].attrib['unmask'] = '1'
                        break
                # preprocessor is going to delete the input file. we don't want to do that so send in a copy
                copy_fasta_id = job.fileStore.writeGlobalFile(fasta_path)                    
                # run the preprocessor on the whole thing
                pp_job = job.addChild(CactusPreprocessor([copy_fasta_id], params, eventNames=[event]))
                pp_id = pp_job.rv(0)
                # mix in the masking for the to-mask contigs with our original file
                fa_merge_job = pp_job.addFollowOnJobFn(merge_fa, event, fasta_id, pp_id, set(unmasked_contigs.keys()),
                                                       disk=fasta_id.size * 4)
                fasta_id = fa_merge_job.rv()
                # delete the pp_id
                fa_merge_job.addFollowOnJobFn(clean_jobstore_files, file_ids=[pp_id])

            else:
                RealtimeLogger.warning('Remasking activated but no maskers active in preprocessor: doing nothing')

    return fasta_id

def merge_fa(job, event, fasta1_id, fasta2_id, contigs):
    """ splice two fastas together, taking contigs list from second one """
    work_dir = job.fileStore.getLocalTempDir()
    
    fa_1 = os.path.join(work_dir, '{}.fa'.format(event))
    fa_2 = os.path.join(work_dir, '{}.fa.pp'.format(event))
    job.fileStore.readGlobalFile(fasta1_id, fa_1)
    job.fileStore.readGlobalFile(fasta2_id, fa_2)
    fa_merge = os.path.join(work_dir, '{}-remasked.fa'.format(event))
    with open(fa_merge, 'w') as out_file:
        for seq_record in SeqIO.parse(fa_1, 'fasta'):
            if seq_record.id not in contigs:
                SeqIO.write(seq_record, out_file, 'fasta')
        for seq_record in SeqIO.parse(fa_2, 'fasta'):
            if seq_record.id in contigs:
                SeqIO.write(seq_record, out_file, 'fasta')

    analysisString1 = cactus_call(parameters=["cactus_analyseAssembly", os.path.basename(fa_1)],
                                 work_dir=work_dir, check_output=True)
    job.fileStore.logToMaster("Assembly stats before remasking for %s: %s" % (event, analysisString1))

    analysisString2 = cactus_call(parameters=["cactus_analyseAssembly", os.path.basename(fa_merge)],
                                 work_dir=work_dir, check_output=True)
    job.fileStore.logToMaster("Assembly stats after remasking for %s: %s" % (event, analysisString2))
                
    return job.fileStore.writeGlobalFile(fa_merge)


    
    
    
    
