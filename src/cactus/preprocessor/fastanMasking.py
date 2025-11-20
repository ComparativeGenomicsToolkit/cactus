#!/usr/bin/env python3
"""Uses FasTAN to mask repeats
"""

import os
import re
import sys
import shutil

from cactus.shared.common import cactus_cpu_count

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import cactusRootPath
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import makeURL
from cactus.shared.common import get_faidx_subpath_rename_cmd
from cactus.shared.common import cactus_clamp_memory

from toil.realtimeLogger import RealtimeLogger


class FasTANMaskJob(RoundedJob):
    def __init__(self, fastaID, fastanOpts, fastanPrefilterOpts, eventName=None, unmask=False):
        memory = cactus_clamp_memory(12*fastaID.size)
        disk = 5*(fastaID.size)
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.fastaID = fastaID
        self.fastanOpts = fastanOpts
        self.fastanPrefilterOpts = fastanPrefilterOpts
        self.eventName = eventName if eventName else 'seq'
        self.unmask = unmask

    def run(self, fileStore):
        """
        mask repeats with FasTAN.  FasTAN ignores existing masking.  So if unmask is false, the old 
        masking will be explicitly merged back in. 
        """
        # download fasta
        work_dir = fileStore.getLocalTempDir()
        raw_fa_path = os.path.join(work_dir, '{}.fa'.format(self.eventName))
        in_fa_path = os.path.join(work_dir, '{}.filter.fa'.format(self.eventName))
        out_1aln_path = os.path.join(work_dir, '{}.mask.1aln'.format(self.eventName))
        fileStore.readGlobalFile(self.fastaID, raw_fa_path)

        # get rid of small or single-base contigs that might crash FasTAN
        filter_cmd = ['cactus_redPrefilter', raw_fa_path]
        if self.fastanPrefilterOpts:
            assert '-x' not in self.fastanPrefilterOpts and '--extract' not in self.fastanPrefilterOpts
            filter_cmd += self.fastanPrefilterOpts.split()
        cactus_call(parameters=filter_cmd, outfile=in_fa_path)

        if os.path.getsize(in_fa_path) > 0:
            # compute stats for existing masking
            pre_mask_size = 0
            if not self.unmask:
                awkres = cactus_call(parameters=[['cactus_softmask2hardmask', '-b', raw_fa_path],
                                                 ['awk', '{sum += $3-$2} END {print sum}']],
                                     check_output=True, rt_log_cmd=False).strip()
                pre_mask_size = int(float(awkres)) if awkres else 0
                
            # run FasTAN
            fastan_cmd = ['FasTAN', in_fa_path, out_1aln_path]
            if self.fastanOpts:
                fastan_cmd += self.fastanOpts.split()
            cactus_call(parameters=fastan_cmd)

            # extract the BED
            out_bed = out_1aln_path + '.bed'
            cactus_call(parameters=['tanbed', out_1aln_path], outfile=out_bed)

            # merge the FasTAN masking back in
            mask_fa_path = os.path.join(work_dir, '{}.mask.fa'.format(self.eventName))
            softmask_cmd = ['cactus_fasta_softmask_intervals.py', '--origin=zero']
            if self.unmask:
                softmask_cmd += ['--unmask']
            softmask_cmd += [out_bed]            
            cactus_call(parameters=softmask_cmd, infile=raw_fa_path, outfile=mask_fa_path)
            
            awkres = cactus_call(parameters=[['cactus_softmask2hardmask', '-b', mask_fa_path],
                                             ['awk', '{sum += $3-$2} END {print sum}']],
                                 check_output=True, rt_log_cmd=False).strip()
            
            post_mask_size = int(float(awkres)) if awkres else 0
            RealtimeLogger.info('FasTAN masked {} bp of {}, increasing masking from {} to {}'.format(
                post_mask_size - pre_mask_size, self.eventName, pre_mask_size, post_mask_size))
        else:
            RealtimeLogger.info('Skipping FasTAN for {} because contigs are too small'.format(self.eventName))
            mask_fa_path = raw_fa_path

        return fileStore.writeGlobalFile(mask_fa_path)
