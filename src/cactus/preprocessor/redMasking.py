#!/usr/bin/env python3
"""Uses RED to mask repeats
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


class RedMaskJob(RoundedJob):
    def __init__(self, fastaID, redOpts, redPrefilterOpts, eventName=None, unmask=False):
        memory = cactus_clamp_memory(6*fastaID.size)
        disk = 5*(fastaID.size)
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.fastaID = fastaID
        self.redOpts = redOpts
        self.redPrefilterOpts = redPrefilterOpts
        self.eventName = eventName if eventName else 'seq'
        self.unmask = unmask

    def run(self, fileStore):
        """
        mask repeats with RED.  RED ignores existing masking.  So if unmask is false, the old 
        masking will be explicitly merged back in. 
        """
        # download fasta
        work_dir = fileStore.getLocalTempDir()
        red_in_dir = os.path.join(work_dir, 'red-in-{}'.format(self.eventName))
        red_out_dir = os.path.join(work_dir, 'red-out-{}'.format(self.eventName))
        os.makedirs(red_in_dir)
        os.makedirs(red_out_dir)
        raw_fa_path = os.path.join(work_dir, '{}.fa'.format(self.eventName))
        in_fa_path = os.path.join(red_in_dir, '{}.filter.fa'.format(self.eventName))
        out_fa_path = os.path.join(red_out_dir, '{}.filter.msk'.format(self.eventName))
        fileStore.readGlobalFile(self.fastaID, raw_fa_path)

        # get rid of small or single-base contigs that might crash Red
        filter_cmd = ['cactus_redPrefilter', raw_fa_path]
        if self.redPrefilterOpts:
            assert '-x' not in self.redPrefilterOpts and '--extract' not in self.redPrefilterOpts
            filter_cmd += self.redPrefilterOpts.split()
        cactus_call(parameters=filter_cmd, outfile=in_fa_path)

        if os.path.getsize(in_fa_path) > 0:
            # preserve existing masking
            pre_mask_size = 0
            if not self.unmask:
                bed_path = os.path.join(work_dir, '{}.input.masking.bed'.format(self.eventName))
                cactus_call(parameters=['cactus_softmask2hardmask', '-b', in_fa_path], outfile=bed_path)
                awkres = cactus_call(parameters=['awk', '{sum += $3-$2} END {print sum}', bed_path],
                                                check_output=True, rt_log_cmd=False).strip()
                try:
                    pre_mask_size = int(float(awkres)) if awkres else 0
                except ValueError as e:
                    print(f"Error converting awkres to int: {e}")
                    pre_mask_size = 0
                
            # run red
            red_cmd = ['Red', '-gnm', red_in_dir, '-msk', red_out_dir]
            if self.redOpts:
                red_cmd += self.redOpts.split()
            cactus_call(parameters=red_cmd)

            # merge the exsiting masking back in
            if not self.unmask:
                if pre_mask_size:
                    cactus_call(infile=out_fa_path, outfile=out_fa_path + '.remask',
                                parameters=['cactus_fasta_softmask_intervals.py', '--origin=zero', bed_path])
                    out_fa_path = out_fa_path + '.remask'

                awkres = cactus_call(parameters=[['cactus_softmask2hardmask', '-b', out_fa_path],
                                                 ['awk', '{sum += $3-$2} END {print sum}']],
                                     check_output=True, rt_log_cmd=False).strip()
                post_mask_size = int(awkres) if awkres else 0
                RealtimeLogger.info('Red masked {} bp of {}, increasing masking from {} to {}'.format(
                    post_mask_size - pre_mask_size, self.eventName, pre_mask_size, post_mask_size))
        else:
            RealtimeLogger.info('Skipping Red for {} because contigs are too small'.format(self.eventName))

        # put the filtered contigs back
        cactus_call(parameters=filter_cmd + ['-x'], outfile=out_fa_path, outappend=True)

        return fileStore.writeGlobalFile(out_fa_path)
