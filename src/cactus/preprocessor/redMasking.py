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

from toil.realtimeLogger import RealtimeLogger


class RedMaskJob(RoundedJob):
    def __init__(self, fastaID, redOpts, eventName=None, unmask=False):
        memory = max(4*1024*1024*1024, 8*fastaID.size)
        disk = 3*(fastaID.size)
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.fastaID = fastaID
        self.redOpts = redOpts
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
        in_fa_path = os.path.join(red_in_dir, '{}.fa'.format(self.eventName))
        out_fa_path = os.path.join(red_out_dir, '{}.msk'.format(self.eventName))
        fileStore.readGlobalFile(self.fastaID, in_fa_path)

        # preserve existing masking
        pre_mask_size = 0
        if not self.unmask:
            bed_path = os.path.join(work_dir, '{}.input.masking.bed'.format(self.eventName))
            cactus_call(parameters=['cactus_softmask2hardmask', '-b', in_fa_path], outfile=bed_path)
            pre_mask_size = int(cactus_call(parameters=['awk', '{sum += $3-$2} END {print sum}', bed_path], check_output=True).strip())
            
        # run red
        red_cmd = ['Red', '-gnm', red_in_dir, '-msk', red_out_dir]
        if self.redOpts:
            red_cmd += self.redOpts.split()
        cactus_call(parameters=red_cmd)

        # merge the exsiting masking back in
        if not self.unmask:
            cactus_call(infile=out_fa_path, outfile=out_fa_path + '.remask',
                        parameters=['cactus_fasta_softmask_intervals.py', '--origin=zero', bed_path])
            out_fa_path = out_fa_path + '.remask'

        post_mask_size = int(cactus_call(parameters=[['cactus_softmask2hardmask', '-b', out_fa_path],
                                                     ['awk', '{sum += $3-$2} END {print sum}']], check_output=True).strip())
        RealtimeLogger.info('Red masked {} bp of {}, increasing masking from {} to {}'.format(
            post_mask_size - pre_mask_size, self.eventName, pre_mask_size, post_mask_size))

        return fileStore.writeGlobalFile(out_fa_path)
