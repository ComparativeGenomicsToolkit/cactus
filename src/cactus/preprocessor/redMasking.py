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
from cactus.preprocessor.checkPreprocessedSequence import check_sequence_preserved
from cactus.preprocessor.maskingCommon import prefilter_cmd, masked_base_count
from cactus.preprocessor.maskingCommon import extract_masking_bed, soft_mask_intervals
from cactus.preprocessor.maskingCommon import log_masking_delta

from toil.realtimeLogger import RealtimeLogger


class RedMaskJob(RoundedJob):
    def __init__(self, fastaID, redOpts, redPrefilterOpts, eventName=None, unmask=False):
        memory = cactus_clamp_memory(12*fastaID.size)
        disk = 5*(fastaID.size)
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True)
        self.fastaID = fastaID
        self.redOpts = redOpts
        self.redPrefilterOpts = redPrefilterOpts
        self.eventName = eventName if eventName else 'seq'
        self.unmask = unmask

    def run(self, fileStore):
        """
        mask repeats with RED.  RED ignores existing masking, so the intervals it reports
        are merged with the ones the input already carried.

        Only the intervals are taken from RED's output; the sequence itself always comes
        from our own copy of the input.  RED writes its output with no error checking of
        any kind, so a full disk gives a silently truncated file and a zero exit status,
        and applying one file's coordinates to another file's sequence would then drop
        masking off the end without a word.  Working from the input instead means a bad
        RED can cost us masking but never sequence.  It is also checked for outright, so
        it should not get that far.
        """
        # download fasta
        work_dir = fileStore.getLocalTempDir()
        red_in_dir = os.path.join(work_dir, 'red-in-{}'.format(self.eventName))
        red_out_dir = os.path.join(work_dir, 'red-out-{}'.format(self.eventName))
        os.makedirs(red_in_dir)
        os.makedirs(red_out_dir)
        raw_fa_path = os.path.join(work_dir, '{}.fa'.format(self.eventName))
        in_fa_path = os.path.join(red_in_dir, '{}.filter.fa'.format(self.eventName))
        red_msk_path = os.path.join(red_out_dir, '{}.filter.msk'.format(self.eventName))
        out_fa_path = os.path.join(work_dir, '{}.mask.fa'.format(self.eventName))
        fileStore.readGlobalFile(self.fastaID, raw_fa_path)

        # get rid of small or single-base contigs that might crash Red
        filter_cmd = prefilter_cmd(raw_fa_path, self.redPrefilterOpts)
        cactus_call(parameters=filter_cmd, outfile=in_fa_path)

        if os.path.getsize(in_fa_path) > 0:
            # measure the masking the input already carried, for the log line below
            pre_mask_size = masked_base_count(in_fa_path)

            # run red
            red_cmd = ['Red', '-gnm', red_in_dir, '-msk', red_out_dir]
            if self.redOpts:
                red_cmd += self.redOpts.split()
            cactus_call(parameters=red_cmd)

            # RED has been seen returning less sequence than it was given, so make sure
            # its output still describes the same bases before believing its intervals
            check_sequence_preserved(in_fa_path, red_msk_path, event_name=self.eventName,
                                     step_name='Red')

            # take the intervals RED masked.  this picks up the N runs as well, which is
            # harmless: they end up soft-masked rather than left as upper case Ns
            red_bed = os.path.join(work_dir, '{}.red.masking.bed'.format(self.eventName))
            extract_masking_bed(red_msk_path, red_bed)

            # and apply them to the input, not to RED's output.  see maskingCommon
            soft_mask_intervals(in_fa_path, red_bed, out_fa_path, unmask=self.unmask)

            log_masking_delta('Red', self.eventName, pre_mask_size,
                              masked_base_count(out_fa_path))
        else:
            RealtimeLogger.info('Skipping Red for {} because contigs are too small'.format(self.eventName))

        # put the filtered contigs back
        cactus_call(parameters=filter_cmd + ['-x'], outfile=out_fa_path, outappend=True)

        return fileStore.writeGlobalFile(out_fa_path)
