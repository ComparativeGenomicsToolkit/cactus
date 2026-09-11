#!/usr/bin/env python3
"""Uses RED to mask repeats
"""

import math
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
from cactus.shared.common import cactus_walltime
from cactus.preprocessor.checkPreprocessedSequence import check_sequence_preserved
from cactus.preprocessor.maskingCommon import prefilter_cmd, masked_base_count
from cactus.preprocessor.maskingCommon import extract_masking_bed, soft_mask_intervals
from cactus.preprocessor.maskingCommon import log_masking_delta

from toil.realtimeLogger import RealtimeLogger


def red_memory_estimate(fasta_size, longest_record_bytes):
    """Peak memory Red needs for a fasta of this shape, in bytes.

    Red's footprint is the k-mer table plus a fixed cost per base of the *longest
    single sequence*, not of the genome.  It reads one sequence at a time and,
    while scanning one, holds it three times over -- the raw record, the encoded
    copy it scans, and the original bases it writes back out -- plus a four-byte
    score per base, so about seven bytes per base of that one sequence.  The table
    is four bytes per k-mer with k = floor(log4(genome size)) clamped to [12, 15],
    so it tops out at 4 GiB however large the genome gets.

    longest_record_bytes is a record of the fasta including its header and
    newlines, so it already runs a little above the base count; the multiplier
    below adds the rest of the margin.  Checked against Red at 27e4480, where
    table + 7 * longest lands within 4% of the real peak at every scale and the
    multiplier leaves about 30% on top of that:

        genome                    k   longest   table+7L    estimate   measured
        chr14, 101 Mbp, 1 seq    13    101 Mb     0.99 GB     1.37 GB    1.03 GB
        chr15-19, 422 Mbp, 5     14    100 Mb     1.78 GB     2.35 GB    1.82 GB
        chr1-6, 1.24 Gbp, 6      15    248 Mb     6.06 GB     7.88 GB    6.10 GB

    The estimate used to be 12 * the fasta size, from before Red was rewritten to
    stream sequences and to stop building a Viterbi matrix it never read.  On a
    chromosome-scale assembly that overshoots badly, because the whole genome
    stands in for what is really the longest chromosome: the salamander fastas in
    the VGP runs asked for 231-331 GiB, used 57-73 GiB, and should now want tens of
    GB.  Over about 2.5 GB of input this estimate is below the old one even in the
    worst case, a genome delivered as a single sequence, where the two agree that
    the longest sequence is the whole thing.
    """
    # Red picks k from the non-N genome size; using the file size can only round it
    # up, which errs towards a bigger table than Red will really allocate.
    k = min(15, max(12, int(math.log(max(fasta_size, 4), 4))))
    table_bytes = 4 * (4 ** k)
    return int(1.25 * (table_bytes + 8 * longest_record_bytes))


# How much faster Red is than it was when the VGP 577-way logs were made.  Named because
# FasTAN borrows Red's rate and has to multiply this back out -- it got no such speedup.
RED_SPEEDUP = 3.0

# Seconds of Red per GB of input fasta.  Across the 625 Red runs of the VGP 577-way (0.13 to
# 10 Gb of genome) the p99 was 2799 s/Gb and the worst 5147 s/Gb.  Everything else this job
# runs -- the prefilter, the softmask/hardmask conversions, extracting and applying the
# intervals -- came to well under 100 s each even on the largest genome, and is covered by
# cactus_walltime()'s safety factor.
RED_SECS_PER_GB = 2799 / RED_SPEEDUP


class RedMaskJob(RoundedJob):
    def __init__(self, fastaID, redOpts, redPrefilterOpts, eventName=None, unmask=False,
                 longestRecordSize=None):
        # Without a measurement, fall back to assuming the whole input is one
        # sequence, which is the worst case for Red's memory.
        if longestRecordSize is None:
            longestRecordSize = fastaID.size
        memory = cactus_clamp_memory(red_memory_estimate(fastaID.size, longestRecordSize))
        disk = 5*(fastaID.size)
        RoundedJob.__init__(self, memory=memory, disk=disk, preemptable=True,
                            walltime=cactus_walltime(RED_SECS_PER_GB * fastaID.size / 1e9,
                                                     io_bytes=2 * fastaID.size))
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
