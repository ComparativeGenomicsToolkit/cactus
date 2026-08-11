#!/usr/bin/env python3
"""Helpers shared by the tool-driven repeat maskers (Red, FasTAN).

Both jobs have the same shape: drop the contigs the tool cannot cope with, run the
tool, take the intervals it reports, and soft-mask a fasta at those intervals.

The one rule worth stating out loud is that the intervals are always applied to our
own copy of the input, never to whatever fasta the tool wrote out.  Applying one
file's coordinates to another file's sequence assumes the tool returned the bases it
was given, and cactus_fasta_softmask_intervals.py used to clip an interval that ran
off the end rather than complain.  A tool that half-wrote its output therefore cost
sequence silently.  Working from the input means such a tool can cost masking, which
the checks in the callers will catch, but never sequence.
"""

from cactus.shared.common import cactus_call

from toil.realtimeLogger import RealtimeLogger


def prefilter_cmd(raw_fa_path, prefilter_opts):
    """Command that drops the tiny and low-entropy contigs that upset the maskers.

    The same command with -x appended yields exactly the contigs this one removes.
    """
    cmd = ['cactus_redPrefilter', raw_fa_path]
    if prefilter_opts:
        assert '-x' not in prefilter_opts and '--extract' not in prefilter_opts
        cmd += prefilter_opts.split()
    return cmd


def masked_base_count(fa_path):
    """Number of masked bases in a fasta.

    Counts soft-masked bases and hard-masked Ns alike, since that is what
    cactus_softmask2hardmask -b reports.
    """
    awkres = cactus_call(parameters=[['cactus_softmask2hardmask', '-b', fa_path],
                                     ['awk', '{sum += $3-$2} END {print sum}']],
                         check_output=True, rt_log_cmd=False).strip()
    return int(float(awkres)) if awkres else 0


def extract_masking_bed(fa_path, bed_path):
    """Write the masked intervals of a fasta to bed_path."""
    cactus_call(parameters=['cactus_softmask2hardmask', '-b', fa_path], outfile=bed_path)


def soft_mask_intervals(in_fa_path, bed_path, out_fa_path, unmask=False):
    """Copy in_fa_path to out_fa_path, soft-masked at the intervals in bed_path.

    Masking already present in the input survives by never being overwritten, unless
    unmask is set, which strips it and leaves only the intervals in the bed.
    """
    cmd = ['cactus_fasta_softmask_intervals.py', '--origin=zero']
    if unmask:
        cmd += ['--unmask']
    cmd += [bed_path]
    cactus_call(parameters=cmd, infile=in_fa_path, outfile=out_fa_path)


def log_masking_delta(tool_name, event_name, pre_mask_size, post_mask_size):
    """Report what the masker changed.

    A decrease should not be possible now that the intervals are applied to the input,
    so it is worded rather than signed, to keep a negative from reading as normal.
    """
    RealtimeLogger.info('{} masked {} bp of {}, {} masking from {} to {}'.format(
        tool_name, abs(post_mask_size - pre_mask_size), event_name,
        'decreasing' if post_mask_size < pre_mask_size else 'increasing',
        pre_mask_size, post_mask_size))
