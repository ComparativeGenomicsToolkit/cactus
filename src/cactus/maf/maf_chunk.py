#!/usr/bin/env python3

"""
cactus-phast helpers for chunking a MAF (.maf / .maf.gz / .taf / .taf.gz)
into reference-coordinate chunks using a taffy index (.tai).

The Toil job functions here read options.refGenome and options.chunkCores off
the cactus-phast options namespace, so they're not standalone-reusable without
passing a compatible options object. The pure helpers (parse_bed_ranges,
_slice_region, plan_chunks, plan_chunks_in_regions) take plain arguments;
parse_bed_ranges in particular is also imported by cactus-hal2maf so the two
tools interpret --bedRanges identically.

The model:
  - A single sequential job builds a .tai if the user didn't supply one.
  - A planning step turns the reference-genome sequence lengths into a flat
    list of chunk specs (contig, start, end), each at most chunk_size bp on
    a single contig.
  - A single multi-core 'chunker' job localizes the source MAF once, runs
    one `taffy view -m -r contig:start-end -c -o chunk.maf.gz` per chunk
    via GNU parallel, and writes each chunk to the Toil file store.
  - Downstream stages take a chunk file ID (or a small list of consecutive
    chunk file IDs for joint 4d-site extraction) and run a linear, single-
    process pipeline against it. No further parallel orchestration needed.
"""

import os
import shlex

from toil.realtimeLogger import RealtimeLogger

from cactus.shared.common import cactus_call


def taffy_index_job(job, maf_id, maf_basename):
    """ Build a .tai for the given source MAF. Single sequential pass. """
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, maf_basename)
    RealtimeLogger.info("Reading MAF file from job store to {}".format(maf_path))
    job.fileStore.readGlobalFile(maf_id, maf_path)
    RealtimeLogger.info('Building taffy index for {}'.format(maf_basename))
    cactus_call(parameters=['taffy', 'index', '-i', maf_path])
    return job.fileStore.writeGlobalFile(maf_path + '.tai')


def get_ref_sequence_lengths(hal_path, ref_genome):
    """ Run halStats --sequenceStats and return {contig: length}. """
    res = cactus_call(parameters=['halStats', hal_path, '--sequenceStats', ref_genome],
                      check_output=True)
    lengths = {}
    for line in res.strip().split('\n')[1:]:
        toks = line.strip().split(',')
        if len(toks) == 4:
            lengths[toks[0].strip()] = int(toks[1].strip())
    return lengths


def get_aligned_ref_contigs(hal_path, ref_genome):
    """ Return the set of `ref_genome` sequence names that have at least one base
    aligned to the parent genome, via halAlignedExtract. Sequences absent from
    the result are completely unaligned — in the MAF they only ever appear as
    reference-only columns, which phyloP cannot score. halAlignedExtract is a
    single linear scan of the reference's top-segment array (it touches no other
    genome and no bottom segments). The aligned-region BED can be huge, so it's
    streamed through a single-pass awk that emits each contig name once — no
    O(n log n) sort. Because halAlignedExtract reports only alignment to the
    PARENT, this is equivalent to "appears in a multi-row MAF column" only for a
    LEAF reference; callers must restrict it to leaf references (an internal node
    can align to its descendants, which this does not see). For a parentless
    genome (the HAL root) halAlignedExtract exits 0 with empty output — it does
    not error — so this would return an empty set; callers must guard against
    that too. """
    res = cactus_call(parameters=[['halAlignedExtract', hal_path, ref_genome],
                                  ['awk', '!seen[$1]++{print $1}']],
                      check_output=True)
    return set(name for name in res.split('\n') if name.strip())


def _slice_region(contig, start, end, chunk_size):
    """ Slice [start, end) on `contig` into chunks of ~chunk_size bp. The last
    piece may absorb a small remainder rather than create a tiny tail chunk
    (only if doing so keeps it under 1.5x chunk_size). Returns a list of
    (contig, s, e) tuples in ascending order. """
    L = end - start
    if L <= 0:
        return []
    if L <= chunk_size:
        return [(contig, start, end)]
    # L > chunk_size, so n is always >= 1 here.
    n = L // chunk_size
    last = L - (n - 1) * chunk_size
    # Avoid a final chunk much larger than chunk_size; absorbing was for tiny
    # remainders, but if `last` would exceed 1.5x, split off another chunk.
    if last > chunk_size * 1.5:
        n += 1
        last = L - (n - 1) * chunk_size
    out = []
    for i in range(n):
        s = start + i * chunk_size
        size = chunk_size if i < n - 1 else last
        out.append((contig, s, s + size))
    return out


def plan_chunks(ref_sequence_lengths, chunk_size):
    """ Slice each reference contig into chunks of ~chunk_size bp. Big contigs
    get split; small contigs become their own single-chunk piece (no cross-
    contig packing — each chunk is one contiguous range on one contig).

    Returns a list of (contig, start, end) tuples ordered by (contig, start). """
    chunks = []
    for contig in sorted(ref_sequence_lengths.keys()):
        chunks += _slice_region(contig, 0, ref_sequence_lengths[contig], chunk_size)
    return chunks


def parse_bed_ranges(bed_path):
    """ Parse a BED file into a list of (sequence, start, end) tuples from the
    first three whitespace-separated columns (0-based, half-open). Blank lines,
    '#' comments, and 'track'/'browser' header lines (matched on the first token,
    so a contig literally named e.g. 'track1' is NOT skipped) are ignored, as are
    lines with fewer than three columns. A line with >=3 columns whose start/end
    aren't integers is a hard error naming the file and the offending line text
    (rather than a bare ValueError escaping). Shared by cactus-hal2maf and
    cactus-phast so both interpret --bedRanges identically. """
    ranges = []
    with open(bed_path, 'r') as bed_file:
        for line in bed_file:
            if not line.strip() or line.startswith('#'):
                continue
            toks = line.split()
            if toks[0] in ('track', 'browser') or len(toks) < 3:
                continue
            try:
                ranges.append((toks[0], int(toks[1]), int(toks[2])))
            except ValueError:
                raise RuntimeError('unparseable BED line in {}: {!r}'.format(
                    bed_path, line.rstrip('\n')))
    return ranges


def plan_chunks_in_regions(bed_regions, ref_sequence_lengths, chunk_size):
    """ Like plan_chunks, but restricted to `bed_regions` (a list of
    (contig, start, end) parsed from a BED). Each region is clamped to its
    contig's length; regions on the same contig that overlap or touch are merged
    (so chunks never produce duplicate reference positions); the merged spans are
    then sliced into chunks of ~chunk_size bp. Regions whose contig is absent
    from `ref_sequence_lengths` are dropped.

    Returns (chunks, dropped_contigs, clamped) where chunks is the
    (contig, start)-ordered chunk list, dropped_contigs is the sorted list of BED
    contigs not present in the reference, and clamped is the list of input
    (contig, start, end) intervals that extended past [0, contig_length] (so they
    were clamped, or dropped if nothing remained). Note: in-bounds intervals that
    are empty (end <= start) are silently dropped and appear in neither list. """
    by_contig = {}
    dropped = set()
    clamped = []
    for (contig, start, end) in bed_regions:
        if contig not in ref_sequence_lengths:
            dropped.add(contig)
            continue
        L = ref_sequence_lengths[contig]
        s = max(0, start)
        e = min(end, L)
        if start < 0 or end > L:
            clamped.append((contig, start, end))
        if e > s:
            by_contig.setdefault(contig, []).append((s, e))
    chunks = []
    for contig in sorted(by_contig.keys()):
        spans = sorted(by_contig[contig])
        cs, ce = spans[0]
        merged = []
        for (s, e) in spans[1:]:
            if s <= ce:            # overlapping or touching -> extend
                ce = max(ce, e)
            else:
                merged.append((cs, ce))
                cs, ce = s, e
        merged.append((cs, ce))
        for (s, e) in merged:
            chunks += _slice_region(contig, s, e, chunk_size)
    return chunks, sorted(dropped), clamped


def filter_chunks_to_indexed(chunks, tai_path, ref_genome):
    """ Drop chunk specs whose contig has no entry in the source .tai. The
    .tai's named-seq lines look like `<refGenome>.<contig>\t...`; chunks for
    contigs that don't appear are unalignable (no MAF blocks) and would
    cause taffy view to fail with 'Region not found in taffy index'. """
    prefix = ref_genome + '.'
    present = set()
    with open(tai_path) as f:
        for line in f:
            tok = line.split('\t', 1)[0]
            if tok.startswith(prefix):
                present.add(tok[len(prefix):])
    return [c for c in chunks if c[0] in present]


def chunker_job(job, options, maf_id, tai_id, source_basename, chunks, filter_cmd=None):
    """ Single multi-core job: localize the source MAF + .tai once, then run
    one `taffy view -m -r contig:start-end` per chunk via GNU parallel.
    Returns a list of (contig, start, end, chunk_id) tuples in plan order.
    This is the only job in the workflow that holds the full source MAF and
    the only job that uses GNU parallel.

    `filter_cmd`, if given, is a shell pipeline fragment spliced between
    taffy view and bgzip — e.g. `<strip-perl> | mafDuplicateFilter -k -m -`.
    With None, the chunks are emitted as taffy view's bgzipped output.

    Chunks for contigs whose ref range has no aligned MAF blocks (sparse
    alt / random / unplaced contigs, centromeric gaps) come back as a
    header-only MAF; downstream consumers' _maf_has_blocks check skips
    them. """
    work_dir = job.fileStore.getLocalTempDir()
    src_path = os.path.join(work_dir, source_basename)
    tai_path = src_path + '.tai'
    RealtimeLogger.info('Reading source MAF from job store to {}'.format(src_path))
    job.fileStore.readGlobalFile(maf_id, src_path)
    job.fileStore.readGlobalFile(tai_id, tai_path)

    cmd_file = os.path.join(work_dir, 'chunker.cmds')
    chunk_paths = []
    with open(cmd_file, 'w') as f:
        for i, (contig, start, end) in enumerate(chunks):
            # Encode region coords in the filename so any downstream error
            # message immediately tells you which slice of the MAF blew up.
            chunk_path = os.path.join(work_dir, 'chunk.{}_{}-{}.maf.gz'.format(contig, start, end))
            chunk_paths.append(chunk_path)
            region = '{}.{}:{}-{}'.format(options.refGenome, contig, start, end)
            if filter_cmd:
                # bash -o pipefail wrapper so any failure in taffy view /
                # filter_cmd / bgzip propagates to parallel's --halt.
                pipeline = 'taffy view -m -r {region} -i {src} | {filt} | bgzip > {out}'.format(
                    region=shlex.quote(region),
                    src=shlex.quote(src_path),
                    filt=filter_cmd,
                    out=shlex.quote(chunk_path))
                f.write('bash -o pipefail -c {}\n'.format(shlex.quote(pipeline)))
            else:
                f.write('taffy view -m -c -r {region} -i {src} -o {out}\n'.format(
                    region=shlex.quote(region),
                    src=shlex.quote(src_path),
                    out=shlex.quote(chunk_path)))

    n_par = options.chunkCores or 1
    RealtimeLogger.info('Chunking source MAF into {} chunks with -j {}'.format(
        len(chunks), n_par))
    cactus_call(parameters=['bash', '-c',
        'set -eo pipefail && cat {} | parallel --halt now,fail=1 -j {} \'{{}}\''.format(
            shlex.quote(cmd_file), n_par)])

    # Free the source MAF + .tai before the upload loop so peak local disk
    # during writeGlobalFile is just the chunks (not src + chunks).
    os.remove(src_path)
    os.remove(tai_path)

    # writeGlobalFile each chunk and pair the file IDs back with their plan
    # entries; remove the local copy as we go to keep disk pressure bounded.
    out = []
    for (contig, start, end), chunk_path in zip(chunks, chunk_paths):
        chunk_id = job.fileStore.writeGlobalFile(chunk_path)
        out.append((contig, start, end, chunk_id))
        os.remove(chunk_path)
    return out


def group_chunks_for_fit(chunks_with_ids, group_size):
    """ Partition the chunker's output into consecutive same-contig groups of
    size <= group_size for joint 4d-site extraction. Each group is a list of
    (contig, start, end, chunk_id) tuples on the same contig. Groups never
    cross contig boundaries. """
    groups = []
    current = []
    current_contig = None
    for spec in chunks_with_ids:
        contig = spec[0]
        if current and (contig != current_contig or len(current) >= group_size):
            groups.append(current)
            current = []
        current_contig = contig
        current.append(spec)
    if current:
        groups.append(current)
    return groups
