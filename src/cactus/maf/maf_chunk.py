#!/usr/bin/env python3

"""
cactus-phast helpers for chunking a MAF (.maf / .maf.gz / .taf / .taf.gz)
into reference-coordinate chunks using a taffy index (.tai).

The job functions here read options.refGenome and options.chunkCores off
the cactus-phast options namespace, so they're not standalone-reusable
without passing a compatible options object.

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


def plan_chunks(ref_sequence_lengths, chunk_size):
    """ Slice each reference contig into chunks of ~chunk_size bp. Big contigs
    get split; the last piece may absorb a small remainder rather than create
    a tiny tail chunk (only if doing so keeps it under 1.5x chunk_size). Small
    contigs become their own single-chunk piece (no cross-contig packing —
    each chunk is one contiguous range on one contig).

    Returns a list of (contig, start, end) tuples ordered by (contig, start). """
    chunks = []
    for contig in sorted(ref_sequence_lengths.keys()):
        L = ref_sequence_lengths[contig]
        if L <= 0:
            continue
        if L <= chunk_size:
            chunks.append((contig, 0, L))
            continue
        # L > chunk_size, so n is always >= 1 here.
        n = L // chunk_size
        last = L - (n - 1) * chunk_size
        # Avoid a final chunk much larger than chunk_size; absorbing was for
        # tiny remainders, but if `last` would exceed 1.5x, split off another
        # chunk instead.
        if last > chunk_size * 1.5:
            n += 1
            last = L - (n - 1) * chunk_size
        for i in range(n):
            start = i * chunk_size
            size = chunk_size if i < n - 1 else last
            chunks.append((contig, start, start + size))
    return chunks


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
