#!/usr/bin/env python3
"""Check that a finished HAL still contains the sequence that went into it.

cactus builds the HAL out of cactus_consolidated's output rather than out of the
input fasta, so nothing downstream ever compares the two.  A fasta that was
silently truncated on the way in -- the failure this exists to catch -- yields a
perfectly valid HAL, and every tool that reads it carries the damage forward.
This walks the alignment back to its input and insists the bases are the same.

What is allowed to differ:

  case         soft-masking is a case change, and every masking preprocessor
               makes them, so both sides are compared upper-cased
  IUPAC codes  cactus_sanitizeFastaHeaders folds R, Y, S ... to N on the way in,
               so the input is folded the same way before it is compared
  names        the sanitizer cuts headers at whitespace, adds an "id=<genome>|"
               prefix (which cactus_consolidated later strips again) and, in
               pangenome mode, drops everything through the last '#' and
               rewrites "chr1:10-15" as "chr1_sub_9_15".  hal_sequence_name()
               replays all of that, and which convention a given HAL used is
               worked out from the HAL itself
  order        sequences are paired by name, never by position
  subpaths     input that was cut into pieces before the run -- by hand, as the
               HPRC recipe does for a misjoined contig, or by a clipping
               preprocessor -- reaches the HAL as "contig_sub_START_END".  The
               piece is checked against exactly that range of the input contig.
               Neither pipeline clips anything itself: cactus-graphmap-split
               passes whole contigs through, and cactus-graphmap-join clips the
               graphs, not the HAL

Anything else -- a base that changed, a sequence that is short, missing or
unaccounted for -- is reported.

Hard-masking is deliberately not tolerated.  The in-pipeline --validate compares
the HAL against the sequence the preprocessor produced, so any hard-masking has
already been applied to both sides; checkPreprocessedSequence is what covers the
step that does the masking.  Running this by hand against a pre-preprocessor
seqfile with a hard-masking preprocessor configured will therefore report base
differences, which is why --validate uses the post-preprocessor sequence.
"""

import gzip
import os
import re
import sys
import timeit
import xml.etree.ElementTree as ET

from cactus.preprocessor.checkPreprocessedSequence import check_digests_match
from cactus.preprocessor.checkPreprocessedSequence import digest_stream
from cactus.preprocessor.checkPreprocessedSequence import iter_fasta_stream
from cactus.preprocessor.checkPreprocessedSequence import preprocessed_fasta_id
from cactus.preprocessor.checkPreprocessedSequence import sequence_checksum
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import cactus_call
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import cactus_walltime
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactusRootPath
from cactus.shared.common import catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import findRequiredNode
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import importSingularityImage
from cactus.shared.common import makeURL
from cactus.shared.common import setupBinaries
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.version import cactus_commit

from toil.common import Toil
from toil.fileStores import FileID
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options

from sonLib.bioio import getTempFile
from sonLib.nxnewick import NXNewick

# a piece cut out of a contig is named after the range it came from: 0-based
# start, end open.  Cutting can nest (a chromosome split out of an already cut
# assembly), so this is applied until it stops matching.
_SUBPATH = re.compile(r'^(.+)_sub_(\d+)_(\d+)$')

# hal2fasta writes 80-column fasta by default, which costs it a getSubString
# call per line.  Wide lines are ~30% faster and the digest does not care.
_LINE_WIDTH = '1000000'

def _convert_faidx_range(name):
    """Rewrite a trailing samtools range the way cactus_sanitizeFastaHeaders -p does.

    "chr1:10-15" (1-based, end inclusive) becomes "chr1_sub_9_15" (0-based, end
    open).  Any other colon becomes an underscore, because a bare ':' crashes
    parts of the HAL toolchain.
    """
    if ':' not in name:
        return name
    colon = name.rfind(':')
    dash = name.rfind('-')
    if dash > colon + 1:
        start_text, end_text = name[colon + 1:dash], name[dash + 1:]
        if start_text.isdigit() and end_text.isdigit():
            start, end = int(start_text), int(end_text)
            # samtools treats 0 like 1, and so does the sanitizer
            if start > 0:
                start -= 1
            name = '{}_sub_{}_{}'.format(name[:colon], start, end)
    return name.replace(':', '_')


def hal_sequence_name(header, pangenome=False):
    """The name the first word of a fasta header ends up with inside the HAL.

    Mirrors cactus_sanitizeFastaHeaders' addUniqueFastaPrefix() followed by
    cactus_consolidated's stripUniqueIdsFromLeafSequences().  Keep the three in
    step; cactus_validateTest pins this one against the real binary.
    """
    fields = header.split()
    name = fields[0] if fields else ''
    if pangenome:
        pound = name.rfind('#')
        if pound >= 0:
            name = name[pound + 1:]
        name = _convert_faidx_range(name)
    # the sanitizer adds "id=<genome>|" unless the header already carries one,
    # and cactus_consolidated drops that first '|'-separated token again
    if name.startswith('id='):
        bar = name.find('|')
        if bar >= 0:
            name = name[bar + 1:]
    return name


def _fasta_parts(path):
    """The files making up one genome's sequence; cactus accepts a directory."""
    if os.path.isdir(path):
        return [os.path.join(path, part) for part in sorted(os.listdir(path))]
    return [path]


def _open_fasta(path):
    """Open a fasta for reading, transparently un-gzipping it."""
    handle = open(path, 'rb')
    if handle.read(2) == b'\x1f\x8b':
        handle.close()
        return gzip.open(path, 'rb')
    handle.seek(0)
    return handle


def input_digest(fasta_path):
    """Digest an input genome, under its own header names.

    Bases are folded the way cactus_sanitizeFastaHeaders folds them.  Zero-length
    records are dropped: the sanitizer refuses to pass those on, so the HAL never
    has them.
    """
    records = []
    for part in _fasta_parts(fasta_path):
        with _open_fasta(part) as handle:
            records += [record for record in digest_stream(handle, part, fold_ambiguity=True)
                        if record[1] > 0]
    return records


def map_hal_names(records, pangenome):
    """Rename a digest to the names the HAL should be holding it under.

    Returns (records, raw name by HAL name).  Renaming is pure string work, so
    both conventions can be tried without reading the fasta twice.
    """
    renamed = []
    raw_by_hal = {}
    for name, length, crc in records:
        hal_name = hal_sequence_name(name, pangenome=pangenome)
        renamed.append((hal_name, length, crc))
        raw_by_hal[hal_name] = name
    return renamed, raw_by_hal


def expected_records(fasta_path, pangenome=False):
    """Digest an input genome as the HAL should hold it."""
    return map_hal_names(input_digest(fasta_path), pangenome)


def read_read_only(job, file_id, path):
    """Put a job store file where a job can read it, without copying it if possible.

    Everything this tool reads is read-only and read *sequentially*: halStats
    wants metadata, hal2fasta streams one genome start to finish, and the input
    fastas are digested in one pass.  That is an access pattern a network
    filesystem is happy to serve directly, so toil is told it may hand back a
    symlink into the job store rather than a copy.

    cactus-hal2maf goes the other way and copies the HAL into every batch on
    purpose: it queries the file aggressively and in parallel, and pointing that
    at network storage would hurt the whole cluster.  The trade is worth making
    there and not here -- do not "fix" this by copying.  Copying would also cost
    far more: a copy is the whole HAL per batch, where reading through the
    symlink costs only the genomes that batch actually looks at.  The total
    number of concurrent readers is batchCount x batchParallelGenomes, which is
    the knob to turn down if the filesystem does complain.

    Note toil's --symlinkJobStoreReads only grants permission (it is True by
    default); read_file symlinks only when the caller asks as well, which is
    what this exists to do.  Without it toil hardlinks when the job store and
    the work directory share a filesystem and otherwise copies.  A job store
    that cannot symlink at all (S3) falls back to a real download, which is why
    the jobs still ask for disk enough to hold one.
    """
    RealtimeLogger.info('Reading {} from the job store to {}'.format(
        os.path.basename(path), path))
    job.fileStore.readGlobalFile(file_id, path, symlink=True)
    return path


def hal_genomes(hal_path):
    """Every genome in the alignment, ancestors included."""
    return cactus_call(parameters=['halStats', '--genomes', hal_path],
                       check_output=True).split()


def hal_summary(hal_path):
    """Everything worth knowing about an alignment, from one halStats call.

    Returns (tree, {genome: (length in bases, is a leaf)}).  Bare halStats
    prints the newick and then a row per genome, and it reads only metadata, so
    this stays cheap however large the file is.  The lengths matter because they
    are the only per-genome size available before anything is fetched: with a
    seqfile full of URLs there is no local file to measure.
    """
    output = cactus_call(parameters=['halStats', hal_path], check_output=True)
    tree_string, info = None, {}
    for line in output.splitlines():
        line = line.strip()
        if not line or line.startswith('hal v'):
            continue
        if tree_string is None and line.endswith(';'):
            tree_string = line
            continue
        fields = [field.strip() for field in line.split(',')]
        if len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit():
            info[fields[0]] = (int(fields[2]), int(fields[1]) == 0)
    if tree_string is None or not info:
        raise RuntimeError('could not read the genome list out of "halStats {}"'.format(hal_path))
    return MultiCactusTree(NXNewick().parseString(tree_string, addImpliedRoots=False)), info


def looks_like_pangenome(mc_tree, graph_event):
    """True if the alignment has the shape cactus-pangenome makes.

    Those alignments are stars: every genome is a leaf hanging straight off a
    single ancestor.  The minigraph leaf is the other giveaway, but only the
    per-chromosome HALs still have it -- cactus-graphmap-join strips it out of
    the merged one -- so the shape has to be able to carry the decision alone.
    A progressive alignment of a flat seqfile is a star too -- and so is every
    two-genome subproblem cactus-prepare emits, which is why a pair on its own is
    not taken as evidence -- so validate_genome() checks the answer against the
    names in the HAL rather than trusting it.
    """
    root = mc_tree.rootId
    names = set()
    for node in mc_tree.breadthFirstTraversal():
        if node == root:
            continue
        if not mc_tree.isLeaf(node) or mc_tree.getParent(node) != root:
            return False
        if mc_tree.hasName(node):
            names.add(mc_tree.getName(node))
    return graph_event in names or len(names) > 2


def hal_records(hal_path, genome, work_dir, lengths_only=False, job_memory=None):
    """Digest one genome as the HAL actually holds it."""
    if lengths_only:
        records = []
        sizes = cactus_call(parameters=['halStats', '--chromSizes', genome, hal_path],
                            check_output=True)
        for line in sizes.splitlines():
            fields = line.split('\t')
            if len(fields) == 2:
                records.append((fields[0], int(fields[1]), 0))
        return records

    # hal2fasta only streams, so this goes through a file rather than into the
    # digest directly: cactus_call has no way to hand back a pipe, and buffering
    # a genome in memory is worse than spending the scratch space
    fa_path = os.path.join(work_dir, 'hal_genome.fa')
    try:
        cactus_call(parameters=['hal2fasta', '--lineWidth', _LINE_WIDTH, hal_path, genome],
                    outfile=fa_path, job_memory=job_memory)
        with open(fa_path, 'rb') as handle:
            return digest_stream(handle, fa_path)
    finally:
        if os.path.exists(fa_path):
            os.remove(fa_path)


def resolve_subpath(name, lengths_by_name):
    """Map a cut-out piece of a contig back to (input name, start, end).

    None if the name is not a subpath of anything in lengths_by_name, or if the
    range it claims runs off the end of that sequence.
    """
    start, end = 0, None
    current = name
    while current not in lengths_by_name:
        match = _SUBPATH.match(current)
        if not match:
            return None
        piece_start, piece_end = int(match.group(2)), int(match.group(3))
        # the window we are tracking is relative to [piece_start, piece_end) of
        # the name one level up, so shift it into those coordinates
        end = piece_end if end is None else piece_start + end
        start = piece_start + start
        current = match.group(1)
    if end is None:
        # the name was never a subpath at all; the caller pairs those by name
        end = lengths_by_name[current]
    if start >= end or end > lengths_by_name[current]:
        return None
    return current, start, end


def _subpath_records(fasta_path, requests, raw_by_hal, lengths_only=False):
    """Digest the slices of the input that the HAL's cut-out pieces claim.

    requests maps a HAL sequence name to (base name, start, end).
    """
    if lengths_only:
        return {hal_name: (hal_name, end - start, 0)
                for hal_name, (_, start, end) in requests.items()}

    wanted = {}
    for hal_name, (base, start, end) in requests.items():
        wanted.setdefault(raw_by_hal[base], []).append((hal_name, start, end))

    found = {}
    for part in _fasta_parts(fasta_path):
        with _open_fasta(part) as handle:
            for raw_name, sequence in iter_fasta_stream(handle, fold_ambiguity=True):
                for hal_name, start, end in wanted.get(raw_name, []):
                    piece = sequence[start:end]
                    found[hal_name] = (hal_name, len(piece), sequence_checksum(piece))
    for hal_name in requests:
        # only reachable if the fasta changed under us between the two passes
        found.setdefault(hal_name, (hal_name, 0, 0))
    return found


def _summarise(records, limit=5):
    head = ', '.join('{} ({} bp)'.format(name, length) for name, length, _ in records[:limit])
    if len(records) > limit:
        head += ', ...'
    return head


def validate_genome(hal_path, genome, fasta_path, work_dir, pangenome=False,
                    force_convention=False, lengths_only=False, allow_missing=False,
                    source=None, job_memory=None):
    """Compare one genome in the HAL against its input fasta.

    fasta_path is the copy on local disk, which under toil is a scratch file
    nobody will recognise; source is what to call it in a message, and defaults
    to naming no file at all.

    Returns a list of problems, empty when the genome checks out.
    """
    context = 'the HAL copy of genome "{}"'.format(genome)
    named_source = ' in {}'.format(source) if source else ''
    problems = []

    raw = input_digest(fasta_path)
    if lengths_only:
        # halStats gives no bases to checksum, so drop the input's checksums too
        raw = [(name, length, 0) for name, length, _ in raw]
    actual = hal_records(hal_path, genome, work_dir, lengths_only=lengths_only,
                         job_memory=job_memory)
    actual_by_name = {record[0]: record for record in actual}

    expected, raw_by_hal = map_hal_names(raw, pangenome)
    if not force_convention:
        # the naming convention was guessed from the shape of the alignment, so
        # check the guess against the names actually in the HAL.  Renaming is
        # pure string work, so trying the other convention costs nothing.
        other, other_raw_by_hal = map_hal_names(raw, not pangenome)
        hits = sum(1 for name, _, _ in expected if name in actual_by_name)
        other_hits = sum(1 for name, _, _ in other if name in actual_by_name)
        if other_hits > hits:
            expected, raw_by_hal = other, other_raw_by_hal

    if len(raw_by_hal) != len(expected):
        # the sanitizer rejects this outright, so it means the fasta being
        # compared against is not the one the alignment was built from
        counts = {}
        for name, _, _ in expected:
            counts[name] = counts.get(name, 0) + 1
        repeats = sorted(name for name, count in counts.items() if count > 1)
        problems.append('{}: the input fasta{} has {} duplicated sequence name(s) after '
                        'sanitisation: {}'.format(context, named_source, len(repeats),
                                                  ', '.join(repeats[:10])))
        return problems

    lengths_by_name = {name: length for name, length, _ in expected}
    expected_by_name = {record[0]: record for record in expected}

    paired_expected, paired_actual = [], []
    subpath_requests, unresolved = {}, []
    matched_names, subpath_bases = set(), set()
    for record in actual:
        name = record[0]
        if name in expected_by_name:
            paired_expected.append(expected_by_name[name])
            paired_actual.append(record)
            matched_names.add(name)
            continue
        resolved = resolve_subpath(name, lengths_by_name)
        if resolved is None:
            unresolved.append(record)
        else:
            subpath_requests[name] = resolved
            subpath_bases.add(resolved[0])

    if unresolved:
        problems.append(
            '{} contains {} sequence(s) that are not in the input fasta{}: {}'.format(
                context, len(unresolved), named_source, _summarise(unresolved)))

    if subpath_requests:
        found = _subpath_records(fasta_path, subpath_requests, raw_by_hal,
                                 lengths_only=lengths_only)
        for name in sorted(subpath_requests):
            paired_expected.append(found[name])
            paired_actual.append(actual_by_name[name])

    absent = [expected_by_name[name] for name in lengths_by_name
              if name not in matched_names and name not in subpath_bases]
    if absent:
        absent.sort(key=lambda record: -record[1])
        message = ('{} is missing {} of the {} input sequence(s){}, totalling {} bp: {}'
                   .format(context, len(absent), len(expected), named_source,
                           sum(length for _, length, _ in absent), _summarise(absent)))
        if allow_missing:
            problems.append('WARNING: ' + message)
        else:
            problems.append(message + '. Sequences are legitimately left out by pangenome '
                            'chromosome binning and by aligning a subtree, so pass '
                            '--allowMissing if that is what happened here.')

    try:
        check_digests_match(paired_expected, paired_actual, context)
    except RuntimeError as e:
        problems.append(str(e))

    return problems


def validate_hal_export(job, hal_path, work_dir, mc_tree, root_node, seq_id_map):
    """Check a freshly exported HAL against the sequence that was aligned.

    This is what --validate runs, from inside export_hal, so the HAL is still on
    local disk and the comparison costs one pass over it.  The sequences here
    have already been through cactus_sanitizeFastaHeaders, so no header mangling
    has to be replayed, and nothing that went in is allowed to be missing.

    Returns the problems rather than raising them.  export_hal is expensive and
    a bad alignment fails this check every time it is tried, so the caller fails
    from a cheap follow-on: raising here would have toil rebuild the whole HAL
    and check it again, five times over on a cluster, before giving up -- and
    with --doubleMem doubling the request each round until nothing can schedule
    it.
    """
    genomes = [name for name in (mc_tree.getName(node)
                                 for node in mc_tree.breadthFirstTraversal(root_node))
               if seq_id_map.get(name)]
    in_hal = set(hal_genomes(hal_path))

    problems = []
    missing = [name for name in genomes if name not in in_hal]
    if missing:
        problems.append('the HAL is missing {} genome(s) that were aligned: {}'.format(
            len(missing), ', '.join(missing)))

    fa_path = os.path.join(work_dir, 'validate_input.fa')
    for genome in genomes:
        if genome not in in_hal:
            continue
        seq_id = preprocessed_fasta_id(seq_id_map[genome])
        read_read_only(job, seq_id, fa_path)
        try:
            problems += validate_genome(hal_path, genome, fa_path, work_dir,
                                        force_convention=True, job_memory=job.memory)
        finally:
            # one genome at a time, so the job only needs room for the biggest
            job.fileStore.deleteLocalFile(seq_id)

    if problems:
        for problem in problems:
            job.fileStore.logToMaster('HAL VALIDATION: {}'.format(problem))
        return problems

    job.fileStore.logToMaster('Validated {} genome(s) of {} against their input sequence'.format(
        len(genomes), os.path.basename(hal_path)))
    return []


VALIDATE_SECS_PER_GB = 150

def fail_hal_validation(job, hal_name, problems):
    """Fail the workflow because the HAL did not match its input.

    A job of its own so that the retries toil will do here are free, instead of
    rebuilding the alignment each time.
    """
    raise RuntimeError('HAL validation of {} failed with {} problem(s): {}'.format(
        hal_name, len(problems), ' | '.join(problems)))


def validate_workflow(job, hal_id, seq_id_map, seq_paths, unreadable, options, config):
    """Plan the check, then run it in as many batches as were asked for."""
    # the plan job only reads the alignment's tree, which is metadata: no page-cache
    # allowance needed, unlike the batches below
    plan_job = job.addChildJobFn(validate_plan, hal_id, sorted(seq_id_map), unreadable,
                                 options, config, cores=1, disk=int(hal_id.size * 1.1),
                                 walltime=cactus_walltime())
    return plan_job.addFollowOnJobFn(validate_all, hal_id, seq_id_map, seq_paths,
                                     plan_job.rv(), options, walltime=cactus_walltime()).rv()


def validate_plan(job, hal_id, seqfile_genomes, unreadable, options, config):
    """Work out which genomes to check, and under which naming convention.

    Returns (problems so far, genomes to check, whether names are pangenome-style,
    each of those genomes' length in the HAL).
    """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = read_read_only(job, hal_id, os.path.join(
        work_dir, os.path.basename(options.halFile.replace(' ', '.'))))

    mc_tree, hal_info = hal_summary(hal_path)
    in_hal = set(hal_info)
    leaves = {name for name, (_, is_leaf) in hal_info.items() if is_leaf}
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"),
                                    "assemblyName", default="_MINIGRAPH_")
    pangenome = True if options.pangenome else looks_like_pangenome(mc_tree, graph_event)
    RealtimeLogger.info('{} has {} genomes ({} of them leaves) and looks like a {} '
                        'alignment{}'.format(
                            os.path.basename(hal_path), len(in_hal), len(leaves),
                            'pangenome' if pangenome else 'progressive',
                            ' (forced with --pangenome)' if options.pangenome else ''))

    problems = []
    missing_genomes = sorted(set(seqfile_genomes) - in_hal)
    if missing_genomes:
        message = 'the HAL is missing {} genome(s) from the seqfile: {}'.format(
            len(missing_genomes), ', '.join(missing_genomes))
        if options.allowMissing:
            problems.append('WARNING: ' + message)
        else:
            problems.append(message + '. Outgroups, --root and pangenome chromosome '
                            'binning all leave genomes out, so pass --allowMissing if '
                            'that is what happened here.')

    # a leaf came from an input fasta, so one the seqfile does not account for
    # means we are validating the wrong pair of files.  A genome the seqfile does
    # name but whose file could not be read is a different complaint, below.
    unaccounted = sorted(leaves - set(seqfile_genomes) - set(unreadable))
    if unaccounted:
        problems.append('the HAL contains {} leaf genome(s) that the seqfile does not '
                        'mention: {}'.format(len(unaccounted), ', '.join(unaccounted)))

    # a seqfile can name sequence files that were never made -- cactus-prepare
    # writes one entry per ancestor -- so this only matters for genomes that are
    # actually going to be looked at
    blocked = sorted(name for name in unreadable if name in in_hal)
    if blocked:
        problems.append('cannot check {} genome(s), their sequence file is missing: {}'.format(
            len(blocked), ', '.join('{} ({})'.format(name, unreadable[name]) for name in blocked)))

    genomes = sorted((set(seqfile_genomes) & in_hal) - set(blocked))
    lengths = {name: hal_info[name][0] for name in genomes}
    return problems, genomes, pangenome, lengths


def validate_all(job, hal_id, seq_id_map, seq_paths, plan, options):
    """Fan out one job per genome and gather what they find.

    One job per genome rather than batches of them, because read_read_only()
    means a job no longer pays to get at the HAL.  Batching exists in
    cactus-hal2maf to amortise a copy of the file over as much work as
    possible; with nothing to amortise, the small jobs win on every count.
    They retry one genome instead of a fiftieth of the run, they finish inside
    any sane walltime, the cluster's own scheduler decides how many go at once,
    and -- the big one -- each asks for memory the size of its own genome
    rather than a share of an 800 GB file.
    """
    problems, genomes, pangenome, lengths = plan
    if not genomes:
        return problems

    RealtimeLogger.info('Checking {} genome(s), one job each'.format(len(genomes)))
    results = []
    for genome in genomes:
        # the HAL's own length for the genome, which is all we have to go on
        # when the seqfile points at URLs that nothing has fetched yet.  The
        # input fasta can be bigger than the HAL's copy of it -- that is what
        # clipping means -- hence the headroom.
        length = lengths[genome]
        # room for the input fasta and for the copy extracted back out of the
        # hal; the hal itself is symlinked and costs nothing, except on a job
        # store that cannot symlink, where it has to be downloaded
        disk = int(hal_id.size * 1.1) + 3 * length
        # this job reads one genome out of the HAL, start to finish, so the
        # pages HDF5 mmaps -- which SLURM's cgroup accounting charges to the
        # job -- are bounded by that genome, not by the size of the file.  That
        # is why the whole-file allowance cactus-hal2maf needs is not needed
        # here: its batches touch the whole HAL, and one of these does not.
        #
        # measured on a 577-genome, 800 GB alignment: hal2fasta peaked at 42 MiB
        # across all 577, and no job came near its allowance.  So this does not
        # scale with the genome either -- a flat figure covers the streaming and
        # leaves room for the subpath check, which holds one contig of the input
        # in memory and so is bounded by the largest chromosome, not the genome.
        memory = options.validateMemory if options.validateMemory else \
            cactus_clamp_memory(4 * 1024**3)
        # two passes over the genome: hal2fasta streams it out of the HAL and the check reads it
        # back against the input fasta.  hal2fasta alone measured 60 s/GB at the p99 over the 576
        # per-ancestor exports of the VGP 577-way (HAL2FASTA_SECS_PER_GB in cactus_hal2seqfile),
        # and the comparison is the cheaper half.  The HAL arrives by symlink, so only the input
        # fasta and the extracted copy are I/O.
        results.append(job.addChildJobFn(validate_one, hal_id, seq_id_map[genome],
                                         seq_paths.get(genome), genome, pangenome, options,
                                         cores=1, disk=disk, memory=memory,
                                         walltime=cactus_walltime(VALIDATE_SECS_PER_GB * length / 1e9,
                                                                  io_bytes=2 * length)).rv())
    return job.addFollowOnJobFn(gather_problems, problems, results, walltime=cactus_walltime()).rv()


def gather_problems(job, problems, results):
    """Flatten what the per-genome jobs found onto what the plan already knew."""
    for result in results:
        problems += result
    return problems


def validate_one(job, hal_id, seq_source, source, genome, pangenome, options):
    """Check a single genome against its input fasta.

    seq_source is a URL when the seqfile pointed at one.  Fetching it here
    rather than on the leader is the whole point: a seqfile of 500 URLs would
    otherwise be downloaded one after another before any checking could start,
    and would leave every one of them sitting in the job store.  This way they
    come down in parallel, and each is thrown away as soon as it has been read.
    """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = read_read_only(job, hal_id, os.path.join(
        work_dir, os.path.basename(options.halFile.replace(' ', '.'))))

    fetched = None
    # a FileID subclasses str, so "is this a URL?" has to be asked the other way
    # round: anything that is not already in the job store is one to fetch
    if isinstance(seq_source, FileID):
        seq_id = preprocessed_fasta_id(seq_source)
    else:
        RealtimeLogger.info('Fetching {} for {}'.format(seq_source, genome))
        fetched = job.fileStore.import_file(seq_source)
        seq_id = fetched

    try:
        fa_path = read_read_only(job, seq_id, os.path.join(work_dir, 'input.fa'))
        RealtimeLogger.info('Checking {}'.format(genome))
        problems = validate_genome(hal_path, genome, fa_path, work_dir,
                                   pangenome=pangenome,
                                   force_convention=options.pangenome,
                                   lengths_only=options.lengthsOnly,
                                   allow_missing=options.allowMissing,
                                   source=source,
                                   job_memory=job.memory)
        # report as we go: the leader only prints its summary at the very end,
        # and on a run of hundreds of genomes you want to know at minute three
        for problem in problems:
            RealtimeLogger.warning(problem)
        return problems
    finally:
        if fetched is not None:
            # nothing else will ever want it, and 500 of these would otherwise
            # leave the whole input sitting in the job store
            job.fileStore.deleteGlobalFile(fetched)


def main():
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("seqFile", help="Seq file the alignment was built from")
    parser.add_argument("halFile", help="HAL file to check")
    parser.add_argument("--allowMissing", action="store_true",
                        help="Report sequences and genomes that are in the seqfile but not in "
                        "the HAL as warnings rather than errors. Expected when the alignment "
                        "was restricted with --root, has outgroups that were not exported, or "
                        "is a pangenome whose contigs did not all get binned to a chromosome")
    parser.add_argument("--lengthsOnly", action="store_true",
                        help="Compare only sequence lengths, which is much faster because the "
                        "HAL bases are never read. Catches truncation but not a substitution")
    parser.add_argument("--pangenome", action="store_true",
                        help="Force the pangenome header conventions (SAMPLE#HAP#CONTIG and "
                        "CONTIG:START-END). This is worked out from the HAL by default, so it "
                        "is only needed to overrule that")

    parser.add_argument("--validateMemory", type=human2bytesN, default=None,
                        help="Memory in bytes for each genome's job (defaults to an estimate "
                        "from that genome's own size). Standard suffixes like K, Ki, M, Mi, G "
                        "or Gi are supported (default=bytes)")

    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        if options.restart:
            problems = toil.restart()
        else:
            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)

            hal_id = toil.importFile(makeURL(options.halFile))
            seq_id_map, seq_paths, unreadable = resolve_seqfile(toil, options.seqFile)
            if not seq_id_map and not unreadable:
                raise RuntimeError('no sequences found in seqfile {}'.format(options.seqFile))
            problems = toil.start(Job.wrapJobFn(validate_workflow, hal_id, seq_id_map,
                                                seq_paths, unreadable, options, config,
                                                walltime=cactus_walltime()))

    end_time = timeit.default_timer()
    logger.info("cactus-validate has finished after {} seconds".format(end_time - start_time))

    # through the logger rather than stderr: a run of this size is watched
    # through --logFile, and a finding that only ever reached someone's terminal
    # is a finding nobody has a record of
    failures = [problem for problem in problems if not problem.startswith('WARNING: ')]
    for problem in problems:
        if problem.startswith('WARNING: '):
            logger.warning(problem[len('WARNING: '):])
    if failures:
        logger.error('VALIDATION OF {} FAILED: {} problem(s) found'.format(
            options.halFile, len(failures)))
        for problem in failures:
            logger.error('  {}'.format(problem))
        sys.exit(1)
    logger.info('{}: validated against {}'.format(options.halFile, options.seqFile))


def resolve_seqfile(toil, seqfile_path):
    """Work out where each genome's sequence is, importing only what has to be.

    Returns (genome -> file id or URL, genome -> the seqfile's own path for it,
    genome -> path that could not be read).

    Local files are imported here, which on a file job store is a symlink and so
    costs nothing.  A URL is left alone and handed to the job that needs it:
    importing 500 of them on the leader would download the entire input, one
    file at a time, before any checking could begin.

    The second map is only ever used to name a file in a message -- by the time
    a genome is checked, what is on disk is a scratch copy whose name would mean
    nothing to anyone.

    A local path that is not there yet is not fatal: cactus-prepare writes a
    seqfile entry for every ancestor, and only the genomes that turn out to be
    in the HAL have to be readable.
    """
    seq_file = SeqFile(seqfile_path)
    seq_id_map, seq_paths, unreadable = {}, {}, {}
    for genome, path in seq_file.pathMap.items():
        is_url = '://' in path
        if not is_url and not os.path.exists(path):
            unreadable[genome] = path
            continue
        seq_paths[genome] = path
        if is_url:
            seq_id_map[genome] = path
            continue
        if os.path.isdir(path):
            # cactus takes a directory of fastas for one genome, as cactus-align does
            merged = getTempFile()
            catFiles([os.path.join(path, part) for part in sorted(os.listdir(path))], merged)
            path = merged
        seq_id_map[genome] = toil.importFile(makeURL(path))
    return seq_id_map, seq_paths, unreadable


if __name__ == '__main__':
    main()
