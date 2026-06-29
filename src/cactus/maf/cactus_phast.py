#!/usr/bin/env python3

"""
cactus-phast: train a phyloFit neutral model and/or compute phyloP conservation
scores from a cactus-hal2maf MAF or TAF.

Input formats: .maf, .maf.gz, .taf, .taf.gz (auto-detected by taffy). The
single multi-core chunker job materializes each per-chunk shard as bgzipped
MAF via `taffy view -m -c`.

Two modes:
  --mode phyloFit : extract 4d sites from a gene annotation, train a neutral
                    model with phyloFit, save the .mod.
  --mode phyloP   : compute per-base phyloP scores using a neutral model.
                    If --neutralModel is not given but --geneAnnotation is,
                    train the model on the fly first.
"""

import os, sys
import re
import shlex
import shutil
import timeit
import xml.etree.ElementTree as ET

from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import cactus_cpu_count
from cactus.shared.version import cactus_commit
from cactus.progressive.cactus_prepare import human2bytesN

from cactus.maf.maf_chunk import (taffy_index_job, get_ref_sequence_lengths,
                                  get_aligned_ref_contigs,
                                  plan_chunks, plan_chunks_in_regions, parse_bed_ranges,
                                  filter_chunks_to_indexed,
                                  chunker_job, group_chunks_for_fit)

from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.humanize import bytes2human


def main():
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("inMaf", help="Input alignment as produced by cactus-hal2maf. "
                                      "Accepts .maf, .maf.gz, .taf, or .taf.gz (auto-detected "
                                      "by taffy). The chunker job materializes per-chunk "
                                      "files as bgzipped MAF regardless of source format.")
    parser.add_argument("halFile", help="HAL file matching inMaf (used for tree, species list, chrom sizes).")
    parser.add_argument("refGenome", help="Reference genome name (matches cactus-hal2maf --refGenome).")
    parser.add_argument("output",
                        help="Primary output. For --mode phyloFit a .mod file. For --mode phyloP "
                             "a .wig.gz file. Other intermediates are written next to it using the "
                             "same base name.")

    parser.add_argument("--mode", choices=["phyloFit", "phyloP"], default="phyloP",
                        help="phyloFit: train a neutral model only. phyloP: compute conservation "
                             "scores; trains a model first if only --geneAnnotation is given. "
                             "[default=phyloP]")

    parser.add_argument("--neutralModel", help="Pre-trained neutral model (.mod). Use with --mode phyloP "
                                              "to skip training.")
    parser.add_argument("--root", help="Restrict the workflow to species that descend from this "
                                       "internal-node name in the HAL tree (the 'root' of the "
                                       "subtree to use). phyloFit then trains only on the "
                                       "subtree's leaves; phyloP scores are computed against "
                                       "that subtree. The reference genome must be a leaf "
                                       "descendant of the named node.")
    parser.add_argument("--subtree", nargs='+', metavar='NAME',
                        help="Compute one output track per named internal node of the HAL tree. "
                             "For each name: phyloP --subtree <name> emits a per-base wig "
                             "expressing lineage-specific deviation from neutral on that subtree. "
                             "The reference genome must be a leaf descendant of every named "
                             "subtree (so the per-base score at each reference position has "
                             "subtree-side data). As a special case, naming the (effective) "
                             "tree root — the HAL tree root, or --root if set — produces the "
                             "standard whole-tree conservation track for that node (no phyloP "
                             "--subtree flag, no .s tag). So `--subtree squamata sauropsida root` "
                             "produces three tracks: two lineage-specific plus the global, all "
                             "anchored at the reference. With --root, every named subtree must "
                             "be a descendant of --root (or equal it, for the global-within-"
                             "clade track), and the reference must lie inside both. One model is "
                             "trained and reused across all tracks. NOTE: phyloP --subtree fits "
                             "a per-site likelihood for the subtree alone and will crash on "
                             "columns where the subtree has no aligned data; pick subtrees with "
                             "consistent coverage in your alignment.")
    parser.add_argument("--geneAnnotation", help="Gene annotation (genePred .gp, GFF/GFF3, GTF, or "
                                                "UCSC ncbiRefSeq.txt-style genePred with leading bin "
                                                "column). May be a local path or an http(s)/ftp URL. "
                                                "The pipeline downloads, converts to genePred, filters "
                                                "to records with cdsEnd > cdsStart, and runs "
                                                "genePredSingleCover to produce a non-redundant "
                                                "single-cover annotation. Required for --mode phyloFit; "
                                                "required for --mode phyloP if --neutralModel is not "
                                                "given. Coordinates must use the same contig names as "
                                                "--refGenome in the input MAF.")

    parser.add_argument("--mafIndex",
                        help="Pre-built taffy index (.tai) for the input MAF. Defaults to "
                             "`<inMaf>.tai` if it exists. If neither this option nor the default "
                             "sibling exists, an index is built as a sequential pre-pass and "
                             "saved next to the output for re-use.")

    parser.add_argument("--bedRanges",
                        help="Restrict the analysis to the reference ranges in this BED file "
                             "(coordinates in --refGenome, as with cactus-hal2maf --bedRanges). "
                             "Only these ranges are chunked, scored, and (when training) used for "
                             "4d-site extraction. Use it to skip contigs you don't want to process "
                             "— e.g. select only the large chromosomes and leave out tens of "
                             "thousands of tiny alt/unplaced contigs. BED sequence names must be "
                             "the reference contig names (the same names used in --geneAnnotation), "
                             "without the genome prefix; coordinates are 0-based half-open and "
                             "overlapping/touching intervals on a contig are merged. "
                             "[default: the whole reference]")

    parser.add_argument("--keepUnalignedContigs", action="store_true",
                        help="By default, when the reference is a leaf genome, cactus-phast queries "
                             "the HAL (via halAlignedExtract, a single scan of the reference's top "
                             "segments) for reference contigs that have no alignment to anything — "
                             "these only ever produce reference-only MAF columns that phyloP cannot "
                             "score — and excludes them so the chunker doesn't waste time "
                             "extracting them. This flag disables that check and processes every "
                             "reference contig. The check is skipped automatically for an "
                             "internal/ancestral reference (halAlignedExtract reports alignment to "
                             "the parent only, so it cannot see contigs that align solely to "
                             "descendants).")

    # phyloFit-related options
    parser.add_argument("--substMod", default="REV",
                        help="Substitution model passed to phyloFit --subst-mod. REV is the "
                             "model used by every published UCSC phyloP track (100-way, 470-way, "
                             "447-way Zoonomia, 200-mammals, etc.). SSREV is REV with a strand-"
                             "symmetric constraint and is the better choice for splice-site / "
                             "intronic-conservation analyses (it makes scores independent of "
                             "transcription strand) but produces a wider score distribution "
                             "and isn't comparable to the published tracks. [default=REV]")
    parser.add_argument("--precision", default="HIGH",
                        help="phyloFit --precision (HIGH|MED|LOW). Tighter EM convergence "
                             "yields more accurate parameter estimates at the cost of more "
                             "iterations. phyloFit is one job per workflow but its model is "
                             "reused by every phyloP chunk, so HIGH typically pays for itself. "
                             "[default=HIGH]")
    parser.add_argument("--modFreqs", action="store_true",
                        help="After phyloFit, run modFreqs to symmetrize background base frequencies.")

    # phyloP-related options
    parser.add_argument("--method", default="LRT", choices=["SPH", "LRT", "SCORE", "GERP"],
                        help="phyloP --method. LRT is recommended by phast for --wig-scores. "
                             "[default=LRT]")
    parser.add_argument("--phyloPMode", default="CONACC", choices=["CON", "ACC", "NNEUT", "CONACC"],
                        help="phyloP --mode. CONACC produces signed scores (positive=conservation, "
                             "negative=acceleration); CON / ACC are one-sided p-values; NNEUT is "
                             "two-sided. [default=CONACC]")
    parser.add_argument("--bigwig", action="store_true",
                        help="Also export a .bw file (run wigToBigWig on the merged wig).")

    # chunking
    parser.add_argument("--chunkSize", type=int, default=1_000_000,
                        help="Reference-coordinate size of each chunk, in bp. The source MAF is "
                             "split into chunks of this size by a single multi-core 'chunker' "
                             "job; each chunk is then a single Toil job for phyloP scoring (and, "
                             "in groups of --fitChunkGroup, for 4d-site extraction). "
                             "[default=1000000 (1Mb)]")
    parser.add_argument("--fitChunkGroup", type=int, default=10,
                        help="Number of consecutive same-contig chunks fed to one 4d-site "
                             "extraction job. Lets msa_view --4d amortize startup over a larger "
                             "region while keeping phyloP scoring at single-chunk granularity. "
                             "[default=10, i.e. 10 chunks of --chunkSize per 4d-extraction job]")

    parser.add_argument("--chunkCores", type=int, default=None,
                        help="Cores for the chunker job that splits the source MAF into "
                             "per-chunk shards via parallel `taffy view -r`. Required for "
                             "batch systems other than singleMachine; defaults to all "
                             "available cores on singleMachine.")

    # phyloFit job sizing
    parser.add_argument("--phyloFitMemory", type=human2bytesN, default=None,
                        help="Memory for the (single) phyloFit job. If unset, Toil's default "
                             "applies (~2 Gi); phyloFit on 1.4 GB SS / 577 species peaks at "
                             "~300 MiB so the default is fine for typical runs. Override only "
                             "if you observe an OOM, e.g. with much larger SS files or many "
                             "more species.")
    parser.add_argument("--phyloFitCores", type=int, default=None,
                        help="Cores for the phyloFit job (passed to phyloFit's --threads, which "
                             "uses OpenMP-parallelized tree likelihoods in phast v1.9.7+). "
                             "Required for batch systems other than singleMachine; defaults to "
                             "all available cores on singleMachine.")

    # Cactus boilerplate
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the specified pre-built container image rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    # validate inputs
    if not any(options.inMaf.endswith(s) for s in ('.maf', '.maf.gz', '.taf', '.taf.gz')):
        raise RuntimeError('inMaf must end with .maf, .maf.gz, .taf, or .taf.gz')

    if options.mode == 'phyloFit':
        if not options.geneAnnotation:
            raise RuntimeError('--geneAnnotation is required for --mode phyloFit')
        if options.neutralModel:
            logger.warning('--neutralModel is ignored when --mode phyloFit')
            options.neutralModel = None
    else:  # phyloP
        if not options.neutralModel and not options.geneAnnotation:
            raise RuntimeError('--mode phyloP requires either --neutralModel (skip training) or '
                               '--geneAnnotation (train first).')
        if options.neutralModel and options.geneAnnotation:
            logger.warning('--geneAnnotation is ignored when --neutralModel is given '
                           '(no training; the model is reused as-is)')
            options.geneAnnotation = None

    if options.chunkSize <= 0:
        raise RuntimeError('--chunkSize must be positive (got {})'.format(options.chunkSize))
    if options.fitChunkGroup < 1:
        raise RuntimeError('--fitChunkGroup must be >= 1 (got {})'.format(options.fitChunkGroup))
    if options.bedRanges and not options.bedRanges.endswith('.bed'):
        raise RuntimeError('file passed to --bedRanges must end with .bed')

    # Context-dependent substitution models (U2/U2S/R2/U3/...) need adjacent
    # column context, but our 4d-site extraction collapses each per-chunk SS
    # to --tuple-size 1 before aggregating, which destroys that context. The
    # resulting model would be silently wrong. Guard.
    context_dependent_subst_mods = {'U2', 'U2S', 'R2', 'R2S', 'U3', 'U3S', 'R3', 'R3S'}
    if options.substMod in context_dependent_subst_mods:
        raise RuntimeError(
            '--substMod {} is context-dependent (needs adjacent-column tuples) and is not '
            'compatible with the per-chunk --tuple-size 1 SS aggregation used here. Use a '
            'single-column model: REV (default), SSREV, JC69, F81, HKY85, etc.'.format(
                options.substMod))

    # Validate the annotation extension up front so unrecognized formats
    # (e.g. .gtf, .bed12) don't silently fall through to the genePred branch
    # in localize_annotation, awk on the wrong column, and emit zero 4d sites.
    if options.geneAnnotation:
        ann_low = options.geneAnnotation.lower()
        recognized = ('.gff', '.gff.gz', '.gff3', '.gff3.gz',
                      '.gtf', '.gtf.gz',
                      '.gp', '.gp.gz', '.genepred', '.genepred.gz',
                      '.txt', '.txt.gz')
        if not any(ann_low.endswith(s) for s in recognized):
            raise RuntimeError(
                '--geneAnnotation extension not recognized. Expected one of: '
                '{}. Got: {}'.format(', '.join(recognized), options.geneAnnotation))

    if options.subtree and options.mode == 'phyloFit':
        raise RuntimeError('--subtree is a phyloP-only option (it controls the LRT partition at '
                           'scoring time); it has no effect in --mode phyloFit.')
    if options.subtree:
        # nargs='+' guarantees at least one element; reject duplicates so users
        # don't accidentally request the same track twice.
        seen = set()
        dups = []
        for s in options.subtree:
            if s in seen:
                dups.append(s)
            else:
                seen.add(s)
        if dups:
            raise RuntimeError('--subtree values must be unique; duplicates: {}'.format(
                sorted(set(dups))))

    # If the HAL is a local file, parse the tree once now and reuse it for
    # every --root / --subtree validator + the HAL-root output-naming check
    # below. Remote HAL URLs are deferred to phast_setup (which validates
    # again as defense-in-depth).
    cli_mc_tree = None
    if (options.root or options.subtree) and os.path.isfile(options.halFile):
        cli_mc_tree, _ = _parse_hal_tree(options.halFile)
    if options.root and cli_mc_tree is not None:
        _validate_root_in_hal(options.halFile, options.refGenome, options.root,
                              mc_tree=cli_mc_tree)
    if options.subtree and cli_mc_tree is not None:
        for sub in options.subtree:
            _validate_subtree_in_hal(options.halFile, options.refGenome, sub, options.root,
                                     mc_tree=cli_mc_tree)

    # output sanity
    if options.mode == 'phyloFit':
        if not options.output.endswith('.mod'):
            logger.warning('--mode phyloFit output should end with .mod (got {})'.format(options.output))
    else:
        if not (options.output.endswith('.wig') or options.output.endswith('.wig.gz')):
            logger.warning('--mode phyloP output should end with .wig or .wig.gz (got {})'.format(options.output))

    # Encode --root into the output basename so different runs of the same MAF
    # don't clobber each other. The .s<subtree> token is added per-track at
    # export time (a single invocation may produce N tracks via --subtree).
    # Skipped if the user already included .r<root> in the output, or if
    # --root names the HAL tree's actual root (i.e. the whole tree —
    # equivalent to no restriction). The HAL-root check only runs when the
    # HAL is a local file; for remote URLs the redundant .r tag is accepted.
    name_root = options.root
    if name_root and cli_mc_tree is not None:
        # cli_mc_tree was parsed once above for validation; reuse it. (When
        # name_root is set, cli_mc_tree is always populated by the parse-once
        # block above, since options.root is part of that block's predicate.)
        if name_root == cli_mc_tree.getName(cli_mc_tree.rootId):
            logger.info('--root {!r} is the HAL tree root (whole tree); not '
                        'adding .r tag to output name.'.format(name_root))
            name_root = None
    options.output = _apply_root_subtree_suffix(options.output, name_root, None)

    # absolute paths for any local outputs
    if not options.output.startswith('s3://'):
        options.output = os.path.abspath(options.output)
        # make sure the output's parent directory exists — Toil's
        # exportFile to a file:// URL won't create parents, and a missing
        # dir would only fail at the very end of a possibly-multi-hour
        # run on the per-track export job. Mirrors cactus-hal2chains.
        out_dir = os.path.dirname(options.output)
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir)

    # apply cpu overrides (same idiom as cactus-hal2maf)
    is_singlemachine = options.batchSystem.lower() in ['single_machine', 'singlemachine']
    if options.chunkCores is None:
        if is_singlemachine:
            options.chunkCores = cactus_cpu_count()
            if options.maxCores:
                options.chunkCores = min(options.chunkCores, int(options.maxCores))
            logger.info('Setting chunkCores to {}'.format(options.chunkCores))
        else:
            raise RuntimeError('--chunkCores must be specified for batch systems other than singleMachine')
    # phyloFit only runs when training is needed (i.e. no --neutralModel was
    # supplied for phyloP, or --mode phyloFit). Don't demand --phyloFitCores
    # in score-only runs.
    will_train = (options.mode == 'phyloFit') or (options.mode == 'phyloP' and not options.neutralModel)
    if will_train and options.phyloFitCores is None:
        if is_singlemachine:
            options.phyloFitCores = cactus_cpu_count()
            if options.maxCores:
                options.phyloFitCores = min(options.phyloFitCores, int(options.maxCores))
            logger.info('Setting phyloFitCores to {}'.format(options.phyloFitCores))
        else:
            raise RuntimeError('--phyloFitCores must be specified for batch systems other than singleMachine')

    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        if options.restart:
            toil.restart()
        else:
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)

            maf_id = toil.importFile(makeURL(options.inMaf))
            hal_id = toil.importFile(makeURL(options.halFile))

            # Try to import a .tai. --mafIndex overrides the default sibling path.
            # If the import fails and the user didn't ask for a specific one, the
            # workflow builds it as a pre-pass and exports it next to the output.
            # We narrow the auto-discover catch to "file simply isn't there" — any
            # other failure (S3 auth, network, permission denied, etc.) re-raises
            # so we don't silently fall back to a multi-hour rebuild on what was
            # actually a transient or configurable error.
            tai_path = options.mafIndex or (options.inMaf + '.tai')
            tai_id = None
            tai_built = False
            try:
                tai_id = toil.importFile(makeURL(tai_path))
            except FileNotFoundError as e:
                if options.mafIndex:
                    raise
                logger.warning('No taffy index found at {} ({}). Will build one as a sequential '
                               'pre-pass and export it next to the output for re-use.'
                               .format(tai_path, type(e).__name__))
                tai_built = True

            ann_id = None
            if options.geneAnnotation:
                ann_id = toil.importFile(makeURL(options.geneAnnotation))
            model_id = None
            if options.neutralModel:
                model_id = toil.importFile(makeURL(options.neutralModel))
            bed_id = None
            if options.bedRanges:
                bed_id = toil.importFile(makeURL(options.bedRanges))

            toil.start(Job.wrapJobFn(phast_workflow, config, options,
                                     maf_id, tai_id, tai_built, hal_id, ann_id, model_id,
                                     bed_id))

    end_time = timeit.default_timer()
    logger.info("cactus-phast finished in {} seconds".format(end_time - start_time))


def export_file(job, file_id, out_path):
    """ Export a file ID to its final URL. """
    job.fileStore.exportFile(file_id, makeURL(out_path))
    return True


def export_name_map_job(job, options, species_list):
    """ Write a sed script mapping phast canonical names (everything before
    the first '.' in each HAL genome name) back to the full HAL genome name
    and export it alongside the other outputs. Skipped silently if every
    HAL genome is already a single-token name (no rename happens, so no
    map is needed). The user can apply it with
    `sed -f <out>.name_map.sed phyloFit.mod` to restore HAL names in the
    trained tree. """
    pairs = [(phast_name(g), g) for g in species_list if phast_name(g) != g]
    if not pairs:
        RealtimeLogger.info('No HAL genome names contain ".", so phast names already '
                            'match HAL names; skipping name-map export.')
        return
    work_dir = job.fileStore.getLocalTempDir()
    out_path = os.path.join(work_dir, 'name_map.sed')
    with open(out_path, 'w') as f:
        f.write('# Map phast canonical names -> HAL genome names. Apply with:\n')
        f.write('#   sed -f name_map.sed <phast-output>\n')
        for phast, original in pairs:
            f.write('s/\\b{}\\b/{}/g\n'.format(phast, original))
    job.fileStore.exportFile(job.fileStore.writeGlobalFile(out_path),
                             makeURL(name_map_export_path(options)))
    RealtimeLogger.info('Exported name-map sed script ({} entries) to {}'.format(
        len(pairs), name_map_export_path(options)))


# -----------------------------------------------------------------------------
# Output naming helpers
# -----------------------------------------------------------------------------

def _apply_root_subtree_suffix(output, root, subtree):
    """ Inject .r<root> and/or .s<subtree> tokens into the output filename
    just before the final extension (.wig.gz / .wig / .mod), so different
    runs of the same MAF with different --root / --subtree settings don't
    overwrite each other. Tokens already present in the user-supplied path
    (e.g. they typed `.rAnc0` themselves) are not duplicated. """
    if not (root or subtree):
        return output

    ext = ''
    stem = output
    for e in ('.wig.gz', '.wig', '.mod'):
        if output.endswith(e):
            ext = e
            stem = output[:-len(e)]
            break

    stem_tokens = set(stem.split('.'))
    additions = []
    if root and ('r' + root) not in stem_tokens:
        additions.append('r' + root)
    if subtree and ('s' + subtree) not in stem_tokens:
        additions.append('s' + subtree)

    if not additions:
        return output

    new_output = stem + '.' + '.'.join(additions) + ext
    logger.info('Output renamed to encode --root / --subtree: {} -> {}'.format(
        output, new_output))
    return new_output


def output_basename(options):
    """ Strip .wig.gz / .wig / .mod from the output path to get the shared base
    used for sibling exports (model file, .4d.ss, .bw, etc.). """
    out = options.output
    for ext in ('.wig.gz', '.wig', '.mod'):
        if out.endswith(ext):
            return out[:-len(ext)]
    return out


def model_export_path(options):
    if options.mode == 'phyloFit':
        # phyloFit-mode output IS the model; user picked the path.
        return options.output
    # Model is full-tree (or --root clade); shared across all phyloP tracks
    # so no .s token. options.output already carries .r<root> if applicable.
    return output_basename(options) + '.mod'


def ss_export_path(options):
    # 4d.ss is computed pre-scoring and shared across all phyloP tracks.
    return output_basename(options) + '.4d.ss'


def wig_export_path(options, subtree):
    """ Per-track wig path. `subtree` is the per-track --subtree value, or
    None for the default global track. options.output is the global-track
    path (already carrying .r<root> if --root was set); for an LRT track we
    append .s<subtree> just before the .wig.gz / .wig extension. """
    if not subtree:
        return options.output
    return _apply_root_subtree_suffix(options.output, None, subtree)


def bigwig_export_path(options, subtree):
    """ Per-track bigwig path: same .s<subtree> handling as the wig. """
    wig = wig_export_path(options, subtree)
    for ext in ('.wig.gz', '.wig'):
        if wig.endswith(ext):
            return wig[:-len(ext)] + '.bw'
    return wig + '.bw'


def tai_export_path(options):
    """ The .tai is exported next to <output> with the input MAF's basename. """
    out_dir = os.path.dirname(options.output) or '.'
    return os.path.join(out_dir, os.path.basename(options.inMaf) + '.tai')


def name_map_export_path(options):
    """ Convenience sed script mapping phast canonical names back to the
    original HAL genome names. Exported alongside the .mod / .4d.ss. """
    return output_basename(options) + '.name_map.sed'


# -----------------------------------------------------------------------------
# Helpers for species-name normalization and annotation chrom rewriting
# -----------------------------------------------------------------------------

def phast_name(genome):
    """ phast's msa_view internally splits species source names at the first
    '.' to separate species from chrom. So `GCA_028858775.2` becomes
    `GCA_028858775` in the SS / model files. We must use this canonical form
    everywhere phast needs to match species names. """
    return genome.split('.', 1)[0]


def make_strip_perl_cmd():
    """ Return a perl one-liner that rewrites MAF s-lines from
    `s\\t<genome>.<contig>\\t...` to `s\\t<phast_name(genome)>\\t...`. The
    rewritten MAF lines up with downstream phast tools (msa_view, phyloFit,
    phyloP), which all key on the post-first-dot name.

    Equivalent to a per-species substitution list given the collision check
    at workflow start (which guarantees no two HAL genomes share a
    phast_name), but O(1) per s-line vs O(N_species). At 400-577 species
    the explicit list was ~30 KB of regex per s-line — measurable wall-time
    cost during stream rewriting. s-lines whose species name has no dot
    (bare names like 'hg38' or 'Anc0') are unmatched and pass through
    unchanged, which is correct since phast_name is a no-op there. """
    return r"""perl -wpe 's/^(s\s+)([^.]+)\.\S+/$1$2/'"""


def newick_leaves(tree_str):
    """ Return the leaf labels of a Newick tree string in the order they
    appear. A leaf is a name immediately preceded by '(' or ',' and followed
    by ':' / ',' / ')' (no inner '(...)' between). Ancestor labels are the
    names that immediately follow ')' and are excluded. """
    leaves = []
    s = tree_str
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '(' or c == ',':
            # leaf or subtree starts at i+1
            j = i + 1
            # skip whitespace
            while j < n and s[j] in ' \t\n':
                j += 1
            if j < n and s[j] != '(':
                # leaf label
                k = j
                while k < n and s[k] not in ':,();':
                    k += 1
                label = s[j:k].strip()
                if label:
                    leaves.append(label)
                i = k
                continue
        i += 1
    return leaves


def rewrite_tree_to_phast_names(tree_str, species_list):
    """ Rewrite a Newick tree's leaf labels so dotted genome names become their
    phast-canonical form (split at first dot). Names without dots are
    untouched. """
    out = tree_str
    sorted_species = sorted({g for g in species_list if '.' in g}, key=len, reverse=True)
    for g in sorted_species:
        out = re.sub(r'(?<![\w.])' + re.escape(g) + r'(?![\w.])', phast_name(g), out)
    return out


def _parse_hal_tree(hal_path):
    """ Parse the HAL tree into a MultiCactusTree. Returns the parsed tree. """
    from sonLib.nxnewick import NXNewick
    from cactus.progressive.multiCactusTree import MultiCactusTree
    tree_str = cactus_call(parameters=['halStats', hal_path, '--tree'],
                           check_output=True).strip()
    return MultiCactusTree(NXNewick().parseString(tree_str, addImpliedRoots=False)), tree_str


def _hal_root_name(hal_path):
    """ Return the name of the root node of the HAL tree. """
    mc_tree, _ = _parse_hal_tree(hal_path)
    return mc_tree.getName(mc_tree.rootId)


def _root_descendant_leaves(mc_tree, root, ref_genome):
    """ Validate that `root` is an internal-node name in the HAL tree and
    that `ref_genome` is one of its leaf descendants. Returns the list of
    leaf descendants (including ref_genome). Raises RuntimeError on any
    failure with a clear, immediate message. """
    if root not in mc_tree.nameToId:
        raise RuntimeError(
            "--root {!r} is not a node name in the HAL tree. Available "
            "internal-node names: {}".format(
                root,
                sorted(n for n, i in mc_tree.nameToId.items() if not mc_tree.isLeaf(i))))
    root_id = mc_tree.nameToId[root]
    if mc_tree.isLeaf(root_id):
        raise RuntimeError(
            "--root {!r} is a leaf in the HAL tree, not an internal node; "
            "it has no descendant species.".format(root))
    descendant_names = mc_tree.getChildNames(root)  # leaves + internal
    descendant_leaves = [n for n in descendant_names
                         if mc_tree.isLeaf(mc_tree.nameToId[n])]
    if ref_genome not in descendant_leaves:
        # Truncate the leaf dump for trees with hundreds of leaves so the error
        # message stays readable.
        preview = descendant_leaves[:20]
        more = len(descendant_leaves) - len(preview)
        raise RuntimeError(
            "Reference genome {!r} is not a leaf descendant of --root {!r}. "
            "Descendant leaves: {}{}".format(
                ref_genome, root, preview,
                ' (...and {} more)'.format(more) if more > 0 else ''))
    return descendant_leaves


def _validate_root_in_hal(hal_path, ref_genome, root, mc_tree=None):
    """ Quick sanity check used in main() when the HAL is a local file: parse
    the tree and verify the --root name is valid and is an ancestor of the
    reference. Defense-in-depth — the full workflow does this again in
    phast_setup, but doing it here lets the user discover misnamed nodes at
    invocation rather than after Toil job-store init.

    `mc_tree` may be passed in by callers that have already parsed the tree
    (e.g. main() validating multiple --subtree values), to avoid a redundant
    halStats --tree spawn per validation. """
    if mc_tree is None:
        mc_tree, _ = _parse_hal_tree(hal_path)
    _root_descendant_leaves(mc_tree, root, ref_genome)


def _subtree_leaves(mc_tree, subtree, root=None):
    """ Validate that `subtree` is an internal-node name in the HAL tree.
    If `root` is given, also require that `subtree` is a strict descendant
    of `root` (the pipeline first prunes to `root`'s clade, and the LRT
    partitions that pruned tree at `subtree`). Returns the leaf descendants
    of `subtree`. """
    if subtree not in mc_tree.nameToId:
        raise RuntimeError(
            "--subtree {!r} is not a node name in the HAL tree. Available "
            "internal-node names: {}".format(
                subtree,
                sorted(n for n, i in mc_tree.nameToId.items() if not mc_tree.isLeaf(i))))
    sub_id = mc_tree.nameToId[subtree]
    if mc_tree.isLeaf(sub_id):
        raise RuntimeError(
            "--subtree {!r} is a leaf in the HAL tree, not an internal node; "
            "it has no descendant species.".format(subtree))
    if root is not None:
        if subtree == root:
            raise RuntimeError(
                "--subtree {!r} is the same node as --root; the LRT then has "
                "no 'rest of the tree' to compare against. Pick a strict descendant "
                "of --root, or drop --root.".format(subtree))
        root_descendants = set(mc_tree.getChildNames(root))
        if subtree not in root_descendants:
            raise RuntimeError(
                "--subtree {!r} is not a descendant of --root {!r}. When both "
                "options are set, the pipeline first prunes to --root's clade "
                "and then partitions that clade at --subtree, so --subtree must be "
                "a strict internal-node descendant of --root.".format(
                    subtree, root))
    descendant_names = mc_tree.getChildNames(subtree)
    return [n for n in descendant_names if mc_tree.isLeaf(mc_tree.nameToId[n])]


def _validate_subtree_in_hal(hal_path, ref_genome, subtree, root=None, mc_tree=None):
    """ Quick sanity check used in main() when the HAL is a local file.
    Validates --subtree, its containment relationship with --root (when set),
    and that the reference is a leaf of --subtree (so the per-base wig
    anchored at the reference has subtree-side data at every column).

    Special case: if `subtree` equals the effective tree root (the HAL root,
    or `root` if set), this is the standard global-conservation track on
    that (sub)tree — no phyloP --subtree partition. Only the existence /
    internal-node checks apply.

    `mc_tree` may be passed in by callers that have already parsed the tree
    (e.g. main() validating multiple --subtree values), to avoid a redundant
    halStats --tree spawn per call. """
    if mc_tree is None:
        mc_tree, _ = _parse_hal_tree(hal_path)
    if root is not None:
        # Always validate --root first so the user sees consistent error
        # ordering (root → subtree → ref).
        _root_descendant_leaves(mc_tree, root, ref_genome)

    effective_root = root if root else mc_tree.getName(mc_tree.rootId)
    if subtree == effective_root:
        # Global track on the (sub)tree: only verify the name exists and is internal.
        if subtree not in mc_tree.nameToId:
            raise RuntimeError("--subtree {!r} is not a node name in the HAL tree.".format(subtree))
        if mc_tree.isLeaf(mc_tree.nameToId[subtree]):
            raise RuntimeError("--subtree {!r} is a leaf, not an internal node.".format(subtree))
        return

    leaves = _subtree_leaves(mc_tree, subtree, root)
    if ref_genome not in leaves:
        raise RuntimeError(
            "Reference genome {!r} must be a leaf descendant of --subtree {!r}. "
            "phyloP --subtree partitions the tree at the named node and produces "
            "a per-base score anchored at the reference; for that score to have "
            "subtree-side data at each reference position, the reference must "
            "lie within the subtree. Pick a subtree that contains {!r}.".format(
                ref_genome, subtree, ref_genome))
    if len(leaves) == 1:
        raise RuntimeError(
            "--subtree {!r}'s only leaf descendant is the reference {!r}; the "
            "LRT cannot estimate a subtree rate from a single species. Pick a "
            "larger subtree that contains the reference plus at least one "
            "other species.".format(subtree, ref_genome))


def localize_annotation(work_dir, ann_path, contig, ref_genome, idx):
    """ Filter the prepped genePred to records on `contig` and rewrite the
    chrom column from `<contig>` to `ref_genome` so it matches the stripped
    MAF reference name `<refGenome>` (rather than the original
    `<refGenome>.<contig>`). The input is always the single-cover genePred
    produced by prep_annotation; column layout is name=0, chrom=1, strand=2,
    txStart=3, txEnd=4, cdsStart=5, cdsEnd=6 (no leading bin). """
    out_path = os.path.join(work_dir, 'ann_batch_{}.local.gp'.format(idx))
    chrom_col = 1
    with open(ann_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            cols = line.rstrip('\n').split('\t')
            if len(cols) > chrom_col and cols[chrom_col] == contig:
                cols[chrom_col] = ref_genome
                fout.write('\t'.join(cols) + '\n')
    return out_path


# -----------------------------------------------------------------------------
# Top-level workflow
# -----------------------------------------------------------------------------

def phast_workflow(job, config, options, maf_id, tai_id, tai_built, hal_id, ann_id, model_id,
                   bed_id=None):
    """ Orchestrate setup -> index -> plan -> shard -> (phyloFit?) -> (phyloP?) """

    # fail fast if any required external binary is missing rather than crashing
    # mid-workflow on a remote worker after hours of localization
    check_job = job.addChildJobFn(phast_check_tools, options)

    # Normalize the annotation up front: convert to genePred, filter to
    # CDS-bearing transcripts, single-cover. Only matters when we'll train.
    # Must sit on the linear chain leading to train_workflow (which consumes
    # prep_job.rv()) — Toil rejects RV promises across non-descendant jobs.
    need_train = (options.mode == 'phyloFit') or (options.mode == 'phyloP' and model_id is None)
    if need_train and ann_id is not None:
        # Annotations are tiny relative to MAF/HAL; default memory is plenty
        # for a human-scale ncbiRefSeq.txt.gz (~10 MB compressed).
        prep_job = check_job.addFollowOnJobFn(prep_annotation, options, ann_id,
                                              disk=max(2 * 1024**3, int(ann_id.size * 20)))
        ann_id = prep_job.rv()
        setup_parent = prep_job
    else:
        setup_parent = check_job

    setup_job = setup_parent.addFollowOnJobFn(phast_setup, options, hal_id,
                                              disk=int(hal_id.size * 1.1))
    species_list = setup_job.rv(0)
    tree_str = setup_job.rv(1)
    ref_seq_lengths = setup_job.rv(2)
    effective_root_name = setup_job.rv(3)
    aligned_contigs = setup_job.rv(4)

    # Convenience sed script: maps phast canonical names back to HAL genome
    # names. Only relevant when at least one HAL genome name contains '.';
    # the job no-ops otherwise.
    setup_job.addFollowOnJobFn(export_name_map_job, options, species_list)

    # build .tai if missing, else just chain through the existing one
    if tai_built:
        idx_job = setup_job.addFollowOnJobFn(taffy_index_job, maf_id, os.path.basename(options.inMaf),
                                             disk=int(maf_id.size * 1.1))
        tai_id = idx_job.rv()
        idx_job.addFollowOnJobFn(export_file, tai_id, tai_export_path(options))
        plan_parent = idx_job
    else:
        plan_parent = setup_job

    plan_job = plan_parent.addFollowOnJobFn(plan_chunks_job, options, ref_seq_lengths, tai_id,
                                            bed_id, aligned_contigs)
    chunk_specs = plan_job.rv()

    # The single multi-core chunker job: localizes the source MAF once,
    # parallelizes per-chunk taffy view via GNU parallel, returns a list of
    # (contig, start, end, chunk_id). All downstream processing is single-
    # process per chunk (or per chunk-group for 4d extraction).
    # Disk: source MAF (~S) + chunks accumulated as parallel runs (~S total
    # across all chunks, since they cover the same data). The source is
    # removed before the writeGlobalFile loop and chunks are removed as
    # they're uploaded, so peak ≈ 2S at the moment parallel finishes. Add
    # 4 GiB slack for taffy scratch and Toil worker overhead.
    chunker_disk = int(maf_id.size * 2.0) + 4 * 1024**3
    # taffy view is small (~50 MiB peak at 32 cores on a 577-way 300 GB MAF);
    # default memory is fine. Splice strip + mafDuplicateFilter into the
    # chunker pipeline so each chunk lands pre-filtered for downstream phast
    # tools (msa_view --4d and phyloP both need the same filtered MAF, no
    # point doing it twice per chunk).
    filter_cmd = '{strip} | mafDuplicateFilter -k -m -'.format(strip=make_strip_perl_cmd())
    chunk_job = plan_job.addFollowOnJobFn(chunker_job, options, maf_id, tai_id,
                                          os.path.basename(options.inMaf), chunk_specs,
                                          filter_cmd,
                                          disk=chunker_disk,
                                          cores=options.chunkCores)
    chunks = chunk_job.rv()

    # ----- phyloFit branch (also runs in phyloP mode if no model was given) -----
    if need_train:
        train_job = chunk_job.addFollowOnJobFn(train_workflow, options, chunks,
                                               species_list, tree_str, ann_id)
        # train_workflow returns the model file id
        trained_model_id = train_job.rv()
        if options.mode == 'phyloFit':
            # done after training
            return
        else:
            model_id = trained_model_id
            score_parent = train_job
    else:
        score_parent = chunk_job

    # ----- phyloP branch (one track per --subtree value, or one default) -----
    # Default (no --subtree): one global track, signaled by passing None down.
    # When --subtree X is given, X==effective_root_name resolves at runtime in
    # phyloP_workflow as "no phyloP --subtree flag" (the standard global track
    # on the (sub)tree). Each track has its own merge/wigToBigWig/export
    # follow-on chain so failures and outputs are independent.
    track_inputs = options.subtree if options.subtree else [None]
    track_rvs = []
    for sub in track_inputs:
        score_job = score_parent.addFollowOnJobFn(phyloP_workflow, options, chunks,
                                                  model_id, ref_seq_lengths, species_list,
                                                  sub, effective_root_name)
        track_rvs.append(score_job.rv())
    return track_rvs


def phast_check_tools(job, options):
    """ Verify every required external binary is callable so a missing tool
    fails the workflow early with an informative error, rather than crashing
    mid-workflow on a remote worker after the source MAF + HAL have been
    localized. Mirrors hal2chains_check_tools' approach: invoke each tool
    via cactus_call (which routes through the configured binariesMode, so
    docker / singularity / local all probe the right place) and treat
    FileNotFoundError or "exited 127" as "missing". A non-zero exit with a
    usage banner (the typical no-args response) means the tool is present. """
    install_hints = {
        'phyloFit': 'phast (build-tools/downloadPhast)',
        'phyloP':   'phast (build-tools/downloadPhast)',
        'msa_view': 'phast (build-tools/downloadPhast)',
        'modFreqs': 'phast (build-tools/downloadPhast)',
        'taffy':    'taffy (submodules/taffy or build-tools)',
        'halStats': 'hal (submodules/hal)',
        'mafDuplicateFilter': 'maftools (build-tools/downloadMafTools)',
        'wigToBigWig': 'UCSC tools '
                       '(https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig)',
        'genePredSingleCover': 'UCSC tools '
                       '(https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredSingleCover)',
        'gff3ToGenePred': 'UCSC tools '
                       '(https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred)',
        'gtfToGenePred': 'UCSC tools '
                       '(https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred)',
        'parallel': 'GNU parallel (apt-get install parallel / brew install parallel)',
        'bgzip':    'htslib (apt-get install tabix / brew install htslib)',
    }

    required = ['phyloFit', 'phyloP', 'msa_view', 'taffy', 'halStats',
                'mafDuplicateFilter', 'parallel', 'bgzip']
    if options.modFreqs:
        required.append('modFreqs')
    if options.bigwig:
        required.append('wigToBigWig')
    if options.geneAnnotation:
        # genePredSingleCover is always run on the prepped annotation; the
        # GFF/GTF converters are only needed if the user supplied one of
        # those formats.
        required.append('genePredSingleCover')
        ann_low = options.geneAnnotation.lower()
        if ann_low.endswith('.gff') or ann_low.endswith('.gff.gz') \
                or ann_low.endswith('.gff3') or ann_low.endswith('.gff3.gz'):
            required.append('gff3ToGenePred')
        if ann_low.endswith('.gtf') or ann_low.endswith('.gtf.gz'):
            required.append('gtfToGenePred')

    missing = []
    for tool in required:
        try:
            cactus_call(parameters=[tool])
        except FileNotFoundError:
            missing.append(tool)
        except RuntimeError as e:
            # 127 is the standard exit code for "command not found" both from
            # /bin/sh and from `docker run`. Anything else means the tool
            # actually ran (printing usage / an error and exiting non-zero),
            # which is fine — it's present.
            if 'exited 127' in str(e):
                missing.append(tool)

    if missing:
        lines = ['Required tool(s) not found:']
        for t in missing:
            lines.append('  {} — install: {}'.format(
                t, install_hints.get(t, '(see cactus-phast docs)')))
        raise RuntimeError('\n'.join(lines))


def phast_setup(job, options, hal_id):
    """ Read tree + species list + ref chrom sizes from HAL. """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))
    job.fileStore.readGlobalFile(hal_id, hal_path)

    ref_seq_lengths = get_ref_sequence_lengths(hal_path, options.refGenome)
    if not ref_seq_lengths:
        raise RuntimeError('No sequences found for refGenome {} in {}'.format(
            options.refGenome, options.halFile))

    species_list = cactus_call(parameters=['halStats', '--genomes', hal_path],
                               check_output=True).strip().split()
    mc_tree, tree_str = _parse_hal_tree(hal_path)

    # If the user picked a subtree, validate (again — main() does this too when
    # the HAL is local) and restrict species_list + tree_str to that subtree.
    # Everything downstream (msa_view --aggregate, phyloFit, phyloP) keys on
    # these two values, so the rest of the pipeline transparently scopes to the
    # clade. phyloP simply ignores MAF rows whose species isn't in the model.
    if options.root:
        descendant_leaves = _root_descendant_leaves(mc_tree, options.root,
                                                        options.refGenome)
        if options.root == mc_tree.getName(mc_tree.rootId):
            # --root names the HAL tree root: identity prune. Skip the
            # extractSubTree / species_list filter / Newick rewrite — every
            # genome is already in the (sub)tree by definition.
            RealtimeLogger.info('--root {!r} is the HAL tree root; using full tree.'.format(
                options.root))
        else:
            from sonLib.nxnewick import NXNewick
            descendants = set([options.root] + mc_tree.getChildNames(options.root))
            species_list = [g for g in species_list if g in descendants]
            sub_tree = mc_tree.extractSubTree(options.root)
            tree_str = NXNewick().writeString(sub_tree)
            RealtimeLogger.info('Restricted to clade rooted at {!r}: {} leaves '
                                '({} HAL genomes total in subtree)'.format(
                                    options.root, len(descendant_leaves), len(species_list)))

    # Defense-in-depth: --subtree validation may not have run in main() if the
    # HAL was a remote URL. Check each subtree value here too. The containment
    # check operates on the FULL HAL tree (mc_tree) which is still the
    # original parse, so the root/subtree relationship is checked correctly
    # regardless of whether the --root block above subset tree_str.
    if options.subtree:
        effective_root = options.root if options.root else mc_tree.getName(mc_tree.rootId)
        for sub in options.subtree:
            if sub == effective_root:
                # Global track within (sub)tree — only existence + internal-node check.
                if sub not in mc_tree.nameToId:
                    raise RuntimeError("--subtree {!r} not in HAL tree.".format(sub))
                if mc_tree.isLeaf(mc_tree.nameToId[sub]):
                    raise RuntimeError("--subtree {!r} is a leaf, not internal.".format(sub))
                RealtimeLogger.info('Track {!r}: global conservation on the (sub)tree '
                                    '(no phyloP --subtree flag).'.format(sub))
                continue
            sub_leaves = _subtree_leaves(mc_tree, sub, options.root)
            if options.refGenome not in sub_leaves:
                raise RuntimeError(
                    "Reference genome {!r} must be a leaf descendant of --subtree {!r}. Pick "
                    "a subtree that contains the reference.".format(options.refGenome, sub))
            if len(sub_leaves) == 1:
                raise RuntimeError(
                    "--subtree {!r}'s only leaf descendant is the reference {!r}; the LRT "
                    "needs at least one other species in the subtree.".format(
                        sub, options.refGenome))
            RealtimeLogger.info('Track {!r}: lineage-specific LRT over {} leaves{}.'.format(
                sub, len(sub_leaves),
                ' (within --root clade)' if options.root else ''))

    # phast splits source names at the first '.' to derive the species label.
    # If two HAL genomes collapse to the same phast-canonical name (e.g.
    # 'hg38.1' and 'hg38.2' both -> 'hg38'), the SS file, trained model, and
    # tree would silently conflate them. Refuse to run.
    collisions = {}
    for g in species_list:
        collisions.setdefault(phast_name(g), []).append(g)
    bad = {k: v for k, v in collisions.items() if len(v) > 1}
    if bad:
        msg = ['Genome name collisions detected: phast splits source names at the first "."',
               'so the following HAL genomes would all be treated as one species:']
        for canonical, originals in sorted(bad.items()):
            msg.append('  {!r}  <-  {}'.format(canonical, ', '.join(sorted(originals))))
        msg.append('Rename the offending genomes in the HAL (e.g. via halRenameGenomes) so they '
                   'differ before the first ".", then re-run.')
        raise RuntimeError('\n'.join(msg))

    ref_contigs_sorted = sorted(ref_seq_lengths.keys())
    RealtimeLogger.info('Found {} reference sequences for {} ({}), total {} bp; {} genomes in HAL.'.format(
        len(ref_seq_lengths), options.refGenome,
        ', '.join(ref_contigs_sorted[:10]) + ('...' if len(ref_contigs_sorted) > 10 else ''),
        sum(ref_seq_lengths.values()), len(species_list)))

    # The effective root name — the (sub)tree root for downstream "is this
    # a global track?" comparisons. If --root is set, that's what we pruned
    # to; otherwise it's the HAL tree root.
    effective_root_name = options.root if options.root else mc_tree.getName(mc_tree.rootId)

    # While the HAL is already localized in this job (it is the only job that
    # reads it — copying it again on a cluster would be expensive), find which
    # reference contigs are actually aligned to anything. Fully-unaligned contigs
    # only ever yield reference-only MAF columns that phyloP can't score, so we
    # drop them at planning time unless --keepUnalignedContigs. aligned_contigs
    # stays None (meaning "don't filter") whenever the check is skipped or fails.
    #
    # halAlignedExtract reports only alignment to the PARENT (top segments with
    # hasParent()). That equals "appears in a scoreable multi-row MAF column" only
    # for a LEAF reference. An internal/ancestral reference can align to its
    # descendants (via its bottom segments / the column iterator's parse-down),
    # which halAlignedExtract does not see, so the filter would wrongly drop
    # scoreable contigs there. So only run it for a leaf reference (this also
    # covers the HAL root, which is internal and additionally has no parent).
    aligned_contigs = None
    if not options.keepUnalignedContigs:
        ref_node = mc_tree.nameToId.get(options.refGenome)
        if ref_node is None or not mc_tree.isLeaf(ref_node):
            RealtimeLogger.info('Reference {!r} is not a leaf genome (internal/ancestral, or not a '
                                'named node in the HAL tree); skipping the unaligned-contig check '
                                '(halAlignedExtract only reports alignment to the parent, not to '
                                'descendants). Pass --bedRanges to restrict manually if '
                                'needed.'.format(options.refGenome))
        else:
            try:
                aligned = get_aligned_ref_contigs(hal_path, options.refGenome)
            except RuntimeError as e:
                # Don't let an auxiliary optimization break an otherwise-fine run:
                # degrade to processing every contig (old behavior).
                aligned = None
                RealtimeLogger.warning('halAlignedExtract failed ({}); processing all reference '
                                       'contigs without the unaligned-contig filter.'.format(e))
            if aligned is not None and not aligned:
                # Empty result (e.g. a pathological all-unaligned genome, or an
                # unexpected name mismatch) — skip the filter rather than drop
                # every contig and hard-fail.
                RealtimeLogger.warning('halAlignedExtract reported no aligned contigs for {!r}; '
                                       'processing all reference contigs without the '
                                       'unaligned-contig filter.'.format(options.refGenome))
            elif aligned:
                aligned_contigs = aligned
                n_unaligned = len(set(ref_seq_lengths) - aligned_contigs)
                RealtimeLogger.info('halAlignedExtract: {} of {} reference contigs are aligned; '
                                    '{} unaligned will be skipped (use --keepUnalignedContigs to '
                                    'keep them).'.format(
                                        len(aligned_contigs), len(ref_seq_lengths), n_unaligned))

    return species_list, tree_str, ref_seq_lengths, effective_root_name, aligned_contigs


def plan_chunks_job(job, options, ref_seq_lengths, tai_id, bed_id=None, aligned_contigs=None):
    """ Slice the reference into chunks of --chunkSize and filter to contigs
    present in the source .tai. With --bedRanges, only the BED-selected
    reference ranges are sliced (rather than every contig). When `aligned_contigs`
    is given (the default, computed in phast_setup), contigs with no alignment in
    the HAL are dropped. Returns the flat chunk list consumed by chunker_job. """
    work_dir = job.fileStore.getLocalTempDir()
    tai_path = os.path.join(work_dir, 'in.maf.tai')
    job.fileStore.readGlobalFile(tai_id, tai_path)

    if bed_id is not None:
        bed_path = os.path.join(work_dir, 'ranges.bed')
        job.fileStore.readGlobalFile(bed_id, bed_path)
        bed_ranges = parse_bed_ranges(bed_path)
        if not bed_ranges:
            raise RuntimeError('--bedRanges contained no usable intervals (each line needs at '
                               'least 3 columns: sequence start end).')
        chunks_all, dropped, clamped = plan_chunks_in_regions(bed_ranges, ref_seq_lengths,
                                                              options.chunkSize)
        if dropped:
            shown = ', '.join(dropped[:20]) + (' ...' if len(dropped) > 20 else '')
            RealtimeLogger.info('--bedRanges: {} BED sequence(s) are not sequences of reference '
                                '{!r}; skipping: {}'.format(len(dropped), options.refGenome, shown))
        if clamped:
            shown = ', '.join('{}:{}-{}'.format(c, s, e) for (c, s, e) in clamped[:10]) \
                + (' ...' if len(clamped) > 10 else '')
            RealtimeLogger.warning('--bedRanges: {} interval(s) extend past [0, contig length] and '
                                   'were clamped or dropped — check the BED matches the {!r} '
                                   'assembly: {}'.format(len(clamped), options.refGenome, shown))
        if not chunks_all:
            raise RuntimeError(
                '--bedRanges selected no reference sequence: none of its sequences are sequences '
                'of {!r} (or all intervals are empty after clamping to contig lengths). BED '
                'sequence names must match the reference contig names without the genome '
                'prefix.'.format(options.refGenome))
    else:
        chunks_all = plan_chunks(ref_seq_lengths, options.chunkSize)

    # Drop contigs the HAL says are unaligned (phast_setup computed this from the
    # already-localized HAL; None means the check was skipped / --keepUnalignedContigs).
    if aligned_contigs is not None:
        before = len(chunks_all)
        unaligned = sorted({c[0] for c in chunks_all} - aligned_contigs)
        chunks_all = [c for c in chunks_all if c[0] in aligned_contigs]
        if unaligned:
            shown = ', '.join(unaligned[:20]) + (' ...' if len(unaligned) > 20 else '')
            RealtimeLogger.info('Excluding {} unaligned reference contig(s) ({} chunks) with no '
                                'alignment in the HAL: {}'.format(
                                    len(unaligned), before - len(chunks_all), shown))
        if not chunks_all:
            raise RuntimeError(
                'All selected reference contigs are unaligned in the HAL — nothing for phyloP to '
                'score. Use --keepUnalignedContigs to process them anyway.')

    chunks = filter_chunks_to_indexed(chunks_all, tai_path, options.refGenome)
    if not chunks:
        raise RuntimeError(
            'No {} appear in the source MAF .tai. Either the MAF is empty or the reference name '
            'does not match the MAF s-line prefixes.'.format(
                'selected --bedRanges sequences' if bed_id is not None
                else 'contigs from {!r}'.format(options.refGenome)))
    skipped = len(chunks_all) - len(chunks)
    if skipped:
        RealtimeLogger.info('Skipping {} chunks for contigs not present in MAF index '
                            '(unaligned / out-of-MAF reference contigs).'.format(skipped))
    n_contigs = len({c[0] for c in chunks})
    RealtimeLogger.info('Planned {} chunks of ~{} bp each (covering {} contigs){}'.format(
        len(chunks), options.chunkSize, n_contigs,
        ' from --bedRanges' if bed_id is not None else ''))
    return chunks


# -----------------------------------------------------------------------------
# Annotation preprocessing
# -----------------------------------------------------------------------------

def prep_annotation(job, options, ann_id):
    """ Convert the user's gene annotation into a single-cover, non-redundant
    genePred suitable for msa_view --4d.

    Format is detected from the input filename extension:
      .gff / .gff3 (.gz)        : converted via gff3ToGenePred
      .gtf (.gz)                : converted via gtfToGenePred
      .txt (.gz)                : UCSC db dump (e.g. ncbiRefSeq.txt.gz) — leading
                                  'bin' column is stripped via cut -f2-
      .gp / .genePred (.gz)     : used as-is (just decompressed if .gz)

    The genePred is then filtered to records with cdsEnd > cdsStart (drops
    non-coding transcripts and entries with no CDS), run through
    genePredSingleCover to produce one transcript per genomic position
    (required so 4d-site extraction doesn't double-count overlapping
    transcripts), and sorted for determinism. """
    work_dir = job.fileStore.getLocalTempDir()
    ann_basename = os.path.basename(options.geneAnnotation)
    src_path = os.path.join(work_dir, ann_basename)
    job.fileStore.readGlobalFile(ann_id, src_path)

    ann_low = ann_basename.lower()
    is_gz = ann_low.endswith('.gz')
    base = ann_low[:-3] if is_gz else ann_low

    # If the input is gzipped, decompress to a sibling temp file first so we
    # can hand a real path to gff3ToGenePred / gtfToGenePred (no shell pipes).
    if is_gz:
        import gzip
        decompressed = os.path.join(work_dir, 'ann.decompressed')
        with gzip.open(src_path, 'rb') as fin, open(decompressed, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        input_path = decompressed
    else:
        input_path = src_path

    gp_path = os.path.join(work_dir, 'ann.raw.gp')
    if base.endswith('.gff') or base.endswith('.gff3'):
        cactus_call(parameters=['gff3ToGenePred', input_path, gp_path])
    elif base.endswith('.gtf'):
        cactus_call(parameters=['gtfToGenePred', input_path, gp_path])
    elif base.endswith('.txt'):
        # UCSC db dump: leading 'bin' column to strip.
        with open(input_path) as fin, open(gp_path, 'w') as fout:
            for line in fin:
                _, _, rest = line.partition('\t')
                fout.write(rest)
    else:
        # Plain genePred (.gp / .genePred); just copy to the standard name.
        shutil.copyfile(input_path, gp_path)

    if os.path.getsize(gp_path) == 0:
        raise RuntimeError('Annotation conversion of {} produced an empty genePred.'
                           .format(ann_basename))

    # genePred (no leading bin) cols 0-indexed: name=0, chrom=1, strand=2,
    # txStart=3, txEnd=4, cdsStart=5, cdsEnd=6. Drop records with no CDS.
    cdsfilt_path = os.path.join(work_dir, 'ann.cdsfilt.gp')
    with open(gp_path) as fin, open(cdsfilt_path, 'w') as fout:
        for line in fin:
            cols = line.rstrip('\n').split('\t')
            if len(cols) > 6 and int(cols[6]) > int(cols[5]):
                fout.write(line)

    if os.path.getsize(cdsfilt_path) == 0:
        raise RuntimeError(
            'After filtering for cdsEnd > cdsStart, no transcripts remained in '
            '--geneAnnotation {}. Check that the file contains coding transcripts.'
            .format(options.geneAnnotation))

    nr_path = os.path.join(work_dir, 'ann.prepped.gp')
    # LC_ALL=C: byte-stable sort regardless of locale, so the prepped
    # annotation is deterministic across machines.
    cactus_call(parameters=['bash', '-c',
        "set -eo pipefail && genePredSingleCover {inp} /dev/stdout | LC_ALL=C sort > {out}".format(
            inp=shlex.quote(cdsfilt_path), out=shlex.quote(nr_path))])

    if os.path.getsize(nr_path) == 0:
        raise RuntimeError('genePredSingleCover produced an empty annotation from {}.'
                           .format(options.geneAnnotation))

    with open(cdsfilt_path) as f:
        n_in = sum(1 for _ in f)
    chroms = set()
    with open(nr_path) as f:
        n_out = 0
        for line in f:
            n_out += 1
            cols = line.split('\t', 2)
            if len(cols) >= 2:
                chroms.add(cols[1])
    RealtimeLogger.info('Prepped annotation {}: {} cds-filtered transcripts -> '
                        '{} single-cover records on {} contigs ({})'.format(
                            ann_basename, n_in, n_out, len(chroms),
                            ', '.join(sorted(chroms)[:10]) + ('...' if len(chroms) > 10 else '')))
    return job.fileStore.writeGlobalFile(nr_path)


# -----------------------------------------------------------------------------
# phyloFit subgraph
# -----------------------------------------------------------------------------

def train_workflow(job, options, chunks, species_list, tree_str, ann_id):
    """ Per-chunk-group 4d extraction -> aggregate -> phyloFit -> optional
    modFreqs. Returns the trained model file ID. """
    # Per-chunk SS files must be projected to the leaves of tree_str so that
    # msa_view --aggregate (which requires its species CSV to cover every row
    # in every input file) doesn't reject leaves outside the tree. Without
    # --root, this list is just every leaf in the HAL tree (a no-op
    # projection); with --root, it's restricted to clade leaves.
    leaves_csv = ','.join(phast_name(g) for g in newick_leaves(tree_str))
    extract_job = job.addChildJobFn(extract_4d_all, options, chunks, ann_id,
                                    species_list, leaves_csv)
    ss_results = extract_job.rv()  # flat list of ss file ids (one per chunk-group)

    # 4d-site SS files are tiny even at 447-way (~300 MB aggregated). 4 GiB
    # disk covers concat scratch comfortably.
    aggregate_job = extract_job.addFollowOnJobFn(aggregate_4d, options, ss_results,
                                                 species_list, tree_str,
                                                 disk=4 * 1024**3)
    aggregate_id = aggregate_job.rv()
    aggregate_job.addFollowOnJobFn(export_file, aggregate_id, ss_export_path(options))

    # phyloFit on 577-way × 1.4 GB SS peaked at ~300 MiB; default memory is
    # fine unless the user overrides via --phyloFitMemory.
    fit_kwargs = dict(disk=4 * 1024**3, cores=options.phyloFitCores)
    if options.phyloFitMemory:
        fit_kwargs['memory'] = cactus_clamp_memory(options.phyloFitMemory)
    fit_job = aggregate_job.addFollowOnJobFn(phyloFit_job, options, aggregate_id,
                                             tree_str, species_list, **fit_kwargs)
    model_id = fit_job.rv()

    if options.modFreqs:
        mf_job = fit_job.addFollowOnJobFn(mod_freqs_job, options, model_id,
                                          disk=2 * 1024**3)
        model_id = mf_job.rv()
        export_parent = mf_job
    else:
        export_parent = fit_job

    export_parent.addFollowOnJobFn(export_file, model_id, model_export_path(options))
    return model_id


def extract_4d_all(job, options, chunks, ann_id, species_list, leaves_csv):
    """ Group consecutive same-contig chunks into groups of --fitChunkGroup
    and fan out one extract_4d_chunk_group job per group. Each group job is
    1-core, single-process. """
    groups = group_chunks_for_fit(chunks, options.fitChunkGroup)
    RealtimeLogger.info('4d extract: planning {} chunk-group jobs across {} chunks'.format(
        len(groups), len(chunks)))

    # Per-job resources: a chunk-group decompresses K chunks into one combined
    # MAF, then runs msa_view --4d on it. Worst-case ~K * (10x compressed
    # chunk size) of disk + a similar amount in memory.
    rvs = []
    for group in groups:
        # Sum compressed sizes of chunks in this group (each chunk_id is a Toil
        # FileID with .size). Decompressed MAF can be ~10x; msa_view scratch
        # adds another factor.
        group_compressed = sum(c[3].size for c in group)
        per_disk = max(512 * 1024**2,
                       int(group_compressed * 12) + int(ann_id.size * 2))
        # msa_view --4d peak doesn't scale with compressed input size — it
        # tracks the aligned-column count over the (small) CDS subset, which
        # is bounded by the reference-coordinate window per group. Observed
        # ~700 MiB on a 577-way ~10 Mb chunk-group; default memory is fine.
        rvs.append(job.addChildJobFn(extract_4d_chunk_group, options, group, ann_id,
                                     leaves_csv,
                                     disk=per_disk, cores=1).rv())
    return rvs


def extract_4d_chunk_group(job, options, group, ann_id, leaves_csv):
    """ Single-process 4d-site extraction for one group of consecutive same-
    contig chunks. Cats the chunks (Python-side, dedupe MAF headers), then
    runs msa_view --4d | msa_view --tuple-size 1. Chunks come pre-filtered
    (strip + mafDuplicateFilter) from the chunker.

    Returns the file ID of the per-group .tup.ss, or None if the group had
    no MAF blocks / no overlapping CDS. """
    import gzip  # shutil is imported at module scope
    work_dir = job.fileStore.getLocalTempDir()
    contig = group[0][0]
    summary = '{}.{}:{}-{}'.format(options.refGenome, contig, group[0][1], group[-1][2])
    # Coord-encoded label used for every per-group temp file so any error
    # message points back at the exact slice of the MAF.
    region_tag = '{}_{}-{}'.format(contig, group[0][1], group[-1][2])
    RealtimeLogger.info('4d extract group ({} chunks, {})'.format(len(group), summary))

    # localize the K chunk files
    chunk_paths = []
    for (c, s, e, chunk_id) in group:
        p = os.path.join(work_dir, 'chunk.{}_{}-{}.maf.gz'.format(c, s, e))
        job.fileStore.readGlobalFile(chunk_id, p)
        chunk_paths.append(p)

    # localize the prepped annotation (small) and project to this contig
    ann_path = os.path.join(work_dir, 'ann.prepped.gp')
    job.fileStore.readGlobalFile(ann_id, ann_path)
    local_ann = localize_annotation(work_dir, ann_path, contig, options.refGenome, '0')
    if not os.path.isfile(local_ann) or os.path.getsize(local_ann) == 0:
        RealtimeLogger.info('4d extract group: no annotation features on {}; skipping'.format(contig))
        return None

    # Decompress + concatenate the chunks into one MAF (drop comment lines
    # `##maf` etc. on chunks 1+). Chunks are already strip+dedup filtered.
    combined = os.path.join(work_dir, 'combined.{}.maf'.format(region_tag))
    with open(combined, 'wb') as out:
        for i, p in enumerate(chunk_paths):
            with gzip.open(p, 'rb') as f:
                if i == 0:
                    shutil.copyfileobj(f, out)
                else:
                    for line in f:
                        if not line.startswith(b'#'):
                            out.write(line)
            os.remove(p)

    if not _maf_has_blocks(combined):
        RealtimeLogger.info('4d extract group ({}): no MAF blocks; skipping'.format(summary))
        return None

    ss_path = os.path.join(work_dir, 'group.{}.4d.ss'.format(region_tag))
    cactus_call(parameters=['msa_view', '--4d', '--features', local_ann,
                            '-i', 'MAF', combined, '-o', 'SS'],
                outfile=ss_path)

    if not _ss_has_length(ss_path):
        RealtimeLogger.info('4d extract group ({}): SS LENGTH=0; skipping'.format(summary))
        return None

    # Project to leaves only so msa_view --aggregate downstream sees a
    # consistent species set across all per-group SS files.
    tup_path = os.path.join(work_dir, 'group.{}.tup.ss'.format(region_tag))
    cactus_call(parameters=['msa_view', '-i', 'SS', '-o', 'SS', '--tuple-size', '1',
                            '-l', leaves_csv, ss_path],
                outfile=tup_path)
    return job.fileStore.writeGlobalFile(tup_path)


def _maf_has_blocks(maf_path):
    """ Quick check for at least one alignment block (`a` line) in a MAF. """
    with open(maf_path) as f:
        for line in f:
            if line.startswith('a'):
                return True
    return False


def _ss_has_length(ss_path):
    """ Quick check that an SS file has LENGTH > 0 (msa_view --4d sometimes
    emits a header-only SS when no CDS columns overlap). """
    with open(ss_path) as f:
        for line in f:
            if line.startswith('LENGTH'):
                # `LENGTH = 12345`
                parts = line.split('=')
                if len(parts) == 2:
                    try:
                        return int(parts[1].strip()) > 0
                    except ValueError:
                        return False
    return False


def aggregate_4d(job, options, ss_results, species_list, tree_str):
    """ Collect all per-chunk-group SS files and run msa_view --aggregate.
    `ss_results` is a flat list of file IDs (or None for groups that yielded
    no 4d sites). """
    work_dir = job.fileStore.getLocalTempDir()
    in_paths = []
    for ss_id in ss_results:
        if ss_id is None:
            continue
        p = job.fileStore.readGlobalFile(ss_id)
        in_paths.append(p)

    if not in_paths:
        raise RuntimeError('No 4d sites were extracted; check that --geneAnnotation contig names '
                           'match the reference genome contigs in the input MAF.')

    out_path = os.path.join(work_dir, 'all.4d.ss')
    # Aggregate over the *leaves* of the tree (phast-canonical names), not all
    # genomes — the cactus447 doc builds this list from MAF s-line names which
    # only contain leaves. Including ancestors here would either pad rows with
    # all-Ns or silently collide with a leaf if any phast_name() coincides.
    leaves = [phast_name(g) for g in newick_leaves(tree_str)]
    if not leaves:
        raise RuntimeError(
            'newick_leaves returned no leaves from the HAL tree — refusing to silently fall '
            'back to the full species list (which would include ancestors and pad rows with '
            'all-Ns). Tree string was: {!r}'.format(tree_str))
    species_csv = ','.join(leaves)
    # --unordered-ss: our per-chunk SS files are tuple-count summaries (no
    # column order). Without this flag msa_view --aggregate refuses with
    # "msa_concat_from_files requires an ordered alignment".
    cactus_call(parameters=['msa_view', '--aggregate', species_csv,
                            '--in-format', 'SS', '--out-format', 'SS',
                            '--unordered-ss'] + in_paths,
                outfile=out_path)
    RealtimeLogger.info('Aggregated {} 4d-site SS files over {} leaf species into {}'.format(
        len(in_paths), len(leaves), out_path))
    return job.fileStore.writeGlobalFile(out_path)


def phyloFit_job(job, options, ss_id, tree_str, species_list):
    """ Single big serial job: run phyloFit. """
    work_dir = job.fileStore.getLocalTempDir()
    ss_path = os.path.join(work_dir, 'all.4d.ss')
    job.fileStore.readGlobalFile(ss_id, ss_path)
    tree_path = os.path.join(work_dir, 'tree.nh')
    rewritten = rewrite_tree_to_phast_names(tree_str, species_list)
    with open(tree_path, 'w') as tf:
        tf.write(rewritten + '\n')

    out_root = os.path.join(work_dir, 'phyloFit')
    threads = max(1, int(options.phyloFitCores or 1))
    RealtimeLogger.info('Running phyloFit (subst-mod={}, precision={}, threads={}); '
                        'this may take many hours.'.format(
                            options.substMod, options.precision, threads))
    # Stream phyloFit's per-iteration EM progress to the leader log so a
    # multi-hour fit shows up live (otherwise nothing prints until the
    # process finishes).
    cactus_call(parameters=['phyloFit', '--EM',
                            '--precision', options.precision,
                            '--msa-format', 'SS',
                            '--subst-mod', options.substMod,
                            '--threads', str(threads),
                            '--tree', tree_path,
                            '--out-root', out_root,
                            ss_path],
                realtimeStderrPrefix='phyloFit')
    mod_path = out_root + '.mod'
    if not os.path.isfile(mod_path):
        raise RuntimeError('phyloFit did not produce expected output {}'.format(mod_path))
    return job.fileStore.writeGlobalFile(mod_path)


def mod_freqs_job(job, options, model_id):
    """ Symmetrize background base frequencies via modFreqs. """
    work_dir = job.fileStore.getLocalTempDir()
    in_path = os.path.join(work_dir, 'in.mod')
    out_path = os.path.join(work_dir, 'out.mod')
    job.fileStore.readGlobalFile(model_id, in_path)

    # parse BACKGROUND line and pass GC fraction
    gc = None
    with open(in_path) as f:
        for line in f:
            if line.startswith('BACKGROUND:'):
                toks = line.split()
                if len(toks) >= 5:
                    a, c, g, t = float(toks[1]), float(toks[2]), float(toks[3]), float(toks[4])
                    gc = c + g
                break
    if gc is None:
        raise RuntimeError('modFreqs: could not parse BACKGROUND line in model file')
    cactus_call(parameters=['modFreqs', in_path, '{:.6f}'.format(gc)], outfile=out_path)
    return job.fileStore.writeGlobalFile(out_path)


# -----------------------------------------------------------------------------
# phyloP subgraph
# -----------------------------------------------------------------------------

def phyloP_workflow(job, options, chunks, model_id, ref_seq_lengths, species_list,
                    subtree, effective_root_name):
    """ Per-chunk phyloP scoring -> merge -> export wig (+ optional bigwig)
    for one output track.

    `subtree` is the user-given --subtree value for this track (or None when
    no --subtree was supplied at all). `effective_root_name` is the name of
    the (sub)tree root — the HAL tree root, or --root if set. If
    `subtree == effective_root_name`, this track is the standard global
    conservation track (no phyloP --subtree flag, no .s tag in output). """
    track_subtree = None if (subtree is None or subtree == effective_root_name) else subtree
    if track_subtree:
        RealtimeLogger.info('phyloP track {!r}: lineage-specific LRT'.format(track_subtree))
    else:
        RealtimeLogger.info('phyloP track: global conservation (no --subtree)')
    score_job = job.addChildJobFn(phyloP_all, options, chunks, model_id, species_list,
                                  track_subtree)
    per_chunk_wigs = score_job.rv()  # flat list of (contig, start, wig_id_or_None) per chunk

    # disk estimates: per-base wig text is ~10 bytes/ref_bp uncompressed; bgzip
    # shrinks ~10x; account for input + concat + bgzip output. 30 bytes/ref_bp
    # is a safe upper bound. For human-sized 3 Gb refs this is ~90 GiB.
    total_ref_bp = sum(ref_seq_lengths.values())
    merge_disk = max(2 * 1024**3, 30 * total_ref_bp)
    # bgzip on the concatenated wig is CPU-bound and parallelizable. 8 cores is
    # enough to be much faster than single-threaded without monopolizing a node.
    merge_cores = min(options.chunkCores or 8, 8)
    merge_job = score_job.addFollowOnJobFn(phyloP_merge, options, per_chunk_wigs,
                                           disk=merge_disk,
                                           cores=merge_cores)
    wig_id = merge_job.rv()
    merge_job.addFollowOnJobFn(export_file, wig_id, wig_export_path(options, track_subtree))

    if options.bigwig:
        # decompressed wig + .bw output + bbiFile scratch
        bw_disk = max(2 * 1024**3, 15 * total_ref_bp)
        # wigToBigWig memory is dominated by an in-memory index of every
        # chrom's per-base values; empirically ~10-12 bytes/ref_bp on the
        # 8-way hg38 run (which OOM'd at 31 GiB / 32 GiB cap). 24 bytes/bp
        # leaves comfortable headroom; floor at 8 GiB so smaller tests don't
        # over-request.
        bw_mem = max(8 * 1024**3, 24 * total_ref_bp)
        bw_job = merge_job.addFollowOnJobFn(wig_to_bigwig_job, options, wig_id, ref_seq_lengths,
                                            disk=bw_disk,
                                            memory=cactus_clamp_memory(bw_mem),
                                            cores=merge_cores)
        bw_job.addFollowOnJobFn(export_file, bw_job.rv(),
                                bigwig_export_path(options, track_subtree))

    return wig_id


def phyloP_all(job, options, chunks, model_id, species_list, subtree):
    """ Fan out one 1-core phyloP_chunk job per chunk. """
    rvs = []
    for chunk_spec in chunks:
        contig, start, end, chunk_id = chunk_spec
        # decompressed chunk + raw wig output + a small per-job scratch
        per_disk = max(256 * 1024**2,
                       int(chunk_id.size * 12) + 30 * (end - start))
        # phyloP streams MAF columns; working set is one column's state vector
        # (n_species × small). Observed peak ~400 MiB on 577-way; default
        # memory is fine.
        rvs.append(job.addChildJobFn(phyloP_chunk, options, chunk_spec, model_id,
                                     subtree,
                                     disk=per_disk, cores=1).rv())
    return rvs


def phyloP_chunk(job, options, chunk_spec, model_id, subtree):
    """ Single-process phyloP scoring for one chunk. Linear pipeline:
    decompress chunk -> phyloP -> awk-clip the per-base wig to the chunk's
    ref range. Chunks come pre-filtered (strip + mafDuplicateFilter) from
    the chunker. Returns (contig, start, wig_id), or (contig, start, None)
    for empty chunks. """
    contig, c_start, c_end, chunk_id = chunk_spec
    work_dir = job.fileStore.getLocalTempDir()
    # Coord-encoded suffix on every per-chunk temp file so any error message
    # (or a stuck process inspected via lsof) points back at the exact slice.
    region_tag = '{}_{}-{}'.format(contig, c_start, c_end)
    chunk_path = os.path.join(work_dir, 'chunk.{}.maf.gz'.format(region_tag))
    job.fileStore.readGlobalFile(chunk_id, chunk_path)
    model_path = os.path.join(work_dir, 'neutral.mod')
    job.fileStore.readGlobalFile(model_id, model_path)

    # Decompress chunk for phyloP (which reads MAF, not gzipped MAF).
    sub_maf = os.path.join(work_dir, 'sub.{}.maf'.format(region_tag))
    cactus_call(parameters=['bash', '-c',
        'set -eo pipefail && zcat {chunk} > {out}'.format(
            chunk=shlex.quote(chunk_path), out=shlex.quote(sub_maf))])
    os.remove(chunk_path)

    if not _maf_has_blocks(sub_maf):
        RealtimeLogger.info('phyloP: chunk {}.{}:{}-{} has no MAF blocks'.format(
            options.refGenome, contig, c_start, c_end))
        return (contig, c_start, None)

    raw_wig = os.path.join(work_dir, 'raw.{}.wig'.format(region_tag))
    phyloP_cmd = ['phyloP', '--method', options.method, '--mode', options.phyloPMode,
                  '--wig-scores', '--msa-format', 'MAF', '--chrom', contig]
    if subtree:
        phyloP_cmd += ['--subtree', subtree]
    phyloP_cmd += [model_path, sub_maf]
    try:
        cactus_call(parameters=phyloP_cmd, outfile=raw_wig)
    except RuntimeError as e:
        # Known per-chunk degeneracies in phast that aren't workflow failures
        # — skip the chunk for this track and let other chunks fill in. All
        # three cases stem from per-chunk tree pruning (phyloP first prunes
        # the model tree to leaves with data in the chunk's MAF; what remains
        # can be too-small or numerically pathological for the LRT):
        #   - "no node named '<subtree>'": the named internal node collapsed
        #     to a single leaf after pruning (only the reference left from
        #     its clade in this chunk).
        #   - "NaN detected in matrix exponentiation parameters": the per-
        #     branch rescaling fit produced a NaN branch_scale, typically on
        #     a chunk where the pruned subtree has so little data that the
        #     log-likelihood gradient is undefined. v1.9.7's rescaling fix
        #     (PR #106) covers many but not all of these cases.
        #   - "col_lrts_sub: delta_lnl = X <= -0.1": phast's column-wise LRT
        #     sanity check (phast_fit_column.c:894). The alt model should
        #     fit at least as well as the null, and small negative slack
        #     (between 0 and -0.1) is silently clamped to 0; below -0.1
        #     phast die()s. Triggered on data-sparse subtrees where the
        #     column optimization converges to a slightly-worse-than-null
        #     fit (catshark VertebratesAnc1 hits this on hundreds of chunks).
        msg = str(e)
        skip = None
        if subtree and "no node named '{}'".format(subtree) in msg:
            skip = 'subtree {!r} collapsed to a leaf after tree pruning'.format(subtree)
        elif 'NaN detected in matrix exponentiation parameters' in msg:
            skip = 'phyloP per-branch rescaling produced NaN'
        elif 'col_lrts_sub: delta_lnl' in msg:
            skip = 'phyloP column-LRT sanity check (delta_lnl <= -0.1) tripped on a sparse-data column'
        if skip is not None:
            RealtimeLogger.warning(
                'phyloP: skipping chunk {}.{}:{}-{} for track {!r}: {}'.format(
                    options.refGenome, contig, c_start, c_end,
                    subtree if subtree else '(global)', skip))
            os.remove(sub_maf)
            return (contig, c_start, None)
        # Unhandled phyloP failure: surface the chunk's identity in the leader
        # log before re-raising so the failure isn't only visible in the
        # worker stderr.
        RealtimeLogger.error('phyloP failed on chunk {}.{}:{}-{} for track {!r}'.format(
            options.refGenome, contig, c_start, c_end,
            subtree if subtree else '(global)'))
        raise
    os.remove(sub_maf)

    # Step 3: clip the wig to [c_start+1, c_end] (wig is 1-based; chunk coords
    # are 0-based half-open). phyloP scores entire MAF blocks even when those
    # blocks extend past the requested region, so adjacent chunks would have
    # overlapping positions that wigToBigWig refuses.
    out_wig = os.path.join(work_dir, 'out.{}.wig'.format(region_tag))
    clip_fixed_step_wig(raw_wig, out_wig, c_start + 1, c_end)

    if os.path.getsize(out_wig) == 0:
        return (contig, c_start, None)
    return (contig, c_start, job.fileStore.writeGlobalFile(out_wig))


_FIXED_STEP_KV_RE = re.compile(r'(chrom|start|step)=(\S+)')


def clip_fixed_step_wig(in_path, out_path, clip_start, clip_end):
    """ Clip a fixedStep WIG to a [clip_start, clip_end] (1-based, inclusive)
    reference range. Each fixedStep block tracks position via its `start=`
    and `step=` header fields; we drop lines whose position falls outside
    the range and rewrite the block header's `start=` if the first kept
    position differs. variableStep blocks (which we don't expect from
    phyloP --wig-scores but handle defensively) are skipped entirely. """
    with open(in_path) as inp, open(out_path, 'w') as out:
        in_block = False
        printed_header = False
        chrom = None
        pos = 0
        step = 1
        for line in inp:
            if line.startswith('fixedStep'):
                in_block = True
                printed_header = False
                chrom = None
                pos = 0
                step = 1
                for k, v in _FIXED_STEP_KV_RE.findall(line):
                    if k == 'chrom':
                        chrom = v
                    elif k == 'start':
                        pos = int(v)
                    elif k == 'step':
                        step = int(v)
                continue
            if line.startswith('variableStep'):
                in_block = False
                continue
            if in_block:
                if clip_start <= pos <= clip_end:
                    if not printed_header:
                        out.write('fixedStep chrom={} start={} step={}\n'.format(
                            chrom, pos, step))
                        printed_header = True
                    out.write(line)
                pos += step


def phyloP_merge(job, options, per_chunk_wigs):
    """ Concatenate all per-chunk wigs sorted globally by (contig, start) so
    the merged wig is in the order wigToBigWig requires. `per_chunk_wigs` is
    a flat list of (contig, start, wig_id_or_None) — None entries are chunks
    that legitimately produced no scores (no MAF blocks / no aligned data). """
    work_dir = job.fileStore.getLocalTempDir()
    flat = [(c, s, wid) for c, s, wid in per_chunk_wigs if wid is not None]
    flat.sort(key=lambda x: (x[0], x[1]))
    in_paths = [job.fileStore.readGlobalFile(wig_id) for _, _, wig_id in flat]

    if not in_paths:
        raise RuntimeError('No phyloP scores were produced.')

    cat_path = os.path.join(work_dir, 'merged.wig')
    catFiles(in_paths, cat_path)

    if options.output.endswith('.gz'):
        gz_path = cat_path + '.gz'
        # bgzip parallelizes well; use whatever cores Toil allocated to this job.
        threads = max(1, int(job.cores or 1))
        cactus_call(parameters=['bgzip', '-f', '--threads', str(threads), cat_path])
        return job.fileStore.writeGlobalFile(gz_path)
    return job.fileStore.writeGlobalFile(cat_path)


def wig_to_bigwig_job(job, options, wig_id, ref_seq_lengths):
    """ Build a BigWig from the merged wig. """
    work_dir = job.fileStore.getLocalTempDir()
    wig_path = os.path.join(work_dir, 'merged.wig')
    if options.output.endswith('.gz'):
        gz_path = wig_path + '.gz'
        job.fileStore.readGlobalFile(wig_id, gz_path)
        threads = max(1, int(job.cores or 1))
        cactus_call(parameters=['bgzip', '-d', '--threads', str(threads), gz_path])
    else:
        job.fileStore.readGlobalFile(wig_id, wig_path)

    sizes_path = os.path.join(work_dir, 'chrom.sizes')
    with open(sizes_path, 'w') as f:
        for c, l in ref_seq_lengths.items():
            f.write('{}\t{}\n'.format(c, l))

    bw_path = os.path.join(work_dir, 'out.bw')
    cactus_call(parameters=['wigToBigWig', wig_path, sizes_path, bw_path])
    return job.fileStore.writeGlobalFile(bw_path)


if __name__ == '__main__':
    main()
