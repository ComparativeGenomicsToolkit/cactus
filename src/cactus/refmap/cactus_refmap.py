#!/usr/bin/env python3
"""
cactus-refmap: map assemblies to a single reference with minimap2 and emit a
cactus-ready PAF (query and target columns stamped as 'id=<event>|<contig>').

This is a *rescue* / fallback alignment source for the minigraph-cactus pangenome
pipeline.  It is NOT a graph aligner and builds no graph: it produces direct
assembly->reference pairwise alignments that cactus-align can consume alongside
(or in place of) the minigraph graphmap PAF.  The intended use is to replace
minigraph mappings where the graphmap produces conflicting / overlapping records
(e.g. in segmental duplications), since a direct alignment to the single-copy
reference backbone is far less ambiguous.

Name compatibility with the graphmap PAF matters: cactus only links records that
share sequence names.  So we parse the seqfile and sanitize headers in pangenome
mode (the same way cactus-graphmap does) and stamp 'id=<event>|<contig>' onto the
minimap2 output, so the emitted names line up with the graphmap PAF.

This file was a dead, unused early-pangenome prototype; it has been rewritten from
scratch for the fallback-rescue use case.
"""

import os
import xml.etree.ElementTree as ET

from toil.common import Toil
from toil.job import Job

from cactus.shared.common import (setupBinaries, importSingularityImage, makeURL,
                                   cactus_call, cactusRootPath, cactus_override_toil_options,
                                   catFiles, cactus_clamp_memory)
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.progressive_decomposition import parse_seqfile, get_event_set
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers


def rescue_refmap_workflow(job, seq_id_map, reference, preset):
    """Sanitize headers (pangenome-style, to match graphmap names), then map every
    non-reference event against the reference and consolidate to a single PAF."""
    sanitize_job = job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True)
    return sanitize_job.addFollowOnJobFn(map_all_to_ref, sanitize_job.rv(), reference, preset).rv()


def index_reference(job, ref_id, preset):
    """Build the minimap2 index of the reference ONCE (with the preset's -k/-w), so each map_one
    loads it instead of re-indexing the whole reference per assembly -- indexing dominates a
    reference-vs-assembly run and is identical every time.  This index job also serves as the lead:
    its children (the map_one jobs) run after it and read its return value (the .mmi), and its
    follow-on (consolidate) runs after them, so map_all_to_ref's return resolves before the rescue
    follow-on -- the same ordering the old _empty lead gave."""
    ref_path = job.fileStore.readGlobalFile(ref_id)
    mmi = job.fileStore.getLocalTempFile()
    cactus_call(parameters=["minimap2", "-x", preset, "-t", str(job.cores), "-d", mmi, ref_path])
    return job.fileStore.writeGlobalFile(mmi)


def map_all_to_ref(job, seq_id_map, reference, preset):
    """Index the reference once, then minimap2 each non-reference event against that index and
    consolidate one PAF."""
    if reference not in seq_id_map:
        raise RuntimeError("reference event '{}' not found in seqfile events: {}".format(
            reference, list(seq_id_map.keys())))
    ref_id = seq_id_map[reference]
    # index generation (mm_idx_gen) is the memory-hungry step, so the index job gets the big
    # allocation; the map jobs only LOAD the index, so they get less.
    index_job = job.addChildJobFn(index_reference, ref_id, preset, cores=8,
                                  disk=4 * ref_id.size,
                                  memory=cactus_clamp_memory(max(32 * 2**30, 16 * ref_id.size)))
    mmi_id = index_job.rv()
    paf_ids = {}
    for event, seq_id in seq_id_map.items():
        if event == reference:
            continue
        paf_ids[event] = index_job.addChildJobFn(
            map_one, event, seq_id, mmi_id, preset, cores=8,
            disk=4 * (seq_id.size + ref_id.size),
            memory=cactus_clamp_memory(max(16 * 2**30, 8 * ref_id.size))).rv()
    return index_job.addFollowOnJobFn(consolidate, paf_ids).rv()


def map_one(job, event, asm_id, mmi_id, preset):
    """minimap2 <ref.mmi> <asm> -> PAF, loading the prebuilt reference index.  Names already match
    the graphmap PAF: the fastas were sanitized with pangenome=True, which stamps id=<event>|<contig>
    onto the headers, so minimap2's query (assembly) and target (reference) columns are already
    id=<event>|<contig>.  The .mmi carries the reference sequence, so -c base-level alignment works
    exactly as against the fasta (minimap2 warns that -x's indexing params are overridden by the
    index -- harmless, it was built with this same preset)."""
    asm_path = job.fileStore.readGlobalFile(asm_id)
    mmi_path = job.fileStore.readGlobalFile(mmi_id)
    raw_paf = job.fileStore.getLocalTempFile()
    # -c emits the base-level CIGAR (cg:Z:) that cactus-align consumes; secondary mappings
    # (tp:A:S) are kept so the gap-fill can tell whether the reference maps a gap uniquely.
    cactus_call(parameters=["minimap2", "-c", "-x", preset, "-t", str(job.cores),
                            "-o", raw_paf, mmi_path, asm_path])
    return job.fileStore.writeGlobalFile(raw_paf)


def consolidate(job, paf_ids):
    """Concatenate the per-assembly PAFs into a single output PAF."""
    out_paf = job.fileStore.getLocalTempFile()
    catFiles([job.fileStore.readGlobalFile(pid) for pid in paf_ids.values() if pid is not None], out_paf)
    return job.fileStore.writeGlobalFile(out_paf)


def get_options():
    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("seqFile", type=str,
                        help="Cactus seqfile (the newick tree is ignored; events must include the reference).")
    parser.add_argument("reference", type=str,
                        help="Event name (from seqFile) to use as the single reference target.")
    parser.add_argument("outputFile", type=str,
                        help="Output PAF: assembly->reference, cactus-prefixed (id=<event>|<contig>).")
    parser.add_argument("--minimapPreset", default="asm5",
                        help="minimap2 -x preset (default asm5; use asm20 for more divergent input).")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"),
                        help="Cactus configuration file.")
    parser.add_argument("--latest", action="store_true",
                        help="Use the latest version of the docker container.")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"], default=None,
                        help="The way to run the Cactus binaries.")
    parser.add_argument("--containerImage", default=None,
                        help="Use the specified prebuilt container image rather than pulling one from quay.io.")
    return parser.parse_args()


def main():
    options = get_options()
    cactus_override_toil_options(options)

    with Toil(options) as toil:
        setupBinaries(options)
        importSingularityImage(options)

        config_node = ET.parse(options.configFile).getroot()
        config_wrapper = ConfigWrapper(config_node)
        config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)

        # pangenome=True so event/contig naming matches the graphmap PAF
        mc_tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile, config_wrapper, pangenome=True)
        event_set = get_event_set(mc_tree, config_wrapper, {}, mc_tree.getRootName())

        leaves = [mc_tree.getName(leaf) for leaf in mc_tree.getLeaves()]
        if options.reference not in leaves:
            raise RuntimeError("--reference '{}' not found in seqfile events: {}".format(options.reference, leaves))

        seq_id_map = {}
        for genome, seq in input_seq_map.items():
            if genome in event_set:
                seq_id_map[genome] = toil.importFile(makeURL(seq))

        if not toil.options.restart:
            paf_id = toil.start(Job.wrapJobFn(rescue_refmap_workflow, seq_id_map,
                                              options.reference, options.minimapPreset))
        else:
            paf_id = toil.restart()

        toil.exportFile(paf_id, makeURL(options.outputFile))


if __name__ == "__main__":
    main()
