# Relase 2.5.2 2023-05-15

This Release patches some bugs in the pangenome pipeline and makes it a bit more user-friendly

- Fix support for multiple referenes via `--reference` and `--vcfReference`
- Fix bug where certain combinations of options (ie returning filtered but not clipped index) could lead to crash
- Fix crash when handling non-ascii characters in vg crash reports
- Fix the `--chrom-vg` option in `cactus-pangenome`
- New option `--mgCores` to specify number of cores for minigraph construction (rather than lumping in with `--mapCores` which is also used for mapping)
- Better defaults for number of cores used in pangenome pipeline on singlemachine.
- Fix bug where small contigs in the reference sample could lead to crashes if they couldn't map to themselves (and `--refContigs` was not used to specify chromosomes). `--refContigs` is now automatically set if not specifed. 
- Update to vg 1.48.0
- Update pangenome paper citation from preprint to published version.


# Release 2.5.1 2023-04-19

This Release mostly patches some bugs in the pangenome pipeline

- `cactus-pangenome` now saves PAF file in the output directory
- Ship version of `vg` that is patched to not make too-slow `giraffe` indexes for some complex graphs
- Fix bug where `.` characters in reference sample name could lead to strange errors at end of pipeline
- Better sample name validation for all pangenome tools to prevent confusion around `.`s.
- Update `taffy` to fix issues where `--filterGapCausingDupes` could lead to crashes in `cactus-hal2maf`
- Strip defaults of `taffy` commands from being specified in `cactus-hal2maf` -- they are now taken from `taffy`

# Release 2.5.0 2023-04-03

This Release greatly simplifies the interface for building pangenomes

- Introduction of `cactus-pangenome` command that, like `cactus` for the progressive aligner, runs the whole pipeline in one shot. Its inputs are a list of fasta sequences and sample names and it outputs the pangenome graph and various indexes, vcf, etc. Intermediate outputs are exported at the end of each stage, so low-level commands can be used to repeat or continue the workflow.
- `#` characters in fasta contig names no longer need to be cleaned out with special invocation of `cactus-preprocess` to avoid conflicts with `vg` indexes.  This now happens automatically within the pipeline.
- The lower level pangenome commands are all still supported and documented.  But `cactus-align-batch` is now deprecated.  This means it is removed from the documentation (except when describing older results) and from continuous integration tests, and it will give a warning when used. Users should now just run `cactus-align` instead of `cactus-align-batch` where applicable, using the former's `--batch` option. `cactus-align-batch` was a hack required to scale up before the cactus v2.0 rewrite but used nested Toil workflows and needed to go. 
- Documentation updated to focus on the simpler interface
- `GFAffix` updated to fix a rare crash
- `cactus-graphmap-join` (and `cactus-pangenome`) will not fail in the event a VCF indexing error (ex from chromosomes that are too long for `tbi`). Instead it will give a warning and produce no index.
- New `cactus-hal2maf` option `--keepGapCausingDupes` is changed to `--filterGapCausingDupes` and turned off by default. The underlying code has bugs that cause problems on certain datasets, and is not ready to be activated by default.

# Release 2.4.4 2023-03-16

This release includes some new export tools for the UCSC Genome Browser

- `cactus-hal2chains` created in order to convert HAL output from Cactus into sets of pairwise alignment chains, using either `halLiftover` or `halSynteny`
- `cactus-maf2bigmaf` created to convert `.maf` output from `cactus-hal2maf` to BigMaf and BigMaf Summary files for display on the Genome Browser
- `cactus-hal2maf` typo fixed where 3 (instead of 30) was set for the default value of `--maximumGapLength`
- Boost TAFFY normalization defaults in `cactus-hal2maf`, bringing `--maxmimumGapLength` to 100, and `--maximumBlockLengthToMerge` to 1000, and adding the heuristic block-breaking dupe filter from `taffy norm`. The latter is on by default to prevent over fragmentation, but can be disabled with `--kepGapCausingDupes`
- Remove `--onlyOrthologs` and `--noDupes` options from `cactus-hal2maf` and replace with the `--dupeMode` option. `--dupeMode single` is now the recommended way of getting at most one row / species.  More information about this added to the documentation.
- `--maxRefNFrac` option added to `cactus-hal2maf` to filter out blocks where the reference sequence is mostly Ns (default to filter out >95% Ns).
- Change abPOA scoring matrix to be more consistent with lastz parameters used by cactus, where `N` bases are penalized when aligned with other characters. Before, they could be aligned to anything. This will hopefully make the above filter less necessary.
- Fix bug where `cactus-blast --restart` would not work.

# Release 2.4.3 2023-03-07

This release patches a critical pangenome indexing bug introduced in v2.3.0, where a typo in the refactor of `cactus-graphmap-join` effectively caused *all* variation to be removed from the allele-frequency-filtered (ie .d2) graphs.

Changes
- Fix typo in `cactus-graphmap-join` where minimum fragment length (default 1000) was passed as minimum depth to `vg clip`, overriding the correct value.
- Introduce stub filtering in `cactus-graphmap-join` that cleans out all dangling nodes. Resulting graphs will have exactly two stubs (tips) per reference chromosome (just like minigraph). This can be toggled off via the `removeStubs` config parameter.
- Update HAL to fix a bug in `halRenameSequences`

# Release 2.4.2 2023-02-14

Changes include:

- Update HAL to fix `halAppendSubtree` crash in certain cases (for real this time)
- Upgrade to Toil 5.9.2 (ditto)

# Release 2.4.1 2023-02-08

Changes include:

- Fix bug that could cause `cactus-graphmap-join` to crash when running `vg clip`
- Upgrade to Toil 5.9.0
- Include `taffy` binary with bgzip support enabled
- Better error when trying to install on Python version < 3.7

# Release 2.4.0 2023-01-09

This release drastically increases Cactus's default chaining parameters, resulting in much cleaner and more linear alignments.

Other changes include

- Upgrade to Toil 5.8.0, which allows GPU counts to be assigned to jobs (which can be passed from cactus via the `--gpu` option). Toil only currently supports this functionality in single machine mode.
- Fix bug introduced in v2.3.1 that broke GPU-enabled lastz preprocessing
- Include latest vg release.
- Include latest version of taffy (fka taf).

# Relase 2.3.1 2022-12-23

This release contains bugfixes changes that should result in cleaner alignments both in pangenome mode in some cases.

- Re-activate pangenome-specific paf filters (these were deactivated by accident back in v2.2.0, and were important in difficult regions in acrocentric chromosomes in HPRC graph).
- Pangenome pipeline now included in legacy binaries. 
- Use mash distance (by default) to determine the minigraph construction order
- Update HAL to fix `halAppendSubtree` crash in certain cases
- Log `blossom5` stderr to hopefully help debug failurs associated with it.


# Release 2.3.0 2022-11-21

This release contains substantial changes to the Minigraph-Cactus pangenome pipeline, namely in `cactus-graphmap-join`. 

- The SAMPLE.HAPLOTYPE naming convention now only needed for non-haploid samples: you no longer need to add `.0` to non-reference samples
- All paths, including reference, are stored as W-lines in the GFA output (previously the reference was stored as a P-line)
- XG/GBWT/GG output replaced with GBZ
- Latest vg now shipped with Cactus (was stuck at v.1.40.0 for a while).  **In most cases vg version >= 1.44.0 is now required to use the output.**
- Giraffe indexing should now be more efficient
- `cactus-graphmap-join` now only needs to be run once. Previously, up to 3 times was suggested. The single invocation of `cactus-graphmap-join` can still create up to three graphs and sets of indexes.  Sensible defaults are provided to encourage users to, for example, index the filtered graph.
- Many options from `cactus-graphmap-join` have been removed or moved into the config xml
- `--delFilter` now enabled by default in `cactus-graphmap`.  In addition, it is broadened to disallow any contig from inducing a deletion greater than (half) its length. (this ratio can be controlled using the `delFilterQuerySizeThreshold` config parameter).
- `cactus-graphmap-join` will no longer ever produce a graphs with edges that aren't covered by any paths (was rare but possible before).
- `minigraph` version included update to fix assertion failure bug.
- [mafTools](https://github.com/ComparativeGenomicsToolkit/maftools) binaries now included in Cactus.
- Toil caching turned off by default when using the slurm batch system.

# Release 2.2.4 2022-11-03

This release contains some small bugfixes in addition to improved MAF output support

- New tool `cactus-hal2maf` added to speed up HAL->MAF conversion (replaces the old hal2mafMP.py from HAL), and includes TAF-based normalization and reference gaps.
- Big documentation refactor
- `--includeRoot` option added to `cactus-prepare`
- Phast binaries now included in Cactus release (previously only halPhyloP was included)
- `cactus --root` option regression from v2.2.0 is fixed
- Better error-handling in case of degree-2 ancestral nodes (1-parent, 1-descendant) in input tree
- Update HAL to version with improved memory utilization for `hal2maf`, `halExport` and `halAppendSubtree`, and new tool `halRemoveSubtree`

# Release 2.2.3 2022-10-04

This release contains yet another patch regarding the `minimumBlockDegreeToCheckSupport` filter option, this time in progressive alignment. It was supposed to have been toggled off entirely in v2.2.0 but instead was applied to all blocks. 

- `minimumBlockDegreeToCheckSupport` of <= 0 now interpreted as disabling the filter, rather than applying it to all blocks with > 0 support.  Before v2.2.0 it was applied to blocks with support > 10, since v2.2.0 it was applied to blocks with support > 0, and now it is disabled. 

# Release 2.2.2 2022-09-28

This release fixes a critcal bug in `cactus-align --pangenome` that was introduced in v2.2.0 and causes massive underalignments in pangenome graphs with more than 10 or so samples. Huge thanks for @minglibio for catching this!

- `cactus-align --pangenome` fixed to properly set pangenome overrides in the config.
- `gfaffix` updated to fix `cactus-graphmap-join` crash while normalizing some types of hairpins in the graph.  also, un-covered nodes left by `gfaffix` now filtered out before they can cause errors in `deconstruct`.
- fix `cactus-graphmap` regression where it woudn't run with docker binaries due to directory issues.
- `cactus-minigraph` now adds (non-reference) genomes in order of decreasing size by default.
- `--realTimeLogging` enabled by default
- better error message when invalid `--root` supplied to cactus
- print cactus commit at the top of each log

# Release 2.2.1 2022-09-07

This release patches recent regressions in the "blast" phase:

- Fix `cactus-blast` crash when using GPU on mammalian-sized genomes 
- Use tree distance (rather than 1) for determining lastz parameters for outgroups.

Other changes include:

- Fix regression in 2.2.0 where "legacy" binary release was built with same compiler options as normal release (and therefore no more portable).
- Make WDL output from `cactus-prepare` a bit cleaner.  Revise Terra best practices in the documentation to be much more efficient. 
- Update to Toil 5.7.1
- Update to minigraph 0.19


# Release 2.2.0 2022-08-19

This release contains a major update to the "blast" phase, where chaining logic is introduced to select lastz anchors, replacing the old quality-based heuristic. It also uses 1 fewer outgroup (2 instead of 3) by default, and no longer explicitly computes self-alignments, which should result in faster runtimes. 

Other changes include:

- Complete rewrite and drastic simplification of all code used to genereate lastz anchors
- PAF format now used natively throughout Cactus (replacing lastz cigars)
- Major refactor and cleanup of the "progressive" python module, removing vestiges of old Progressive Cactus repo
- Rewrite and simplifcation of the "Cactus Workflow" Python code.
- Intermediate files (project, multicactus project, experiment XML) all done away with.
- More explicit error message for "illegal instruction" signal (which commonly confused people trying to run on older CPUs)
- Fasta contig name checking and prefixing done at beginning of each tool (this should prevent cryptic `halAppendSubtree` errors in the pangenome pipeline)
- Update to newest SegAlign, which should fix an overflow bug that occurs when repeatmasking some data.
- Increase binary compatibility by linking with newer libxml2
- Add `cactus-terra-helper` tool to force-resume Terra workflows (when its own call caching fails)
- Small cleanup of `cactus-preprocess` interface

# Release 2.1.1 2022-06-15

This release includes:

- Update Segalign to fix crash while lastz-repeatmasking certain (fragmented?) assemblies using GPUs.
- Add [cactus-update-prepare](doc/cactus-update-prepare.md) which generates scripts for updating HAL alignments (Thanks @thiagogenez)
- Upgrade release (CPU) Docker images from Ubuntu 18.04 to Ubuntu 22.04.
- Upgrade release GPU Docker image from Ubuntu 18.04 / Cuda 10.2 to Ubuntu 20.04 / Cuda 11.4.3 (the most recent Cuda currently supported by Terra)

# Release 2.1.0  2022-06-02

This release introduces a major overhaul to the Minigraph-Cactus Pangenome Pipeline, including:

- Total documentation rewrite (doc/pangenome.md) with more explanations, a new example (data included) for yeast, and detailed instructions that exactly reproduce a HPRC pangenome.
- Incorporation of latest minigraph version that can write base alignments.  These alignments, via GAF cigars, are now used by the Minigraph-Cactus pipeline rather than the raw minimizers.
- Masking with dna-brnn is no longer needed or recommended (but it is still supported).  Instead, a graph with the full sequences is constructed and any trimming is done based on the alignment in postprocessing (via `cactus-graphmap-join`).
- Better Continuous Integration testing for the entire pangenome pipeline.

Graphs constructed with the new, simpler pipeline should be slightly more accurate and much cleaner.

Other changes include:

- Fix bug in Cactus (since v2.0) that sometimes caused spurious tiny self-alignments.
- Update to newer version of abPOA (improves stability, and some corner case accuracy)
- Fix Dockerfile so that Cactus Docker images are now much (5X) smaller.
- Fix HAL support for remote files in Cactus Docker images.
- Update HAL library to patched version that works for alignment updates (as described in doc/updating-alignments.md)

The `--gpu` option still doesn't always work.  When using the GPU outside the gpu Docker Release, it is still advised to set gpuLastz="true" in src/cactus/cactus_progressive_config.xml (and rerun `pip install -U`).  


# Release 2.0.5  2022-01-25

This release fixes fixes a major (though rare) bug where the reference phase could take forever on some inputs.  It includes a newer version of lastz which seems to fix some crashes as well.

- Debug symbols no longer stripped from `cactus_consolidated` binary in Release.
- Update to Toil 3.5.6
- Update examples to use Python 3.8 specifically (was previously just python3, but this is often python3.6, support for which was dropped in Toil 3.5.6)
- cactus-prepare WDL output can now batch up `hal_append_subtree` jobs
- Fix bug where "reference" phase within cactus_consolidated could take ages on some input
- Fix bug where --realTimeLogging flag would cause infinite loop after cactus_consolidated within some Docker invocations.
- Upgrade to more recent version of lastz

# Release 2.0.4   2021-11-12

This release fixes a bug introduced in 2.0.0 where ancestral sequences could not be specified in the input, which prevented the recommended producedure for updating existing alignments from working. 

- Fix `Assertion `cap_getSequence(cap) == sequence' failed` error when ancestral fasta provided in input seqfile.
- Several minor pangenome updates, mostly in `cactus-graphmap-join`

# Release 2.0.3   2021-07-22

This release fixes some issues in pangenome normalization and CAF running time.

- Fix new regression that caused CAF's secondary filter to sometimes take forever.  This code has been causing occaisional slowdowns for some time, but should finally be fixed once and for all.
- Fix cactus-preprocess to work on zipped fasta inputs even when not running dna-brnn.
- Fix normalization in cactus-graphmap-join
- Update to abPOA v1.2.5


# Release 2.0.2   2021-07-07

This release primarily addresses stability issues during pangenome construction.

Changelog:
- Use latest abpoa, which fixes bug where aligning >1024 sequences would lead to a segfault
- Update to Toil 5.4.0
- More consistently apply filters to minimap2 output in the fallback stage of graphmap-split
- Build abpoa with AVX2 SIMD extensions instead of SSE4.1 in order to work around instability when building pangenomes.  This ups the hardware requirements for releases, unfortunately, as AVX2 is slightly newer.
- Clean up CAF config parameters
- Fix CAF secondary filter worst-case runtime issue.  It was very rare but could add days to runtime.
- Slightly tune minimap2 thresholds used for chromosome splitting
- Normalization option added to cactus-graphmap-join (should be used to work around soon-to-be addressed underalignment bug)

# Release 2.0.1   2021-06-19

This a patch release that fixes an issue where the new `--consCores` option could not be used with `--maxCores` (Thanks @RenzoTale88).  It also reverts some last minute CAF parameter changes to something more tested (though known to be slow in some cases with large numbers of secondary alignments)

Changelog:
- Fix bug where `cactus` doesn't work when both `--maxCores` and `--consCores` are specified.
- Static binaries script more portable.
- Revert CAF trimming parameters to their previous defaults. 

# Release 2.0.0   2021-06-18

This release includes a major update to the Cactus workflow which should dramatically improve both speed and robustness. Previously, Cactus used a multiprocess architecture for all cactus graph operations (everything after the "blast" phase).  Each process was run in its own Toil job, and they would communicate via the CactusDisk database that ran as its own separate service process (ktserver by default). Writing to and from the database was often a bottleneck, and it would fail sporadically on larger inputs with frustrating "network errors". This has all now been changed to run as a single multithreaded executable, `cactus_consolidated`.  Apart from saving on database I/O, `cactus_consolidated` now uses the much-faster, SIMD-accelerated abPOA by default instead of cPecan for performing multiple sequence alignments within the BAR phase.

Cactus was originally designed for a heterogeneous compute environment where a handful of large memory machines ran a small number of jobs, and much of the compute could be farmed off to a large number of smaller machines.  While lastz jobs from the "preprocess" and "blast" phases (or `cactus-preprocess` and `cactus-blast`) can still be farmed out to smaller machines, the rest of cactus (`cactus-align`) can now only be run on more powerful systems.  The exact requirements depend as usual on genome size and divergence, but roughly 64 cores / 512G RAM are required for distant mammals. 

This release also contains several fixes and usability improvements for the pangenome pipeline, and finally includes `halPhyloP`.

Changlelog:
- Fold all post-blast processing into single binary executable,`cactus_consolidated`
- New option, `--consCores`, to control the number of threads for each `cactus_consolidated` process.
- Cactus database (ktserver) no longer used.
- abPOA now default base aligner, replacing cPecan
- cPecan updated to include multithreading support via MUM anchors (as opposed to spawning lastz processes), and can be toggled on in the config
- Fix bug in how `cactus-prepare` transmits Toil size parameters
- `cactus-prepare-join` tool added to combine and index chromosome output from `cactus-align-batch`
- `cactus-graphmap-split` fixes
- Update to latest Segalign
- Update to Toil 5.3
- Update HAL
- Add `halPhyloP` to binary release and docker images

# Release 1.3.0   2021-02-11

This release introduces the [Cactus Pangenome Pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md), which can be used to align together samples from the same species in order to create a pangenome graph:

- `cactus_bar` now has a POA-mode via the abpoa aligner, which scales better than Pecan for large numbers of sequences and is nearly as accurate if the sequences are highly similar
- `cactus-refmap` tool added to produce cactus alignment anchors with all-to-reference minimap2 alignments instead of all-to-all lastz
- `cactus-graphmap` tool added to produce cactus alignment anchors with all-to-reference-graph minigraph alignments instead of all-to-all lastsz
- `--maskAlpha` option added to `cactus-preprocess` to softmask (or clip out) satellite sequence using `dna-brnn.
- `cactus_bar` now has an option to ignore masked sequence with a given length threshold.
- `cactus-graphmap-split` tool added to split input fasta sequences into chromosomes using minigraph alignments in order to create alignment subproblems and improve scaling.
- `cactus-align-batch` tool added to align several chromsomes at once using one job per chromosome. (`--batch` option added to `cactus-align` to achieve the same using many jobs per chromosome)
- `--outVG` and `outGFA` options added to `cactus-align` to output pangenome graphs in addtion to hal.

Other changes:
- `cactus-prepare` scheduling bug fix
- `--database redis` option added to use Redis instead of Kyoto Tycoon
- `cactus-blast --restart` bug fix
- "Legacy" binary release provided for those whose hardware is too old to run the normal release. 

# Release 1.2.3   2020-10-05

- Fix bug where `cactus_fasta_softmask_intervals.py` was expecting 1-based intervals from GPU lastz repeatmasker.

hal2vg version included: [v1.0.1](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.1/hal2vg)
GPU Lastz version used in GPU-enabled Docker image: [8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7)

# Release 1.2.2   2020-10-02

- hal2vg updated to version 1.0.1
- GPU lastz updated for more disk-efficient repeat masking and better error handling
- Fixed memory in `cactus_convertAlignmentsToInternalNames`
- CAF fixes targeted towards pangenome construction

hal2vg version included: [v1.0.1](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.1/hal2vg)
GPU Lastz version used in GPU-enabled Docker image: [8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/8b63a0fe1c06b3511dfc4660bd0f9fb7ad7176e7)

# Release 1.2.1   2020-08-31

This release fixes bugs related to GPU lastz

- Cactus fixed to correctly handle GPU lastz repeatmasking output, as well to keep temporary files inside Toil's working directory.
- `cactus-prepare --wdl` updated to support `--preprocessorBatchSize > 1`
- `cactus_covered_intervals` bug fix and speedup
- GPU lastz updated to fix crash

GPU Lastz version used in GPU-enabled Docker image: [12af3c295da7e1ca87e01186ddf5b0088cb29685](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/12af3c295da7e1ca87e01186ddf5b0088cb29685)
hal2vg version included: [v1.0.0](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.0/hal2vg)


# Release 1.2.0   2020-08-21

This release adds GPU lastz repeatmasking support in the preprocessor, and includes hal2vg

Notable Changes:
 - GPU lastz repeat masking is about an order of magnitude faster than CPU masking, and should provide better results.  It's toggled on in the config file or by using the GPU-enabled Docker image.
 - hal2vg (pangenome graph export) included in the binary release as well as docker images.
 - update hal to [f8f3fa2dada4751b642f0089b2bf30769967e68a](https://github.com/ComparativeGenomicsToolkit/hal/commit/f8f3fa2dada4751b642f0089b2bf30769967e68a)

GPU Lastz version used in GPU-enabled Docker image: [f84a94663bbd6c42543f63b50c5843b0b5025dda](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/f84a94663bbd6c42543f63b50c5843b0b5025dda)
hal2vg version included: [v1.0.0](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.0.0/hal2vg)

# Release 1.1.1   2020-07-31

This release fixes how Kent tools required for `hal2assemblyHub.py` were packaged in 1.1.0 (thanks @nathanweeks).  

Notable Changes:
 - The required shared libaries to run the Kent tools are added to the Docker Image
 - The same Kent tools are removed from the binary release.  They were included under the assumption that they were statically built and fully standalone, but they are not.  Instead, instrucitons are provided to guide interested users to installing them on their own. 

GPU Lastz version used in GPU-enabled Docker image: [3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa)

# Release 1.1.0   2020-07-30

This release contains some important bug fixes, as well as major improvements to `cactus-prepare` functionality.

Notable Changes:
 - `cactus-prepare` improvements including:
    - WDL / Terra support
    - GPU lastz support
    - bug fixes
- Upgrade to Toil 4.1.0
- Fix bug causing `cactus-reference` to run forever in presence of 0-length branches
- Major speed improvement for `cactus-caf` by fixing secondary alignment filter.
- Include HAL python tools in Docker images and binary release
- Fix static binary for `cPecanRealign`
- Provide GPU-enabled Docker image for release
- Included HAL tools contains several crash fixes

GPU Lastz version used in GPU-enabled Docker image: [3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa](https://github.com/ComparativeGenomicsToolkit/SegAlign/commit/3e14c3b8ceeb139b3b929b5993d96d8e5d3ef9fa)

# Release 1.0.0   2020-04-19

This is the first official release of Cactus.  The goal is to provide a
stable, track-able version of Cactus to users.  The releases is provided in
source, static-compiled binaries, and Docker formats.

Notable Changes:
 - Kyoto Cabinet and Typhoon are now included as a sub-module.  This is due to
   the lack of consistent, stable releases and the difficulty in compiling it.
 - Cactus is now available as a tar file of static-linked binary executables,
   along with a wheel of the Cactus Python packages.  This avoids compilation and dependency problems. These should work on most Linux platforms when Docker is not used.
 - Added support to run Cactus in individual steps. This works around problems with using Toil in some distributed environments by dividing up alignment into tasks that can be run manually on separate machines.
   See *Running step by step (experimental)* in `README.md`.
 - Conversion to Python 3, allowing Toil to drop Python 2 support.


