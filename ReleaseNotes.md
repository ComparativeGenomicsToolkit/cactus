# Release 2.1.0

Changelog:
- Option to use stable minigraph coordinates.  This obviates need for phony mingraph genome, whose small contigs can lead to stitching woes in BAR.  
- Fix cactus-preprocess to work on zipped fasta inputs even when not running dna-brnn.
- Do centromere clipping using softmasked regions in pangenome pipeline (which means using hash graph instead of packed graph before clipping)

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


