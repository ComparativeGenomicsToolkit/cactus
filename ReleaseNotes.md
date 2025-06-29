# Release 2.9.9 2025-09-11

This release adds FastGA support, along with some bug fixes

- (prototype) [FastGA](https://github.com/thegenemyers/FASTGA) support added to Progressive Cactus (via `--fastga`). 
- Fix gpu functionality with singularity binaries
- Fix `cactus-update-prepare` to give a valid `halRemoveGenomes` command
- Raise an error as soon as empty input FASTA is detected (as opposed to producing empty alignment and eventually crashing)
- Update to Toil 8.2.0
- Fix `--restart` flag bug in `cactus-minigraph`
- Fix `make` error that could (quite rarely) bug out by making a `bin` file instead of `bin/` directory
- Fix VCF construction crash when multiple VCF references are given, but not all present in every reference contig
- Update vg to v1.66

Note: as in previous releases (since v2.9.4), you must specify a patched config file when running `cactus` or `cactus-align` along with GPUs in order not to generate ancestral contigs that are too big for `2bit` as used in KegAlign (download the patch from the [releases](https://github.com/ComparativeGenomicsToolkit/cactus/releases) page).

# Release 2.9.8 2025-04-15

This release resolves some issues in Progressive Cactus

- `--queryhspbest=100000` added to default `lastz` (but not kegalign) parameters in order to prevent tool from taking forever in some worst-case (ie not fully masked) situations
- `paffy chunk` fixed to split FASTA contigs across multiple lines (instead of 1 / line).  This prevents kegalign crashes on large contigs (which can happen now that ancestors are big).
- Remove `--gpu` option from `cactus-align`.  Not only was it inconsistent and unneeded, it prevented `cactus-align` from even being run from the gpu docker image if no gpus were available. 
- Specify memory to Toil jobs that read in FASTA files (`faffy chunk`, `cacatus_santize_fasta_headers`). These tools read entire contigs into memory, and can therefore take non-trivial amounts for large contigs.

# Release 2.9.7 2025-03-17

- Resolve `vcfwave` normalization crash on reference chromosomes with zero variation, ex `chrEBV`. 
- Fix bug introduced in v2.9.6 where per-chromosome outputs could be misnamed.

# Release 2.9.6 2025-03-15

This release contains yet another `gfaffix` patch (hopefully the last for a long time), as well as a fix for `deconstruct`

- Update `gfaffix` to 0.2.1. Tests so far show this version to be much more stable on big graphs.
- Resolve `vg deconstruct` crash on reference chromosomes with zero variation, ex `chrEBV`. 

# Release 2.9.5 2025-03-12

This release patches a few issues that arose in the previous release

- Toil version updated from 8.0.0 to 8.1.0b1. This fixes slurm partition support by adding the `--slurmPartition` option.
- Fix (via interface update) how environment variables are passed through to Toil (which was constributing factor to partition selection issue above).
- Update Cactus's minimum Python requirement to 3.9. This will prevent cryptic Toil install errors. 
- Fix `--snarlStats` option, which previously returned an invalid gzip file
- Fix `mash` binary to have same user/group IDs as other files in the Docker image. The fact that it didn't apparently caused singularity issues for some users.
- Revert from `gfaffix` 0.2.0 back to 0.1.5b, as the decollapsing process of the former still crashes in some cases. 

# Release 2.9.4 2025-02-27

This release fixes compatibility with newer Python versions and improves the `--vcfwave` option, along with a few other patches

- `vcfwave` more parallelized (on chunks rather than whole-chromosome)
- vcf normalization turned on by default after `vcfwave`, including allele merging from `collapse_bubble` package
- fix `--defautCores` option
- fix compatibility with Python 3.13
- update to vg 1.6.3
- update Toil to 8.0.0
- add `--snarl-stats` option to print a table of large bubbles and their reference coordinates
- default ancestral genome construction parameters changed to give larger, more contiguous ancestors
- update to gfaffix 0.2.0
- fix recent regression where non-agctn IUPAC fasta characters would cause errors in `cactus_sanitizeFastaHeaders`, they are now read in as Ns (once again)
- minigraph graph (`.sv.gfa.gz`) now output with PanSN names, as opposed to Cactus's internal id=event| prefix.
- use more sensitive seeds for abPOA's progressive ordering
- reference genomes now added by mash distance order (instead of always first) by default.  this behaviour can be adjusted in the config with `minigraphSortInput` and `minigraphSortReference`.

# Releaes 2.9.3 2024-11-18

This release adds some new options to the pangenome pipeline, and hopefully improves robustness overall

- Faster path normalization (`vg paths -n`) for pangenomes via vg upgrade to v1.61.0
- Sanity checks added to better detect corrupted intermediate FASTA files
- Switch off abPOA's progressive mode unless input sequences have same length (otherwise sort by length)
- `--lastTrain / --scoresFile` options added to learn and/or use custom scoring models for multiple alignment using `last-train`.
- Update to latest `vcflib`.  Also add `vcflib` installation command as option to `BIN-INSTALL` instructions
- Make `--maxLen` default value consistent between `cactus-align --pangenome` and `cactus-pangenome`.  Previously it was 100X bigger in the former, which made it very easy to have wildly different performance between the all-at-once and step-by-step versions of the pipeline
- Fix bug where `--binariesMode singularity` could potentially attempt to write temporary files outside specified workDir
- Tighten disk usage estimate for `tile_alignments` job
- Patch `mafTools` to fix a bug where `taffy` normalization in `cactus-hal2maf` would crash if 1-character genome names were present in the input

# Release 2.9.2 2024-10-14

This release patches a couple bugs

- fix broken `--collapse` option 
- give Toil exact disk requirement for merge_aligments job, fixing a potential over-estimate
- update abpoa to latest release (v1.5.3)

# Release 2.9.1 2024-09-25

This release updates the pangenome pipeline, and adds `KegAlign` to progressive cactus. 

- GPU lastz implementation changed from `SegAlign` to `KegAlign`, and should be more robust and better supported as a result. 
- Path normalization added to pangenome pipeline to make sure no two equivalent alleles through any site have different paths. `AT` fields in VCF should now always be consistent with the graph as a result. 
- Always chop nodes to 1024bp by default in pangenome pipeline. This ensures that all outputs (gfa, gbz, vcf etc) have compatible node ids.  Before, only GBZ and downstream graphs were chopped which was too confusing.  Old logic can be re-activated using the config XML though.  
- Fix recent bug where using the `--mgMemory` option would crash `cactus-pangenome`
- (Experimental) `--collapse` option added to pangenome pipeline to incorporate nearby self-alignments, including on the reference path.
- Left shifting VCFs (`bcftools norm -f`) no longer run by default (except on `vcfwave`-normalized outputs), since it can cause conflicts with PanGenie by writing multiple variants at the same location.

# Release-2.9.0 2024-07-29

This release addresses two important scaling issues in the pangenome pipeline.

- The haplotype sampling index (`--haplo`) can now be built without giraffe indexes (`--giraffe`).  This significantly reduces peak memory consumption when using `--haplo`, especially for big diverse pangenomes.
- Previously, you could not align more than roughly 500 samples with Minigraph-Cactus, no matter how small the input genomes were. This bottleneck has been removed: you can now align as many genomes as your system resources allow. For very small genomes, this could be well into the tens of thousands.
- Two bugs were recently found in `vcfwave`, which can be run with the `--vcfwave` option since v2.8.2. First the `AT` field is wrong in the output. Second, and more seriously, genotypes can be incorrect. The latter seems specific to multiallelic sites (but I'm not sure).  This release now strips `AT` fields (they are not relevant after re-alignment anyway). It also splits multiallelic sites before running `vcfwave` in an attempt to work around the genotyping bug.  


# Release-2.8.4 2024-06-21

This release updates `vcfbub` in order to fix a longstanding issue where this tool can produce invalid VCFs. 

- `vcfbub` updated to `v0.1.1` which resolves a bub where records could be missing columns in the presence of `.` genotypes
- run `bcftools view` as sanity check on generated VCFs to prevent various normalization steps from ever silently producing invalid output.

# Release-2.8.3 2024-06-12

This release fixes some bugs and updates to the latest Toil.

- Fix broken `--restart` option in `cactus-graphmap`
- Raise Toil job memory requirement for `filter-paf-deletions`
- Update to `vg` v1.57.0
- Update `Toil` to v7.0
- Fix bug where trim-outgroups job could requeset way too little memory when there are no outgroups
- Fix typo that broke `cactus-maf2bigmaf` on uncompressed inputs
- More robust implementation of `vcfwave`
- Fix bug where RED preprocessing crashed `awk` returned a number in scientific notation

# Release 2.8.2 2024-05-09

This release fixes some bugs and adds a (docker-only) `vcfwave` normalization option for pangenomes.

- Use correct `bigChain.as` that allows chain scores to be huge (instead of capping them)
- Update `odgi`, `vg`, `abPOA` and `taffy` to their latest releases
- Fix `cactus-hal2maf` and `cactus-pangenome --odgi` to work with `--binariesMode docker`
- `bcftools norm -f` now run by default on all non-raw VCF outputs (toggle off in the config)
- Minigraph fasta file renamed from `.gfa.fa` to `.sv.gfa.fa` to be less confusing
- Gap and empty MAF block filtering moved from `cactus-hal2maf` to `cactus-maf2bigmaf`.  So MAF output will now have a reference base for every position.
- Fix `cactus-preprocess` to do only RED masking by default (there was previously a bug where it ran RED then lastz after).  The `--maskMode` option is also fixed to work properly.
- Update to Toil 6.1.0

# Release 2.8.1 2024-04-04

This release patches some major recent bugs, including a major bug introduced in `cactus-hal2maf` in v2.7.2 that could produce negative-stranded (and out of order) reference rows. 

- Do not apply RED masker to contigs that are likely to crash it (tiny contigs and extremely low information contigs)
- Add `--coverage` option to `cactus-hal2maf` to include table of coverage statistics in the output
- Fix bug where `:start-end` contig suffixes caused the pangenome pipeline to crash.  They are now correctly handled as subranges
- Turn off `abPOA` seeding by default, after finding (what must be a fairly rare) case where it doesn't work.
- Improve `cactus-hal2chains` interface
- Add range support to `cactus-hal2maf` via `--start/--length` or `--bedRanges`
- Deprecate `cactus-maf2bigmaf --chromSizes` (use `--halFile` instead, as it handles "."s in genome names properly)
- Fix bug where reference row could be lost in `cactus-hal2maf` MAF due to sorting error. 

# Release 2.8.0 2024-03-13

This release significantly changes the preprocessor step of Progressive Cactus in order to be more robust and efficient in the presence of unmasked repeats, something that seems more prevalent with newer, T2T assemblies. 

- Replace lastz repeatmasking with REepeat Detector ([RED](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5#Sec41)) in the Progressive Cactus preprocessor.  RED is more sensitive and orders of magnitude faster than the old lastz masking pipeline. Crucially, it is able to mask regions that would slip by RepeatMasker/WindowMasker/lastz in new T2T ape genomes that would otherwise break Cactus downstream. Tests so far show this change to make Cactus much faster and more robust.  The old lastz pipeline can still be toggled back on in the config.
- Delete many unneeded files that previously collected in the jobstore directory until the end of execution. This was a particular issue in large `cactus-pangenome` runs where the jobstore would creep up to several terabytes for HPRC-sized inputs.
- No longer require manually editing the blast chunksize in the config when running on Slurm (to reduce the number of jobs).  It's now scaled up automatically on slurm environments (by a factor controlled in the config).
- Fix bug introduced in last release where Cactus would not work on AWS/MESOS clusters unless `--defaultMemory` and `--maxMemory` options were specified (and in bytes). 
- Update to the latest `taffy` and `vg`

# Release 2.7.2 2024-02-23

This release improves MAF output, along with some other fixes

- `--maxMemory` option given more teeth. It is now used to clamp most large Toil jobs. On single-machine it defaults to system memory. This should prevent errors where Toil requrests more memory than available, halting the pipeline in an un-resumable state.
- Update to latest `taffy` and use newer MAF normalization. This should result in larger blocks and fewer gaps. MAF rows will now be sorted phylogenetically rather than alphabetically
- Better handle `.` characters in genome names during MAF processing. Previously neither duplicate filtering nor bigmaf summary creation could handle dots, but that should be fixed now.
- Duplicating filtering now done automatically in `cactus-maf2bigmaf`. 
- Disable support for multifurcations (aka polytomies or internal nodes with more than 2 children) in Progressive Cactus. I'm doing this because I got spooked by a drop in coverage I noticed recently in a 4-child alignment. This regression appears to be linked to the new PAF chaining logic that's been added over the past several months. Until that's resolved, Cactus will exit with an error if it sees degree > 2 in the tree. This behaviour can, however, be overridden in the XML configuration file.  

# Release 2.7.1 2024-01-19

This release adds some options to tune outgroup selection, as well as updates many included dependencies and tools

- Add `--chromInfo` option to specify sex chromosomes of input genomes, in order to make sure outgroups are selected accordingly
- Add `--maxOutgroups` option so that the number of outgroups can be toggled via the command line (previously required using a modified configuration file).
- Update to Toil 6.0.0
- Update to vg 1.54.0
- Update to odgi 0.8.4
- Update to latest taffy (fixing bug in paf export)
- Update to abPOA 1.5.1
- Fix Dockerfile so that Phast binaries are included
- `--indexMemory` now acts as upper limit on chromosome-level jobs.

# Release 2.7.0 2023-12-05

This release changes how outgroups are used during chaining during progressive alignment, and adds some pangenome options 

- Add `--xg` option for `xg` pangenome output.
- Add (experimental) `cactus-pagnenome --noSplit` option in order to bypass reference chromosome splitting. This was previously only possible by running step-by-step and not using `cactus-graphmap-join`.
- Add pangenome tutorial (developed for recent hackathon) to documentation.
- Update to `vg` version [1.53.0](https://github.com/vgteam/vg/releases/tag/v1.53.0). 
- Updated local alignment selection criteria. At each internal node of the guide tree Cactus picks a set of pairwise local alignments between the genomes being aligned to construct an initial sequence graph representing the whole genome alignment. This sequence graph is then refined and an ancestral sequence inferred to complete the alignment process for the internal node. The pairwise local alignments are generated with LASTZ (or SegAlign if using the GPU mode). To create a reliable subset of local alignments Cactus employs a chaining process that organizes the pairwise local alignments into pairwise chains of syntenic alignments, using a process akin to the chains and nets procedure used by the UCSC Browser. Previously, each genome being aligned, including both ingroup and outgroup genomes, was used to select a set of primary chains. That is, for each genome sequence non-overlapping chains of pairwise alignments were chosen, each of which could be to any of the other genomes in the set. Only these primary chains were then fed into the Cactus process to construct the sequence graph. This heuristic works reasonably well, in effect it allows each subsequence to choose a sequence in another genome with which it shares a most recent common ancestor. In the new, updated version we tweak this process slightly to avoid rare edge cases. Now each sequence in each ingroup genome picks primary chains only to other ingroup genomes. Thus the set of primary chains for ingroup genomes does not include any outgroup alignments. The outgroup genomes then get to pick primary chains to the ingroups, effectively voting on which parts of the ingroups they are syntenic too. The result of this change is that the outgroups are effectively only used to determine ancestral orderings and do not ever prevent the syntenic portions of two ingroups from aligning together.

# Release 2.6.13 2023-11-15

This release fixes an issue where Toil can ask for way too much memory for minigraph construction
- Cut default minigraph construction memory estimate by half
- Add `--mgMemory` option to override minigraph construction memory estimate no matter what
- Exit with a clear error message (instead of more cryptic crash) when user tries to run container binaries in a container
- Fix double Toil delete that seems to cause fatal error in some environments
- Fix `gfaffix` regular expression bug that could cause paths other than the reference to be protoected from collapse.

# Release 2.6.12 2023-11-07

The release contains fixes some recent regressions:

- Include more portable (at least on Ubuntu) `gfaffix` binary.
- Fix error where gpu support on singularity is completely broken.
- Fix `export_hal` and `export_vg` job memory estimates when `--consMemory` not provided.

# Release 2.6.11 2023-10-31

This release fixes a bug introduced in v2.6.10 that prevents diploid samples from working with `cactus-pangenome`

- Remove stray `assert False` from diploid mash distance that was accidentally included in previous release

# Release 2.6.10 2023-10-30

This release contains bug fixes for MAF export and the pangenome pipeline

- Patch `taffy` to fix bug where sometimes length fields in output MAF can be wrong when using `cactus-hal2maf --filterGapCausingDupes`
- Fix regression `cactus-graphmap-split / cactus-pangenome` so that small poorly aligned reference contigs (like some tiny unplaced GRCh38 contigs) do not get unintentionally filtered out. These contigs do not help the graph in any way, but the tool should do what it says and make a component for every single reference contig no matter what, which it is now fixed to do.
- Cactus will now terminate with a clear error message if any other `--batchSystem` than `single_machine` is attempted from *inside* a docker container.
- Mash distance order into `minigraph` construction fixed so that haplotypes from the same sample are always added contiguously in the presence of ties.
- CI fixed to run all `hal` tests, and not just a subset.
- `pip install wheel` added to installation instructions, as apparently that's needed to install Cactus with some (newer?) Pythons.

# Release 2.6.9 2023-10-20

This release contains some bug fixes and changes to docker image uploads

- GFAffix updated to latest release
- CI no longer pushes a docker image to quay.io for every single commit.
- CPU docker release now made locally as done for GPU
- `--binariesMode docker` will automatically point to release image (using GPU one as appropriate)
- `--consMemory` overrides `hal2vg` memory as well
- `--defaulMemory` defaults to `4Gi` when using docker binaries
- SegAlign job memory specification increased to something more realistic
- `--lastzMemory` option added to override SegAlign memory -- highly recommended on SLURM
- chromosome (.vg / .og) outputs from pangenome pipeline will have ref paths of form `GRCh38#0#chr1` instead of `GRCh38#chr1` to be more consistent with full-genome indexes (and PanSN in general)

# Release 2.6.8 2023-09-28

This release includes several bug fixes for the pangenome pipeline

- Fix bug where mash distances used to determine minigraph construction order could be wrong when input sample names differ by only their last character
- Fix some job memory specifications to better support slurm environments
- Guarantee that pangenome components have exactly two tips (one for each reference path endpoint). This is required for vg's now haplotype sampling logic.
- Add warning message when genomes that are too diverse are input to `cactus-pangenome`
- Update hal
- `--consMemory` option now overrides memory for hal export, in addition to `cactus_consolidated`
- Update `vg` to v1.51.0
- Fix bug where samples passed in to `--reference` (except the first) could be dropped as reference in the final output if they are missing from the first chromosome
- Port CI to OpenStack

# Release 2.6.7 2023-08-16

This release includes a patched vg and gfaffix

- Update to `vg` version `1.50.1` which patches the path name incompatibility bug in `1.50.0`
- Revert minigraph-cactus Reference path name convention introduced in v2.6.6 (ie haplotypes can be left unspecified)
- Upgrade to GFAffix 0.1.5 which fixes a crash among other things (thanks @danydoerr!)
- Fix bug where large input contig sizes (>2Gb) would break the pangenome pipeline with some versions of `awk`. 

# Release 2.6.6 2023-08-03

This release fixes a compatibility problem between Cactus and the newly released `vg` version `1.50.0`.

- Patch `hal2vg` to never write reference path names of the form `SAMPLE#CONTIG`, as `vg` now fails with an exception when reading them with `convert`. Instead, `SAMPLE#HAPLOTYPE#CONTIG` is always used, even if there is no haplotype specified (in which case it's set to 0).
- Add (prototype) `--haplo` option to build new subsampling-compatible giraffe indexes. 

# Release 2.6.5 2023-07-27

This release patches a Toil bug that broke GPU support on single-machine.

- Update to Toil v5.12, which fixes issue where trying to use GPUs on single machine batch systems would lead to a crash
- Make cactus more robust to numeric and duplicate internal node labels on input tree (ie ignore them instead of crashing with a cryptic scheduling error)
- Fix `hal2chains --targetGenomes` option.

# Release 2.6.4 2023-07-02

This is another minor patch release to fix support for multiple values to `--reference`.

- Don't fail if a given reference contig contains no sequences from the 2nd reference. This issue prevented completing of HPRC graphs with `--reference GRCh38 CHM13` because CHM13 wasn't present in the non-chromosome components.
- Fix `--dupeMode consensus` to output MAF with rows sorted (and, importantly, leaving the first row as reference). 

# Release 2.6.3 2023-06-29

This release contains a single, minor patch that only applies when passing multiple values to `--reference`.  

- Pangenome stub filtering changed to apply only to the first genome passed via `--reference` in order to be consistent with gap filter.ing (and enforce maximum of two stubs per graph component).

# Release 2.6.2 2023-06-28

This release patches a few bugs introduced or found in v2.6.1

- Docker container fixed to include runtime libatomic dependency for odgi
- `--refContigs` option fixed in `cactus-graphmap-split`
- `cactus-pangenome` fixed to properly output intermediate GAF and unfiltered PAF alignment files

# Release 2.6.1 2023-06-27

This Release adds SLURM cluster support for Cactus (both progressive and pangenome). It also adds some new visualization features to the pangenome pipeline, along with several bugfixes.

- SLURM support added (via Toil v5.11).  It's important to read the documentation about `--consMemory` and `--indexMemory` when running large aligments.
- ODGI now integrated into Cactus.  See `cactus-pangenome` options `--viz, --odgi, --chrom-og` and `--draw` for incorporating it into your output.
- Add `--chop` option to `cactus-pangenome` make sure all output graphs have nodes chopped down to at most `1024bp`. By default, the `.gfa.gz` is unchopped and the `.gbz` is, which can lead to annoying confusion when trying to compare node IDs across different files.
- Fix bug where `_MINIGRAPH_` paths ended up in the `.full` output graphs.
- `vg clip` crash fixed.
- `mash` distance for minigraph order now computed by sample (and not by haplotype).  Fixes issue where, for example, diploid ordering would be dependent on whether assembly has chrX.
- If `--refContigs` is not specified, minigraph-cactus now uses naming in addition to size to try to guess reference contigs.
- `--dupeMode consensus` option added to `cactus-hal2maf` in order to use `maf_stream` to merge multiple rows into consensus rows, which may be the best compromise to get the data into PhyloP or the Browser.
- `halReplaceGenome` patched to fix a regression from late 2022 where updating large alignments could lead to a crash.
- Fix `cactus-update-prepare replace` to print the `halRemoveGenome` command rather than quietly running it on the input.
- `--consMemory` and `--indexMemory` options bugs fixed.
- Fix bug with multuple `--reference` samples; they are now all promoted properly to REFERENCE-sense paths. 

# Release 2.6.0 2023-06-09

This Release significantly updates Cactus's chaining logic which, in early tests at least, allows alignment of T2T-quality assemblies as well as WindowMasked (and not RepeatMasked) genomes. 

- Improved chaining of lastz's PAF output in order to support alignment of T2T-quality genomes
- Minimum chain length in the Cactus graph now determined by branch length, so that more closely-related genomes can be chained more aggressively while not losing sensitivity along more distant branches where the alignment is more fragmented.
- Early experiments show that the above changes make Cactus much less sensitive to the input repeat masking. Genomes that previously required masking with RepeatMasker were able to align with the WindowMasking-based fastas directly from NCBI. 
- Update to Toil 5.10.0
- Update to latest Taffy
- Specify memory requirements for all Toil jobs (in Progressive Cactus). Cactus Consolidated memory is estimated conservatively, but can be overridden with `--consMemory`. 

# Release 2.5.2 2023-05-15

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


