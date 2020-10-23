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


