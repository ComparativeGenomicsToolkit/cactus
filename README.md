# Cactus
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/cactus.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/cactus)

Cactus is a reference-free whole-genome alignment program, as well as a pagenome graph construction toolkit.  

## Getting Cactus

* Use the precompiled binaries (Linux X86) or Docker image from [the latest release](https://github.com/ComparativeGenomicsToolkit/cactus/releases)
* See [below](#installing-manually-from-source) for details on building from source. 

## Align Genomes from Different Species

* See the [Progressive Cactus documenation](doc/progressive.md)
* Please cite the [Progressive Cactus paper](https://doi.org/10.1038/s41586-020-2871-y) when using Cactus.  Additional descriptions of the core algorithms can be found [here](https://doi.org/10.1101/gr.123356.111) and [here](https://doi.org/10.1089/cmb.2010.0252).
* Please cite the [HAL paper](https://doi.org/10.1093/bioinformatics/btt128) when using HAL tools. 

## Align Genomes from the Same Species and Build Pangenome Graphs 

* See the [Minigraph-Cactus Pangenome Pipeline documenatation](doc/pangenome.md)
* Please cite the [Minigraph-Cactus paper](https://doi.org/10.1038/s41587-023-01793-w).

## Acknowledgements

Cactus uses many different algorithms and individual code contributions, principally from Joel Armstrong, Glenn Hickey, Mark Diekhans and Benedict Paten. We are particularly grateful to:

- Yung H. Tsin and Nima Norouzi for contributing their 3-edge connected components program code, which is crucial in constructing the cactus graph structure, see: Tsin,Y.H., "A simple 3-edge-connected component algorithm," Theory of Computing Systems, vol.40, No.2, 2007, pp.125-142.
- Bob Harris for providing endless support for his [LastZ](https://github.com/lastz/lastz) pairwise, blast-like genome alignment tool.
- Melissa Jane Hubiz and Adam Siepel for halPhyloP and [Phast](http://compgen.cshl.edu/phast/).
- Sneha Goenka and Yatish Turakhia for [SegAlign](https://github.com/gsneha26/SegAlign), the GPU-accelerated version of LastZ.
- Yan Gao et al. for [abPOA](https://github.com/yangao07/abPOA)
- Heng Li for [minigraph](https://github.com/lh3/minigraph), [minimap2](https://github.com/lh3/minimap2), [gfatools](https://github.com/lh3/gfatools) and [dna-brnn](https://github.com/lh3/dna-rnn)
- Dany Doerr for [GFAffix](https://github.com/marschall-lab/GFAffix), used to optionally clean pangenome graphs.
- The vg team for [vg](https://github.com/vgteam/vg), used to process pangenome graphs.
- The authors of [Mash](https://github.com/marbl/Mash)

## Mailing List

Please subscribe to the [cactus-announce](https://groups.google.com/d/forum/cactus-announce) low-volume mailing list so we may reach about releases and other announcements.

## Installing Manually From Source

**Cactus requires Python >= 3.7 along with Python development headers and libraries**

Clone cactus and submodules
```
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
```

Create the Python virtual environment.  Install virtualenv first if needed with `python3 -m pip install virtualenv`.
```
cd cactus
virtualenv -p python3 cactus_env
echo "export PATH=$(pwd)/bin:\$PATH" >> cactus_env/bin/activate
echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus_env/bin/activate
source cactus_env/bin/activate
python3 -m pip install -U setuptools pip
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
```

If you have Docker installed, you can now run Cactus.  All binaries, such as `lastz` and `cactus-consolidated` will be run via Docker.  Singularity binaries can be used in place of docker binaries with the `--binariesMode singularity` flag.  Note, you must use Singularity 2.3 - 2.6 or Singularity 3.1.0+. Singularity 3 versions below 3.1.0 are incompatible with cactus (see [issue #55](https://github.com/ComparativeGenomicsToolkit/cactus/issues/55) and [issue #60](https://github.com/ComparativeGenomicsToolkit/cactus/issues/60)).

By default, cactus will use the image, `quay.io/comparative-genomics-toolkit/cactus:<CACTUS_COMMIT>` when running binaries. This is usually okay, but can be overridden with the `CACTUS_DOCKER_ORG` and `CACTUS_DOCKER_TAG` environment variables.  For example, to use GPU release 2.4.4, run `export CACTUS_DOCKER_TAG=v2.4.4-gpu` before running cactus.

### Compiling Binaries Locally
In order to compile the binaries locally and not use a Docker image, you need some dependencies installed.  On Ubuntu (we've tested on 20.04 and 22.04), you can look at the [Cactus Dockerfile](./Dockerfile) for guidance. To obtain the `apt-get` command:
```
grep apt-get Dockerfile | head -1 | sed -e 's/RUN //g' -e 's/apt-get/sudo apt-get/g'
```

Progressive Cactus can be built on ARM cpus including on Mac (with packages installed via Brew), but Minigraph-Cactus is currently X86-only.

To build Cactus, run
```
make -j 8
```
In order to run the Minigraph-Cactus pipeline, you must also run
```
build-tools/downloadPangenomeTools
```

In order to toggle between local and Docker binaries, use the `--binariesMode` command line option. If `--binariesMode` is not specified, local binaries will be used if found in `PATH`, otherwise a Docker image will be used.
