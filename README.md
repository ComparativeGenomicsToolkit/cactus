# Cactus
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/cactus.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/cactus)

Cactus is a reference-free whole-genome alignment program, as well as a pangenome graph construction toolkit.  

## Getting Cactus

* Use the precompiled binaries (Linux X86) or Docker image from [the latest release](https://github.com/ComparativeGenomicsToolkit/cactus/releases)
* See [below](#installing-manually-from-source) for details on building from source. 

## Getting help

Please subscribe to the [cactus-announce](https://groups.google.com/d/forum/cactus-announce) low-volume mailing list so we may reach out about releases and other announcements.

To ask questions or request help, please use the [Cactus GitHub Discussions](https://github.com/ComparativeGenomicsToolkit/cactus/discussions).

To file a bug report or enhancement request against the code or documentation, create a [GitHub Issue](https://github.com/ComparativeGenomicsToolkit/cactus/issues).

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
- B Gulhan, R Burhans, R Harris, M Kandemir, M Haeussler, A Nekrutenko for [KegAlign](https://github.com/galaxyproject/KegAlign), the GPU-accelerated version of LastZ.
- Yan Gao et al. for [abPOA](https://github.com/yangao07/abPOA)
- Heng Li for [minigraph](https://github.com/lh3/minigraph), [minimap2](https://github.com/lh3/minimap2), [gfatools](https://github.com/lh3/gfatools) and [dna-brnn](https://github.com/lh3/dna-rnn)
- Dany Doerr for [GFAffix](https://github.com/marschall-lab/GFAffix), used to optionally clean pangenome graphs.
- The vg team for [vg](https://github.com/vgteam/vg), used to process pangenome graphs.
- The authors of [Mash](https://github.com/marbl/Mash)
- Andrea Guarracino, Erik Garrison and co-authors for [odgi](https://github.com/pangenome/odgi). Make sure to [cite odgi](https://doi.org/10.1093/bioinformatics/btac308) when using it or its visualizations.
- Hani Z. Girgis for [RED](http://toolsmith.ens.utulsa.edu/)
- Erik Garrison and co-authors for [vcfwave](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md). [vcflib citation](https://doi.org/10.1371/journal.pcbi.1009123)
- Martin Frith and collaborators for `last-train` from [last](https://gitlab.com/mcfrith/last). [last-train citation](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btw742)
- Han Cao for [merge_duplicates](https://github.com/Han-Cao/collapse-bubble) vcf cleanup script
- Gene Myers, Richard Durbin and Chenxi Zhou for [FastGA](https://github.com/thegenemyers/FASTGA/)

## Installing Manually From Source

The instructions below are meant primarily for developers.  Everyone else should try to use the precompiled binaries (Linux X86) or Docker image from [the latest release](https://github.com/ComparativeGenomicsToolkit/cactus/releases) instead. 

The top-level cactus interface (`cactus`, `cactus-pangenome`, `cactus-hal2maf`, etc) is awlays a Python package that is `pip install`ed into a Python virtualenv.  This package runs several other tools as subprocesses, which are compiled into binaries.

### Contents

* [Cloning Cactus](#cloning-cactus)
* [Creating the Python virtualenv](#creating-the-python-virtualenv)
* [Compiling the binaries](#compiling-the-binaries)
* [Building on Mac (Read first if you're on a Mac!)](#building-on-mac) 

### Cloning Cactus

Cactus contains many submodules, so it is necessary to clone with `--recursive` or to run `git submodule update --init --recursive` inside the cactus directory after cloning.
```
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
```

### Creating the Python virtualenv

**Cactus requires Python >= 3.9 along with Python development headers and libraries**

Install virtualenv first if needed with `python3 -m pip install virtualenv`.

Create the Python virtual environment run (from inside `cactus/`):
```
python3 -m virtualenv  cactus_env
echo "export PATH=$(pwd)/bin:\$PATH" >> cactus_env/bin/activate
echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus_env/bin/activate
echo "export LD_LIBRARY_PATH=$(pwd)/lib:$LD_LIBRARY_PATH" >> cactus_env/bin/activate
source cactus_env/bin/activate
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
```

If you have Docker installed, you can now run Cactus.  All binaries, such as `lastz` and `cactus-consolidated` will be run via Docker using the latest release.  Singularity binaries can be used in place of docker binaries with the `--binariesMode singularity` flag.  Note, you must use Singularity 2.3 - 2.6 or Singularity 3.1.0+. Singularity 3 versions below 3.1.0 are incompatible with cactus (see [issue #55](https://github.com/ComparativeGenomicsToolkit/cactus/issues/55) and [issue #60](https://github.com/ComparativeGenomicsToolkit/cactus/issues/60)).

By default, cactus will use the image corresponding to the latest release when running docker binaries. This is usually okay, but can be overridden with the `CACTUS_DOCKER_ORG` and `CACTUS_DOCKER_TAG` environment variables.  For example, to use GPU release 2.4.4, run `export CACTUS_DOCKER_TAG=v2.4.4-gpu` before running cactus.

### Compiling the binaries

In order to compile the binaries locally and not use a Docker image, you need some dependencies installed.  On Ubuntu (we've tested on 20.04 and 22.04), you can look at the [Cactus Dockerfile](./Dockerfile) for guidance. To obtain the `apt-get` command:
```
grep apt-get Dockerfile | head -1 | sed -e 's/RUN //g' -e 's/apt-get/sudo apt-get/g'
```

Progressive Cactus can be built on ARM cpus including on [Mac](#building-on-mac), but Minigraph-Cactus is currently X86-only.

To build Cactus, run (from inside `cactus/`):
```
make -j 8
```
In order to run the Minigraph-Cactus pipeline, you must also run
```
build-tools/downloadPangenomeTools
```
If you want to work with MAF, including running `cactus-hal2maf`, you must also run
```
build-tools/downloadMafTools
```

In order to toggle between local and Docker binaries, use the `--binariesMode` command line option. If `--binariesMode` is not specified, local binaries will be used if found in `PATH`, otherwise a Docker image will be used.

### Building on Mac

These are the steps I used to build Cactus on a new M4 Mac Mini with MacOS Sequoia 15.5:

#### Developer Tools

Install command-line developer tools.  I did this by typing `make` on the command line (in Terminal), and accepting the prompt in the pop-up window to install them.  The version installed, as obtained from `pkgutil --pkg-info=com.apple.pkg.CLTools_Executables` was

```
package-id: com.apple.pkg.CLTools_Executables
version: 16.4.0.0.1.1747106510
volume: /
location: /
install-time: 1751461503
```

#### Homebrew

I pasted the install commanid from the [Homebrew homepage](https://brew.sh/) into the Terminal and ran it.  For me this command was the following, but you're probably better off to get it from the webpage

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then I installed the brew dependencies for Cactus

```
brew install coreutils libomp hdf5 libdeflate parallel wget samtools bcftools 
```

When installing `libomp` above (use `brew reinstall libomp` if you missed it), it printed some messages about setting:

```
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
```

So I added some flags to this effect to ~/.zprofile.  Doing so, and including `CFLAGS`, is crucial for the build to work.

```
printf "export LDFLAGS=\"-L/opt/homebrew/opt/libomp/lib\" >> ~/.zprofile
printf "export CPPFLAGS="-I/opt/homebrew/opt/libomp/include" >> ~/.zprofile
printf "export CFLAGS="-I/opt/homebrew/opt/libomp/include" >> ~/.zprofile
```

Make sure to reload the profile to apply the changes immediately:
```
source ~/.zprofile
```

#### Virtualenv

Python3 seemed to be installed by default.  In order to install virtualenv I ran

```
python3 -m pip install virtualenv
```

#### That's it

You should now be able to run the [steps](#installing-manually-from-source) to clone, setup the virtualenv and compile the binaries exactly as they are described above. 

**Note that Minigraph-Cactus is not (yet) supported on Mac.**

I recommend running
```
make evolver_test_local
```
in order to test your installation (you will need to have run `build-tools/downloadMafTools` when setting up the binaries for this to work)
