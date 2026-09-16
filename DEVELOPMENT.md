# Notes on developing and debugging cactus

## Overriding make settings
A file include.local.mk can be created in the root directory
to override make variables, including setting environment variables.
This should not be committed.

It is included before any build decision is made, so it is the tidiest place for a setting
that belongs to one checkout rather than to your shell.  For instance, put
`CACTUS_NATIVE_BUILD = 1` there on a workstation whose binaries never leave it, and there is
nothing to remember or unset later.  It cannot leak into anything we distribute: the
Dockerfile deletes it before building, and `build-tools/makeBinRelease` builds a fresh clone.

## Environment variables controlling how cactus is built
- CACTUS_NATIVE_BUILD - compile for the machine doing the compiling (`-march=native`)?
  - unset <default> - portable: the baseline in include.mk, which any reasonably modern
    CPU of the same architecture can run
  - 1 - native.  Only for binaries that stay on the build machine.  A native binary
    SIGILLs on any older CPU, which is the normal outcome of building on a cluster head
    node and running on a worker.  No effect on ARM or with CACTUS_LEGACY_ARCH.

- CACTUS_PORTABLE_BUILD - force the portable baseline, overriding CACTUS_NATIVE_BUILD.
  Set by the Dockerfile and by build-tools/makeBinRelease, since those builds get shipped.
  Redundant otherwise, as portable is the default.

- CACTUS_LEGACY_ARCH - build for pre-AVX2 CPUs (`-msse2`), for the legacy binary release.

- CACTUS_ARCH_FLAGS - the baseline as literal compiler flags, e.g.
  `CACTUS_ARCH_FLAGS="-march=x86-64-v3 -mtune=znver3"`.  Overrides all of the above.

Run `make arch-flags` to print what the baseline resolves to without building anything.
The reasoning behind each value is in the comments in include.mk.

## Environment variables controlling how cactus is run
- CACTUS_BINARIES_MODE - how are cactus programs found?
  - docker <default>
  - singularity
  - local
- CACTUS_DOCKER_MODE - is Docker being used?
  - 1 <default>
  - 0
- CACTUS_USE_LOCAL_IMAGE - is Docker image on local server?
  - 0 <default>
  - 1

## Environment variables controlling tests
- SON_TRACE_DATASETS location of test data set, currently available with
    git clone https://github.com/ComparativeGenomicsToolkit/cactusTestData

- SONLIB_TEST_LENGTH  filters tests by maximum run time length category (case-insensitive)
  - SHORT - tests taking less than ~10 seconds, with some exceptions <default>
  - MEDIUM - tests taking less than ~100 seconds
  - LONG - test taking less than ~1000 seconds
  - VERG_LONG - test taking even longer

- CACTUS_TEST_LOG_LEVEL - Set log-level used for the test, may not set it for all test, but very useful for Toil.
  


## Running tests with docker in single machine mode
    make docker
    export CACTUS_USE_LOCAL_IMAGE=1
    make test

## Debugging hints
   - The main Cactus Python process will print out a stack trace of all of the Python
     threads if sent a SIGUSR1 signal.  They will then continue execution.  This
     maybe useful in determining the state of cactus.
   
   

