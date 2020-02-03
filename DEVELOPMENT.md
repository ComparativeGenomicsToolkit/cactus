# Notes on developing and debugging cactus

## Overriding make settings
A file include.local.mk can be created in the root directory
to override make variables, including setting environment variables.
This should not be committed.

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
    git clone https://github.com/UCSantaCruzComputationalGenomicsLab/cactusTestData

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
