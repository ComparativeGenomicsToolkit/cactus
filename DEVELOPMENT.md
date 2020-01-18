# Notes on developing and debugging cactus


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
- SON_TRACE_DATASETS location of test data set, currently available in
    git clone https://github.com/UCSantaCruzComputationalGenomicsLab/cactusTestData

## Running tests with docker in single machine mode
    make docker
    export CACTUS_USE_LOCAL_IMAGE=1
    make test
