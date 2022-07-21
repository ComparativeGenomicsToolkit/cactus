#!/usr/bin/env bash

# Script to make a test alignment

# exit when any command fails
set -e

# log the commands its running
set -x

# Cleanup job store when done
trap "rm -rf ./jobStore" EXIT

# Parameters for the the alignment is passed in as an argument
params_file=$1
echo "Parameters file" $params_file

# Example file to use with cactus_prepare
example_file=$2
echo "Example file" $example_file

# Directory to run compute in
working_dir=$3
echo "Working directory" $working_dir

# Run cactus prepare
cactus-prepare ${example_file}  --outDir ${working_dir}

# Switch example file to relative path
base_example_file=`basename ${example_file}`

## Preprocessor
time cactus-preprocess ./jobstore/0 ./${example_file} ${working_dir}/${base_example_file} --realTimeLogging --logInfo --retryCount 0

## Generate the paf file
time cactus-blast ./jobstore/1 ${working_dir}/${base_example_file} ${working_dir}/Anc0.paf --root Anc0 --realTimeLogging --logInfo --retryCount 0

# Run the cactus step
time cactus-align ./jobstore/2 ${working_dir}/${base_example_file} ${working_dir}/Anc0.paf ${working_dir}/Anc0.hal --root Anc0 --realTimeLogging --logInfo --retryCount 0 --maxCores 2

# Now generate the maf
time hal2maf ${working_dir}/Anc0.hal ${working_dir}/output.maf --onlySequenceNames

# Report number of PAFs
wc -l ${working_dir}/Anc0.paf
