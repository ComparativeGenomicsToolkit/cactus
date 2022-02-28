#!/usr/bin/env bash

# Script to make a test alignment

# exit when any command fails
set -e

# log the commands its running
set -x

# Directory to run compute in
working_dir=$(mktemp -d -t temp_chains-XXXXXXXXXX)
# Make sure we cleanup the temp dir
trap "rm -rf ${working_dir}" EXIT
# and the jobStore
trap "rm -rf ./jobStore" EXIT

# Parameters for the the alignment is passed in as an argument
params_file=$1
echo "Parameters file" $params_file

# Example file to use with cactus_prepare
example_file=$2
echo "Example file" $example_file

# Ingroup sequence names (hard coded)
ingroupSequenceNames="simMouse.chr6 simHuman.chr6"

# Run cactus prepare
cactus-prepare ${example_file}  --outDir ${working_dir}

# Switch example file to relative path
base_example_file=`basename ${example_file}`

## Preprocessor
cactus-preprocess ./jobstore/0 ./${example_file} ${working_dir}/${base_example_file} --inputNames c a b --realTimeLogging --logInfo --retryCount 0

## Generate the paf file
cactus-blast ./jobstore/1 ${working_dir}/${base_example_file} ${working_dir}/anc1.paf --root anc1 --realTimeLogging --logDebug --retryCount 0

# Run the cactus step
cactus-align ./jobstore/2 ${working_dir}/${base_example_file} ${working_dir}/anc1.paf ${working_dir}/anc1.hal --root anc1 --realTimeLogging --logInfo --retryCount 0 --maxCores 2

# Now generate the maf
hal2maf ${working_dir}/anc1.hal ${working_dir}/output.maf --onlySequenceNames

# Get the true alignment
wget -q https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/all.maf -O ${working_dir}/mammals-truth.maf

# Run mafComparator
mafComparator --maf1 ${working_dir}/mammals-truth.maf --maf2 ${working_dir}/output.maf --out ${working_dir}/maf_comparison.xml

# Report the stats for the relevant comparison
for i in ${ingroupSequenceNames}
do
grep "${i}" -A 2 ${working_dir}/maf_comparison.xml > ${working_dir}/maf_comparison.xml.2
mv ${working_dir}/maf_comparison.xml.2 ${working_dir}/maf_comparison.xml
done
grep "singleHomologyTest" -A 2 ${working_dir}/maf_comparison.xml | grep "all"

# Report number of PAFs
wc -l ${working_dir}/anc1.paf
