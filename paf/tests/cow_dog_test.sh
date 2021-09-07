#!/usr/bin/env bash

# Script to make a test alignment

# exit when any command fails
set -e

# log the commands its running
set -x

# Temporary chain files
temp_dir=$(mktemp -d -t temp_chains-XXXXXXXXXX)
# Make sure we cleanup the temp dir
trap "rm -rf ${temp_dir}" EXIT

# Get the sequence files
for i in simCow.chr6 simDog.chr6 simHuman.chr6
do
wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/${i} -O ${temp_dir}/${i}
done

# The tree
tree='((simCow:0.18908,simDog:0.16303)anc1:0.3,simHuman:0.3)anc0:0.1;'

# Generate the paf alignments
python3 ../cactus_local_alignment.py ./jobStore "simCow ${temp_dir}/simCow.chr6 simDog ${temp_dir}/simDog.chr6 simHuman ${temp_dir}/simHuman.chr6" --logLevel DEBUG --outputFile ${temp_dir}/output.paf --params ../../src/cactus/cactus_progressive_config.xml --speciesTree "${tree}" --referenceEvent anc1

# Split the alignment into primary and secondary
grep "tl:i:1" ${temp_dir}/output.paf > ${temp_dir}/primary.paf
grep -v "tl:i:1" ${temp_dir}/output.paf > ${temp_dir}/secondary.paf

# Run consolidate
cactus_consolidated --logLevel INFO --outputFile ${temp_dir}/output.c2h --params ../../src/cactus/cactus_progressive_config.xml --alignments ${temp_dir}/primary.paf --secondaryAlignments ${temp_dir}/secondary.paf --sequences "simCow ${temp_dir}/simCow.chr6 simDog ${temp_dir}/simDog.chr6 simHuman ${temp_dir}/simHuman.chr6" --speciesTree "${tree}" --outgroupEvents simHuman --referenceEvent anc1 --outputHalFastaFile ${temp_dir}/output.fasta

# Generate the hal file
halAppendCactusSubtree ${temp_dir}/output.c2h ${temp_dir}/output.fasta '(simCow:0.18908,simDog:0.16303)anc1;' ${temp_dir}/output.hal --outgroups simHuman

# Now generate the maf
hal2maf ${temp_dir}/output.hal ./output.maf --onlySequenceNames

# Copy back the alignment file
mv ${temp_dir}/output.paf ./output.paf
