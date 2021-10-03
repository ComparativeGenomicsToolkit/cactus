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

# The tree
tree='((a:0.5,b:0.5)anc1:0.3,c:0.5)anc0:0.1;'
ingroup_tree='(a:0.5,b:0.5)anc1;'
# The sequences
sequences="a ${temp_dir}/simCow.chr6 b ${temp_dir}/simDog.chr6 c ${temp_dir}/simHuman.chr6"
sequenceNames="simCow.chr6 simDog.chr6 simHuman.chr6"
# Outgroup
outgroups="c"

## Uncomment these next two lines to instead compare simHuman-simMouse using simDog as the outgroup
#sequences="a ${temp_dir}/simMouse.chr6 b ${temp_dir}/simHuman.chr6 c ${temp_dir}/simDog.chr6"
#sequenceNames="simMouse.chr6 simHuman.chr6 simDog.chr6"

# Get the sequence files
for i in ${sequenceNames}
do
wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/${i} -O ${temp_dir}/${i}
done

# Generate the paf alignments
python3 ../cactus_local_alignment.py ./jobStore  --logLevel DEBUG --outputFile ${temp_dir}/output.paf "${sequences}" --params ../../src/cactus/cactus_progressive_config.xml --speciesTree "${tree}" --referenceEvent anc1

# Split the alignment into primary and secondary
grep "tl:i:1" ${temp_dir}/output.paf > ${temp_dir}/primary.paf
grep -v "tl:i:1" ${temp_dir}/output.paf > ${temp_dir}/secondary.paf

# Run consolidate
cactus_consolidated --logLevel INFO --outputFile ${temp_dir}/output.c2h --params ../../src/cactus/cactus_progressive_config.xml --alignments ${temp_dir}/primary.paf --secondaryAlignments ${temp_dir}/secondary.paf  --sequences "${sequences}" --speciesTree "${tree}" --outgroupEvents "${outgroups}" --referenceEvent anc1 --outputHalFastaFile ${temp_dir}/output.fasta

# Generate the hal file
halAppendCactusSubtree ${temp_dir}/output.c2h ${temp_dir}/output.fasta "${ingroup_tree}" ${temp_dir}/output.hal --outgroups "${outgroups}"

# Now generate the maf
hal2maf ${temp_dir}/output.hal ./output.maf --onlySequenceNames

# Copy back the alignment file
mv ${temp_dir}/output.paf ./output.paf
