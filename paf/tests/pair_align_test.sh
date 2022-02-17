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

# Parameters for the the alignment is passed in as an argument
params_file=$1
echo "Parameters file" $params_file

# The tree
tree='((a:0.5,b:0.5)anc1:0.3,c:0.5)anc0:0.1;'
ingroup_tree='(a:0.5,b:0.5)anc1;'
# Outgroup
outgroups="c"

if [ "$2" = "cow-dog" ]; then
  # Cow dow
  sequences="a ${working_dir}/simCow.chr6 b ${working_dir}/simDog.chr6 c ${working_dir}/simHuman.chr6"
  sequence_files="${working_dir}/simCow.chr6 ${working_dir}/simDog.chr6 ${working_dir}/simHuman.chr6"
  sequenceNames="simCow.chr6 simDog.chr6 simHuman.chr6"
  ingroupSequenceNames="simCow.chr6 simDog.chr6"
else
  # Human mouse
  sequences="a ${working_dir}/simMouse.chr6 b ${working_dir}/simHuman.chr6 c ${working_dir}/simDog.chr6"
  sequence_files="${working_dir}/simMouse.chr6 ${working_dir}/simHuman.chr6 ${working_dir}/simDog.chr6"
  sequenceNames="simMouse.chr6 simHuman.chr6 simDog.chr6"
  ingroupSequenceNames="simMouse.chr6 simHuman.chr6"
fi

# Get the sequence files
for i in ${sequenceNames}
do
wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/${i} -O ${working_dir}/${i}
done

# Generate the paf alignments
cactus_local_alignment.py ./jobStore --logLevel INFO --outputFile ${working_dir}/output.paf "${sequences}" --params $params_file --speciesTree "${tree}" --referenceEvent anc1

# Split the alignment into primary and secondary
cp ${working_dir}/output.paf ./output_ugly.paf
grep "tl:i:1" ${working_dir}/output.paf > ${working_dir}/primary.paf || true
grep -v "tl:i:1" ${working_dir}/output.paf > ${working_dir}/secondary.paf || true

# Run consolidate
cactus_consolidated --logLevel INFO --outputFile ${working_dir}/output.c2h --params $params_file --alignments ${working_dir}/primary.paf --secondaryAlignments ${working_dir}/secondary.paf  --sequences "${sequences}" --speciesTree "${tree}" --outgroupEvents "${outgroups}" --referenceEvent anc1 --outputHalFastaFile ${working_dir}/output.fasta

# Generate the hal file
halAppendCactusSubtree ${working_dir}/output.c2h ${working_dir}/output.fasta "${ingroup_tree}" ${working_dir}/output.hal --outgroups "${outgroups}"

# Now generate the maf
hal2maf ${working_dir}/output.hal ${working_dir}/output.maf --onlySequenceNames

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
wc -l ${working_dir}/output.paf