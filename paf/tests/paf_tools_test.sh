#!/usr/bin/env bash

# Script to test pafs

# exit when any command fails
set -e

# log the commands its running
#set -x

working_dir=$(mktemp -d -t temp_chains-XXXXXXXXXX)
# Make sure we cleanup the temp dir
trap "rm -rf ${working_dir}" EXIT

sequenceNames="simCow.chr6 simDog.chr6"
for i in ${sequenceNames}
do
wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/${i} -O ${working_dir}/${i}.fa
done

lastz ${working_dir}/*.fa --format=paf > ${working_dir}/output.paf

# Run paf_view
echo "minimum local alignment identity"
paf_view -i ${working_dir}/output.paf ${working_dir}/*.fa | cut -f9 | sort | head -n1

# Run paf_view with trim
echo "paf_trim minimum local alignment identity"
paf_trim -i ${working_dir}/output.paf | paf_view ${working_dir}/*.fa | cut -f9 | sort | head -n1

# Run paf_view with invert
echo "paf_invert minimum local alignment identity"
paf_invert -i ${working_dir}/output.paf | paf_view ${working_dir}/*.fa | cut -f9 | sort | head -n1

# Run paf_view with chain
echo "paf_chain minimum local alignment identity"
paf_chain -i ${working_dir}/output.paf | paf_view ${working_dir}/*.fa | cut -f9 | sort | head -n1

# Run paf_view with shatter
echo "paf_shatter minimum local alignment identity (will be low as equal to worst run of matches)"
paf_shatter -i ${working_dir}/output.paf | paf_view ${working_dir}/*.fa | cut -f9 | sort | head -n1

# Run paf_view with tile
echo "paf_tile minimum local alignment identity"
paf_tile -i ${working_dir}/output.paf | paf_view ${working_dir}/*.fa | cut -f9 | sort | head -n1