#! /bin/sh

# This script demonstrates using LAST and maf-join to construct a
# multiple alignment of the human, mouse, chicken, and fugu
# mitochondrial genomes.

# You have to run this from inside the examples directory, else it
# won't find the files.

# Put the LAST programs and scripts into the command search path:
PATH=$PATH:../src:../scripts

# Make a LAST database of the human sequence:
lastdb -c humanMito humanMito.fa

# Align the mouse sequence to the human sequence:
lastal -e25 -j4 humanMito mouseMito.fa | last-split | maf-sort > hm.maf

# Align the chicken sequence to the human sequence:
lastal -e25 -j4 humanMito chickenMito.fa | last-split | maf-sort > hc.maf

# Align the fugu sequence to the human sequence:
lastal -e25 -j4 humanMito fuguMito.fa | last-split | maf-sort > hf.maf

# Join the pairwise alignments into a multiple alignment:
maf-join hm.maf hc.maf hf.maf

# Clean up the intermediate files that we made:
rm humanMito.??? h?.maf
