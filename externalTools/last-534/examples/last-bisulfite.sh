#! /bin/sh

# Align bisulfite-converted DNA reads to a genome.

# This assumes that the reads are all from the converted strand
# (i.e. they have C->T conversions, not G->A conversions).

[ $# -gt 1 ] || {
    cat <<EOF
Typical usage:

  lastdb -uBISF my_f mygenome.fa
  lastdb -uBISR my_r mygenome.fa

  $(basename $0) my_f my_r reads.fastq > results.maf

EOF
    exit 2
}
my_f=$1
my_r=$2
shift 2

# Try to get the LAST programs into the PATH, if they aren't already:
PATH=$PATH:$(dirname $0)/../src

tmp=${TMPDIR-/tmp}/$$
trap 'rm -f $tmp.*' EXIT

# Convert C to t, and all other letters to uppercase:
perl -pe 'y/Cca-z/ttA-Z/ if $. % 4 == 2' "$@" > "$tmp".q

lastal -pBISF -s1 -Q1 -e120 "$my_f" "$tmp".q > "$tmp".f
lastal -pBISR -s0 -Q1 -e120 "$my_r" "$tmp".q > "$tmp".r

last-merge-batches "$tmp".f "$tmp".r | last-split -m0.1 |
perl -F'(\s+)' -ane '$F[12] =~ y/ta/CG/ if /^s/ and $s++ % 2; print @F'
