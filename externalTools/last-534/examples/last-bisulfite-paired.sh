#! /bin/sh

# Align paired bisulfite-converted DNA reads to a genome.

# This assumes that reads1.fastq are all from the converted strand
# (i.e. they have C->T conversions) and reads2.fastq are all from the
# reverse-complement (i.e. they have G->A conversions).

# "GNU parallel" needs to be installed.

[ $# -eq 4 ] || {
    cat <<EOF
Typical usage:

  lastdb -uBISF my_f mygenome.fa
  lastdb -uBISR my_r mygenome.fa

  $(basename $0) my_f my_r reads1.fastq reads2.fastq > results.maf

EOF
    exit 2
}

# Try to get the LAST programs into the PATH, if they aren't already:
PATH=$PATH:$(dirname $0)/../src

tmp=${TMPDIR-/tmp}/$$
trap 'rm -f $tmp.*' EXIT

cat > $tmp.script << 'EOF'
t=$1.$$

lastal -pBISF -s1 -Q1 -e120 -i1 "$2" "$4" > $t.t1f
lastal -pBISR -s0 -Q1 -e120 -i1 "$3" "$4" > $t.t1r
last-merge-batches $t.t1f $t.t1r > $t.t1
rm $t.t1f $t.t1r

lastal -pBISF -s0 -Q1 -e120 -i1 "$2" "$5" > $t.t2f
lastal -pBISR -s1 -Q1 -e120 -i1 "$3" "$5" > $t.t2r
last-merge-batches $t.t2f $t.t2r > $t.t2
rm $t.t2f $t.t2r

last-pair-probs -m0.1 $t.t1 $t.t2 |
perl -F'(\s+)' -ane '$F[12] =~ y/ta/CG/ if /^s/ and $s++ % 2; print @F'
rm $t.t1 $t.t2
EOF

# Convert C to t, and all other letters to uppercase:
perl -pe 'y/Cca-z/ttA-Z/ if $. % 4 == 2' "$3" | split -l400000 -a5 - $tmp.1

# Convert G to a, and all other letters to uppercase:
perl -pe 'y/Gga-z/aaA-Z/ if $. % 4 == 2' "$4" | split -l400000 -a5 - $tmp.2

parallel --gnu --xapply sh $tmp.script $tmp "$1" "$2" ::: $tmp.1* ::: $tmp.2*
