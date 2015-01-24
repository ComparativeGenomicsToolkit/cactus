#! /bin/sh

# Sort MAF-format alignments by sequence name, then strand, then start
# position, then end position, of the top sequence.  Also, merge
# identical alignments.  Comment lines starting with "#" are written
# at the top, in unchanged order.  If option "-d" is specified, then
# alignments that appear only once are omitted (like uniq -d).

# Minor flaws, that do not matter for typical MAF input:
# 1) It might not work if the input includes TABs.
# 2) Preceding whitespace is considered part of the sequence name.  I
# want to use sort -b, but it seems to be broken in different ways for
# different versions of sort!
# 3) Alignments with differences in whitespace are considered
# non-identical.

# This script uses perl instead of specialized commands like uniq.
# The reason is that, on some systems (e.g. Mac OS X), uniq doesn't
# work with long lines.

# Make "sort" use a standard ordering:
LC_ALL=C
export LC_ALL

uniqOpt=1
whichSequence=1
while getopts hdn: opt
do
    case $opt in
	h)  cat <<EOF
Usage: $(basename $0) [options] my-alignments.maf

Options:
  -h  show this help message and exit
  -d  only print duplicate alignments
  -n  sort by the n-th sequence (default: 1)
EOF
	    exit
	    ;;
	d)  uniqOpt=2
            ;;
	n)  whichSequence="$OPTARG"
	    ;;
    esac
done
shift $((OPTIND - 1))

baseField=$((6 * $whichSequence))
a=$(($baseField - 4))
a=$a,$a
b=$(($baseField - 1))
b=$b,$b
c=$(($baseField - 3))
c=$c,$c
d=$(($baseField - 2))
d=$d,$d

# 1) Add digits to "#" lines, so that sorting won't change their order.
# 2) Replace spaces, except in "s" lines.
# 3) Join each alignment into one big line.
perl -pe '
s/^#/sprintf("#%.9d",$c++)/e;
y/ /\a/ unless /^s/;
y/\n/\b/ if /^\w/;
' "$@" |

sort -k$a -k$b -k${c}n -k${d}n |  # sort the lines

# Print only the first (or second) of each run of identical lines:
perl -ne '$c = 0 if $x ne $_; $x = $_; print if ++$c == '$uniqOpt |

# 1) Remove the digits from "#" lines.
# 2) Restore spaces and newlines.
perl -pe '
s/^#.{9}/#/;
y/\a\b/ \n/;
'
