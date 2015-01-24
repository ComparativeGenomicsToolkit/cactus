#! /usr/bin/env python

# Read MAF-format alignments, and write them, after moving the Nth
# sequence to the top in each alignment.

# Before writing, if the top sequence would be on the - strand, then
# flip all the strands.  But don't do this if the top sequence is
# translated DNA.

# Seems to work with Python 2.x, x>=4

import fileinput, itertools, optparse, os, signal, string, sys

def filterComments(lines):
    for i in lines:
        if i.startswith("#"): print i,
        else: yield i

def mafInput(lines):
    for k, v in itertools.groupby(lines, str.isspace):
        if not k: yield list(v)

def indexOfNthSequence(mafLines, n):
    for i, line in enumerate(mafLines):
        if line.startswith("s"):
            if n == 1: return i
            n -= 1
    raise Exception("encountered an alignment with too few sequences")

def rangeOfNthSequence(mafLines, n):
    """Get the range of lines associated with the Nth sequence."""
    start = indexOfNthSequence(mafLines, n)
    stop = start + 1
    while stop < len(mafLines):
        line = mafLines[stop]
        if not (line.startswith("q") or line.startswith("i")): break
        stop += 1
    return start, stop

complement = string.maketrans('ACGTNSWRYKMBDHVacgtnswrykmbdhv',
                              'TGCANSWYRMKVHDBtgcanswyrmkvhdb')
# doesn't handle "U" in RNA sequences
def revcomp(seq):
    return seq[::-1].translate(complement)

def flippedMafS(words):
    alnStart = int(words[2])
    alnSize = int(words[3])
    strand = words[4]
    seqSize = int(words[5])
    alnString = words[6]
    newStart = seqSize - alnStart - alnSize
    if strand == "-": newStrand = "+"
    else:             newStrand = "-"
    newString = revcomp(alnString)
    out = words[0], words[1], newStart, alnSize, newStrand, seqSize, newString
    return map(str, out)

def flippedMafP(words):
    flippedString = words[1][::-1]
    return words[:1] + [flippedString]

def flippedMafQ(words):
    qualityString = words[2]
    flippedString = qualityString[::-1]
    return words[:2] + [flippedString]

def flippedMafLine(mafLine):
    words = mafLine.split()
    if   words[0] == "s": return flippedMafS(words)
    elif words[0] == "p": return flippedMafP(words)
    elif words[0] == "q": return flippedMafQ(words)
    else: return words

def maxlen(s):
    return max(map(len, s))

def sLineFieldWidths(mafLines):
    sLines = (i for i in mafLines if i[0] == "s")
    sColumns = zip(*sLines)
    return map(maxlen, sColumns)

def joinedMafS(words, fieldWidths):
    formatParams = itertools.chain(*zip(fieldWidths, words))
    return "%*s %-*s %*s %*s %*s %*s %*s\n" % tuple(formatParams)

def joinedMafLine(words, fieldWidths):
    if words[0] == "s":
        return joinedMafS(words, fieldWidths)
    elif words[0] == "q":
        words = words[:2] + [""] * 4 + words[2:]
        return joinedMafS(words, fieldWidths)
    elif words[0] == "p":
        words = words[:1] + [""] * 5 + words[1:]
        return joinedMafS(words, fieldWidths)
    else:
        return " ".join(words) + "\n"

def flippedMaf(mafLines):
    flippedLines = map(flippedMafLine, mafLines)
    fieldWidths = sLineFieldWidths(flippedLines)
    return (joinedMafLine(i, fieldWidths) for i in flippedLines)

def isCanonicalStrand(mafLine):
    words = mafLine.split()
    strand = words[4]
    if strand == "+": return True
    alnString = words[6]
    if "/" in alnString or "\\" in alnString: return True  # frameshifts
    alnSize = int(words[3])
    gapCount = alnString.count("-")
    if len(alnString) - gapCount < alnSize: return True  # translated DNA
    return False

def mafSwap(opts, args):
    inputLines = fileinput.input(args)
    for mafLines in mafInput(filterComments(inputLines)):
        start, stop = rangeOfNthSequence(mafLines, opts.n)
        mafLines[1:stop] = mafLines[start:stop] + mafLines[1:start]
        if not isCanonicalStrand(mafLines[1]):
            mafLines = flippedMaf(mafLines)
        for i in mafLines: print i,
        print  # blank line after each alignment

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = "%prog [options] my-alignments.maf"
    description = "Change the order of sequences in MAF-format alignments."
    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-n", type="int", default=2,
                  help="move the Nth sequence to the top (default: %default)")
    (opts, args) = op.parse_args()
    if opts.n < 1: op.error("option -n: should be >= 1")

    try: mafSwap(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
