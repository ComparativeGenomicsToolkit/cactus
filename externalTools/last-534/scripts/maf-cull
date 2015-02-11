#! /usr/bin/env python

# Read MAF-format alignments.  Write them, omitting alignments whose
# coordinates in the top-most sequence are contained in those of >=
# cullingLimit higher-scoring alignments.

# Alignments on opposite strands are not considered to contain each
# other.

# The alignments must be sorted by sequence name, then strand, then
# start coordinate.

# This algorithm is not theoretically optimal, but it is simple and
# probably fast in practice.  Optimal algorithms are described in:
# Winnowing sequences from a database search.
# Berman P, Zhang Z, Wolf YI, Koonin EV, Miller W.
# J Comput Biol. 2000 Feb-Apr;7(1-2):293-302.
# (Use a "priority search tree" or an "interval tree".)

# Seems to work with Python 2.x, x>=4

import fileinput, itertools, operator, optparse, os, signal, sys

# The intervals must have size > 0.

def isFresh(oldInterval, newInterval):
    return oldInterval.end > newInterval.start

def freshIntervals(storedIntervals, newInterval):
    """Yields those storedIntervals that overlap newInterval."""
    return [i for i in storedIntervals if isFresh(i, newInterval)]

def isDominated(dog, queen):
    return dog.score < queen.score and dog.end <= queen.end

def isWanted(newInterval, storedIntervals, cullingLimit):
    """Is newInterval dominated by < cullingLimit storedIntervals?"""
    dominators = (i for i in storedIntervals if isDominated(newInterval, i))
    return len(list(dominators)) < cullingLimit

# Check that the intervals are sorted by start position, and further
# sort them in descending order of score.
def sortedIntervals(intervals):
    oldStart = ()
    for k, v in itertools.groupby(intervals, operator.attrgetter("start")):
        if k < oldStart: raise Exception("the input is not sorted properly")
        oldStart = k
        for i in sorted(v, key=operator.attrgetter("score"), reverse=True):
            yield i

def culledIntervals(intervals, cullingLimit):
    """Yield intervals contained in < cullingLimit higher-scoring intervals."""
    storedIntervals = []
    for i in sortedIntervals(intervals):
        storedIntervals = freshIntervals(storedIntervals, i)
        if isWanted(i, storedIntervals, cullingLimit):
            yield i
            storedIntervals.append(i)

class Maf:
    def __init__(self, lines):
        self.lines = lines
        try:
            aLine = lines[0]
            aWords = aLine.split()
            scoreGen = (i for i in aWords if i.startswith("score="))
            scoreWord = scoreGen.next()
            self.score = float(scoreWord.split("=")[1])
        except: raise Exception("missing score")
        try:
            sLine = lines[1]
            sWords = sLine.split()
            seqName = sWords[1]
            alnStart = int(sWords[2])
            alnSize = int(sWords[3])
            strand = sWords[4]
            self.start = seqName, strand, alnStart
            self.end = seqName, strand, alnStart + alnSize
        except: raise Exception("can't interpret the MAF data")

def mafInput(lines):
    for k, v in itertools.groupby(lines, str.isspace):
        if not k: yield Maf(list(v))

def filterComments(lines):
    for i in lines:
        if i.startswith("#"): print i,
        else: yield i

def mafCull(opts, args):
    inputLines = fileinput.input(args)
    inputMafs = mafInput(filterComments(inputLines))
    for maf in culledIntervals(inputMafs, opts.limit):
        for i in maf.lines: print i,
        print  # blank line after each alignment

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = "%prog [options] my-alignments.maf"
    description = "Cull alignments whose top-sequence coordinates are contained in LIMIT or more higher-scoring alignments."
    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-l", "--limit", type="int", default=2,
                  help="culling limit (default: %default)")
    (opts, args) = op.parse_args()
    if opts.limit < 1: op.error("option -l: should be >= 1")

    try: mafCull(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
