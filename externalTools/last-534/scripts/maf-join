#! /usr/bin/env python

# Copyright 2009, 2010, 2011 Martin C. Frith

# Join two or more sets of MAF-format multiple alignments into bigger
# multiple alignments.  The 'join field' is the top genome, which
# should be the same for each input.  Each input should be sorted by
# position in the top genome.

# WARNING: Alignment columns with a gap in the top genome are joined
# arbitrarily!!!

import sys, os, fileinput, optparse, signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # stop spurious error message

class peekable:  # Adapted from Python Cookbook 2nd edition
    """An iterator that supports a peek operation."""
    def __init__(self, iterable):
        self.it = iter(iterable)
        self.cache = []
    def __iter__(self):
        return self
    def next(self):
        if self.cache: return self.cache.pop()
        else: return self.it.next()
    def peek(self):
        if not self.cache: self.cache.append(self.it.next())
        return self.cache[0]

def maxLen(things): return max(map(len, things))

class MafBlock:
    def __init__(self, chr, beg, end, strand, chrSize, seq, prob):
        self.chr = chr  # chromosome names
        self.beg = beg  # alignment begin coordinates
        self.end = end  # alignment end coordinates
        self.strand = strand
        self.chrSize = chrSize  # chromosome sizes
        self.seq = seq  # aligned sequences, including gaps
        self.prob = prob  # probabilities (may be empty)

    def __nonzero__(self):
        return len(self.seq) > 0

    def __cmp__(self, other):
        return cmp(self.chr[:1] + self.beg[:1], other.chr[:1] + other.beg[:1])

    def before(self, other):
        return (self.chr[0], self.end[0]) <= (other.chr[0], other.beg[0])

    def after(self, other):
        return (self.chr[0], self.beg[0]) >= (other.chr[0], other.end[0])

    def addLine(self, line):
        words = line.split()
        if line.startswith('s'):
            self.chr.append(words[1])
            self.beg.append(int(words[2]))
            self.end.append(int(words[2]) + int(words[3]))
            self.strand.append(words[4])
            self.chrSize.append(words[5])
            self.seq.append(list(words[6]))
        elif line.startswith('p'):
            self.prob.append(words[1])

    def write(self):
        beg = map(str, self.beg)
        size = [str(e-b) for b, e in zip(self.beg, self.end)]
        seq = [''.join(i) for i in self.seq]
        columns = self.chr, beg, size, self.strand, self.chrSize, seq
        widths = map(maxLen, columns)
        print 'a'
        for row in zip(*columns):
            widthsAndFields = zip(widths, row)
            field0 = "%-*s" % widthsAndFields[0]  # left-justify
            fields = ["%*s" % i for i in widthsAndFields[1:]]  # right-justify
            print 's', field0, ' '.join(fields)
        pad = ' '.join(' ' * i for i in widths[:-1])
        for i in self.prob:
            print 'p', pad, i
        print  # blank line afterwards

def topSeqBeg(maf): return maf.beg[0]
def topSeqEnd(maf): return maf.end[0]
def emptyMaf(): return MafBlock([], [], [], [], [], [], [])

def joinOnFirstItem(x, y):
    if x[0] != y[0]:
        raise ValueError('join fields not equal:\n'+str(x[0])+'\n'+str(y[0]))
    return x + y[1:]

def mafEasyJoin(x, y):
    '''Join two MAF blocks on the top sequence.'''
    xJoin = zip(x.chr, x.beg, x.end, x.strand, x.chrSize, x.seq)
    yJoin = zip(y.chr, y.beg, y.end, y.strand, y.chrSize, y.seq)
    joined = joinOnFirstItem(xJoin, yJoin)
    chr, beg, end, strand, chrSize, seq = zip(*joined)
    prob = x.prob + y.prob
    return MafBlock(chr, beg, end, strand, chrSize, seq, prob)

def countNonGaps(s): return len(s) - s.count('-')

def nthNonGap(s, n):
    '''Get the start position of the n-th non-gap.'''
    for i, x in enumerate(s):
        if x != '-':
            if n == 0: return i
            n -= 1
    raise ValueError('non-gap not found')

def nthLastNonGap(s, n):
    '''Get the end position of the n-th last non-gap.'''
    return len(s) - nthNonGap(s[::-1], n)

def mafSlice(maf, alnBeg, alnEnd):
    '''Return a slice of a MAF block, using coordinates in the alignment.'''
    beg = [b + countNonGaps(s[:alnBeg]) for b, s in zip(maf.beg, maf.seq)]
    end = [e - countNonGaps(s[alnEnd:]) for e, s in zip(maf.end, maf.seq)]
    seq = [i[alnBeg:alnEnd] for i in maf.seq]
    prob = [i[alnBeg:alnEnd] for i in maf.prob]
    return MafBlock(maf.chr, beg, end, maf.strand, maf.chrSize, seq, prob)

def mafSliceTopSeq(maf, newTopSeqBeg, newTopSeqEnd):
    '''Return a slice of a MAF block, using coordinates in the top sequence.'''
    lettersFromBeg = newTopSeqBeg - topSeqBeg(maf)
    lettersFromEnd = topSeqEnd(maf) - newTopSeqEnd
    alnBeg = nthNonGap(maf.seq[0], lettersFromBeg)
    alnEnd = nthLastNonGap(maf.seq[0], lettersFromEnd)
    return mafSlice(maf, alnBeg, alnEnd)

def jumpGaps(sequence, index):
    '''Return the next index of the sequence where there is a non-gap.'''
    nextIndex = index
    while sequence[nextIndex] == '-': nextIndex += 1
    return nextIndex

def gapsToAdd(sequences):
    '''Return new gaps and their positions, needed to align the non-gaps.'''
    gapInfo = [[] for i in sequences]
    gapBeg = [0 for i in sequences]
    try:
        while True:
            gapEnd = [jumpGaps(s, p) for s, p in zip(sequences, gapBeg)]
            gapSize = [e-b for b, e in zip(gapBeg, gapEnd)]
            maxGapSize = max(gapSize)
            for s, e, i in zip(gapSize, gapEnd, gapInfo):
                if s < maxGapSize:
                    newGap = maxGapSize - s
                    i.append((newGap, e))
            gapBeg = [e+1 for e in gapEnd]
    except IndexError: return gapInfo

def chunksAndGaps(s, gapsAndPositions, oneGap):
    '''Yield chunks of "s" interspersed with gaps at given positions.'''
    oldPosition = 0
    for gapLen, position in gapsAndPositions:
        yield s[oldPosition:position]
        yield oneGap * gapLen
        oldPosition = position
    yield s[oldPosition:]

def mafAddGaps(maf, gapsAndPositions):
    '''Add the given gaps at the given positions to a MAF block.'''
    maf.seq = [sum(chunksAndGaps(i, gapsAndPositions, ['-']), [])
               for i in maf.seq]
    maf.prob = [''.join(chunksAndGaps(i, gapsAndPositions, '~'))
                for i in maf.prob]

def mafJoin(mafs):
    '''Intersect and join overlapping MAF blocks.'''
    newTopSeqBeg = max(map(topSeqBeg, mafs))
    newTopSeqEnd = min(map(topSeqEnd, mafs))
    mafs = [mafSliceTopSeq(i, newTopSeqBeg, newTopSeqEnd) for i in mafs]
    topSeqs = [i.seq[0] for i in mafs]
    gapInfo = gapsToAdd(topSeqs)
    for maf, gapsAndPositions in zip(mafs, gapInfo):
        mafAddGaps(maf, gapsAndPositions)
    return reduce(mafEasyJoin, mafs)

def mafInput(lines):
    '''Read lines and yield MAF blocks.'''
    maf = emptyMaf()
    for line in lines:
        if line.isspace():
            if maf: yield maf
            maf = emptyMaf()
        else:
            maf.addLine(line)
    if maf: yield maf

def sortedMafInput(lines):
    '''Read lines and yield MAF blocks, checking that they are in order.'''
    old = emptyMaf()
    for maf in mafInput(lines):
        if maf < old: sys.exit(progName + ": MAF blocks not sorted properly")
        yield maf
        old = maf

def allOverlaps(sequences, beg, end):
    '''Yield all combinations of MAF blocks that overlap in the top genome.'''
    assert beg < end
    if not sequences:
        yield ()
        return
    for i in sequences[0]:
        if topSeqEnd(i) <= beg: continue
        if topSeqBeg(i) >= end: break  # assumes they're sorted by position
        newBeg = max(beg, topSeqBeg(i))
        newEnd = min(end, topSeqEnd(i))
        for j in allOverlaps(sequences[1:], newBeg, newEnd):
            yield (i,) + j

def nextWindow(window, input, referenceMaf):
    '''Yield "relevant" MAF blocks, based on overlap with referenceMaf.'''
    for maf in window:
        if not maf.before(referenceMaf): yield maf
    try:
        while True:
            maf = input.peek()
            if maf.after(referenceMaf): break
            maf = input.next()
            if not maf.before(referenceMaf): yield maf
    except StopIteration: pass

def overlappingMafs(sortedMafInputs):
    '''Yield all combinations of MAF blocks that overlap in the top genome.'''
    if not sortedMafInputs: return
    head, tail = sortedMafInputs[0], sortedMafInputs[1:]
    windows = [[] for t in tail]
    for h in head:  # iterate over MAF blocks in the first input
        windows = [list(nextWindow(w, t, h)) for w, t in zip(windows, tail)]
        for i in allOverlaps(windows, topSeqBeg(h), topSeqEnd(h)):
            yield (h,) + i

op = optparse.OptionParser(usage="%prog sorted-file1.maf sorted-file2.maf ...")
(opts, args) = op.parse_args()

progName = os.path.basename(sys.argv[0])

inputs = [peekable(sortedMafInput(fileinput.input(i))) for i in args]

for mafs in overlappingMafs(inputs):
    mafJoin(mafs).write()
