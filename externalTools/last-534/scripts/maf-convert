#! /usr/bin/env python
# Copyright 2010, 2011, 2013, 2014 Martin C. Frith
# Read MAF-format alignments: write them in other formats.
# Seems to work with Python 2.x, x>=4

# By "MAF" we mean "multiple alignment format" described in the UCSC
# Genome FAQ, not e.g. "MIRA assembly format".

from itertools import *
import sys, os, fileinput, math, operator, optparse, signal, string

def maxlen(s):
    return max(map(len, s))

def joined(things, delimiter):
    return delimiter.join(map(str, things))

identityTranslation = string.maketrans("", "")
def deleted(myString, deleteChars):
    return myString.translate(identityTranslation, deleteChars)

def quantify(iterable, pred=bool):
    """Count how many times the predicate is true."""
    return sum(map(pred, iterable))

class Maf:
    def __init__(self, aLine, sLines, qLines, pLines):
        self.namesAndValues = dict(i.split("=") for i in aLine[1:])

        if not sLines: raise Exception("empty alignment")
        cols = zip(*sLines)
        self.seqNames = cols[1]
        self.alnStarts = map(int, cols[2])
        self.alnSizes = map(int, cols[3])
        self.strands = cols[4]
        self.seqSizes = map(int, cols[5])
        self.alnStrings = cols[6]

        self.qLines = qLines
        self.pLines = pLines

def dieUnlessPairwise(maf):
    if len(maf.alnStrings) != 2:
        raise Exception("pairwise alignments only, please")

def insertionCounts(alnString):
    gaps = alnString.count("-")
    forwardFrameshifts = alnString.count("\\")
    reverseFrameshifts = alnString.count("/")
    letters = len(alnString) - gaps - forwardFrameshifts - reverseFrameshifts
    return letters, forwardFrameshifts, reverseFrameshifts

def coordinatesPerLetter(alnString, alnSize):
    letters, forwardShifts, reverseShifts = insertionCounts(alnString)
    if forwardShifts or reverseShifts or letters < alnSize: return 3
    else: return 1

def mafLetterSizes(maf):
    return map(coordinatesPerLetter, maf.alnStrings, maf.alnSizes)

def alnLetterSizes(sLines):
    return [coordinatesPerLetter(i[6], int(i[3])) for i in sLines]

def isTranslated(sLines):
    for i in sLines:
        alnString = i[6]
        if "\\" in alnString or "/" in alnString: return True
        if len(alnString) - alnString.count("-") < int(i[3]): return True
    return False

def insertionSize(alnString, letterSize):
    """Get the length of sequence included in the alnString."""
    letters, forwardShifts, reverseShifts = insertionCounts(alnString)
    return letters * letterSize + forwardShifts - reverseShifts

def symbolSize(symbol, letterSize):
    if symbol == "-": return 0
    if symbol == "\\": return 1
    if symbol == "/": return -1
    return letterSize

def isMatch(alignmentColumn):
    # No special treatment of ambiguous bases/residues: same as NCBI BLAST.
    first = alignmentColumn[0].upper()
    for i in alignmentColumn[1:]:
        if i.upper() != first: return False
    return True

def isGapless(alignmentColumn):
    return "-" not in alignmentColumn

def matchAndInsertSizes(alignmentColumns, letterSizes):
    """Get sizes of gapless blocks, and of the inserts between them."""
    for k, v in groupby(alignmentColumns, isGapless):
        if k:
            sizeOfGaplessBlock = len(list(v))
            yield str(sizeOfGaplessBlock)
        else:
            blockRows = izip(*v)
            blockRows = imap(''.join, blockRows)
            insertSizes = imap(insertionSize, blockRows, letterSizes)
            yield joined(insertSizes, ":")

##### Routines for converting to AXT format: #####

axtCounter = count()

def writeAxt(maf):
    if maf.strands[0] != "+":
        raise Exception("for AXT, the 1st strand in each alignment must be +")

    # Convert to AXT's 1-based coordinates:
    alnStarts = imap(operator.add, maf.alnStarts, repeat(1))
    alnEnds = imap(operator.add, maf.alnStarts, maf.alnSizes)

    rows = zip(maf.seqNames, alnStarts, alnEnds, maf.strands)
    head, tail = rows[0], rows[1:]

    outWords = []
    outWords.append(axtCounter.next())
    outWords.extend(head[:3])
    outWords.extend(chain(*tail))
    outWords.append(maf.namesAndValues["score"])

    print joined(outWords, " ")
    for i in maf.alnStrings: print i
    print  # print a blank line at the end

##### Routines for converting to tabular format: #####

def writeTab(maf):
    aLine, sLines, qLines, pLines = maf

    score = "0"
    expect = None
    mismap = None
    for i in aLine:
        if   i.startswith("score="): score = i[6:]
        elif i.startswith("expect="): expect = i[7:]
        elif i.startswith("mismap="): mismap = i[7:]

    outWords = []
    outWords.append(score)

    for i in sLines: outWords.extend(i[1:6])

    alnStrings = [i[6] for i in sLines]
    alignmentColumns = zip(*alnStrings)
    letterSizes = alnLetterSizes(sLines)
    gapStrings = matchAndInsertSizes(alignmentColumns, letterSizes)
    gapWord = ",".join(gapStrings)
    outWords.append(gapWord)

    if expect: outWords.append(expect)
    if mismap: outWords.append(mismap)

    print "\t".join(outWords)

##### Routines for converting to PSL format: #####

def pslBlocks(alignmentColumns, alnStarts, letterSizes):
    """Get sizes and start coordinates of gapless blocks in an alignment."""
    start1, start2 = alnStarts
    letterSize1, letterSize2 = letterSizes
    size = 0
    for x, y in alignmentColumns:
        if x == "-":
            if size:
                yield size, start1, start2
                start1 += size * letterSize1
                start2 += size * letterSize2
                size = 0
            start2 += symbolSize(y, letterSize2)
        elif y == "-":
            if size:
                yield size, start1, start2
                start1 += size * letterSize1
                start2 += size * letterSize2
                size = 0
            start1 += symbolSize(x, letterSize1)
        else:
            size += 1
    if size: yield size, start1, start2

def pslCommaString(things):
    # UCSC software seems to prefer a trailing comma
    return joined(things, ",") + ","

def gapRunCount(letters):
    """Get the number of runs of gap characters."""
    uniqLetters = map(operator.itemgetter(0), groupby(letters))
    return uniqLetters.count("-")

def pslEndpoints(alnStart, alnSize, strand, seqSize):
    alnEnd = alnStart + alnSize
    if strand == "+": return alnStart, alnEnd
    else: return seqSize - alnEnd, seqSize - alnStart

def caseSensitivePairwiseMatchCounts(columns, isProtein):
    # repMatches is always zero
    # for proteins, nCount is always zero, because that's what BLATv34 does
    standardBases = "ACGTU"
    matches = mismatches = repMatches = nCount = 0
    for i in columns:
        if "-" in i: continue
        x, y = i
        if x in standardBases and y in standardBases or isProtein:
            if x == y: matches += 1
            else: mismatches += 1
        else: nCount += 1
    return matches, mismatches, repMatches, nCount

def writePsl(maf, isProtein):
    dieUnlessPairwise(maf)

    alnStrings = map(str.upper, maf.alnStrings)
    alignmentColumns = zip(*alnStrings)
    letterSizes = mafLetterSizes(maf)

    matchCounts = caseSensitivePairwiseMatchCounts(alignmentColumns, isProtein)
    matches, mismatches, repMatches, nCount = matchCounts
    numGaplessColumns = sum(matchCounts)

    qNumInsert = gapRunCount(maf.alnStrings[0])
    qBaseInsert = maf.alnSizes[1] - numGaplessColumns * letterSizes[1]

    tNumInsert = gapRunCount(maf.alnStrings[1])
    tBaseInsert = maf.alnSizes[0] - numGaplessColumns * letterSizes[0]

    strand = maf.strands[1]
    if max(letterSizes) > 1:
        strand += maf.strands[0]
    elif maf.strands[0] != "+":
        raise Exception("for non-translated PSL, the 1st strand in each alignment must be +")

    tName, qName = maf.seqNames
    tSeqSize, qSeqSize = maf.seqSizes

    rows = zip(maf.alnStarts, maf.alnSizes, maf.strands, maf.seqSizes)
    tStart, tEnd = pslEndpoints(*rows[0])
    qStart, qEnd = pslEndpoints(*rows[1])

    blocks = list(pslBlocks(alignmentColumns, maf.alnStarts, letterSizes))
    blockCount = len(blocks)
    blockSizes, tStarts, qStarts = map(pslCommaString, zip(*blocks))

    outWords = (matches, mismatches, repMatches, nCount,
                qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand,
                qName, qSeqSize, qStart, qEnd, tName, tSeqSize, tStart, tEnd,
                blockCount, blockSizes, qStarts, tStarts)
    print joined(outWords, "\t")

##### Routines for converting to SAM format: #####

def readGroupId(readGroupItems):
    for i in readGroupItems:
        if i.startswith("ID:"):
            return i[3:]
    raise Exception("readgroup must include ID")

def readSequenceLengths(lines):
    """Read name & length of topmost sequence in each maf block."""
    sequenceLengths = {}  # an OrderedDict might be nice
    isSearching = True
    for line in lines:
        if line.isspace(): isSearching = True
        elif isSearching:
            w = line.split()
            if w[0] == "s":
                seqName, seqSize = w[1], w[5]
                sequenceLengths[seqName] = seqSize
                isSearching = False
    return sequenceLengths

def naturalSortKey(s):
    """Return a key that sorts strings in "natural" order."""
    return [(str, int)[k]("".join(v)) for k, v in groupby(s, str.isdigit)]

def karyotypicSortKey(s):
    """Attempt to sort chromosomes in GATK's ridiculous order."""
    if s == "chrM": return []
    if s == "MT": return ["~"]
    return naturalSortKey(s)

def writeSamHeader(sequenceLengths, dictFile, readGroupItems):
    print "@HD\tVN:1.3\tSO:unsorted"
    for k in sorted(sequenceLengths, key=karyotypicSortKey):
        print "@SQ\tSN:%s\tLN:%s" % (k, sequenceLengths[k])
    if dictFile:
        for i in fileinput.input(dictFile):
            if i.startswith("@SQ"): print i,
            elif not i.startswith("@"): break
    if readGroupItems:
        print "@RG\t" + "\t".join(readGroupItems)

mapqMissing = "255"
mapqMaximum = "254"
mapqMaximumNum = float(mapqMaximum)

def mapqFromProb(probString):
    try: p = float(probString)
    except ValueError: raise Exception("bad probability: " + probString)
    if p < 0 or p > 1: raise Exception("bad probability: " + probString)
    if p == 0: return mapqMaximum
    phred = -10 * math.log(p, 10)
    if phred >= mapqMaximumNum: return mapqMaximum
    return str(int(round(phred)))

def cigarParts(beg, alignmentColumns, end):
    if beg: yield str(beg) + "H"

    # (doesn't handle translated alignments)
    isActive = True
    for x, y in alignmentColumns: break
    else: isActive = False
    while isActive:
        size = 1
        if x == "-":
            for x, y in alignmentColumns:
                if x != "-": break
                size += 1
            else: isActive = False
            yield str(size) + "I"
        elif y == "-":
            for x, y in alignmentColumns:
                if y != "-": break
                size += 1
            else: isActive = False
            yield str(size) + "D"
        else:
            for x, y in alignmentColumns:
                if x == "-" or y == "-": break
                size += 1
            else: isActive = False
            yield str(size) + "M"

    if end: yield str(end) + "H"

def writeSam(maf, rg):
    aLine, sLines, qLines, pLines = maf

    if isTranslated(sLines):
        raise Exception("this looks like translated DNA - can't convert to SAM format")

    score = None
    evalue = None
    mapq = mapqMissing
    for i in aLine:
        if i.startswith("score="):
            v = i[6:]
            if v.isdigit(): score = "AS:i:" + v  # it must be an integer
        elif i.startswith("expect="):
            evalue = "EV:Z:" + i[7:]
        elif i.startswith("mismap="):
            mapq = mapqFromProb(i[7:])

    head, tail = sLines[0], sLines[1:]

    s, rName, rStart, rAlnSize, rStrand, rSeqSize, rAlnString = head
    if rStrand != "+":
        raise Exception("for SAM, the 1st strand in each alignment must be +")
    pos = str(int(rStart) + 1)  # convert to 1-based coordinate

    for s, qName, qStart, qAlnSize, qStrand, qSeqSize, qAlnString in tail:
        alignmentColumns = zip(rAlnString.upper(), qAlnString.upper())

        qStart = int(qStart)
        qRevStart = int(qSeqSize) - qStart - int(qAlnSize)
        cigar = "".join(cigarParts(qStart, iter(alignmentColumns), qRevStart))

        seq = deleted(qAlnString, "-")

        qual = "*"
        if qLines:
            q, qualityName, qualityString = qLines[0]
            # don't try to handle multiple alignments for now (YAGNI):
            if len(qLines) > 1 or len(tail) > 1 or qualityName != qName:
                raise Exception("can't interpret the quality data")
            qual = ''.join(j for i, j in izip(qAlnString, qualityString)
                           if i != "-")

        # It's hard to get all the pair info, so this is very
        # incomplete, but hopefully good enough.
        # I'm not sure whether to add 2 and/or 8 to flag.
        if qName.endswith("/1"):
            qName = qName[:-2]
            if qStrand == "+": flag = "99"  # 1 + 2 + 32 + 64
            else:              flag = "83"  # 1 + 2 + 16 + 64
        elif qName.endswith("/2"):
            qName = qName[:-2]
            if qStrand == "+": flag = "163"  # 1 + 2 + 32 + 128
            else:              flag = "147"  # 1 + 2 + 16 + 128
        else:
            if qStrand == "+": flag = "0"
            else:              flag = "16"

        editDistance = sum(1 for x, y in alignmentColumns if x != y)
        # no special treatment of ambiguous bases: might be a minor bug
        editDistance = "NM:i:" + str(editDistance)

        outWords = [qName, flag, rName, pos, mapq, cigar, "*\t0\t0", seq, qual]
        outWords.append(editDistance)
        if score: outWords.append(score)
        if evalue: outWords.append(evalue)
        if rg: outWords.append(rg)
        print "\t".join(outWords)

##### Routines for converting to BLAST-like format: #####

def pairwiseMatchSymbol(alignmentColumn):
    if isMatch(alignmentColumn): return "|"
    else: return " "

def strandText(strand):
    if strand == "+": return "Plus"
    else: return "Minus"

def blastCoordinate(oneBasedCoordinate, strand, seqSize):
    if strand == "-":
        oneBasedCoordinate = seqSize - oneBasedCoordinate + 1
    return str(oneBasedCoordinate)

def chunker(things, chunkSize):
    for i in range(0, len(things), chunkSize):
        yield things[i:i+chunkSize]

def blastChunker(maf, lineSize, alignmentColumns):
    letterSizes = mafLetterSizes(maf)
    coords = maf.alnStarts
    for chunkCols in chunker(alignmentColumns, lineSize):
        chunkRows = zip(*chunkCols)
        chunkRows = map(''.join, chunkRows)
        starts = [i + 1 for i in coords]  # change to 1-based coordinates
        starts = map(blastCoordinate, starts, maf.strands, maf.seqSizes)
        increments = map(insertionSize, chunkRows, letterSizes)
        coords = map(operator.add, coords, increments)
        ends = map(blastCoordinate, coords, maf.strands, maf.seqSizes)
        yield chunkCols, chunkRows, starts, ends

def writeBlast(maf, lineSize, oldQueryName):
    dieUnlessPairwise(maf)

    if maf.seqNames[1] != oldQueryName:
        print "Query= %s" % maf.seqNames[1]
        print "         (%s letters)" % maf.seqSizes[1]
        print

    print ">%s" % maf.seqNames[0]
    print "          Length = %s" % maf.seqSizes[0]
    print

    scoreLine = " Score = %s" % maf.namesAndValues["score"]
    try: scoreLine += ", Expect = %s" % maf.namesAndValues["expect"]
    except KeyError: pass
    print scoreLine

    alignmentColumns = zip(*maf.alnStrings)

    alnSize = len(alignmentColumns)
    matches = quantify(alignmentColumns, isMatch)
    matchPercent = 100 * matches // alnSize  # round down, like BLAST
    identLine = " Identities = %s/%s (%s%%)" % (matches, alnSize, matchPercent)
    gaps = alnSize - quantify(alignmentColumns, isGapless)
    gapPercent = 100 * gaps // alnSize  # round down, like BLAST
    if gaps: identLine += ", Gaps = %s/%s (%s%%)" % (gaps, alnSize, gapPercent)
    print identLine

    strands = map(strandText, maf.strands)
    print " Strand = %s / %s" % (strands[1], strands[0])
    print

    for chunk in blastChunker(maf, lineSize, alignmentColumns):
        cols, rows, starts, ends = chunk
        startWidth = maxlen(starts)
        matchSymbols = map(pairwiseMatchSymbol, cols)
        matchSymbols = ''.join(matchSymbols)
        print "Query: %-*s %s %s" % (startWidth, starts[1], rows[1], ends[1])
        print "       %-*s %s"    % (startWidth, " ", matchSymbols)
        print "Sbjct: %-*s %s %s" % (startWidth, starts[0], rows[0], ends[0])
        print

##### Routines for converting to HTML format: #####

def writeHtmlHeader():
    print '''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
 "http://www.w3.org/TR/html4/strict.dtd">
<html lang="en"><head>
<meta http-equiv="Content-type" content="text/html; charset=UTF-8">
<title>Reliable Alignments</title>
<style type="text/css">
/* Try to force monospace, working around browser insanity: */
pre {font-family: "Courier New", monospace, serif; font-size: 0.8125em}
.a {background-color: #3333FF}
.b {background-color: #9933FF}
.c {background-color: #FF66CC}
.d {background-color: #FF3333}
.e {background-color: #FF9933}
.f {background-color: #FFFF00}
.key {display:inline; margin-right:2em}
</style>
</head><body>

<div style="line-height:1">
<pre class="key"><span class="a">  </span> prob &gt; 0.999</pre>
<pre class="key"><span class="b">  </span> prob &gt; 0.99 </pre>
<pre class="key"><span class="c">  </span> prob &gt; 0.95 </pre>
<pre class="key"><span class="d">  </span> prob &gt; 0.9  </pre>
<pre class="key"><span class="e">  </span> prob &gt; 0.5  </pre>
<pre class="key"><span class="f">  </span> prob &le; 0.5  </pre>
</div>
'''

def probabilityClass(probabilityColumn):
    p = ord(min(probabilityColumn)) - 33
    if   p >= 30: return 'a'
    elif p >= 20: return 'b'
    elif p >= 13: return 'c'
    elif p >= 10: return 'd'
    elif p >=  3: return 'e'
    else: return 'f'

def identicalRuns(s):
    """Yield (item, start, end) for each run of identical items in s."""
    beg = 0
    for k, v in groupby(s):
        end = beg + len(list(v))
        yield k, beg, end
        beg = end

def htmlSpan(text, classRun):
    key, beg, end = classRun
    textbit = text[beg:end]
    if key: return '<span class="%s">%s</span>' % (key, textbit)
    else: return textbit

def multipleMatchSymbol(alignmentColumn):
    if isMatch(alignmentColumn): return "*"
    else: return " "

def writeHtml(maf, lineSize):
    scoreLine = "Alignment"
    try:
        scoreLine += " score=%s" % maf.namesAndValues["score"]
        scoreLine += ", expect=%s" % maf.namesAndValues["expect"]
    except KeyError: pass
    print "<h3>%s:</h3>" % scoreLine

    if maf.pLines:
        probRows = [i[1] for i in maf.pLines]
        probCols = izip(*probRows)
        classes = imap(probabilityClass, probCols)
    else:
        classes = repeat(None)

    nameWidth = maxlen(maf.seqNames)
    alignmentColumns = zip(*maf.alnStrings)

    print '<pre>'
    for chunk in blastChunker(maf, lineSize, alignmentColumns):
        cols, rows, starts, ends = chunk
        startWidth = maxlen(starts)
        endWidth = maxlen(ends)
        matchSymbols = map(multipleMatchSymbol, cols)
        matchSymbols = ''.join(matchSymbols)
        classChunk = islice(classes, lineSize)
        classRuns = list(identicalRuns(classChunk))
        for n, s, r, e in zip(maf.seqNames, starts, rows, ends):
            spans = [htmlSpan(r, i) for i in classRuns]
            spans = ''.join(spans)
            formatParams = nameWidth, n, startWidth, s, spans, endWidth, e
            print '%-*s %*s %s %*s' % formatParams
        print ' ' * nameWidth, ' ' * startWidth, matchSymbols
        print
    print '</pre>'

##### Routines for reading MAF format: #####

def mafInput(lines, isKeepComments):
    a = []
    s = []
    q = []
    p = []
    for i in lines:
        w = i.split()
        if not w:
            if a: yield a, s, q, p
            a = []
            s = []
            q = []
            p = []
        elif w[0] == "a":
            a = w
        elif w[0] == "s":
            s.append(w)
        elif w[0] == "q":
            q.append(w)
        elif w[0] == "p":
            p.append(w)
        elif i[0] == "#" and isKeepComments:
            print i,
    if a: yield a, s, q, p

def isFormat(myString, myFormat):
    return myFormat.startswith(myString)

def mafConvert(opts, args):
    format = args[0].lower()
    if isFormat(format, "sam"):
        if opts.dictionary: d = readSequenceLengths(fileinput.input(args[1:]))
        else: d = {}
        if opts.readgroup:
            readGroupItems = opts.readgroup.split()
            rg = "RG:Z:" + readGroupId(readGroupItems)
        else:
            readGroupItems = []
            rg = ""
        if not opts.noheader: writeSamHeader(d, opts.dictfile, readGroupItems)
    inputLines = fileinput.input(args[1:])
    if isFormat(format, "html") and not opts.noheader: writeHtmlHeader()
    isKeepCommentLines = isFormat(format, "tabular") and not opts.noheader
    oldQueryName = ""
    for maf in mafInput(inputLines, isKeepCommentLines):
        if   isFormat(format, "axt"): writeAxt(Maf(*maf))
        elif isFormat(format, "blast"):
            maf = Maf(*maf)
            writeBlast(maf, opts.linesize, oldQueryName)
            oldQueryName = maf.seqNames[1]
        elif isFormat(format, "html"): writeHtml(Maf(*maf), opts.linesize)
        elif isFormat(format, "psl"): writePsl(Maf(*maf), opts.protein)
        elif isFormat(format, "sam"): writeSam(maf, rg)
        elif isFormat(format, "tabular"): writeTab(maf)
        else: raise Exception("unknown format: " + format)
    if isFormat(format, "html") and not opts.noheader: print "</body></html>"

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = """
  %prog --help
  %prog axt mafFile(s)
  %prog blast mafFile(s)
  %prog html mafFile(s)
  %prog psl mafFile(s)
  %prog sam mafFile(s)
  %prog tab mafFile(s)"""

    description = "Read MAF-format alignments & write them in another format."

    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-p", "--protein", action="store_true",
                  help="assume protein alignments, for psl match counts")
    op.add_option("-n", "--noheader", action="store_true",
                  help="omit any header lines from the output")
    op.add_option("-d", "--dictionary", action="store_true",
                  help="include dictionary of sequence lengths in sam format")
    op.add_option("-f", "--dictfile",
                  help="get sequence dictionary from DICTFILE")
    op.add_option("-r", "--readgroup",
                  help="read group info for sam format")
    op.add_option("-l", "--linesize", type="int", default=60, #metavar="CHARS",
                  help="line length for blast and html formats (default: %default)")
    (opts, args) = op.parse_args()
    if opts.linesize <= 0: op.error("option -l: should be >= 1")
    if opts.dictionary and opts.dictfile: op.error("can't use both -d and -f")
    if len(args) < 1: op.error("I need a format-name and some MAF alignments")
    if opts.dictionary and (len(args) == 1 or "-" in args[1:]):
        op.error("need file (not pipe) with option -d")

    try: mafConvert(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
