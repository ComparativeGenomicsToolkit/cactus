#!/usr/bin/env python3
from collections import defaultdict
from operator import itemgetter

def windowFilter(windowSize, threshold, blockDict, seqLengths):
    if windowSize == 1 and threshold == 1:
        # Don't need to do expensive window-filtering
        return blockDict
    ret = defaultdict(list)
    for seq, blocks in list(blockDict.items()):
        curBlock = 0
        inRegion = False
        regionStart = 0
        for i in range(seqLengths[seq]):
            score = 0
            while curBlock < len(blocks) and blocks[curBlock][1] < i:
                curBlock += 1
            for blockNum in range(curBlock, len(blocks)):
                block = blocks[blockNum]
                if block[0] > i + windowSize:
                    break
                size = min(block[1], i + windowSize) - max(i, block[0])
                if block[2] >= 1:
                    score += size
            score /= float(windowSize)
            if score >= threshold and not inRegion:
                regionStart = i
                inRegion = True
            elif score < threshold and inRegion:
                ret[seq].append((regionStart, i + windowSize - 1))
                inRegion = False
    return ret

def uniquifyBlocks(blocksDict, mergeDistance):
    """Take list of blocks and return sorted list of non-overlapping and
    blocks (merging blocks that are mergeDistance or less apart)."""
    ret = defaultdict(list)
    for chr, blocks in list(blocksDict.items()):
        blocks = sorted(blocks, key=itemgetter(0))
        newBlocks = []
        prevBlock = None
        for block in blocks:
            if prevBlock is None:
                prevBlock = block
            else:
                if prevBlock[1] < block[0] - mergeDistance:
                    newBlocks.append(prevBlock)
                    prevBlock = block
                else:
                    prevBlock = (prevBlock[0], block[1])
        if prevBlock is not None:
            newBlocks.append(prevBlock)
        ret[chr] = newBlocks
    return ret

def getSeparateBedBlocks(bedFile, depth=1):
    """Get dict of sequence -> (start, stop, score) regions from bed file,
    counting BED12 exons separately, filtering for blocks that have
    score > depth.
    """
    ret = defaultdict(list)
    for line in bedFile:
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        fields = line.split('\t')
        chr = fields[0]
        start = int(fields[1])
        stop = int(fields[2])
        score = int(fields[4])
        if len(fields) <= 9:
            if score >= depth:
                ret[chr].append((start, stop, score))
        else:
            assert(len(fields) == 12)
            blockSizes = list(map(int, [x for x in fields[10].split(',') if x != '']))
            blockStarts = list(map(int, [x for x in fields[11].split(',') if x != '']))
            for blockStart, blockSize in zip(blockSizes, blockStarts):
                nonRelativeBlockStart = start + blockStart
                nonRelativeBlockEnd = nonRelativeBlockStart + blockSize
                if score >= depth:
                    ret[chr].append((nonRelativeBlockStart, nonRelativeBlockEnd,
                                     score))
    return ret

def getSeqLengths(fastaFile):
    """Get a dict which maps header -> sequence size."""
    ret = defaultdict(int)
    chr = ""
    for line in fastaFile:
        line = line.strip()
        if len(line) == 0:
            # Blank line
            continue
        if line[0] == '>':
            chr = line.split()[0][1:]
        else:
            ret[chr] += len(line)
    return ret

def complementBlocks(blocksDict, seqLengths):
    """Complement a sorted block-dict."""
    ret = defaultdict(list)
    for chr, blocks in list(blocksDict.items()):
        len = seqLengths[chr]
        start = 0
        for block in blocks:
            ret[chr].append((start, block[0]))
            start = block[1]
        if start != len:
            ret[chr].append((start, len))
    # Add in blocks for the sequences that aren't covered at all.
    for chr, len in list(seqLengths.items()):
        if chr not in ret: # This still works with defaultdicts
            ret[chr].append((0, len))
    return ret

def printTrimmedSeq(header, seq, blocks, outFile):
    for block in blocks:
        outFile.write(">%s|%d\n" % (header, block[0]))
        outFile.write(seq[block[0]:block[1]])
        outFile.write("\n")

def printTrimmedFasta(fastaFile, toTrim, outFile):
    header = None
    seq = None
    for line in fastaFile:
        line = line.strip()
        if len(line) == 0:
            # Blank line
            continue
        if line[0] == '>':
            if seq is not None:
                printTrimmedSeq(header, seq, toTrim[header], outFile)
            seq = ""
            header = line[1:].split()[0]
            continue
        seq += line
    if seq is not None:
        printTrimmedSeq(header, seq, toTrim[header], outFile)

def trimSequences(fastaPath, bedPath, outputPathOrFile, flanking=0, minSize=0,
                  windowSize=10, threshold=0.8, depth=1, complement=False):
    fastaFile = open(fastaPath)
    seqLengths = getSeqLengths(fastaFile)
    with open(bedPath) as bedFile:
        toTrim = windowFilter(windowSize, threshold,
                              getSeparateBedBlocks(bedFile, depth), seqLengths)
    if complement:
        toTrim = complementBlocks(toTrim, seqLengths)
    toTrim = uniquifyBlocks(toTrim, 2*flanking)
    # filter based on size
    toTrim.update((k, [x for x in v if (x[1] - x[0]) >= minSize])
                  for k, v in list(toTrim.items()))
    # extend blocks to include flanking regions
    toTrim.update((k, [(max(x[0] - flanking, 0),
                                     min(x[1] + flanking, seqLengths[k])) for x in v])
                  for k, v in list(toTrim.items()))

    fastaFile.seek(0)
    try:
        outputPathOrFile.write('')
        outputFile = outputPathOrFile
    except:
        # Not a file
        outputFile = open(outputPathOrFile, 'w')
    printTrimmedFasta(fastaFile, toTrim, outputFile)
