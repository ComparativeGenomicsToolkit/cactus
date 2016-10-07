#!/usr/bin/env python

#Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import sys
import os
import re
import math

#########################################################
#########################################################
#########################################################
#misc functions
#########################################################
#########################################################
#########################################################

def sonTraceRootPath():
    """
    function for finding external location
    """
    import sonLib.bioio
    i = os.path.abspath(sonLib.bioio.__file__)
    return os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def linOriginRegression(points):
    """
    computes a linear regression starting at zero
    """
    j = sum([ i[0] for i in points ])
    k = sum([ i[1] for i in points ])
    if j != 0:
        return k/j, j, k
    return 1, j, k

def close(i, j, tolerance):
    """
    check two float values are within a bound of one another
    """
    return i <= j + tolerance and i >= j - tolerance

def getPositiveCoordinateRangeOverlap(uStart1, uEnd1, uStart2, uEnd2):
    """@return: If the two coordinate ranges overlap on the same strand
    returns the overlap range. If no overlap it returns None.
    """
    if uEnd1 < uStart2 or uEnd2 < uStart1:
        return None
    l = [ uStart1, uEnd1, uStart2, uEnd2 ]
    l.sort()
    return l[1], l[2]

def getCoordinateRangeOverlap(uStart1, uEnd1, uStart2, uEnd2):
    """@return: If the two coordinate ranges overlap (even on opposite strands)
    returns the overlap, on the same strand as the first range of coordinates.
    If no overlap it returns None.
    """
    #deal with negative strands
    if uStart1 <= 0:
        if uStart2 <= 0:
            rangeOverlap = getPositiveCoordinateRangeOverlap(-uEnd1, -uStart1, 
                                                             -uEnd2, -uStart2)
        else:
            rangeOverlap = getPositiveCoordinateRangeOverlap(-uEnd1, -uStart1, 
                                                             uStart2, uEnd2)
        if rangeOverlap is not None:
            x, y = rangeOverlap
            return -y, -x
        return None
    if uStart2 <=0:
        return getPositiveCoordinateRangeOverlap(uStart1, uEnd1, 
                                                 -uEnd2, -uStart2)
    return getPositiveCoordinateRangeOverlap(uStart1, uEnd1, uStart2, uEnd2)

def sortAlignments(alignments):
    def cmpFn(a1, a2): #Sort in descending order
        if a1.score > a2.score:
            return 1
        elif a1.score < a2.score:
            return -1
        return 0
    alignments.sort(cmpFn)
    return alignments

def filterOverlappingAlignments(alignments):
    """Filter alignments to be non-overlapping.
    """
    l = []
    alignments = alignments[:]
    sortAlignments(alignments)
    alignments.reverse()
    for pA1 in alignments:
        for pA2 in l:
            if pA1.contig1 == pA2.contig1 and getPositiveCoordinateRangeOverlap(pA1.start1+1, pA1.end1, pA2.start1+1, pA2.end1) is not None: #One offset, inclusive coordinates
                break
            if pA1.contig2 == pA2.contig2 and getPositiveCoordinateRangeOverlap(pA1.start2+1, pA1.end2, pA2.start2+1, pA2.end2) is not None: #One offset, inclusive coordinates
                break
            if pA1.contig2 == pA2.contig1 and getPositiveCoordinateRangeOverlap(pA1.start2+1, pA1.end2, pA2.start1+1, pA2.end1) is not None: #One offset, inclusive coordinates
                break
            if pA1.contig1 == pA2.contig2 and getPositiveCoordinateRangeOverlap(pA1.start1+1, pA1.end1, pA2.start2+1, pA2.end2) is not None: #One offset, inclusive coordinates
                break
        else:
            l.append(pA1)
    l.reverse()
    return l

#########################################################
#########################################################
#########################################################
#continuous math functions
#########################################################
#########################################################
#########################################################

LOG_ZERO_PROB = -1e30000
LOG_ONE_PROB = 0.0
ZERO_PROB = 0.0

def logAdd(x, y):
    if x < y:
        if x <= LOG_ZERO_PROB:
            return y
        return math.log(math.exp(x - y) + 1) + y
    if y <= LOG_ZERO_PROB:
        return x
    return math.log(math.exp(y - x) + 1) + x

#########################################################
#########################################################
#########################################################
#bio functions
#########################################################
#########################################################
#########################################################

dNAMap_reverseComp_IUPAC = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'n':'n' }

dNAMap_reverseComp_Int = { 0:3, 1:2, 2:1, 3:0, 4:4 }

def reverseComplement(seq, rCM=dNAMap_reverseComp_Int):
    seq.reverse()
    i = [ rCM[j] for j in seq ]
    seq.reverse()
    return i

def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
