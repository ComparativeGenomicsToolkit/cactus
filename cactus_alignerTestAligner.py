#!/usr/bin/env python

"""Script for testing cactus_core recursively with randomly generated alignments and sequences.
"""

import random
from sonLib.bioio import logger

from sonLib.bioio import PairwiseAlignment
from sonLib.bioio import cigarWrite
from sonLib.bioio import getRandomOperationList
from sonLib.bioio import fastaRead

from workflow.jobTree.scriptTree.target import Target

class MakeBlastsLoader:
    """Used to get around strange pickling behavious todo with being unable to pickle functions.
    """
    def makeBlastOptions(self, job, sequences, resultsFile):
        return MakeBlasts(job, sequences, resultsFile)
    
class MakeBlasts(Target):
    """Fills the input file with some random alignments.
    """
    def __init__(self, job, sequences, resultsFile):
        self.sequences = sequences
        self.resultsFile = resultsFile
        Target.__init__(self, job, None)

    def run(self, job):
        ##########################################
        #Stuff to make random alignments
        ##########################################
        sequences = []
        for sequenceFile in self.sequences:
            fileHandle = open(sequenceFile, 'r')
            for fastaHeader, seq in fastaRead(fileHandle):
                sequences.append((fastaHeader, len(seq)))
            fileHandle.close()
        
        randomAlignmentNo = random.choice(xrange(5)) * len(sequences)
        
        def getInterval(length):
            i = random.choice(xrange(length))
            j = random.choice(xrange(length))
            if i <= j:
                return i, j, True
            return i, j, False
        
        ##########################################
        #Loop to make random alignments.
        ##########################################
        
        fileHandle = open(self.resultsFile, 'w')
        for i in xrange(randomAlignmentNo):
            contig1, length1 = random.choice(sequences)
            start1, end1, strand1 = getInterval(length1)
            
            contig2, length2 = random.choice(sequences)
            start2, end2, strand2 = getInterval(length2)
            
            operationList = getRandomOperationList(abs(end1-start1), abs(end2-start2))
            
            pairwiseAlignment = \
                PairwiseAlignment(contig1, start1, end1, strand1, \
                                  contig2, start2, end2, strand2, \
                                  random.random(), operationList)
    
            cigarWrite(fileHandle, pairwiseAlignment)
        fileHandle.close()
        
        logger.info("Made the random alignments")
