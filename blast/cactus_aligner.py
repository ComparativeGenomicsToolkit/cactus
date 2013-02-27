#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script for computing alignments for a reconstruction problem.
"""

import os
from optparse import OptionParser

from sonLib.bioio import TempFileTree
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from cactus.blastAlignment.cactus_batch import makeBlastFromOptions
from cactus.blastAlignment.cactus_batch import makeStandardBlastOptions

from cactus.blastAlignment.cactus_alignerTestAligner import MakeBlastsLoader as MakeBlastsTest
    
class MakeBlasts(Target):
    """Take a reconstruction problem and generate the sequences to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, cactusDisk, flowerName, resultsFile, blastOptions, minimumSequenceLength=1, chunkLength=10000000):
        Target.__init__(self)
        self.cactusDisk = cactusDisk
        self.flowerName = flowerName
        self.resultsFile = resultsFile
        self.blastOptions = blastOptions
        self.minimumSequenceLength = minimumSequenceLength
        self.chunkLength = chunkLength
        self.tempFileTree = None
        
    def run(self):
        if self.tempFileTree == None:
            self.tempFileTree = TempFileTree(os.path.join(self.getGlobalTempDir(), "alignmentFiles"))
            headerFile = os.path.join(self.getGlobalTempDir(), "headerFile.txt")
            chunks = list(popenCatch("cactus_makeBlasts %s '%s' %s %s %s" % (getLogLevelString(), self.cactusDisk, self.flowerName, self.minimumSequenceLength, self.chunkLength, headerFile)))
            fn = lambda chunks : self.addChildTarget(Blaster(seqFile, headerFile, chunks, self.tempFileTree.getTempFile(), self.blastOptions))
            for i in xrange(len(chunks)):
                fn((chunks[i],))
                for j in xrange(i+1, len(chunks)):
                    fn((chunks[i], chunks[j]))
            self.setFollowOnTarget(self)
        else:
            catFiles(self.tempFileTree.listFiles(), self.resultsFile)

def catFiles(filesToCat, catFile):
    """Cats a bunch of files into one file. Ensures a no more than MAX_CAT files
    are concatenated at each step.
    """
    MAX_CAT = 25
    system("cat %s > %s" % (" ".join(filesToCat[:MAX_CAT]), catFile))
    filesToCat = filesToCat[MAX_CAT:]
    while len(filesToCat) > 0:
        system("cat %s >> %s" % (" ".join(filesToCat[:MAX_CAT]), catFile))
        filesToCat = filesToCat[MAX_CAT:]
                   
class Blaster(MakeBlasts):
    """Take a reconstruction problem and generate the sequences to be blasted.
    Then setup the follow on blast targets and collation targets.
    """
    def __init__(self, seqFile, headerFile, chunk, outputFile, blastOptions):
        Target.__init__(self, time=0.0099)
        self.seqFile = seqFile
        self.headerFile = headerFile
        self.chunks = chunks
        self.resultsFile = resultsFile
        self.blastOptions = blastOptions
    
    def run(self):
        tempOutputFile = os.path.join(self.getLocalTempDir(), "alignments.cigar")
        def fn(tokens, name):
            return getSequences((self.seqFile, int(tokens[0]), int(tokens[1])), 
                                (self.headerFile,int(tokens[2]), int(tokens[3])),
                                os.path.join(self.getLocalTempDir(), name))
        if len(self.chunks) == 1: #Self blast
            blastCommand = runSelfBlast(fn(self.chunks[0].split(), "1.fa"), tempOutputFile, self.blastOptions.selfBlastString)
        else:
            assert len(self.chunks) == 2 #Not self blast
            blastCommand = makeBlastString(fn(self.chunks[0].split(), "1.fa"), fn(self.chunks[1].split(), "2.fa"), tempOutputFile, self.blastOptions.blastString)
        #Run blast
        open(self.resultsFile, 'w').close() #Just in case the job restarted
        system(command)
        transformCoordinates(tempOutputFile, self.resultsFile)

def makeBlastString(seqFile1, seqFile2, resultsFile, blastString):
    return blastString.replace("CIGARS_FILE", resultsFile).replace("SEQ_FILE_1", seqFile1).replace("SEQ_FILE_2", seqFile2)
    
def runSelfBlast(seqFile, resultsFile, selfBlastString):
    return selfBlastString.replace("CIGARS_FILE", resultsFile).replace("SEQ_FILE", seqFile)
    
def getSequences(seqFile, headerFile, outFile):
    """Create a sequence file containing the sequence file of interest.
    """
    def fn(file):
        fH = open(file[0], 'r')
        fH.seek(file[1])
        string = fH.read(file[2])
        fH.close()
        return string
    fH = open(outFile, 'w')
    for header, sequence in zip(fn(seqFile).split("\n"), fn(headerFile).split("\n")):
        fHO.write(">%s\n%s\n" % (header, sequence))
    fH.close()

def transformCoordinates(inputFile, outputFile):
    """Transforms the coordinates from the internal coordinates 
    to those understood by cactus"""
    system("cactus_convertBlastCoordinates %s %s" % (inputFile, outputFile))
    
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.blastAlignment.cactus_aligner import *
    _test()
    main()
