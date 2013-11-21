import os
import sys
import random
import time

from sonLib.bioio import system
from sonLib.bioio import logger
from cactus.blast.cactus_blastTest import *

"""Tests if changing parameters of lastz creates results similar to the desired default, checking at different evolutionary distances.
"""

"""
def fn(tree, pairs):
    if tree == None:
        return []
    l = fn(tree.left, pairs)
    l2 = fn(tree.right, pairs)
    for i, d1 in l:
        for j, d2 in l2:
            pairs.append((i, j, d1+d2))
    if tree.iD != None:
        l2.append((tree.iD, 0.0))
    return [ (i, d+tree.distance) for i, d in  l + l2 ]
"""

encodeRegion = "ENm001"
speciesPairs = [('human', 'chimp', 0.016696), ('human', 'baboon', 0.07682800000000001), 
                ('human', 'macaque', 0.082821), ('human', 'marmoset', 0.125919), ('human', 'galago', 0.277323), 
                ('human', 'elephant', 0.347192), ('human', 'rfbat', 0.351039), ('human', 'armadillo', 0.35440400000000005), 
                ('human', 'rabbit', 0.367195), ('human', 'dog', 0.37041599999999997), ('human', 'cow', 0.377763), 
                ('human', 'mouse', 0.474132), ('human', 'rat', 0.48255800000000004), ('human', 'hedgehog', 0.48416000000000003), 
                ('human', 'tenrec', 0.49516499999999997), ('human', 'shrew', 0.521589), ("human", "monodelphis", 0.6) ]

settings = ["--step=2 --ambiguous=iupac,100 --ydrop=3000",
            "--step=4 --ambiguous=iupac,100 --ydrop=3000",
            "--step=8 --ambiguous=iupac,100 --ydrop=3000",
            "--step=16 --ambiguous=iupac,100 --ydrop=3000" ]

#Other species to try "rat", "monodelphis", "macaque", "chimp"
encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
regionPath = os.path.join(encodePath, encodeRegion)
tempDir = getTempDirectory(os.getcwd())
tempOutputFile = os.path.join(tempDir, "results1.txt")
tempOutputFile2 = os.path.join(tempDir, "results2.txt")

for species1, species2, distance in speciesPairs:
    seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
    seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
    
    #Run the random
    ##Record time to run
    baseRuntime = runNaiveBlast([ seqFile1, seqFile2 ], tempOutputFile, tempDir,
                  lastzOptions="--ambiguous=iupac,100 --ydrop=3000")
    results1 = loadResults(tempOutputFile)
    logger.info("Loaded first results")
    
    for setting in settings:
        #Run the blast
        ##Record time to run
        runtime = runNaiveBlast([ seqFile1, seqFile2 ], tempOutputFile2, tempDir, 
                      lastzOptions=setting)
        
        #Now compare the results
        results2 = loadResults(tempOutputFile2)
        logger.info("Loaded second results")
        
        resultsComparator = ResultComparator(results2, results1)
        print "RESULTS", "\t".join([species1, species2, str(distance), "_".join(setting.split()), str(resultsComparator.sensitivity),
                         str(resultsComparator.specificity), str(baseRuntime), str(runtime)  ])