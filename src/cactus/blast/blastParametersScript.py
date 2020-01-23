import os

from sonLib.bioio import system
from sonLib.bioio import logger
from cactus.blast.blastTest import *

"""Tests if changing parameters of lastz creates results similar to the desired default, checking at different evolutionary distances. Outputs comma seperated values.
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
                ('human', 'tenrec', 0.49516499999999997),
                ('rat', 'mouse', 0.15406199999999998), ('chimp', 'mouse', 0.47689000000000004),  
                ('baboon', 'mouse', 0.48644000000000004), ('marmoset', 'mouse', 0.48823300000000003), ('macaque', 'mouse', 0.492433), 
                ('galago', 'mouse', 0.506291), ('mouse', 'rabbit', 0.5155689999999999), ('mouse', 'elephant', 0.53935), 
                ('mouse', 'rfbat', 0.543197), ('mouse', 'armadillo', 0.546562), ('mouse', 'dog', 0.562574), 
                ('mouse', 'cow', 0.569921), ('mouse', 'hedgehog', 0.676318), ('mouse', 'tenrec', 0.6873229999999999), 
                ('dog', 'rfbat', 0.30471099999999995), ('cow', 'dog', 0.322109), 
                ('chimp', 'dog', 0.373174), ('baboon', 'dog', 0.382724), ('marmoset', 'dog', 0.384517), ('macaque', 'dog', 0.388717), 
                ('dog', 'elephant', 0.389114), ('dog', 'armadillo', 0.396326), ('galago', 'dog', 0.402575), 
                ('rabbit', 'dog', 0.45563699999999996), ('dog', 'hedgehog', 0.459598), ('dog', 'tenrec', 0.537087), 
                ('mouse', 'dog', 0.562574), ('rat', 'dog', 0.571)
                 ]
settings = ["--step=2 --ambiguous=iupac,100 --ydrop=3000",
            "--step=2 --ambiguous=iupac,100 --ydrop=3000 --notransition",
            "--step=3 --ambiguous=iupac,100 --ydrop=3000",
            "--step=3 --ambiguous=iupac,100 --ydrop=3000 --notransition",
            "--step=4 --ambiguous=iupac,100 --ydrop=3000",
            "--step=4 --ambiguous=iupac,100 --ydrop=3000 --notransition",
            "--step=5 --ambiguous=iupac,100 --ydrop=3000",
            "--step=5 --ambiguous=iupac,100 --ydrop=3000 --notransition",
            "--step=8 --ambiguous=iupac,100 --ydrop=3000",
            "--step=8 --ambiguous=iupac,100 --notransition --ydrop=3000",
            "--step=8 --ambiguous=iupac,100 --nogapped",
            "--step=8 --ambiguous=iupac,100 --nogapped --notransition",
            "--step=16 --ambiguous=iupac,100 --ydrop=3000",
            "--step=16 --ambiguous=iupac,100 --notransition --ydrop=3000",
            "--step=16 --ambiguous=iupac,100 --nogapped",
            "--step=16 --ambiguous=iupac,100 --nogapped --notransition" ]

#Other species to try "rat", "monodelphis", "macaque", "chimp"
encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
regionPath = os.path.join(encodePath, encodeRegion)
tempDir = getTempDirectory(os.getcwd())
tempOutputFile = os.path.join(tempDir, "results1.txt")
tempOutputFile2 = os.path.join(tempDir, "results2.txt")

print((",".join(["species1", "species2", "setting", "distance", "sensitivity", "specificity", "intersectionSize", "unionSize", "totalPairsMissed", "totalPairsGained", "totalTrueHits", "totalPredictedHits", "hitsDifference", "base-runtime", "runtime" ])))

for species1, species2, distance in speciesPairs:
    seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
    seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
    
    #Run the random
    ##Record time to run
    baseRuntime = runNaiveBlast(seqFile1, seqFile2, tempOutputFile, 
                  lastzOptions="--ambiguous=iupac,100 --ydrop=3000")
    results1 = loadResults(tempOutputFile)
    logger.info("Loaded first results")
    
    for setting in settings:
        #Run the blast
        ##Record time to run
        runtime = runNaiveBlast(seqFile1, seqFile2, tempOutputFile2,
                      lastzOptions=setting)
        
        #Now compare the results
        results2 = loadResults(tempOutputFile2)
        logger.info("Loaded second results")
        
        def fm(f):
            return "%.5f" % float(f)
        
        def fm2(f):
            return str(int(f))
        
        resultsComparator = ResultComparator(results1, results2)
        print((",".join([ species1, species2, "_".join(("_".join(setting.split())).split(",")), fm(distance), fm(resultsComparator.sensitivity),
                         fm(resultsComparator.specificity),
                         fm2(resultsComparator.intersectionSize), fm2(resultsComparator.unionSize),
                         fm2(resultsComparator.trueDifference), fm2(resultsComparator.predictedDifference),
                         fm2(resultsComparator.trueHits), fm2(resultsComparator.predictedHits), fm2(resultsComparator.trueHits -resultsComparator.predictedHits), fm(baseRuntime), fm(runtime) ])))
        
system("rm -rf %s" % tempDir)
