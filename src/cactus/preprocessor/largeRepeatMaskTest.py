import os
import argparse
from toil.job import Job
from toil.common import Toil

from cactus.shared.common import makeURL
from cactus.preprocessor.cactus_preprocessor import PreprocessorOptions, PreprocessSequence


def preprocessSequences(job, seqIDs, prepOptions):
    for seqID in seqIDs:
        job.addChild(PreprocessSequence(prepOptions=prepOptions, inSequenceID=seqID, chunksToCompute=[4,6,20,23,25]))

def main():
    parser = argparse.ArgumentParser()
    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()
    mouse = "http://hgwdev.cse.ucsc.edu/~adderan/HumanMouseRatDog/mm10.fa"
    human = "http://hgwdev.cse.ucsc.edu/~adderan/HumanMouseRatDog/hg19.fa"
    
    prepOptions = PreprocessorOptions(chunkSize=3000000, preprocessJob="lastzRepeatMask",
                                      memory=2000000000,
                                      cpu=1,
                                      check=1,
                                      proportionToSample=0.2,
                                      unmask=0,
                                      lastzOptions="--step=3 --ambiguous=iupac,100,100 --ungapped --queryhsplimit=keep,nowarn:1500",
                                      minPeriod=50,
                                      checkAssemblyHub=False)
    with Toil(options) as toil:
        humanID = toil.importFile(makeURL(human))
        #mouseID = toil.importFile(makeURL(mouse))
        toil.start(Job.wrapJobFn(preprocessSequences, seqIDs=[humanID], prepOptions=prepOptions))


if __name__ == "__main__":
    main()
