import time

from cactus.preprocessor.preprocessorTest import *
from cactus.preprocessor.preprocessorTest import TestCase as PreprocessorTestCase
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import LastzRepeatMaskJob
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMask import RepeatMaskOptions

from toil.common import Toil
from toil.job import Job

from cactus.shared.common import makeURL

"""This test compares running the lastz repeat masking script to the underlying repeat masking of input sequences,
comparing two settings of lastz.
"""

@pytest.mark.blast
class TestCase(PreprocessorTestCase):
    @TestStatus.mediumLength
    def testLastzRepeatMask(self):
        #Demo sequences
        sequenceFiles = [ os.path.join(self.encodePath, self.encodeRegion, "%s.ENm001.fa" % species) for species in ('human', "hedgehog") ]
        #Max occurrences of a repeat within the sequence
        maxOccurrence = 1

        for sequenceFile in sequenceFiles:
            #Parse sequences into dictionary
            originalSequences = getSequences(sequenceFile)
            #Get the masked bases
            maskedBasesOriginal = getMaskedBases(originalSequences)
            #Total bases
            totalBases = sum([ len(i) for i in list(originalSequences.values()) ])
            #Calculate number of hard masked bases
            totalNBases = len([ (header, i, base) for (header, i, base) in maskedBasesOriginal if base.upper() == "N" ])

            #Run lastz repeat masker
            startTime = time.time()
            with Toil(self.toilOptions) as toil:
                sequenceID = toil.importFile(makeURL(sequenceFile))
                repeatMaskOptions = RepeatMaskOptions(proportionSampled=1.0,
                                         minPeriod=maxOccurrence,
                                         lastzOpts="--step=1 --ambiguous=iupac,100,100 --ydrop=3000",
                                         fragment=200)

                outputID = toil.start(LastzRepeatMaskJob(repeatMaskOptions=repeatMaskOptions, queryID=sequenceID, targetIDs=[sequenceID]))
                toil.exportFile(outputID, makeURL(self.tempOutputFile))
            print(("It took %s seconds to run lastzMasking" % (time.time()-startTime)))

            #Parse lastz masked sequences into dictionary
            lastzSequences = getSequences(self.tempOutputFile)

            #Check the sequences are the same modulo masking
            self.checkSequenceSetsEqualModuloSoftMasking(originalSequences, lastzSequences)

            #Compare the proportion of bases masked by lastz with original repeat masking
            maskedBasesOriginal = getMaskedBases(originalSequences)
            maskedBasesLastzMasked = getMaskedBases(lastzSequences)
            print((" For the sequence file ", sequenceFile, \
             " the total number of sequences is ", len(originalSequences), \
             " the total number of bases ", totalBases, \
             " the number of bases originally masked was: ", len(maskedBasesOriginal),\
             " the number of bases masked after running lastz repeat masking is: ", len(maskedBasesLastzMasked), \
             " the intersection of these masked sets is: ", len(maskedBasesLastzMasked.intersection(maskedBasesOriginal)), \
             " the total number of bases that are Ns ", totalNBases, \
             " lastz was filter for max-occurrences of more than : ", maxOccurrence))
            #self.assertGreater(len(maskedBasesLastzMasked), len(maskedBasesOriginal))

            #Run lastz repeat masker using heuristic settings for comparison with the slower settings
            startTime = time.time()
            with Toil(self.toilOptions) as toil:
                sequenceID = toil.importFile(makeURL(sequenceFile))
                repeatMaskOptions = RepeatMaskOptions(proportionSampled=1.0,
                                                    minPeriod=maxOccurrence,
                                                    lastzOpts="--step=3 --ambiguous=iupac,100,100 --ungapped --queryhsplimit=keep,nowarn:%i" % (int(maxOccurrence)*20),
                                                    fragment=200)
                outputID = toil.start(LastzRepeatMaskJob(repeatMaskOptions=repeatMaskOptions, queryID=sequenceID, targetIDs=[sequenceID]))
                toil.exportFile(outputID, makeURL(self.tempOutputFile))
            print(("It took %s seconds to run lastzMasking fast" % (time.time()-startTime)))
            lastzSequencesFast = getSequences(self.tempOutputFile)
            maskedBasesLastzMaskedFast = getMaskedBases(lastzSequencesFast)

            self.assertGreater(len(maskedBasesLastzMaskedFast), len(maskedBasesOriginal))
            i = float(len(maskedBasesLastzMaskedFast.intersection(maskedBasesLastzMasked)))
            precision = i/len(maskedBasesLastzMasked)
            recall = i/len(maskedBasesLastzMaskedFast)
            self.assertGreater(precision, 0.93)
            self.assertGreater(recall, 0.93)


if __name__ == '__main__':
    unittest.main()
