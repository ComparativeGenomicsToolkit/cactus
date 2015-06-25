from cactus.preprocessor.preprocessorTest import *
from cactus.preprocessor.preprocessorTest import TestCase as PreprocessorTestCase

"""This test compares running the lastz repeat masking script to the underlying repeat masking of input sequences, 
comparing two settings of lastz.
"""

class TestCase(PreprocessorTestCase):
    def testLastzRepeatMask(self):
        #Demo sequences
        sequenceFiles = [ os.path.join(self.encodePath, self.encodeRegion, "%s.ENm001.fa" % species) for species in 'human', "hedgehog" ]
        #Max occurrences
        maxOccurrences = [ 1, 10 ]
        
        for sequenceFile in sequenceFiles:
            #Parse sequences into dictionary
            originalSequences = getSequences(sequenceFile)
            #Get the masked bases
            maskedBasesOriginal = getMaskedBases(originalSequences)
            #Total bases
            totalBases = sum([ len(i) for i in originalSequences.values() ])
            #Calculate number of hard masked bases
            totalNBases = len([ (header, i, base) for (header, i, base) in maskedBasesOriginal if base.upper() == "N" ])
            
            for maxOccurrence in maxOccurrences:
                #Run lastz repeat masker
                startTime = time.time()
                command = "cactus_lastzRepeatMask.py --proportionSampled=1.0 --minPeriod=%s --lastzOpts='--step=1 --ambiguous=iupac,100,100 --ydrop=3000' --fragment=%s %s %s" % \
                       (maxOccurrence, 200, sequenceFile, self.tempOutputFile)
                print "command to run", command
                popenPush(command, sequenceFile)
                print "It took %s seconds to run lastzMasking" % (time.time()-startTime)
            
                #Parse lastz masked sequences into dictionary
                lastzSequences = getSequences(self.tempOutputFile)
            
                #Check the sequences are the same modulo masking
                self.checkSequenceSetsEqualModuloSoftMasking(originalSequences, lastzSequences)
            
                #Compare the proportion of bases masked by lastz with original repeat masking
                maskedBasesOriginal = getMaskedBases(originalSequences)
                maskedBasesLastzMasked = getMaskedBases(lastzSequences)
                
                print " For the sequence file ", sequenceFile, \
                 " the total number of sequences is ", len(originalSequences), \
                 " the total number of bases ", totalBases, \
                 " the number of bases originally masked was: ", len(maskedBasesOriginal),\
                 " the number of bases masked after running lastz repeat masking is: ", len(maskedBasesLastzMasked), \
                 " the intersection of these masked sets is: ", len(maskedBasesLastzMasked.intersection(maskedBasesOriginal)), \
                 " the total number of bases that are Ns ", totalNBases, \
                 " lastz was filter for max-occurrences of more than : ", maxOccurrence
                 
                 #Run lastz repeat masker using heuristic settings for comparison with the slower settings
                command = "cactus_lastzRepeatMask.py --proportionSampled=1.0 --minPeriod=%s --lastzOpts='--step=3 --ambiguous=iupac,100,100 --ungapped --queryhsplimit=keep,nowarn:%s' --fragment=%s %s %s" % \
                       (maxOccurrence, maxOccurrence*20, 200, sequenceFile, self.tempOutputFile)
                print "command to run (fast)", command
                startTime = time.time()
                popenPush(command, sequenceFile)
                print "It took %s seconds to run lastzMasking fast" % (time.time()-startTime)
                lastzSequencesFast = getSequences(self.tempOutputFile)
                maskedBasesLastzMaskedFast = getMaskedBases(lastzSequencesFast)
                
                i = float(len(maskedBasesLastzMaskedFast.intersection(maskedBasesLastzMasked)))
                print " The number of bases masked after running faster lastz repeat masking is: ", len(maskedBasesLastzMaskedFast), \
                 " the intersection of the original and fast lastz masked sets is (bases): ", len(maskedBasesLastzMaskedFast.intersection(maskedBasesOriginal)), \
                 " the recall of the fast vs. the new is: ", i/len(maskedBasesLastzMasked), \
                 " the precision of the fast vs. the new is: ", i/len(maskedBasesLastzMaskedFast)

        
if __name__ == '__main__':
    unittest.main()