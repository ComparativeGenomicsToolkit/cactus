from cactus.repeats.cactus_repeats import *
from cactus.preprocessor.preprocessorTest import TestCase as PreprocessorTestCase

from toil.common import Toil
from toil.job import Job

from cactus.shared.common import makeURL

"""This test compares running the lastz repeat masking script to the underlying repeat masking of input sequences, 
comparing two settings of lastz.
"""


class TestCase(PreprocessorTestCase):
    def testLastzRepeatMask(self):
        #Demo sequences
        sequenceFiles = [ os.path.join(self.encodePath, self.encodeRegion, "%s.ENm001.fa" % species) for species in 'human', "hedgehog" ]
        #Max occurrences
        
        for sequenceFile in sequenceFiles:
            pass
        
if __name__ == '__main__':
    unittest.main()
