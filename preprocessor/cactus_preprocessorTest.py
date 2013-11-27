from cactus.preprocessor.preprocessorTest import *
from cactus.preprocessor.preprocessorTest import TestCase as PreprocessorTestCase
from cactus.shared.common import cactusRootPath
import xml.etree.ElementTree as ET

"""Runs cactus preprocessor using the lastz repeat mask script to show it working.
"""

class TestCase(PreprocessorTestCase):
    def testCactusPreprocessor(self):
        #Demo sequences
        sequenceNames = [ "%s.ENm001.fa" % species for species in 'human', "hedgehog" ]
        for sequenceName in sequenceNames:
            sequenceFile = os.path.join(self.encodePath, self.encodeRegion, sequenceName)
            #Parse sequences into dictionary
            originalSequences = getSequences(sequenceFile)
            #Make config file
            configFile = os.path.join(self.tempDir, "config.xml")
            rootElem =  ET.Element("preprocessor")
            #<preprocessor chunkSize="10000" proportionToSample="0.2" memory="littleMemory" preprocessorString="cactus_lastzRepeatMask.py --proportionSampled=PROPORTION_SAMPLED --minPeriod=1 --lastzOpts='--step=1 --ambiguous=iupac,100 --ungapped' IN_FILE OUT_FILE "/>
            preprocessor = ET.SubElement(rootElem, "preprocessor")
            preprocessor.attrib["chunkSize"] = "10000"
            preprocessor.attrib["proportionToSample"] = "0.2"
            preprocessor.attrib["string"] = "cactus_lastzRepeatMask.py --proportionSampled=PROPORTION_SAMPLED --minPeriod=1 --lastzOpts='--step=1 --ambiguous=iupac,100 --ungapped' IN_FILE OUT_FILE"
            fileHandle = open(configFile, "w")
            fileHandle.write(ET.tostring(preprocessor))
            fileHandle.close()
            #Run preprocessor
            command = "cactus_preprocessor.py %s %s %s --jobTree %s" % (self.tempDir, configFile, sequenceFile, os.path.join(self.tempDir, "jobTree"))
            system(command)
            #Load the new sequences
            processedSequences = getSequences(os.path.join(self.tempDir, sequenceName))
            #Check they are the same module masking
            self.testSequenceSetsEqualModuleSoftMasking(originalSequences, processedSequences)
            
            #Compare the proportion of bases masked by lastz with original repeat masking
            maskedBasesOriginal = getMaskedBases(originalSequences)
            maskedBasesLastzMasked = getMaskedBases(processedSequences)
            
            print " For the sequence file ", sequenceFile, \
             " the total number of sequences is ", len(originalSequences), \
             " the total number of bases ", totalBases, \
             " the number of bases originally masked was: ", len(maskedBasesOriginal),\
             " the number of bases masked after running lastz repeat masking is: ", len(maskedBasesLastzMasked), \
             " the intersection of these masked sets is: ", len(maskedBasesLastzMasked.intersection(maskedBasesOriginal)), \
             " the total number of bases that are Ns ", totalNBases, \
             " lastz was filter for max-occurrences of more than : ", maxOccurrence
        
if __name__ == '__main__':
    unittest.main()