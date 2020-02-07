import os
import pytest
import unittest

from sonLib.bioio import TestStatus
from cactus.preprocessor.preprocessorTest import getSequences, getMaskedBases
from cactus.preprocessor.preprocessorTest import TestCase as PreprocessorTestCase
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
import xml.etree.ElementTree as ET
from cactus.preprocessor.cactus_preprocessor import runCactusPreprocessor

"""Runs cactus preprocessor using the lastz repeat mask script to show it working.
"""

@pytest.mark.blast
class TestCase(PreprocessorTestCase):
    @TestStatus.mediumLength
    def testCactusPreprocessor(self):
        #Demo sequences
        sequenceNames = [ "%s.ENm001.fa" % species for species in ['human', 'hedgehog'] ]
        sequenceFiles = [ os.path.join(self.encodePath, self.encodeRegion, sequenceName) for sequenceName in sequenceNames ]
        #Make config file
        configFile = os.path.join(self.tempDir, "config.xml")
        rootElem =  ET.Element("preprocessor")
        #<preprocessor chunkSize="10000" proportionToSample="0.2" memory="littleMemory" preprocessorString="cactus_lastzRepeatMask.py --proportionSampled=PROPORTION_SAMPLED --minPeriod=1 --lastzOpts='--step=1 --ambiguous=iupac,100 --ungapped' IN_FILE OUT_FILE "/>
        preprocessor = ET.SubElement(rootElem, "preprocessor")
        preprocessor.attrib["chunkSize"] = "100000"
        preprocessor.attrib["proportionToSample"] = "0.2"
        preprocessor.attrib["preprocessJob"] = "lastzRepeatMask"
        preprocessor.attrib["minPeriod"] = "1"
        preprocessor.attrib["lastzOpts"] = "--step=1 --ambiguous=iupac,100 --ungapped"
        preprocessor.attrib["fragment"] = "200"
        fileHandle = open(configFile, "w")
        fileHandle.write(ET.tostring(rootElem, encoding='unicode'))
        fileHandle.close()
        #Run preprocessor
        tmpToil = os.path.join(self.tempDir, "toil")
        runCactusPreprocessor(outputSequenceDir=self.tempDir, configFile=configFile, inputSequences=sequenceFiles, toilDir=tmpToil)

        for sequenceFile, processedSequenceFile in zip(sequenceFiles, CactusPreprocessor.getOutputSequenceFiles(sequenceFiles, self.tempDir)):
            #Parse sequences into dictionary
            originalSequences = getSequences(sequenceFile)
            #Load the new sequences
            processedSequences = getSequences(processedSequenceFile)

            #Check they are the same module masking
            self.checkSequenceSetsEqualModuloSoftMasking(originalSequences, processedSequences)

            #Compare the proportion of bases masked by lastz with original repeat masking
            maskedBasesOriginal = getMaskedBases(originalSequences)
            maskedBasesLastzMasked = getMaskedBases(processedSequences)
            #Total bases
            totalBases = sum([ len(i) for i in list(originalSequences.values()) ])
            #Calculate number of hard masked bases
            totalNBases = len([ (header, i, base) for (header, i, base) in maskedBasesOriginal if base.upper() == "N" ])

            print((" For the sequence file ", sequenceFile, \
             " the total number of sequences is ", len(originalSequences), \
             " the total number of bases ", totalBases, \
             " the number of bases originally masked was: ", len(maskedBasesOriginal),\
             " the number of bases masked after running lastz repeat masking is: ", len(maskedBasesLastzMasked), \
             " the intersection of these masked sets is: ", len(maskedBasesLastzMasked.intersection(maskedBasesOriginal)), \
             " the total number of bases that are Ns ", totalNBases))
            self.assertGreater(maskedBasesLastzMasked, maskedBasesOriginal)


if __name__ == '__main__':
    if "SON_TRACE_DATASETS" in os.environ:
        unittest.main()
