from cactus.preprocessor.preprocessorTest import *
from cactus.preprocessor.preprocessorTest import TestCase as PreprocessorTestCase
from cactus.shared.common import cactusRootPath
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
import xml.etree.ElementTree as ET
from sonLib.bioio import nameValue, logger

"""Runs cactus preprocessor using the lastz repeat mask script to show it working.
"""

class TestCase(PreprocessorTestCase):
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
        preprocessor.attrib["preprocessorString"] = "cactus_lastzRepeatMask.py --proportionSampled=PROPORTION_SAMPLED --fragment=200 --minPeriod=1 --lastzOpts='--step=1 --ambiguous=iupac,100 --ungapped' IN_FILE OUT_FILE"
        fileHandle = open(configFile, "w")
        fileHandle.write(ET.tostring(rootElem))
        fileHandle.close()
        #Run preprocessor
        inputSequencesStr = nameValue("inputSequences", " ".join(sequenceFiles), quotes=True)
        outputSequenceDirStr = nameValue("outputSequenceDir", self.tempDir)
        configFileStr = nameValue("configFile", configFile)
        tmpToil = "./toil" #os.path.join(self.tempDir, "toil")
        command = "cactus_preprocessor.py %s --logLevel=DEBUG %s %s %s" % (tmpToil, outputSequenceDirStr, configFileStr, inputSequencesStr)
        system(command)
        
        for sequenceFile, processedSequenceFile in zip(sequenceFiles, CactusPreprocessor.getOutputSequenceFiles(sequenceFiles, self.tempDir)):
            print "sequenceFile: %s" % sequenceFile
            print "output sequence file: %s" % processedSequenceFile
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
            totalBases = sum([ len(i) for i in originalSequences.values() ])
            #Calculate number of hard masked bases
            totalNBases = len([ (header, i, base) for (header, i, base) in maskedBasesOriginal if base.upper() == "N" ])
            
            print " For the sequence file ", sequenceFile, \
             " the total number of sequences is ", len(originalSequences), \
             " the total number of bases ", totalBases, \
             " the number of bases originally masked was: ", len(maskedBasesOriginal),\
             " the number of bases masked after running lastz repeat masking is: ", len(maskedBasesLastzMasked), \
             " the intersection of these masked sets is: ", len(maskedBasesLastzMasked.intersection(maskedBasesOriginal)), \
             " the total number of bases that are Ns ", totalNBases
             
            #Now compare to running lastz on its own
            command = "cactus_lastzRepeatMask.py --proportionSampled=0.2 --minPeriod=1 --lastzOpts='--step=1 --ambiguous=iupac,100 --ungapped --queryhsplimit=keep,nowarn:30' --fragment=200 %s %s" % \
                       (sequenceFile, self.tempOutputFile) 
            popenPush(command, sequenceFile)
            lastzSequencesFast = getSequences(self.tempOutputFile)
            maskedBasesLastzMaskedFast = getMaskedBases(lastzSequencesFast)
            
            i = float(len(maskedBasesLastzMaskedFast.intersection(maskedBasesLastzMasked)))
            print " The number of bases masked after running lastz repeat masking without the preprocessor is: ", len(maskedBasesLastzMaskedFast), \
             " the recall of the fast vs. the new is: ", i/len(maskedBasesLastzMasked), \
             " the precision of the fast vs. the new is: ", i/len(maskedBasesLastzMaskedFast)
        
if __name__ == '__main__':
    unittest.main()
