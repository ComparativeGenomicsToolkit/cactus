import unittest
import os, sys
from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import logger
from sonLib.bioio import getTempFile, getTempDirectory, system
from cactus.bar.cactus_realignTest import seqFilePairGenerator
from cactus.shared.common import runCactusExpectationMaximisation
from cactus.bar.cactus_expectationMaximisation import Hmm, writeLastzScoringMatrix, makeBlastScoringMatrix
from sonLib.bioio import getLogLevelString
import xml.etree.cElementTree as ET

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testCactusExpectationMaximisation(self):
        """Runs cactus expectation-maximisation. 
        """
        trial = 0
        for modelType in ("fiveState", "fiveStateAsymmetric", "threeState", "threeStateAsymmetric"):
            for seqFile1, seqFile2 in seqFilePairGenerator():
                tempDir = getTempDirectory(rootDir=os.getcwd())
                jobTreeDir = os.path.join(tempDir, "jobTree")
                alignmentsFile = os.path.join(tempDir, "alignments.cigars")
                computeAlignments(seqFile1, seqFile2, alignmentsFile)
                logger.info("Computed alignments for seqs %s and %s" % (seqFile1, seqFile2))
                outputModelFile = os.path.join(tempDir, "outputModel.txt")
                #First run the script to generate a model and do one iteration of EM to 
                #get the likelihood to compare with the final likelihood
                runCactusExpectationMaximisation(sequenceFiles=[ seqFile1, seqFile2 ], 
                                                 alignmentsFile=alignmentsFile, outputModelFile=outputModelFile, 
                                                 modelType=modelType,
                                                 jobTreeDir=jobTreeDir,
                                                 iterations=1, trials=1, randomStart=False, logLevel=getLogLevelString(),
                                                 setJukesCantorStartingEmissions=0.2,
                                                 #useDefaultModelAsStart=,
                                                 trainEmissions=True,
                                                 optionsToRealign="--diagonalExpansion=6 --splitMatrixBiggerThanThis=100")
                hmm = Hmm.loadHmm(outputModelFile)
                system("rm -rf %s" % jobTreeDir) #Cleanup the old jobTree
                logger.info("For trial %s the likelihood after 1 iteration of EM is %s" % (trial, hmm.likelihood))
                iterations = 5
                runCactusExpectationMaximisation(sequenceFiles=[ seqFile1, seqFile2 ], 
                                                 alignmentsFile=alignmentsFile, outputModelFile=outputModelFile, jobTreeDir=jobTreeDir,
                                                 optionsToRealign="--diagonalExpansion=6 --splitMatrixBiggerThanThis=100",
                                                 iterations=iterations, inputModelFile=outputModelFile, logLevel=getLogLevelString(),
                                                 numberOfAlignmentsPerJob=20) #, updateTheBand=True)
                hmm2 = Hmm.loadHmm(outputModelFile)
                logger.info("For trial %s the likelihood after a further %s iterations of EM is %s" % (trial, iterations, hmm2.likelihood))
                self.assertTrue(hmm.likelihood < hmm2.likelihood)
                hmm2.normalise()
                logger.info("Final transitions: %s" % " ".join(map(str, hmm2.transitions)))
                logger.info("Final emissions: %s" % " ".join(map(str, hmm2.emissions)))
                system("rm -rf %s" % tempDir)
                trial += 1
    
    def testCactusExpectationMaximisationMultipleTrials(self):
        """Runs cactus expectation-maximisation with multiple different trials.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            tempDir = getTempDirectory(rootDir=os.getcwd())
            jobTreeDir = os.path.join(tempDir, "jobTree")
            alignmentsFile = os.path.join(tempDir, "alignments.cigars")
            computeAlignments(seqFile1, seqFile2, alignmentsFile)
            logger.info("Computed alignments for seqs %s and %s" % (seqFile1, seqFile2))
            outputModelFile = os.path.join(tempDir, "outputModel.txt")
            outputModelXMLFile = os.path.join(tempDir, "outputModel.xml")
            outputBlastFile = os.path.join(tempDir, "outputBlast.txt")
            #First run the script to generate a model and do one iteration of EM to 
            #get the likelihood to compare with the final likelihood
            trials=3
            runCactusExpectationMaximisation(sequenceFiles=[ seqFile1, seqFile2 ], 
                                             alignmentsFile=alignmentsFile, outputModelFile=outputModelFile, 
                                             jobTreeDir=jobTreeDir,
                                             trials=trials,
                                             outputTrialHmms=True,
                                             iterations=5, randomStart=True, logLevel=getLogLevelString(),
                                             optionsToRealign="--diagonalExpansion=6 --splitMatrixBiggerThanThis=100",
                                             outputXMLModelFile=outputModelXMLFile,
                                             blastScoringMatrixFile=outputBlastFile)
            trialHmms = [ Hmm.loadHmm(outputModelFile + ("_%i" % i)) for i in xrange(trials) ]
            hmm = Hmm.loadHmm(outputModelFile)
            node = ET.parse(outputModelXMLFile).getroot()
            logger.info("After multiple trials and iterations of EM the best likelihood found was %s, the likelihoods of the variants were: %s" % 
                        (hmm.likelihood, " ".join(map(lambda x : str(x.likelihood), trialHmms))))
            
            matchProbs, gapOpen, gapExtend = makeBlastScoringMatrix(hmm, ("ACTG",))
            logger.info("Gap open: %s, Gap extend: %s, Match probs %s" % (gapOpen, gapExtend, " ".join(map(str, matchProbs))))
            
            self.assertTrue(float(node.attrib["maxLikelihood"]) == hmm.likelihood)
            
            #Now use the blast file to compute a new matrix
            computeAlignments(seqFile1, seqFile2, alignmentsFile, lastzArguments=("--ambiguous=iupac --scores=%s" % outputBlastFile))
            
            system("rm -rf %s" % tempDir)
    
    def testHMMToBlast(self):
        hmmFile = getTempFile()
        fH = open(hmmFile, 'w')
        ##This is an HMM trained from some nanopore data.
        fH.write("1 0.769849545837 0.192124461785 0.0373796764444 0.000454412327473 0.000191903606847 0.615582885081 0.384417114919 0.0 0.0 0.0 0.360492263331 0.0 0.639507736669 0.0 0.0 0.000837020116237 0.0 0.0 0.999162979884 0.0 0.00263503613846 0.0 0.0 0.0 0.997364963862 -83964693614.2\n")
        fH.write("0.124467347093 0.0510185372341 0.055395149667 0.0165710761929 0.0424721479638 0.107314387884 0.0326918092026 0.0169710192012 0.0512318749092 0.0282356959149 0.112202089573 0.0215204542575 0.0384239342042 0.0479915657848 0.046228721193 0.207264189725 0.0660896084237 0.0593325150728 0.0603492612177 0.0582600197325 0.0584825157865 0.0522874453523 0.0584394811677 0.0575754515175 0.0527867639602 0.0495513728754 0.0503223877237 0.0558605997644 0.0825830058041 0.076482194432 0.0786344471539 0.0829629300156 0.061340281862 0.0603951822769 0.0697104777685 0.0624365718074 0.0515829309891 0.044694507015 0.0563446095066 0.0500495229974 0.0535285782675 0.0498949344494 0.0552268591741 0.0513874548672 0.0834191759666 0.0807522229048 0.088970651067 0.0802660390805 0.0585941736742 0.0666707645145 0.0671013181763 0.0540169346484 0.0546857600483 0.0614356162075 0.0637011684438 0.0493325269862 0.0497656880321 0.0565185406621 0.0574882564972 0.0453586109249 0.0756527827468 0.0866857198132 0.0833713253471 0.0696208132772 0.0489551609873 0.0485086714065 0.055299084522 0.0480626877811 0.0360480113306 0.0352929139202 0.0404027968092 0.0362560298668 0.0504819795621 0.0505241023095 0.0581098587574 0.0526378897109 0.107006368897 0.106127475919 0.11956635493 0.10672061329")
        fH.close()
        hmm = Hmm.loadHmm(hmmFile)
        matchProbs, gapOpen, gapExtend = makeBlastScoringMatrix(hmm, ("ACTG",))
        writeLastzScoringMatrix(sys.stdout, matchProbs, gapOpen, gapExtend)
        logger.info("Gap open: %s, Gap extend: %s, Match probs %s" % (gapOpen, gapExtend, " ".join(map(str, matchProbs))))

def computeAlignments(seqFile1, seqFile2, alignmentsFile, lastzArguments="--ambiguous=iupac"):  
    system("cactus_lastz --format=cigar %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace] > %s" % (lastzArguments, seqFile1, seqFile2, alignmentsFile))

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
