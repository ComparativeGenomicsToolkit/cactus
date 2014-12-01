#!/usr/bin/env python

"""Script to train pair-HMM of Cactus.
"""

import math
import os
import random
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import setLoggingFromOptions, logger, system, popenCatch, cigarRead, cigarWrite, nameValue, prettyXml, fastaRead
import numpy
import xml.etree.cElementTree as ET
from itertools import product

SYMBOL_NUMBER=4

class Hmm:
    def __init__(self, modelType):
        self.modelType=modelType
        self.stateNumber = { "fiveState":5, "fiveStateAsymmetric":5, "threeState":3, "threeStateAsymmetric":3}[modelType]
        self.transitions = [0.0] * self.stateNumber**2
        self.emissions = [0.0] * (SYMBOL_NUMBER**2 * self.stateNumber)
        self.likelihood = 0.0
        self.runningLikelihoods = []
        
    def _modelTypeInt(self):
        return { "fiveState":0, "fiveStateAsymmetric":1, "threeState":2, "threeStateAsymmetric":3}[self.modelType]

    def write(self, file):
        f = open(file, 'w')
        f.write(("%s " % self._modelTypeInt()) + " ".join(map(str, self.transitions)) + (" %s\n" % self.likelihood))
        f.write(" ".join(map(str, self.emissions)) + "\n")
        f.close()

    def addExpectationsFile(self, file):
        fH = open(file, 'r')
        l = map(float, fH.readline().split())
        assert len(l) == len(self.transitions)+2
        assert int(l[0]) == self._modelTypeInt()
        self.likelihood += l[-1]
        self.transitions = map(lambda x : sum(x), zip(self.transitions, l[1:-1]))
        assert len(self.transitions) == self.stateNumber**2
        l = map(float, fH.readline().split())
        assert len(l) == len(self.emissions)
        self.emissions = map(lambda x : sum(x), zip(self.emissions, l))
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2
        self.runningLikelihoods = map(float, fH.readline().split()) #This allows us to keep track of running likelihoods
        fH.close()
        return self

    @staticmethod
    def loadHmm(file):
        fH = open(file, 'r')
        l = fH.readline().split()
        assert len(l) > 0
        fH.close()
        return Hmm({ 0:"fiveState", 1:"fiveStateAsymmetric", 2:"threeState", 3:"threeStateAsymmetric"}[int(l[0])]).addExpectationsFile(file)

    def normalise(self):
        """Normalises the EM probs.
        """
        for fromState in xrange(self.stateNumber):
            i = self.stateNumber * fromState
            j = sum(self.transitions[i:i+self.stateNumber])
            for toState in xrange(self.stateNumber):
                self.transitions[i + toState] = self.transitions[i + toState] / j
        for state in xrange(self.stateNumber):
            i = state * SYMBOL_NUMBER * SYMBOL_NUMBER
            j = sum(self.emissions[i:i+SYMBOL_NUMBER * SYMBOL_NUMBER])
            for emission in xrange(SYMBOL_NUMBER * SYMBOL_NUMBER):
                self.emissions[i + emission] = self.emissions[i + emission] / j

    def randomise(self):
        """Randomise the values in the HMM to small values.
        """
        self.transitions = map(lambda x : random.random(), range(self.stateNumber*self.stateNumber))
        self.emissions = map(lambda x : random.random(), range(self.stateNumber*SYMBOL_NUMBER*SYMBOL_NUMBER))
        self.normalise()

    def equalise(self):
        """Initialise the hmm with all equal probabilities.
        """
        self.transitions = [1.0/self.stateNumber] * (self.stateNumber**2)
        self.emissions = [1.0/(SYMBOL_NUMBER*SYMBOL_NUMBER)] * (self.stateNumber * SYMBOL_NUMBER**2)
        
    def setEmissionsToJukesCantor(self, divergence):
        i = (0.25 + 0.75*math.exp(-4.0*divergence/3.0))/4.0
        j = (0.25 - 0.25*math.exp(-4.0*divergence/3.0))/4.0
        for state in xrange(self.stateNumber):
            for x in xrange(SYMBOL_NUMBER):
                for y in xrange(SYMBOL_NUMBER):
                    self.emissions[(state * SYMBOL_NUMBER**2) + x * SYMBOL_NUMBER + y] = i if x == y else j
    
    def tieEmissions(self):
        """Sets the emissions to reflect overall divergence, but not to distinguish between different base identity
        """
        for state in xrange(self.stateNumber):
            a = self.emissions[state*SYMBOL_NUMBER**2:(state+1)*SYMBOL_NUMBER**2]
            identityExpectation = sum(map(lambda i : float(a[i]) if (i % SYMBOL_NUMBER) == (i / SYMBOL_NUMBER) else 0.0, xrange(SYMBOL_NUMBER**2)))
            a = map(lambda i : identityExpectation/SYMBOL_NUMBER if (i % SYMBOL_NUMBER) == (i / SYMBOL_NUMBER) else (1.0 - identityExpectation)/(SYMBOL_NUMBER**2 - SYMBOL_NUMBER), xrange(SYMBOL_NUMBER**2))
            assert sum(a) + 0.001 > 1.0 and sum(a) - 0.001 < 1.0
            self.emissions[state*SYMBOL_NUMBER**2:(state+1)*SYMBOL_NUMBER**2] = a
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2

def expectationMaximisation(target, sequences, alignments, outputModel, options):
    #Iteratively run cactus realign to get expectations and load model.
    if options.inputModel != None: #Read in the model
        target.logToMaster("Loading the model from the input file %s" % options.inputModel)
        hmm = Hmm.loadHmm(options.inputModel)
        target.logToMaster("Loaded the model, has type %s" % hmm.modelType)
        hmm.normalise()
    else:
        target.logToMaster("Making model of type %s" % options.modelType)
        hmm = Hmm(options.modelType)
        if options.randomStart: #Make random parameters
            target.logToMaster("Using random starting parameters")
            hmm.randomise()
        else:
            hmm.equalise()
    if options.setJukesCantorStartingEmissions != None:
        hmm.setEmissionsToJukesCantor(float(options.setJukesCantorStartingEmissions))
    
    #Write out the first version of the output model
    hmm.write(outputModel)

    #Make a set of split alignment files
    alignmentsCounter = 0
    splitAlignmentFiles = []
    fH = None
    for cigar in cigarRead(alignments):
        if fH == None:
            splitAlignmentFiles.append(os.path.join(target.getGlobalTempDir(), "alignments_%s.cigar" % len(splitAlignmentFiles)))
            fH = open(splitAlignmentFiles[-1], 'w')
        alignmentsCounter += 1
        cigarWrite(fH, cigar)
        if alignmentsCounter > options.numberOfAlignmentsPerJob:
            fH.close()
            fH = None
            alignmentsCounter = 0
    if fH != None:
        fH.close()

    #Files to store expectations in
    expectationsFiles = map(lambda i : os.path.join(target.getGlobalTempDir(), "expectation_%i.txt" % i), xrange(len(splitAlignmentFiles)))
    assert len(splitAlignmentFiles) == len(expectationsFiles)

    target.setFollowOnTargetFn(expectationMaximisation2, args=(sequences, splitAlignmentFiles, outputModel, expectationsFiles, 0, [], options))

def expectationMaximisation2(target, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, runningLikelihoods, options):
    if iteration < options.iterations:
        map(lambda x : target.addChildTargetFn(calculateExpectations,
                    args=(sequences, x[0], None if (options.useDefaultModelAsStart and iteration == 0) else modelsFile, x[1], options)), 
            zip(splitAlignments, expectationsFiles))
        target.setFollowOnTargetFn(calculateMaximisation, args=(sequences, splitAlignments, modelsFile, expectationsFiles, iteration, runningLikelihoods, options))
    else:
        #Write out the likelihoods to the bottom of the file
        fH = open(modelsFile, 'a')
        fH.write("\t".join(map(str, runningLikelihoods)) + "\n")
        fH.close()

def calculateExpectations(target, sequences, alignments, modelsFile, expectationsFile, options):
    #Run cactus_realign
    system("cat %s | cactus_realign --logLevel DEBUG %s %s --outputExpectations=%s %s" % (alignments, sequences, nameValue("loadHmm", modelsFile, str), expectationsFile, options.optionsToRealign))

def calculateMaximisation(target, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, runningLikelihoods, options):
    #Load and merge the models files
    if len(expectationsFiles) > 0:
        hmm = Hmm.loadHmm(expectationsFiles[0])
        for expectationsFile in expectationsFiles[1:]:
            hmm.addExpectationsFile(expectationsFile)
        hmm.normalise()
        target.logToMaster("On %i iteration got likelihood: %s for model-type: %s, model-file %s" % (iteration, hmm.likelihood, hmm.modelType, modelsFile))
        runningLikelihoods.append(hmm.likelihood)
        target.logToMaster("On %i iteration got transitions: %s for model-type: %s, model-file %s" % (iteration, " ".join(map(str, hmm.transitions)), hmm.modelType, modelsFile))
        #If not train emissions then load up the old emissions and replace
        if options.trainEmissions:
            if options.tieEmissions:
                hmm.tieEmissions()
            target.logToMaster("On %i iteration got emissions: %s for model-type: %s, model-file %s" % (iteration, " ".join(map(str, hmm.emissions)), hmm.modelType, modelsFile))
        else:
            hmm.emissions = Hmm.loadHmm(modelsFile).emissions
            target.logToMaster("On %i using the original emissions" % iteration)

        #Write out 
        hmm.write(modelsFile)
    
    #Generate a new set of alignments, if necessary
    if options.updateTheBand:
        map(lambda alignments : target.addChildTargetFn(calculateAlignments, args=(sequences, alignments, modelsFile, options)), splitAlignments)
    
    #Start the next iteration
    target.setFollowOnTargetFn(expectationMaximisation2, args=(sequences, splitAlignments, modelsFile, expectationsFiles, iteration+1, runningLikelihoods, options))
    #Call em2

def calculateAlignments(target, sequences, alignments, modelsFile, options):
    temporaryAlignmentFile=os.path.join(target.getLocalTempDir(), "realign.cigar")
    system("cat %s | cactus_realign --logLevel DEBUG %s --loadHmm=%s %s > %s" % (alignments, sequences, modelsFile, options.optionsToRealign, temporaryAlignmentFile))
    system("mv %s %s" % (temporaryAlignmentFile, alignments))

def expectationMaximisationTrials(target, sequences, alignments, outputModel, options):
    if options.inputModel != None or not options.randomStart: #No multiple iterations
        target.setFollowOnTargetFn(expectationMaximisation, args=(sequences, alignments, outputModel, options))
    else:
        target.logToMaster("Running %s random restart trials to find best hmm" % options.trials)
        trialModels = map(lambda i : os.path.join(target.getGlobalTempDir(), "trial_%s.hmm" % i), xrange(options.trials))
        map(lambda trialModel : target.addChildTargetFn(expectationMaximisation, args=(sequences, alignments, trialModel, options)), trialModels)
        target.setFollowOnTargetFn(expectationMaximisationTrials2, args=(sequences, trialModels, outputModel, options))

def expectationMaximisationTrials2(target, sequences, trialModels, outputModel, options):
    trialHmms = map(lambda x : Hmm.loadHmm(x), trialModels)
    if options.outputTrialHmms: #Write out the different trial hmms
        for i in xrange(options.trials):
            trialHmms[i].write(outputModel + ("_%i" % i))
    #Pick the trial hmm with highest likelihood
    hmm = max(trialHmms, key=lambda x : x.likelihood)
    hmm.write(outputModel)
    if options.outputXMLModelFile != None:
        open(options.outputXMLModelFile, 'w').write(prettyXml(hmmsXML(trialHmms)))
    if options.blastScoringMatrixFile != None:
        matchProbs, gapOpen, gapExtend = makeBlastScoringMatrix(hmm, map(lambda x : x[1], reduce(lambda x, y : list(x) + list(y), map(fastaRead, sequences.split()))))
        fH = open(options.blastScoringMatrixFile, 'w')
        writeLastzScoringMatrix(fH, matchProbs, gapOpen, gapExtend)
        fH.close()
    target.logToMaster("Summary of trials:" + prettyXml(hmmsXML(trialHmms)))
    target.logToMaster("Hmm with highest likelihood: %s" % hmm.likelihood)

def hmmsXML(hmms):
    """Converts a bunch of hmms into a stupid little XML file that can be used to plot relevant statistics in a sensible way
    """
    #Get distributions on likelihood convergence and report.
    if len(hmms) == 0:
        raise RuntimeError("No hmms to summarise")
    
    #Do some checks that they are all the same time of hmm
    stateNumber = hmms[0].stateNumber
    modelType = hmms[0].modelType
    for hmm in hmms[1:]:
        if hmm.modelType != modelType:
            raise RuntimeError("Hmms not all of the same type")
        if hmm.stateNumber != stateNumber:
            raise RuntimeError("Hmms do not all have the same number of states")
    
    parent = ET.Element("hmms", { "modelType":str(modelType), "stateNumber":str(stateNumber)})
    
    #For each hmm
    for hmm in hmms:
        child = ET.SubElement(parent, "hmm")
        child.attrib["likelihood"] = str(hmm.likelihood)
        child.attrib["runningLikelihoods"] = "\t".join(map(str, hmm.runningLikelihoods))
        child.attrib["transitions"] = "\t".join(map(str, hmm.transitions))
        child.attrib["emissions"] = "\t".join(map(str, hmm.emissions))
    
    #Plot aggregate distributions

    ##Get the distribution on likelihoods.
    likelihoods = map(lambda x : x.likelihood, hmms) 
    parent.attrib["maxLikelihood"] = str(max(likelihoods))
    parent.attrib["likelihoods"] = "\t".join(map(str, likelihoods))
    parent.attrib["likelihoodAvg"] = str(numpy.average(likelihoods))
    parent.attrib["likelihoodStdDev"] = str(numpy.std(likelihoods))

    def statFn(values, node):
        node.attrib["max"] = str(max(values))
        node.attrib["avg"] = str(numpy.average(values))
        node.attrib["std"] = str(numpy.std(values))
        node.attrib["min"] = str(min(values))
        node.attrib["distribution"] = "\t".join(map(str, values))

    #For each transitions report ML estimate, and distribution of parameters + variance.
    for fromState in range(stateNumber):
        for toState in range(stateNumber):
            statFn(map(lambda x : x.transitions[fromState*stateNumber + toState], hmms), 
                   ET.SubElement(parent, "transition", {"from":str(fromState), "to":str(toState)}))

    #For each emission report ML estimate, and distribution of parameters + variance.
    for state in range(stateNumber):
        for x in range(4):
            for y in range(4):
                statFn(map(lambda z : z.emissions[state * 16 + x * 4 + y], hmms), 
                       ET.SubElement(parent, "emission", {"state":str(state), "x":"ACGT"[x], "y":"ACGT"[y]}))

    return parent

def makeBlastScoringMatrix(hmm, sequences):
    """Converts an hmm into a lastz style scoring matrix
    """
    #convert to a three state hmm
    hmm2 = Hmm("threeState")
    hmm2.transitions = hmm.transitions[:3] + hmm.transitions[hmm.stateNumber*1:hmm.stateNumber*1+3] + hmm.transitions[hmm.stateNumber*2:hmm.stateNumber*2+3]
    hmm2.emissions = hmm.emissions[:3 * SYMBOL_NUMBER**2]
    hmm2.normalise()
    hmm = hmm2
    
    #Get gap distribution, assuming we include reverse complement sequences then it's fraction of GCs
    gcFraction = sum(map(lambda x : sum(map(lambda y : 1.0 if y in 'GC' else 0.0, x)), sequences)) / sum(map(len, sequences))
    logger.debug("Got the GC fraction in the sequences for making the scoring matrix: %s" % gcFraction)
    baseProb = lambda x : gcFraction/2.0 if x in (1,2) else (1.0 - gcFraction)/2.0
  
    #Calculate match matrix
    logger.debug("Original match probs: %s" % " ".join(map(str, hmm.emissions[:SYMBOL_NUMBER**2])))
    matchProbs = [ hmm.emissions[x * SYMBOL_NUMBER + y] / (baseProb(x) * baseProb(y)) for x, y in product(range(SYMBOL_NUMBER), range(SYMBOL_NUMBER)) ]
    logger.debug("Blast emission match probs: %s" % " ".join(map(str, matchProbs)))
    matchContinue = hmm.transitions[0]
    #The 6.94 is the 1/100th the sum of the lastz scoring matrix
    nProb = math.sqrt(math.exp((6.94+sum(map(lambda x : math.log(x * matchContinue), matchProbs)))/len(matchProbs)))
    logger.debug("N prob is: %s" % nProb) #Note it may go above 1!
    weight=100
    matchProbs = map(lambda x : weight*math.log((x * matchContinue) / nProb**2), matchProbs)
    logger.debug("Blast match probs, %s: %s" % (sum(matchProbs)/4.0, " ".join(map(str, matchProbs))))
    
    #Calculate gap open
    gapOpen = weight*math.log((0.5 * (hmm.transitions[1]/nProb + hmm.transitions[2]/nProb)) * \
    ((hmm.transitions[hmm.stateNumber*1 + 0] + hmm.transitions[hmm.stateNumber*2 + 0])/(2*nProb**2)) * \
    ((nProb**2)/matchContinue))
    logger.debug("Gap open: %s" % gapOpen)
    
    #Calculate gap extend
    gapContinue = weight*math.log(0.5 * (hmm.transitions[hmm.stateNumber*1 + 1]/nProb + hmm.transitions[hmm.stateNumber*2 + 2]/nProb))
    logger.debug("Gap continue: %s" % gapContinue)

    return matchProbs, gapOpen, gapContinue

def writeLastzScoringMatrix(fileHandle, matchProbs, gapOpen, gapExtend):
    """# This matches the default scoring set for BLASTZ
    
    bad_score          = X:-1000  # used for sub['X'][*] and sub[*]['X']
    fill_score         = -100     # used when sub[*][*] is not defined
    gap_open_penalty   =  400
    gap_extend_penalty =   30

         A     C     G     T
    A   91  -114   -31  -123
    C -114   100  -125   -31
    G  -31  -125   100  -114
    T -123   -31  -114    91
    """
    fileHandle.write("gap_open_penalty = %s\n" % int(round(-gapOpen)))
    fileHandle.write("gap_extend_penalty = %s\n" % int(round(-gapExtend)))
    bases = "ACGT"
    fileHandle.write("\t\t" + "\t".join(bases) + "\n")
    for x in range(4):
        fileHandle.write("\t%s\t%s\n" % (bases[x], "\t".join(map(lambda x : str(int(round(x))), matchProbs[x*SYMBOL_NUMBER:((x+1)*SYMBOL_NUMBER)]))))
    
class Options:
    """Dictionary representing options, can be used for running pipeline from within another jobTree.
    """
    def __init__(self):
        self.modelType="fiveState"
        self.inputModel=None
        self.iterations=10
        self.trials=3
        self.outputTrialHmms = False
        self.randomStart=False
        self.optionsToRealign="--diagonalExpansion=10 --splitMatrixBiggerThanThis=3000"
        self.updateTheBand=False
        self.numberOfAlignmentsPerJob=100
        self.useDefaultModelAsStart = False
        self.setJukesCantorStartingEmissions=None
        self.tieEmissions = False
        self.trainEmissions=False
        self.outputXMLModelFile = None
        self.blastScoringMatrixFile = None

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: workingDir [options]", version="%prog 0.1")
    options = Options()
    parser.add_option("--sequences", dest="sequences", help="Quoted list of fasta files containing sequences")
    parser.add_option("--alignments", dest="alignments", help="Cigar file ")
    parser.add_option("--inputModel", default=options.inputModel, help="Input model")
    parser.add_option("--outputModel", default="hmm.txt", help="File to write the model in")
    parser.add_option("--outputXMLModelFile", default=options.outputXMLModelFile, help="File to write XML representation of model in - useful for stats")
    parser.add_option("--modelType", default=options.modelType, help="Specify the model type, currently either fiveState, threeState, threeStateAsymmetric")
    parser.add_option("--iterations", default=options.iterations, help="Number of iterations of EM", type=int)
    parser.add_option("--trials", default=options.trials, help="Number of independent EM trials. The model with the highest likelihood will be reported. Will only work if randomStart=True", type=int)
    parser.add_option("--outputTrialHmms", default=options.outputTrialHmms, help="Writes out the final trained hmm for each trial, as outputModel + _i", action="store_true")
    parser.add_option("--randomStart", default=options.randomStart, help="Iterate start model with small random values, else all values are equal", action="store_true")
    parser.add_option("--optionsToRealign", default=options.optionsToRealign, help="Further options to pass to cactus_realign when computing expectation values, should be passed as a quoted string")
    parser.add_option("--updateTheBand", default=options.updateTheBand, help="After each iteration of EM update the set of alignments by realigning them, so allowing stochastic updating of the constraints. This does not alter the input alignments file", action="store_true")
    parser.add_option("--numberOfAlignmentsPerJob", default=options.numberOfAlignmentsPerJob, help="Number of alignments to compute in parallel during expectation/updating the band steps", type=int)
    parser.add_option("--useDefaultModelAsStart", default=options.useDefaultModelAsStart, help="Use the default BAR hmm model as the starting point", action="store_true")
    parser.add_option("--setJukesCantorStartingEmissions", default=options.setJukesCantorStartingEmissions, help="[double] Set the starting hmm emissions by jukes cantor expectation, using given subs/site estimate", type=float)
    parser.add_option("--trainEmissions", default=options.trainEmissions, help="Train the emissions as well as the transitions.", action="store_true")
    parser.add_option("--tieEmissions", default=options.tieEmissions, help="Normalise all emissions to reflect overall level of diversity, but be tied to not reflect differences between different bases, other than identity/difference.", action="store_true")
    parser.add_option("--blastScoringMatrixFile", default=options.blastScoringMatrixFile, help="Calculate a BLAST scoring matrix from the HMM, output in Lastz/Blastz format")
    
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Expected no arguments, got %s arguments: %s" % (len(args), " ".join(args)))
    
    #Log the inputs
    logger.info("Got '%s' sequences, '%s' alignments file, '%s' output model and '%s' iterations of training" % (options.sequences, options.alignments, options.outputModel, options.iterations))

    #This line invokes jobTree  
    i = Stack(Target.makeTargetFn(expectationMaximisationTrials, args=(options.sequences, options.alignments, options.outputModel, options))).startJobTree(options) 
    
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from cactus.bar.cactus_expectationMaximisation import *
    main()
