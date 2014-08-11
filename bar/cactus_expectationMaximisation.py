#!/usr/bin/env python

"""Script to train pair-HMM of Cactus.
"""

import os
import random
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import setLoggingFromOptions, logger, system, popenCatch

SYMBOL_NUMBER=4

class Hmm:
    def __init__(self, modelType):
        self.modelType=modelType
        self.stateNumber = { "fiveState":5, "threeState":3, "threeStateAsymmetric":3}[modelType]
        self.transitions = [0.0] * self.stateNumber**2
        self.emissions = [0.0] * (SYMBOL_NUMBER**2 * self.stateNumber)
        self.likelihood = 0.0
        
    def _modelTypeInt(self):
        return { "fiveState":0, "threeState":1, "threeStateAsymmetric":2}[self.modelType]

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
        self.likelihood = l[-1]
        self.transitions = map(lambda x : sum(x), zip(self.transitions, l[1:-1]))
        assert len(self.transitions) == self.stateNumber**2
        l = map(float, fH.readline().split())
        assert len(l) == len(self.emissions)
        self.emissions = map(lambda x : sum(x), zip(self.emissions, l))
        assert len(self.emissions) == self.stateNumber * SYMBOL_NUMBER**2
        fH.close()
        return self

    @staticmethod
    def loadHmm(file):
        fH = open(file, 'r')
        l = fH.readline().split()
        assert len(l) > 0
        fH.close()
        return Hmm({ 0:"fiveState", 1:"threeState", 2:"threeStateAsymmetric"}[int(l[0])]).addExpectationsFile(file)

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
    for iteration in xrange(options.iterations):
        #Temp file to store model
        modelsFile = os.path.join(target.getLocalTempDir(), "model.txt")
        hmm.write(modelsFile)
        #Run cactus_realign
        system("cat %s | cactus_realign --logLevel DEBUG %s --loadHmm=%s --outputExpectations=%s %s" % (alignments, sequences, modelsFile, modelsFile, options.optionsToRealign))
        #Add to the expectations.
        hmm = Hmm.loadHmm(modelsFile)
        hmm.normalise()
        #Do some logging
        target.logToMaster("After %i iteration got likelihood: %s for model-type: %s" % (iteration, hmm.likelihood, hmm.modelType))
    logger.info("The output file %s" % outputModel)
    hmm.write(outputModel)
    
def expectationMaximisationTrials(target, sequences, alignments, outputModel, options):
    if options.inputModel != None or not options.randomStart: #No multiple iterations
        target.setFollowOnTargetFn(expectationMaximisation, args=(sequences, alignments, outputModel, options))
    else:
        target.logToMaster("Running %s random restart trials to find best hmm" % options.trials)
        trialModels = map(lambda i : os.path.join(target.getGlobalTempDir(), "trial_%s.hmm" % i), xrange(options.trials))
        map(lambda trialModel : target.addChildTargetFn(expectationMaximisation, args=(sequences, alignments, trialModel, options)), trialModels)
        target.setFollowOnTargetFn(expectationMaximisationTrials2, args=(trialModels, outputModel, options))

def expectationMaximisationTrials2(target, trialModels, outputModel, options):
    #Pick the trial hmm with highest likelihood
    hmm = max(map(lambda x : Hmm.loadHmm(x), trialModels), key=lambda x : x.likelihood)
    hmm.write(outputModel)
    target.logToMaster("Hmm with highest likelihood: %s" % hmm.likelihood)
    
class Options:
    """Dictionary representing options, can be used for running pipeline from within another jobTree.
    """
    def __init__(self):
        self.modelType="fiveState"
        self.inputModel=None
        self.iterations=10
        self.trials=3
        self.randomStart=False
        self.optionsToRealign="--diagonalExpansion=10 --splitMatrixBiggerThanThis=3000" 

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: workingDir [options]", version="%prog 0.1")
    options = Options()
    parser.add_option("--sequences", dest="sequences", help="Quoted list of fasta files containing sequences")
    parser.add_option("--alignments", dest="alignments", help="Cigar file ")
    parser.add_option("--inputModel", default=options.inputModel, help="Input model")
    parser.add_option("--outputModel", default="hmm.txt", help="File to write the model in")
    parser.add_option("--modelType", default=options.modelType, help="Specify the model type, currently either fiveState, threeState, threeStateAsymmetric")
    parser.add_option("--iterations", default=options.iterations, help="Number of iterations of EM", type=int)
    parser.add_option("--trials", default=options.trials, help="Number of independent EM trials. The model with the highest likelihood will be reported. Will only work if randomStart=True", type=int)
    parser.add_option("--randomStart", default=options.randomStart, help="Iterate start model with small random values, else all values are equal", action="store_true")
    parser.add_option("--optionsToRealign", default=options.optionsToRealign, help="Further options to pass to cactus_realign when computing expectation values, should be passed as a quoted string")
    
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