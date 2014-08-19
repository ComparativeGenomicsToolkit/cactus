#!/usr/bin/env python

"""Script to train pair-HMM of Cactus.
"""

import math
import os
import random
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import setLoggingFromOptions, logger, system, popenCatch, cigarRead, cigarWrite, nameValue

SYMBOL_NUMBER=4

class Hmm:
    def __init__(self, modelType):
        self.modelType=modelType
        self.stateNumber = { "fiveState":5, "fiveStateAsymmetric":5, "threeState":3, "threeStateAsymmetric":3}[modelType]
        self.transitions = [0.0] * self.stateNumber**2
        self.emissions = [0.0] * (SYMBOL_NUMBER**2 * self.stateNumber)
        self.likelihood = 0.0
        
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
    cigars = list(cigarRead(alignments))
    splitAlignments = []
    while len(cigars) > 0:
        splitAlignments.append(os.path.join(target.getGlobalTempDir(), "alignments_%s.cigar" % len(splitAlignments)))
        fH = open(splitAlignments[-1], 'w')
        for cigar in cigars[:options.numberOfAlignmentsPerJob]:
            cigarWrite(fH, cigar)
        fH.close()
        cigars = cigars[options.numberOfAlignmentsPerJob:]

    #Files to store expectations in
    expectationsFiles = map(lambda i : os.path.join(target.getGlobalTempDir(), "expectation_%i.txt" % i), xrange(len(splitAlignments)))
    assert len(splitAlignments) == len(expectationsFiles)

    target.setFollowOnTargetFn(expectationMaximisation2, args=(sequences, splitAlignments, outputModel, expectationsFiles, 0, options))

def expectationMaximisation2(target, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, options):
    if iteration < options.iterations:
        map(lambda x : target.addChildTargetFn(calculateExpectations,
                    args=(sequences, x[0], None if (options.useDefaultModelAsStart and iterations == 0) else modelsFile, x[1], options)), 
            zip(splitAlignments, expectationsFiles))
        target.setFollowOnTargetFn(calculateMaximisation, args=(sequences, splitAlignments, modelsFile, expectationsFiles, iteration, options))

def calculateExpectations(target, sequences, alignments, modelsFile, expectationsFile, options):
    #Run cactus_realign
    system("cat %s | cactus_realign --logLevel DEBUG %s %s --outputExpectations=%s %s" % (alignments, sequences, nameValue("loadHmm", modelsFile, str), expectationsFile, options.optionsToRealign))

def calculateMaximisation(target, sequences, splitAlignments, modelsFile, expectationsFiles, iteration, options):
    #Load and merge the models files
    if len(expectationsFiles) > 0:
        hmm = Hmm.loadHmm(expectationsFiles[0])
        for expectationsFile in expectationsFiles[1:]:
            hmm.addExpectationsFile(expectationsFile)
        hmm.normalise()
        target.logToMaster("On %i iteration got likelihood: %s for model-type: %s" % (iteration, hmm.likelihood, hmm.modelType))
        #If not train emissions then load up the old emissions and replace
        if options.trainEmissions:
            target.logToMaster("Training expectations")
        else:
            hmm.emissions = Hmm.loadHmm(modelsFile).emissions
            target.logToMaster("Using the original emissions")
            
        #Write out 
        hmm.write(modelsFile)
    
    #Generate a new set of alignments, if necessary
    if options.updateTheBand:
        map(lambda alignments : target.addChildTargetFn(calculateAlignments, args=(sequences, alignments, modelsFile, options)), splitAlignments)
    
    #Start the next iteration
    target.setFollowOnTargetFn(expectationMaximisation2, args=(sequences, splitAlignments, modelsFile, expectationsFiles, iteration+1, options))
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
        self.updateTheBand=False
        self.numberOfAlignmentsPerJob=100
        self.useDefaultModelAsStart = False
        self.setJukesCantorStartingEmissions=None
        self.trainEmissions=False

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
    parser.add_option("--updateTheBand", default=options.updateTheBand, help="After each iteration of EM update the set of alignments by realigning them, so allowing stochastic updating of the constraints. This does not alter the input alignments file", action="store_true")
    parser.add_option("--numberOfAlignmentsPerJob", default=options.numberOfAlignmentsPerJob, help="Number of alignments to compute in parallel during expectation/updating the band steps", type=int)
    parser.add_option("--useDefaultModelAsStart", default=options.useDefaultModelAsStart, help="Use the default BAR hmm model as the starting point", action="store_true")
    parser.add_option("--setJukesCantorStartingEmissions", default=options.setJukesCantorStartingEmissions, help="[double] Set the starting hmm emissions by jukes cantor expectation, using given subs/site estimate", type=float)
    parser.add_option("--trainEmissions", default=options.trainEmissions, help="Train the emissions as well as the transitions.", action="store_true")
    
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