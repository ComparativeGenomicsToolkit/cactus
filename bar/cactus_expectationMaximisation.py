"""Script to train pair-HMM of Cactus.
"""

import os

class HmmExpectations:
    def __init__(self, stateNumber):
        self.logLikelihood = 0.0
        self.stateNumber = stateNumber
        self.transitionExpectations = [0.0] * (stateNumber * stateNumber)
    
    def addExpectationsFile(self, file):
        f = open(file, 'r')
        l = map(float, file.readline().split())
        assert len(l) == len(self.transitionExpectations)+1
        self.logLikelihood = l[-1]
        self.transitionExpectations = map(lambda x : sum(x), zip(self.transitionExpectations, l))
        f.close()
    
    def getMaximisedHmm(self):
        """Calculates maximisation step of EM algorithm for transitions from set of expectations 
        for the transitions, stored in an array representing a STATE_NUMBER**2 matrix such that
        expectations[FROM * STATE_NUMBER + TO] is the expectation of the transition from FROM to TO.
        """
        hmm = Hmm(self.stateNumber)
        for fromState in xrange(self.stateNumber):
            i = self.stateNumber * fromState
            j = self.transitionExpectations[i:i+stateNumber]
            for toState in xrange(self.stateNumber):
                self.hmm[i + toState] = self.transitionExpectations[i + toState] / j
        return hmm

class Hmm:
    def __init__(self, stateNumber):
        self.stateNumber = stateNumber
        self.transitions = [0.0] * (stateNumber * stateNumber)
    
    def write(self, file):
        f = open(file, 'w')
        f.write(" ".join(map(str, self.transitions)) + "\n")
        f.close()

def expectationMaximisation(target, sequences, alignments, outputModel, iterations):
    stateNumber = 5
    #Iteratively run cactus realign to get expectations and load model.
    hmm = Hmm(stateNumber)
    for iteration in xrange(iterations):
        #Temp file to store model
        modelsFile = os.path.join(target.getLocalTempDir(), "model.txt")
        hmm.write(modelsFile)
        #Temp output file to store expectations
        expectationsFile = os.path.join(target.getLocalTempDir(), "expectations.txt")
        #Run cactus_realign
        system("cat %s | cactus_realign %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=3000 --loadHmm=%s --outputExpectations=%s" % (alignments, sequences, modelsFile, expectationsFile))
        #Add to the expectations.
        hmmExpectations = HmmExpectations(stateNumber)
        hmmExpectations.addExpectationsFile(expectationsFile)
        #Get newly maximised HMM
        hmm = hmmExpectations.getMaximisedHmm()
        #Do some logging
        target.logToMaster("After %i iteration got likelihood: %s" % (hmmExpectations.logLikelihood))
    hmm.write(outputModel)

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: workingDir [options]", version="%prog 0.1")
    parser.add_option("--sequences", dest="sequences", help="Quoted list of fasta files containing sequences")
    parser.add_option("--alignments", dest="alignments", help="Cigar file ")
    parser.add_option("--outputModel", default="hmm.txt", help="File to write the model in")
    parser.add_option("--iterations", default="10", help="Number of iterations of EM")
    
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Expected no arguments, got %s arguments: %s" % (len(args), " ".join(args)))
    
    #Log the inputs
    logger.info("Got '%s' sequences, '%s' alignments file, '%s' output model and '%s' iterations of training" % (options.sequences, options.alignments, options.outputModel, options.iterations))

    #This line invokes jobTree  
    i = Stack(Target.makeTargetFn(expectationMaximisation, args=(options.sequences, options.alignments, options.outputModel, options.iterations))).startJobTree(options) 
    
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from cactus.bar.cactus_expectationMaximisation import *
    main()