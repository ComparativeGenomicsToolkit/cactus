#!/usr/bin/env python

"""Script to generate a series of random
"""

from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system
from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def runDbTestScript(options, firstKey=0, keyNumber=0, addRecords=False, setRecords=False, create=False):
    def fn(stringId, bool):
        if bool:
            return stringId
        return ""
    addRecords = fn("--addRecords", addRecords)
    setRecords = fn("--setRecords", setRecords)
    create = fn("--create", create)
    command = "dbTestScript --databaseConf '%s' --firstKey %s --keyNumber %s %s %s --minRecordSize %s --maxRecordSize %s %s" %\
    (options.databaseConf, firstKey, ketNumber, addRecords, setRecords, options.minRecordSize, options.maxRecordSize, create)
    system(command)

class SetupTheDatabase(Target):
    def __init__(self, options):
        self.options = options
        
    def run(self):
        runDBTestScript(self.options, create=True)
        
class AddKeysPhase(SetupTheDatabase):
    def run(self):
        keyIndex = 0
        for jobIndex in xrange(int(self.options.keysPerJob)):
            self.addChildTarget(AddKeys(self.options, keyIndex))
            keyIndex += int(self.options.keysPerJob)
        seld.addFollowOnTarget(SetKeysPhase(self.options))
    
class AddKeys(Target):
    def __init__(self, options, firstKey):
        self.options = options
        self.firstKey = firstKey
        
    def run(self):
        runDbTestScript(self.options, self.firstKey, self.options.keysPerJob, addRecords=True)

class SetKeysPhase(SetupTheDatabase):
    def run(self):
        keyIndex = 0
        for jobIndex in xrange(int(self.options.keysPerJob)):
            self.addChildTarget(SetKeys(self.options, keyIndex))
            keyIndex += int(self.options.keysPerJob)

class SetKeys(AddKeys):
    def run(self):
        runDbTestScript(self.options, self.firstKey, self.options.keysPerJob, setRecords=True)

def main():
    ##########################################
    #Construct the arguments.
    ##########################################

    parser = OptionParser()
    
    parser.add_option("--databaseConf", dest="databaseConf")
    parser.add_option("--keysPerJob", dest="keysPerJob")
    parser.add_option("--totalJobs", dest="totalJobs")
    parser.add_option("--minRecordSize", dest="minRecordSize")
    parser.add_option("--maxRecordSize", dest="maxRecordSize")
    
    Stack.addJobTreeOptions(parser)

    options, args = parser.parse_args()
    setLoggingFromOptions(options)

    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))
    
    Stack(SetupTheDatabase(options)).startJobTree(options)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.dbTest.dbTestScript import *
    _test()
    main()