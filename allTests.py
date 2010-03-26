import unittest

import cactus.setup.cactus_setupTest
import cactus.blastAlignment.cactus_alignerTest
import cactus.blastAlignment.cactus_batchTest
import cactus.core.cactus_coreTest
import cactus.pipeline.cactus_workflowTest
import cactus.baseAlignment.cactus_baseAlignerTest
import cactus.phylogeny.cactus_phylogenyTest
import cactus.reference.cactus_referenceTest

from sonLib.bioio import parseSuiteTestOptions

def allSuites(): 
    cactus_setupSuite = unittest.makeSuite(cactus.setup.cactus_setupTest.TestCase, 'test')
    cactus_alignerSuite = unittest.makeSuite(cactus.blastAlignment.cactus_alignerTest.TestCase, 'test')
    cactus_batchSuite = unittest.makeSuite(cactus.blastAlignment.cactus_batchTest.TestCase, 'test')
    cactus_coreSuite = unittest.makeSuite(cactus.core.cactus_coreTest.TestCase, 'test')
    cactus_workflowSuite = unittest.makeSuite(cactus.pipeline.cactus_workflowTest.TestCase, 'test')
    cactus_baseAlignerSuite = unittest.makeSuite(cactus.baseAlignment.cactus_baseAlignerTest.TestCase, 'test')
    cactus_phylogenySuite = unittest.makeSuite(cactus.phylogeny.cactus_phylogenyTest.TestCase, 'test')
    cactus_referenceSuite = unittest.makeSuite(cactus.reference.cactus_referenceTest.TestCase, 'test')
     
    allTests = unittest.TestSuite((cactus_setupSuite, cactus_alignerSuite, 
                                   cactus_coreSuite, 
                                   cactus_workflowSuite, 
                                   cactus_baseAlignerSuite, cactus_phylogenySuite,
                                   cactus_batchSuite, cactus_referenceSuite))
    return allTests
        
def main():
    parseSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()

import re
import string
import os
suites = []
def fn():
    for f in os.listdir(dir):
        if f[-6:] == "Test.c":
            print "doing file", f
            fH = open(os.path.join(dir, f), 'w')
            fH.write('#include "cactusGlobalsPrivate.h"\n')
            fH.write('\ninline void %sTestSetup() {\n\n}\n\nvoid %sTestTeardown() {\n\n}' % (f[:-6], f[:-6]))
            fH.write('\n\n')
            fH2 = open(os.path.join(dir, f[:-6] + ".h"), 'r')
            functions = []
            for line in fH2.readlines():
                print "got line", line
                p = re.compile("(.*) \*?(.)(.*)_([^\( \t]*)\((.*)\);")
                if p.match(line):
                    g = p.match(line)
                    print "using line", line
                    fH.write("void test%s%s_%s(CuTest* testCase) {\n\t%sTestSetup();\n\t%sTestTeardown();\n}\n\n" % (string.upper(g.group(2)), g.group(3), g.group(4), f[:-6], f[:-6]))
                    functions.append("test%s%s_%s" % (string.upper(g.group(2)), g.group(3), g.group(4)))
            suites.append("%sTestSuite" % f[:-6])
            fH.write("CuSuite* %sTestSuite(void) {\n\tCuSuite* suite = CuSuiteNew();\n" % f[:-6])
            for function in functions:
                fH.write("\tSUITE_ADD_TEST(suite, %s);\n" % function)
            fH.write("\treturn suite;\n}\n")
            fH.close()
            fH2.close()
        
    fH = open(os.path.join(dir, "allTests.c"), 'w')
    fH.write('#include "cactusGlobalsPrivate.h"\n\n')
    for suite in suites:
        fH.write('CuSuite *%s();\n' % suite)
    
    fH.write('\nvoid cactusAPIRunAllTests(void) {\n\tCuString *output = CuStringNew();\n\tCuSuite* suite = CuSuiteNew();\n')
    for suite in suites:
        fH.write('\tCuSuiteAddSuite(suite, %s());\n' % suite)
    fH.write('\tCuSuiteRun(suite);\n\tCuSuiteSummary(suite, output);\n\tCuSuiteDetails(suite, output);\n\tprintf("%s\\n", output->buffer);\n')
    fH.write('}\n\n')
    fH.write("int main(void) {\n\tcactusAPIRunAllTests();\n\treturn 0;\n}\n")
    fH.close()
                
