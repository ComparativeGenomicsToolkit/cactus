import unittest, os, random
from sonLib.bioio import getTempFile
from textwrap import dedent
from cactus.shared.common import cactus_call
from cactus.shared.test import getCactusInputs_encode, silentOnSuccess

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # simple test data
        self.simpleInputCigarPath = getTempFile()
        self.inputCigars = [
            "cigar: simpleSeqB1 0 9 + simpleSeqA1 10 0 - 0 M 8 D 1 M 1" ] #,
            #"cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 0 M 3 I 5 M 1",
            #"cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 0 M 1 I 2 M 2 D 2 M 5",
            #"cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 0 M 2",
            #"cigar: simpleSeqB1 30 32 + simpleSeqA2 7 9 + 0 M 2",
            #"cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 0 M 1",
            #"cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 0 M 5",
            #"cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 0 M 5",
            #"cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 0 M 5",
            #"cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0 M 10" ]
        open(self.simpleInputCigarPath, 'w').write("\n".join(self.inputCigars))
        self.simpleOutputCigarPath = getTempFile()
        self.logLevelString = "DEBUG"

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        os.remove(self.simpleInputCigarPath)
        os.remove(self.simpleOutputCigarPath)

    @silentOnSuccess
    def testMirrorAndOrientAlignments(self):
        cactus_call(parameters=["cactus_mirrorAndOrientAlignments", 
                                 self.logLevelString, 
                                 self.simpleInputCigarPath, 
                                 self.simpleOutputCigarPath])
        outputCigars = open(self.simpleOutputCigarPath).readlines()
        
        print "Input cigars", self.inputCigars
        
        print "Got", "\n".join(outputCigars)

        # For each input alignment check that we have the two, oriented alignments
        for inputCigar in self.inputCigars: 
            name1, start1, end1, strand1 = inputCigar.split()[1:5]
            start1, end1 = int(start1), int(end1)
            coordinates1 = name1, start1, end1, strand1
            name2, start2, end2, strand2 = inputCigar.split()[5:9]
            start2, end2 = int(start2), int(end2)
            coordinates2 = name2, start2, end2, strand2
            score = inputCigar.split()[9]
            ops = inputCigar.split()[10:]
            
            def invertStrand(name, start, end, strand):
                # cigar: simpleSeqB1 0 9 + simpleSeqA1 10 0 - 0 M 8 D 1 M 1
                # cigar: simpleSeqB1 9 0 + simpleSeqA1 0 10 - 0 M 1 D 1 M 8
                assert strand in ("+", "-")
                if strand == "+":
                    return name, end-1, start-1, "+" 
                return name, end+1, start+1, "-" 

            def makeCigar(coordinates1, coordinates2, ops):
                cigar = "cigar: %s %s %s %s" % (coordinates1, coordinates2, score, " ".join(ops))
                print "boo", cigar
                return cigar
            
            def reverseOps(ops):
                l = ops[:]
                l.reverse()
                l2 = []
                for i, j in zip(l[1::2], l[::2]):
                    l2 += [ i, j ]
                return l2
            
            def invertOpsStrands(ops):
                l = [ "I" if op == "D" else ("D" if op == "I" else op) for op in ops[::2] ]
                l2 = []
                for op, length in zip(l, ops[1::2]):
                    l2 += [ op, length ]
                return l2
            
            if strand1 == "+":
                print "goo", inputCigar
                self.assertTrue(inputCigar in outputCigars)
                self.assertTrue(makeCigar(coordinates1, coordinates2, ops) in outputCigars)
            else:
                # Invert the strands
                self.assertTrue(makeCigar(invertStrand(coordinates1), 
                                          invertStrand(coordinates2), reverseOps(ops)) in outputCigars)
                
            if strand2 == "+":
                self.assertTrue(makeCigar(coordinates2, coordinates1, 
                                          invertOpsStrands(ops)) in outputCigars)
            else:
                self.assertTrue(makeCigar(invertStrand(coordinates2), invertStrand(coordinates1), 
                                            invertOpStrands(reverseOps(ops))) in outputCigars)
        
    @silentOnSuccess
    def testBlast_sortAlignmentsByQuery(self):
        cactus_call(parameters=["cactus_blast_sortAlignmentsByQuery", 
                                 self.logLevelString, 
                                 self.simpleInputCigarPath, 
                                 self.simpleOutputCigarPath])
        cigars = open(self.simpleOutputCigarPath).readlines()
        correctCigars = [
            'cigar: simpleSeqB1 0 9 + simpleSeqA1 10 0 - 0 M 8 D 1 M 1\n', 
            'cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 0 M 1 I 2 M 2 D 2 M 5\n', 
            'cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 0 M 2\n', 
            'cigar: simpleSeqB1 30 32 + simpleSeqA2 7 9 + 0 M 2\n', 
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 0 M 3 I 5 M 1\n', 
            'cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 0 M 5\n', 
            'cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 0 M 5\n', 
            'cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 0 M 5\n', 
            'cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0 M 10\n', 
            'cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 0 M 1\n']
        self.assertEqual(correctCigars, cigars)
        
    @silentOnSuccess
    def testSplitAlignmentsOverlaps(self):
        return
        cactus_call(parameters=["cactus_splitAlignmentOverlaps", 
                                 self.logLevelString, 
                                 self.simpleInputCigarPath, 
                                 self.simpleOutputCigarPath])
        outputCigars = open(self.simpleOutputCigarPath).readlines()
        
        # Get start and end coordinates of cigars
        ends = set()
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigarString.split()[1:5]
            ends.add((name1, int(start1)))
            ends.add((name1, int(end1)))
            
        # Count of expected number of chopped up cigars 
        totalExpectedCigars = 0
        
        # Function to split a list of ops into a prefix and suffix list
        def splitPrefixOps(ops, cutPoint):
            pOps, sOps = [], []
            j = 0
            for i in range(0, len(ops), 2):
                op, length = ops[i], int(ops[i+1])
                assert op in ("I", "D", "M")
                if op == "D":
                    pOps.append(op)
                    pOps.append(length)
                    continue
                if j + length <= cutPoint:
                    pOps.append(op)
                    pOps.append(length)
                    j += length
                    if j == cutPoint:
                        break
                else:
                    assert j + length > cutPoint
                    assert j < cutPoint
                    pOps.append(op)
                    pOps.append(cutPoint - j)
                    sOps.append(op)
                    sOps.append(length - (cutPoint - j))
                    break
            sOps += ops[i:]  
            
            return pOps, sOps     
        
        # For each cigar:
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigarString.split()[1:5]
            start1, end1 = int(start1), int(end1)
            assert strand == "+"
            name2, start2, end2, strand2 = inputCigarString.split()[5:9]
            start2, end2 = int(start2), int(end2)
            score = inputCigarString.split()[9]
            ops = inputCigarString.split()[10:]
            
            # For each intermediate chop point
            i = start1
            for j in xrange(start1+1, end1):
                if (name1, j) in ends:
                    # Chop up cigar 
                    coordinates1 = name1, i, j, "+"
                    
                    # Get sublist of ops
                    pOps, subOps = splitPrefixOps(ops, i - start1)
                    subOps, sOps = splitPrefixOps(subOps, j - i)
                    
                    k = sum([ l for op, l in pOps if l != "I"])
                    l = k + sum([ l for op, l in subOps if l != "I"])
                    
                    # Get second coordinates
                    if strand2 == "+":
                        coordinates2 = name2, start2 + k, start2 + l, strand2
                    else:
                        assert strand2 == "-"
                        coordinates2 = name2, start2 - k, start2 - l, strand2
                    
                    choppedCigar = "cigar: %s %s %s %s" % (coordinates1, coordinates2, score, " ".join(subOps))
                    
                    # Check each chopped up cigar is in output
                    self.assertTrue(testCase, choppedCigar in outputCigars)
                    
                    # Inc. number of expected cigars
                    totalExpectedCigars += 1
                    
                    # Check previous coordinate
                    i = j
                    
        # Check we have the expected number of cigars  
        self.assertEquals(totalExpectedCigars, len(outputCigars))

if __name__ == '__main__':
    unittest.main()
