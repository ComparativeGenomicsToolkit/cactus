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
            "cigar: simpleSeqBC 0 9 + simpleSeqAC 10 0 - 0 M 8 D 1 M 1" ,
            "cigar: simpleSeqB1 18 9 - simpleSeqA1 6 2 - 0 M 3 I 5 M 1",
            "cigar: simpleSeqB1 32 30 - simpleSeqA2 7 9 + 0 M 2",
            "cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 0 M 3 I 5 M 1",
            "cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 0 M 1 I 2 M 2 D 2 M 5",
            "cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 0 M 2",
            "cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 0 M 1",
            "cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 0 M 5",
            "cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 0 M 5",
            "cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 0 M 5",
            "cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0 M 10" ]
        with open(self.simpleInputCigarPath, 'w') as fH:
            fH.write("\n".join(self.inputCigars) + "\n")
            
        self.simpleOutputCigarPath = getTempFile()
        self.logLevelString = "DEBUG"

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        os.remove(self.simpleInputCigarPath)
        os.remove(self.simpleOutputCigarPath)
    
    @staticmethod 
    def makeCigar(coordinates1, coordinates2, score, ops):
        cigar = "cigar: %s %s %f %s" % (" ".join(map(str, coordinates2)), 
                                        " ".join(map(str, coordinates1)),              
                                        float(score), " ".join(map(str, ops)))
        return cigar

    @silentOnSuccess
    def testMirrorAndOrientAlignments(self):
        cactus_call(parameters=["cactus_mirrorAndOrientAlignments", 
                                 self.logLevelString, 
                                 self.simpleInputCigarPath, 
                                 self.simpleOutputCigarPath])
        with open(self.simpleOutputCigarPath, 'r') as fh:
            outputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines
            
        # For each input alignment check that we have the two, oriented alignments
        for inputCigar in self.inputCigars: 
            name1, start1, end1, strand1 = inputCigar.split()[5:9]
            start1, end1 = int(start1), int(end1)
            coordinates1 = name1, start1, end1, strand1
            
            name2, start2, end2, strand2 = inputCigar.split()[1:5]
            start2, end2 = int(start2), int(end2)
            coordinates2 = name2, start2, end2, strand2
            
            score = inputCigar.split()[9]
            ops = inputCigar.split()[10:]
            
            def invertStrand(coordinates):
                # cigar: simpleSeqB1 0 9 + simpleSeqA1 10 0 - 0 M 8 D 1 M 1
                # cigar: simpleSeqB1 9 0 + simpleSeqA1 0 10 - 0 M 1 D 1 M 8
                name, start, end, strand = coordinates
                assert strand in ("+", "-")
                if strand == "+":
                    return name, end-1, start-1, "-" 
                return name, end+1, start+1, "+" 
            
            def reverseOps(ops):
                l = ops[:]
                l.reverse()
                l2 = []
                for i, j in zip(l[1::2], l[::2]):
                    l2 += [ i, j ]
                return l2
            
            def invertOpStrands(ops):
                l = [ "I" if op == "D" else ("D" if op == "I" else op) for op in ops[::2] ]
                l2 = []
                for op, length in zip(l, ops[1::2]):
                    l2 += [ op, length ]
                return l2
            
            if strand1 == "+":
                self.assertTrue(self.makeCigar(coordinates1, coordinates2, score, ops) in outputCigars)
            else:
                # Invert the strands
                self.assertTrue(self.makeCigar(invertStrand(coordinates1), 
                                          invertStrand(coordinates2), score, reverseOps(ops)) in outputCigars)
                
            if strand2 == "+":
                self.assertTrue(self.makeCigar(coordinates2, coordinates1, score,
                                          invertOpStrands(ops)) in outputCigars)
            else:
                self.assertTrue(self.makeCigar(invertStrand(coordinates2), invertStrand(coordinates1), 
                                            score, invertOpStrands(reverseOps(ops))) in outputCigars)
        
    @silentOnSuccess
    def testBlast_sortAlignmentsByQuery(self):
        cactus_call(parameters=["cactus_blast_sortAlignmentsByQuery", 
                                 self.logLevelString, 
                                 self.simpleInputCigarPath, 
                                 self.simpleOutputCigarPath])
        
        with open(self.simpleOutputCigarPath, 'r') as fh:
            outputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines
        
        correctCigars = [
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 0 M 3 I 5 M 1',
            'cigar: simpleSeqB1 18 9 - simpleSeqA1 6 2 - 0 M 3 I 5 M 1',
            'cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 0 M 1',
            'cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 0 M 1 I 2 M 2 D 2 M 5',
            'cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 0 M 2',
            'cigar: simpleSeqB1 32 30 - simpleSeqA2 7 9 + 0 M 2',
            'cigar: simpleSeqBC 0 9 + simpleSeqAC 10 0 - 0 M 8 D 1 M 1',
            'cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0 M 10',
            'cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 0 M 5',
            'cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 0 M 5',
            'cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 0 M 5'
            ]
        self.assertEqual(correctCigars, outputCigars)
        
    @silentOnSuccess
    def testSplitAlignmentsOverlaps(self):
        self.inputCigars = [ 'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 0.000000 M 3 I 5 M 1',
                        'cigar: simpleSeqB1 10 19 + simpleSeqA1 3 7 + 0.000000 M 1 I 5 M 3',
                        'cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 0.000000 M 1',
                        'cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 0.000000 M 1 I 2 M 2 D 2 M 5',
                        'cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 0.000000 M 2',
                        'cigar: simpleSeqB1 32 30 - simpleSeqA2 7 9 + 0.000000 M 2',
                        'cigar: simpleSeqBC 8 -1 - simpleSeqAC 1 11 + 0.000000 M 1 D 1 M 8',
                        'cigar: simpleSeqA1 2 6 + simpleSeqB1 9 18 + 0.000000 M 3 D 5 M 1',
                        'cigar: simpleSeqA1 3 7 + simpleSeqB1 10 19 + 0.000000 M 1 D 5 M 3',
                        'cigar: simpleSeqA2 0 10 + simpleSeqB1 18 28 + 0.000000 M 1 D 2 M 2 I 2 M 5',
                        'cigar: simpleSeqA2 6 8 + simpleSeqB1 28 30 + 0.000000 M 2',
                        'cigar: simpleSeqA2 8 6 - simpleSeqB1 31 33 + 0.000000 M 2',
                        'cigar: simpleSeqAC 10 0 - simpleSeqBC 0 9 + 0.000000 M 8 I 1 M 1',
                        'cigar: simpleSeqD 0 5 + simpleSeqC1 0 5 + 0.000000 M 5',
                        'cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0.000000 M 10',
                        'cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 0.000000 M 5',
                        'cigar: simpleSeqC1 15 20 + simpleSeqC1 10 15 + 0.000000 M 5',
                        'cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 0.000000 M 5',
                        'cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 0.000000 M 5',
                        'cigar: simpleSeqC1 5 10 + simpleSeqD 5 10 + 0.000000 M 5',
                        'cigar: simpleSeqC1 0 10 + simpleSeqNonExistent 0 10 + 0.000000 M 10',
                        'cigar: simpleSeqA1 6 7 + simpleSeqZ1 0 1 + 0.000000 M 1' ]
        with open(self.simpleInputCigarPath, 'w') as fH:
            fH.write("\n".join(self.inputCigars) + "\n")
            
        cactus_call(parameters=["cactus_splitAlignmentOverlaps", 
                                 self.logLevelString, 
                                 self.simpleInputCigarPath, 
                                 self.simpleOutputCigarPath])
        
        with open(self.simpleOutputCigarPath, 'r') as fh:
            outputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines
        
        # Get start and end coordinates of cigars
        ends = set()
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigar.split()[5:9]
            ends.add((name1, int(start1)))
            ends.add((name1, int(end1)))
            assert strand1 == "+"
        
        # Count of expected number of chopped up cigars 
        totalExpectedCigars = 0
        
        # Function to split a list of ops into a prefix and suffix list
        def splitPrefixOps(ops, cutPoint):
            pOps, sOps = [], []
            j = 0
            for i in range(0, len(ops), 2):
                op, length = ops[i], int(ops[i+1])
                assert op in ("I", "D", "M")
                if op == "I":
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
                    pOps.append(op)
                    pOps.append(cutPoint - j)
                    sOps.append(op)
                    sOps.append(length - (cutPoint - j))
                    break
            sOps += ops[i+2:]  
            
            return pOps, sOps     
        
        # For each cigar:
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigar.split()[5:9]
            start1, end1 = int(start1), int(end1)
            assert strand1 == "+"
            name2, start2, end2, strand2 = inputCigar.split()[1:5]
            start2, end2 = int(start2), int(end2)
            score = float(inputCigar.split()[9])
            ops = inputCigar.split()[10:]
            
            # For each intermediate chop point
            i = start1
            for j in xrange(start1+1, end1+1):
                if (name1, j) in ends:
                    # Chop up cigar 
                    coordinates1 = name1, i, j, "+"
                    
                    # Get sublist of ops
                    pOps, subOps = splitPrefixOps(ops, i - start1)
                    subOps, sOps = splitPrefixOps(subOps, j - i)
                    
                    x = lambda ops : sum([ int(ops[k+1]) for k in range(0, len(ops), 2) if ops[k] != 'D' ])
                    k = x(pOps)
                    l = k + x(subOps)
                    
                    # Get second coordinates
                    if strand2 == "+":
                        coordinates2 = name2, start2 + k, start2 + l, strand2
                    else:
                        assert strand2 == "-"
                        coordinates2 = name2, start2 - k, start2 - l, strand2
                    
                    choppedCigar = self.makeCigar(coordinates1, coordinates2, score, subOps)
                    
                    # Check each chopped up cigar is in output
                    self.assertTrue(choppedCigar in outputCigars)
                    
                    # Inc. number of expected cigars
                    totalExpectedCigars += 1
                    
                    # Check previous coordinate
                    i = j
                    
        # Check we have the expected number of cigars  
        self.assertEquals(totalExpectedCigars, len(outputCigars))

if __name__ == '__main__':
    unittest.main()
