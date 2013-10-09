import unittest

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testCactusRealign(self):
        seqFile1 = ""
        seqFile2 = ""
        system("lastz")
        
if __name__ == '__main__':
    unittest.main()