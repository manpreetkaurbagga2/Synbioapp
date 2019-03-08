import assignment2
import unittest
from Bio import SeqIO

'''PRESS Enter to get the TestResult'''

class TestSequenceAnalyser(unittest.TestCase):
    
    def setUp(self):
        for s in SeqIO.parse("examplefile.seq","fasta"): #parsing the sequence file
            self.Sequence=s.seq
        self.A=assignment2.SequenceAnalyser(self.Sequence)
        self.B=assignment2.Rbs(self.Sequence)
        self.C=assignment2.RestrictionSites(self.Sequence)
    
    def testorf(self):
       orf=self.A.orf()
       self.assertTrue(len(orf)>0)
       
    def testrbs(self):
        self.B.RbsPosition()
        actual=self.B.RbsPosition()[0]
        expected="AAGGAGGTG"
        self.assertEqual(actual,expected)
        
    def test_EcoRI(self):
        self.C.Pos_EcoRI()
        actual=self.C.Pos_EcoRI()[0]
        expected="GAATTC"
        self.assertEqual(actual,expected)
        
    def test_SpeI(self):
        self.C.Pos_SpeI()
        actual=self.C.Pos_SpeI()[0]
        expected="ACTAGT"
        self.assertEqual(actual,expected)
        
    def test_XbaI(self):
        self.C.Pos_XbaI()
        actual=self.C.Pos_XbaI()[0]
        expected="TCTAGA"
        self.assertEqual(actual,expected)
        
    def test_PstI(self):
        self.C.Pos_PstI()
        actual=self.C.Pos_PstI()[0]
        expected="CTGCAG"
        self.assertEqual(actual,expected)
            
if __name__ == "__main__":
        unittest.main()