from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqUtils

class SequenceAnalyser():
    def __init__(self,seq):
        self.seq=seq
            
    def orf(self):
        
        my_sequences=[] 
        for i in range(len(self.seq)): 
            if self.seq.startswith("ATG",i): 
                for j in range(i,len(self.seq),3): 
                #defining a condition to stop the loop if finds either TGA, TAA or TAG
                    if self.seq.endswith("TGA",i,j) or self.seq.endswith("TAA",i,j) or self.seq.endswith("TAG",i,j):
                        orf=self.seq[i:j]
                        my_sequences.append(orf)
                        break 
        longest_sequence=0
        longest_sequence_orf=my_sequences[0]
        for s in my_sequences: 
            if len(s)>longest_sequence: 
                longest_sequence=len(s)
                longest_sequence_orf=s
        return longest_sequence_orf
    
class Primer(SequenceAnalyser):    
    def __init__(self,seq):
        SequenceAnalyser.__init__(self,seq)
        self.seq=seq
        self.A=seq.count("A")
        self.T=seq.count("T")
        self.G=seq.count("G")
        self.C=seq.count("C")
        
    def meltingtemp(self):
        return 2*(self.A+self.T)+4*(self.G+self.C)
        
    def GCcontent(self):
        return (self.G+self.C)/len(self.seq)*100

def buildingprimer(seq):
    lengthRange=[18,19,20,21,22,23,24,25]
    count=0
    for i in lengthRange:
        startPos=0
        endPos=0
        while endPos<len(seq):
            endPos=startPos+i
            pr=seq[startPos:endPos]
            startPos=endPos
            p=Primer(pr)
            if ((p.meltingtemp() >= 55 and p.meltingtemp() <= 65) and
            (p.GCcontent() >=45 and p.GCcontent() <=55) and
            (primer3.calcHairpin(str(pr)).structure_found !=True)):
                print(pr)
                print ("Melting Temperature : {}".format(p.meltingtemp()))
                print ("GC content : {}".format(p.GCcontent()))
                count=count+1
    print(count)

for s in SeqIO.parse("id170383510.seq.txt","fasta"): #parsing the sequence file
    sequence=s.seq
    abc=SequenceAnalyser(sequence)
    buildingprimer(sequence)
    print(abc.orf())
