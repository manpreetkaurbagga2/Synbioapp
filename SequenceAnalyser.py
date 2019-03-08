from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqUtils
import primer3

class SequenceAnalyser():
    def __init__(self,seq):
        self.seq=seq
        self.ORF_Lenght=0
    
    def orf (self):
        my_sequences=[] 
        for i in range(len(self.seq)): 
            if self.seq.startswith("ATG",i):  
                for j in range(i,len(self.seq),3):                 
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
        self.ORF_Lenght=len(longest_sequence_orf)
        return longest_sequence_orf
    
    def __getattribute__(self, ORF_Lenght):
        return super().__getattribute__(ORF_Lenght)
   
class Primer(SequenceAnalyser):    
    def __init__(self,seq,A,T,G,C):
        SequenceAnalyser.__init__(self,seq)
        self.A=seq.count("A")
        self.T=seq.count("T")
        self.G=seq.count("G")
        self.C=seq.count("C")
        
    def meltingtemp(self):
        return 2*(self.A+self.T)+4*(self.G+self.C)
        
    def GCcontent(self):
        return (self.G+self.C)/len(self.seq)*100

def buildPrimer(seq,length):
    generatedPrimers={}
    clippedForwardSeq = seq[0:length]
    clippedReverseSeq = seq[len(seq)-length:len(seq)]
    forward_Complement = clippedForwardSeq.complement()
    reverse_Complement =clippedReverseSeq.reverse_complement()
    Sequences = [forward_Complement,reverse_Complement]
    OutputString = ['Forward Primer\n ','Reverse Primer\n']
    for l in range(2):
        primers=[]
        clippedSeq = Sequences[l]
        lengthRange=[18,19,20,21,22,23,24,25]
        count=0
        for i in lengthRange:
            startPos=0
            endPos=0
            while endPos<len(clippedSeq):
                primerObject={}
                endPos=startPos+i
                pr=clippedSeq[startPos:endPos]
                startPos=endPos
                p=Primer(pr)
                if ((p.meltingtemp() >= 55 and p.meltingtemp() <= 65) and
                (p.GCcontent() >=45 and p.GCcontent() <=55) and
                (primer3.calcHairpin(str(pr)).structure_found !=True)):
                    primerObject['Primer']=str(pr)
                    primerObject['MeltingTemperature'] = p.meltingtemp()
                    primerObject['GCcontent']=p.GCcontent()
                    primers.append(primerObject)
                    count=count+1
        generatedPrimers[OutputString[l]]=primers
    return generatedPrimers 

class Rbs(SequenceAnalyser):    
    def __init__(self,seq):
        SequenceAnalyser.__init__(self,seq)
        self.rbs=Seq("AAGGAGGTG")
   
    def RbsPosition(self):
        rbs_pos=SeqUtils.nt_search(str(self.seq),self.rbs)
        return rbs_pos

class RestrictionSites(SequenceAnalyser):
    def __init__(self,seq):
        SequenceAnalyser.__init__(self,seq)
        self.EcoRI=Seq("GAATTC")
        self.XbaI=Seq("TCTAGA")
        self.SpeI=Seq("ACTAGT")
        self.PstI=Seq("CTGCAG")
        
    def Pos_EcoRI(self):
        EcoRI_pos=SeqUtils.nt_search(str(self.seq),self.EcoRI)
        return EcoRI_pos
    
    def Pos_XbaI(self):
        XbaI_pos=SeqUtils.nt_search(str(self.seq),self.XbaI)
        return XbaI_pos

    def Pos_SpeI(self):
        SpeI_pos=SeqUtils.nt_search(str(self.seq),self.SpeI)
        return SpeI_pos

    def Pos_PstI(self):
        PstI_pos=SeqUtils.nt_search(str(self.seq),self.PstI)
        return PstI_pos    
  
def menu():
    print("Menu :\n")
    print("Enter 1 to Find ORF of the Given Sequence")
    print("Enter 2 to Find RBS")
    print("Enter 3 to Find Biobrick compatible Restriction Sites")
    print("Enter 4 to Find Forward and Reverse Primers")
    print("Enter q to Quit.")
    choice = input("Enter your Choice ?")
    return choice

   
Sequence = Seq("")
for s in SeqIO.parse("examplefile.seq","fasta"): #parsing the sequence file
    Sequence=s.seq

choice = menu()
while (choice != 'q'):
    if choice == '1':
        a=SequenceAnalyser(Sequence)
        print(a.orf(),file=open("output.txt", "a"))
        more = input("Do you want to get the lenght of the ORF generated ?(Y/N)")
        if more =='y' or more =='Y':
            print("The Lenght of Genereted ORF : {}".format(a.ORF_Lenght),file=open("output.txt", "a"))
            choice = menu()
        else:
            choice =menu()
    elif choice == '2':
        a = Rbs(Sequence)
        print("RBS position : {}\n".format(a.RbsPosition()),file=open("output.txt", "a"))
        choice = menu()
    elif choice == '3':
        a=RestrictionSites(Sequence)
        print("EcoRI position : {}\n".format(a.Pos_EcoRI()),file=open("output.txt", "a"))
        print("XbaI position : {}\n".format(a.Pos_PstI()),file=open("output.txt", "a"))
        print("SpeI position : {}\n".format(a.Pos_SpeI()),file=open("output.txt", "a"))
        print("PstI position : {}\n".format(a.Pos_XbaI()),file=open("output.txt", "a"))
        choice = menu ()
    elif choice == '4':
        lengthToQuery = input("What lenght you want to query from forward and reverse?")
        print("Primers",file=open("output.txt", "a"))
        buildPrimer(Sequence,int(lengthToQuery))
        choice = menu()
    else:
        break
      
