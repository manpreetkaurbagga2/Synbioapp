import assignment2
from Bio.Seq import Seq
from Bio import SeqIO

      
Sequence = Seq("") 
for s in SeqIO.parse("examplefile.seq","fasta"): #parsing the sequence file
    Sequence=s.seq #parsing the sequence file
choice = assignment2.menu() #calling function menu
while (choice != 'q'): #specifying the outputs for various options in menu and what output should be generated
    if choice == '1': #choice 1 should give the orf output
        a=assignment2.SequenceAnalyser(Sequence)
        print(a.orf(),file=open("output.txt", "a")) #exporting the orf in an output.txt file
        more = input("Do you want to get the lenght of the ORF generated ?(Y/N)") #specifying sub-option for orf which needs another input to print the length of orf
        if more =='y' or more =='Y':
            print("The Lenght of Genereted ORF : {}".format(a.ORF_Lenght),file=open("output.txt", "a"))#the legth is also appended in the output file
            choice = assignment2.menu()
        else:
            choice =assignment2.menu()
    elif choice == '2': #choice 2 should give the rbs output
        a = assignment2.Rbs(Sequence)
        print("RBS position : {}\n".format(a.RbsPosition()),file=open("output.txt", "a"))#gives the result in output file
        choice = assignment2.menu()
    elif choice == '3': #choice 3 should give the restriction sites output
        a=assignment2.RestrictionSites(Sequence)
        print("EcoRI position : {}\n".format(a.Pos_EcoRI()),file=open("output.txt", "a"))#gives the result in output file
        print("XbaI position : {}\n".format(a.Pos_PstI()),file=open("output.txt", "a"))
        print("SpeI position : {}\n".format(a.Pos_SpeI()),file=open("output.txt", "a"))
        print("PstI position : {}\n".format(a.Pos_XbaI()),file=open("output.txt", "a"))
        choice = assignment2.menu ()
    elif choice == '4': #choice 4 should give the primers
        lengthToQuery = input("What lenght you want to query from forward and reverse?")#a sub-option to specify the length for which the primers should be made
        print("Primers",file=open("output.txt", "a"))#gives the result in output file
        gp = assignment2.buildPrimer(Sequence,int(lengthToQuery))        
        for i in gp:
            print(i,file=open("output.txt","a"))
            for j in gp[i]:
                print(j,file=open("output.txt","a"))
        choice = assignment2.menu()
    else:
        break
