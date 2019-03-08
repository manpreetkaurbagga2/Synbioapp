# SequenceAnalyser

Sequence Analyser is a program to find various DNA sequences based upon our requirements using an option from a pre-defined Menu which we get on running the code on a seq file (examplefile.seq has been used to run the following code). 

**Firstly**, it finds the **ORF** from the fasta file given as an input and it gives the _length_ of the ORF as well if asked.

**Secondly**, it finds the _position_ of the **RBS** in the sequence given.

**Thirdly**, It gives all possible **Forward and Reverse Primers** (between length (18-25) based upon 3 parameters: Melting Temperature, GC content and hairpin loop formation. 
It requires two inputs from the user, a fasta file and the length of the sequence from which he wants to build the Primers.

**Lastly**, it finds the Biobrick compatible **Restriction sites**, if present in the file.

All the outputs are published in a common text file (if you want to find everything from one fasta file step by step) as you give the inputs.

**How to use?**

1.Run the program.

2.Follow the Menu and choose an option and ans to sub-options if asked.

3.An output.txt file is created for the output you have asked for.

**Unit Testing**

It will run 7 Tests. Run the code. Press Enter when the menu is displayed. 
