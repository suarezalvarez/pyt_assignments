#! /usr/bin/env python

import sequence_dictionaries as seqdicts
import sys
import os

class IncorrectSequenceLetter(ValueError):
    def __init__(self, letter, class_type):
        
        self.letter = letter
        
        self.class_type = type(class_type).__name__

    def __str__(self):
        return "The sequence item %s is not found in the alphabet of class %s" %(self.letter , self.class_type)


class Sequence():       # define class sequence

    alphabet = ""
    weights = {}


    def __init__(self, identifier, sequence):                     # define attributes

        self.__identifier = identifier
        self.__sequence = sequence
    

        for char in self.__sequence:
            if char not in self.alphabet:

                raise IncorrectSequenceLetter(char , self) # raise exception if character is not in alphabet
                
    def get_identifier(self):
        return self.__identifier
    
    def get_sequence(self):
        
        return self.__sequence
    
    def get_mw(self):
    
        weight = 0

        for residue in self.get_sequence():         # calculate mw
            weight += self.weights[residue]
        
        return weight
    

    def has_subsequence(self):
        
        return Sequence.get_sequence() in self.sequence
    

    def __len__(self):
        return len(self.get_sequence())
    

    def __eq__(self,other):

        return self.get_sequence() == other.get_sequence()
    
    def __ne__(self,other):

        return self.get_sequence() != other.get_sequence()
    

    def __add__(self,other):

        if type(self) == type(other):

            return type(self)("+".join([self.get_identifier(),other.get_identifier()]),            # create new instance where the identifier is id1+id2, and the sequence is the concatenation of the 2 sequences
                               self.get_sequence()+other.get_sequence())
        
        else:

            raise TypeError("Sequences must be instances of the same class to be added.")

    def __getitem__(self,key):

        return self.get_sequence()[key]         # return sequence at index key

    def __contains__(self , item):

        return item in self.get_sequence()
    
    
    def __str__(self):
        
        return self.__identifier
    
    def __lt__(self,other):

        return self.get_mw() < other.get_mw()
    
    def __le__(self,other):

        return self.get_mw() <= other.get_mw()
    
    def __gt__(self,other):

        return self.get_mw() > other.get_mw()
    
    def __ge__(self,other):

        return self.get_mw() >= other.get_mw()
    

    def __hash__(self):

        self.hash_attr = (self.get_identifier() , self.get_sequence())

        return self.hash_attr.__hash__()
    






class ProteinSequence(Sequence):

    alphabet = seqdicts.protein_letters

    weights = seqdicts.protein_weights

    pass






class NucleotideSequence(Sequence):

    translate_table = {}

    def translate(self):

        translated_seq = ""
        

        for nucleotide in range(len(self.get_sequence())):

            if self.get_sequence()[nucleotide:nucleotide+3] in self.start:  # slide window of 3 nucleotides to detect start codon
                
                coding_start = nucleotide    # detect start codon
                break

            else:
                coding_start = 0



        for codon_start in range(coding_start,len(self.get_sequence()),3):
            

            codon = self.get_sequence()[codon_start:codon_start+3]     # translate sequence
            
            if len(codon) < 3: # if the sequence length is not multiple of 3, discard last incomplete codon
                break
            
            translated_seq += self.translate_table[codon]

            if codon in self.stop:
                break

                        
        translated_seq = ProteinSequence(self.get_identifier() , translated_seq)
        
        return translated_seq
    




class DNASequence(NucleotideSequence):

    translate_table = seqdicts.dna_table
    alphabet = seqdicts.dna_letters
    weights = seqdicts.dna_weights
    start = seqdicts.dna_start_codons
    stop = seqdicts.dna_stop_codons

    def transcribe(self):
        
        transcribed_seq = ""
        for char in self.get_sequence():

            if char == "T":
                transcribed_seq += "U"


            else:
                transcribed_seq += char

        transcribed_seq = RNASequence(self.get_identifier() , transcribed_seq)
        
        return transcribed_seq
        





class RNASequence(NucleotideSequence):

    translate_table = seqdicts.rna_table
    alphabet = seqdicts.rna_letters
    weights = seqdicts.rna_weights
    start = seqdicts.rna_start_codons
    stop = seqdicts.rna_stop_codons

    def reverse_transcribe(self):
        
        rt_seq = ""
        for char in self.get_sequence():

            if char == "U":
                rt_seq += "T"

            else:
                rt_seq += char

        rt_seq = DNASequence(self.get_identifier() , rt_seq)
        
        return rt_seq





### fasta iterator ###
    
def FASTA_iterator(fasta_filename , class_type = ProteinSequence):

    '''generator of instances of the class Protein()'''

    with open(fasta_filename , 'r') as f:

        seq = ""

        for line in f:

            if line.startswith('>'):                                                    # if header
                
                if seq:
                    
                    try:
                        yield class_type(raw_header, seq)
                                                                 # create instance of class Protein
                    
                    except IncorrectSequenceLetter:
                        sys.stderr.write("The sequence %s is not a valid sequence for class %s" %(raw_header, class_type.__name__))
                        sys.stderr.write("\n")
                    
                    seq = ""                                                 # reset sequence string
            
                
                raw_header = line.strip()                                               # save header of current protein
                raw_header = raw_header.replace(">" , "")


            else:                                                                       # if line is protein sequence

                seq += line.strip()                                                     # append lines of the sequence to "seq"

        yield class_type(raw_header , seq)





## find fastas
        
def find_fastas(directory):

    '''find all fasta files in a directory'''

    fasta_files = []


    for file in os.listdir(directory):

        
        if file.endswith(".fasta") or file.endswith(".fa"):

            fasta_files.append(f'{directory}/{file}')

    
        
    sys.stderr.write(f"{len(fasta_files)} FASTA files found")
    sys.stderr.write('\n')
    return fasta_files








if __name__ == "__main__":


    if len(sys.argv) == 1:      # if there's no specified input file, search for fasta files in the current directory
        
        output_list = []
        sequence_number = 1

        for fastafile in find_fastas(os.getcwd()):

            
            for sequence_instance in FASTA_iterator(fastafile, DNASequence):
                
                id = str(sequence_instance.get_identifier()).ljust(20)
                length_seq = str(len(sequence_instance)).rjust(10)


                mw_seq = sequence_instance.translate().get_mw()

                output_tuple = (id, length_seq, mw_seq)
                output_list.append(output_tuple)
                sequence_number += 1
             

            # write to stderr the progress of the file that is being processed
            sys.stderr.write(f'{fastafile} finished.  \n')

        # write to stderr the number of sequences found
        sys.stderr.write(f'{sequence_number} sequences found. \n')

        
        # sort the sequences by molecular weight
        sys.stderr.write("Sorting the sequences... \n")
        
        output_list.sort(key = lambda x: x[2] , reverse = True)
        sys.stderr.write("Sort process finished. \n")

        # write to stderr that the program finished correctly
        sys.stderr.write("Program finished correctly. \n")
        sys.stderr.write("--------------------------------- \n")

        # write to stdout the sorted list of lengths, molecular weights and identifiers

        for sequence_instance in output_list:
            
            for feature in sequence_instance:
                sys.stdout.write(str(feature))
                sys.stdout.write("\t")
            
            sys.stdout.write("\n")

   


##############################################
    
    elif len(sys.argv) == 2:           # if there's no specified output file, print it in stdout
        
        if os.path.isfile(sys.argv[1]): # if the input is a file

            sequence_number = 1
            output_list = []   

            for sequence_instance in FASTA_iterator(sys.argv[1], DNASequence):
                

                id = str(sequence_instance.get_identifier()).ljust(20)
                length_seq = str(len(sequence_instance)).rjust(10)


                mw_seq = sequence_instance.translate().get_mw()

                output_tuple = (id, length_seq, mw_seq)
                output_list.append(output_tuple)
                sequence_number += 1

            # write to stderr the progress of the file that is being processed
            sys.stderr.write(f'{os.getcwd()}/{sys.argv[1]} finished.  \n')

            # write to stderr the number of sequences found
            sys.stderr.write(f'{sequence_number} sequences found. \n')

            
            # sort the sequences by molecular weight
            sys.stderr.write("Sorting the sequences... \n")
            
            output_list.sort(key = lambda x: x[2] , reverse = True)
            sys.stderr.write("Sort process finished. \n")

            # write to stderr that the program finished correctly
            sys.stderr.write("Program finished correctly. \n")
            sys.stderr.write("--------------------------------- \n")

            # write to stdout the sorted list of lengths, molecular weights and identifiers

            for sequence_instance in output_list:
                
                for feature in sequence_instance:
                    sys.stdout.write(str(feature))
                    sys.stdout.write("\t")
                
                sys.stdout.write("\n")





        elif not os.path.isfile(sys.argv[1]):           # if the input is a directory

            output_list = []
            sequence_number = 1

            for fastafile in find_fastas(sys.argv[1]):

                    
                for sequence_instance in FASTA_iterator(fastafile, DNASequence):
                    
                    id = str(sequence_instance.get_identifier()).ljust(20)
                    length_seq = str(len(sequence_instance)).rjust(10)


                    mw_seq = sequence_instance.translate().get_mw()

                    output_tuple = (id, length_seq, mw_seq)
                    output_list.append(output_tuple)
                    sequence_number += 1
                    

                sys.stderr.write(f'{os.getcwd()}/{fastafile} finished \n')

            sys.stderr.write(f'{sequence_number} sequences found. \n')
            
            # sort the sequences by molecular weight
            sys.stderr.write("Sorting the sequences... \n")
            
            output_list.sort(key = lambda x: x[2] , reverse = True)
            sys.stderr.write("Sort process finished. \n")

            # write to stderr that the program finished correctly
            sys.stderr.write("Program finished correctly. \n")
            sys.stderr.write("--------------------------------- \n")        
            
            # write to stdout the sorted list of lengths, molecular weights and identifiers
            for sequence_instance in output_list:
                
                for feature in sequence_instance:
                    sys.stdout.write(str(feature))
                    sys.stdout.write("\t")
                
                sys.stdout.write("\n")


    
    ##############################################
        

    elif len(sys.argv) == 3:        # if there's input file and output file 
        
        if os.path.isfile(sys.argv[1]): # if the input is a file

            sequence_number = 1
            output_list = []   

            for sequence_instance in FASTA_iterator(sys.argv[1], DNASequence):
                

                id = str(sequence_instance.get_identifier()).ljust(20)
                length_seq = str(len(sequence_instance)).rjust(10)


                mw_seq = sequence_instance.translate().get_mw()

                output_tuple = (id, length_seq, mw_seq)
                output_list.append(output_tuple)
                sequence_number += 1

            # write to stderr the progress of the file that is being processed
            sys.stderr.write(f'{os.getcwd()}/{sys.argv[1]} finished.  \n')

            # write to stderr the number of sequences found
            sys.stderr.write(f'{sequence_number} sequences found. \n')

            
            # sort the sequences by molecular weight
            sys.stderr.write("Sorting the sequences... \n")
            
            output_list.sort(key = lambda x: x[2] , reverse = True)
            sys.stderr.write("Sort process finished. \n")

            # write to stderr that the program finished correctly
            sys.stderr.write("Program finished correctly. \n")
            sys.stderr.write("--------------------------------- \n")

            # write to file the sorted list of lengths, molecular weights and identifiers
            with open(sys.argv[2], 'w') as f:

                for sequence_instance in output_list:
                    
                    for feature in sequence_instance:

                        
                        f.write(str(feature))
                        f.write("\t")

                    f.write('\n')
                        
                
                





        elif not os.path.isfile(sys.argv[1]):           # if the input is a directory
            # create a list of tuples with the identifier, length and molecular weight of each sequence
            output_list = []
            sequence_number = 1

            for fastafile in find_fastas(sys.argv[1]):

                # iterate over the sequences in the fasta file
                for sequence_instance in FASTA_iterator(fastafile, DNASequence):
                    
                    id = str(sequence_instance.get_identifier()).ljust(20)
                    length_seq = str(len(sequence_instance)).rjust(10)


                    mw_seq = sequence_instance.translate().get_mw()

                    output_tuple = (id, length_seq, mw_seq)
                    output_list.append(output_tuple)
                    sequence_number += 1
                    

                sys.stderr.write(f'{os.getcwd()}/{fastafile} finished \n')

            sys.stderr.write(f'{sequence_number} sequences found. \n')
            
            # sort the sequences by molecular weight
            sys.stderr.write("Sorting the sequences... \n")
            
            output_list.sort(key = lambda x: x[2] , reverse = True)
            sys.stderr.write("Sort process finished. \n")

            # write to stderr that the program finished correctly
            sys.stderr.write("Program finished correctly. \n")
            sys.stderr.write("--------------------------------- \n")        
            

            # write to file the sorted list of lengths, molecular weights and identifiers
            with open(sys.argv[2], 'w') as f:

                for sequence_instance in output_list:
                    
                    for feature in sequence_instance:
                        f.write(str(feature))
                        f.write("\t")

                    f.write("\n")




    else:
        print("Provide at least one input file, and at most one input file and one output file.")



