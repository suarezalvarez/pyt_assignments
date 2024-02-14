import sequence_dictionaries as seqdicts



class Sequence():       # define class sequence

    alphabet = ""
    weights = {}


    def __init__(self, identifier, sequence):                     # define attributes

        self.__identifier = identifier
        self.__sequence = sequence
    

        for char in self.__sequence:
            if char not in self.alphabet:

                raise ValueError(f"Impossible to create instance: {char} not possible") # raise exception if character is not in alphabet
                
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
        
        for codon_start in range(0,len(self.get_sequence()),3):

            codon = self.get_sequence()[codon_start:codon_start+3]     # translate sequence
            translated_seq += self.translate_table[codon]

        translated_seq = ProteinSequence(self.get_identifier() , translated_seq)
        
        return translated_seq
    




class DNASequence(NucleotideSequence):

    translate_table = seqdicts.dna_table
    alphabet = seqdicts.dna_letters
    weights = seqdicts.dna_weights

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

    def reverse_transcribe(self):
        
        rt_seq = ""
        for char in self.get_sequence():

            if char == "U":
                rt_seq += "T"

            else:
                rt_seq += char

        rt_seq = DNASequence(self.get_identifier() , rt_seq)
        
        return rt_seq
