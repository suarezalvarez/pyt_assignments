##########################
## Define class Protein ##
##########################


class Protein():

    def __init__(self , identifier , sequence):

        self.identifier = identifier
        self.sequence = sequence

    def get_identifier(self):

        return self.identifier

    def get_sequence(self):
        
        return self.sequence
    
    
    def get_mw(self):

        amino_acid_weights = {
    'A': 89.09,
    'R': 174.2,
    'N': 132.12,
    'D': 133.1,
    'C': 121.16,
    'Q': 146.14,
    'E': 147.13,
    'G': 75.07,
    'H': 155.15,
    'I': 131.17,
    'L': 131.17,
    'K': 146.19,
    'M': 149.21,
    'F': 165.19,
    'P': 115.13,
    'S': 105.09,
    'T': 119.12,
    'W': 204.23,
    'Y': 181.19,
    'V': 117.15
    }
        
        protein_weight = 0
        for residue in self.sequence:
            protein_weight += amino_acid_weights[residue]
        
        return protein_weight
    

    def has_subsequence(self, Protein):

        return Protein.get_sequence() in self.sequence
    
    
    def get_length(self):

        return len(self.sequence)
    


###########################
## Modify FASTA iterator ##
###########################
    
def FASTA_iterator(fasta_filename):

    '''generator of instances of the class Protein()'''

    with open(fasta_filename , 'r') as f:

        seq = ""

        for line in f:

            if line.startswith('>'):                                                    # if header
                
                if seq:
                    
                    raw_header = Protein(identifier = raw_header,       # create instance called with the 
                                         sequence=seq)

                    yield(raw_header)
                    
                    seq = ""                                                 # reset sequence string
            
                
                raw_header = line.strip()                                               # save header of current protein
                raw_header = raw_header.replace(">" , "")


            else:                                                                       # if line is protein sequence

                seq += line.strip()                                                     # append lines of the sequence to "seq"

        raw_header = Protein(identifier = raw_header , sequence = seq)
        yield(raw_header)
