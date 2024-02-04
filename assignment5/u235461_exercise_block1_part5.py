def FASTA_iterator(fasta_filename):
    '''Builds a generator that allows the user to iterate over the proteins of a FASTA field, returning a tuple with the protein ID as first element, and the sequence as second element'''
    with open(fasta_filename , 'r') as f:

        seq = ""

        for line in f:

            if line.startswith('>'):                                                    # if header
                
                if seq:
                    
                    yield (raw_header , seq)
                    seq = ""                                                 # reset sequence string
            
                
                raw_header = line.strip()                                               # save header of current protein
                raw_header = raw_header.replace(">" , "")


            else:                                                                       # if line is protein sequence

                seq += line.strip()                                                     # append lines of the sequence to "seq"

        yield(raw_header , seq)







#### Proteins ratio by residue ####
        
def get_proteins_ratio_by_residue_threshold(filename,
                                            residue,
                                            relative_threshold = 0.03,
                                            absolute_threshold=10):
    '''Given a multi-line protein FASTA file (stored in a file with path defined filename), returns a
float corresponding to the ratio of proteins in the fasta file having a relative frequency higher
or equal than a given threshold provided as an argument named “relative_threshold” and
having an absolute frequency of the same residue higher or equal than a given threshold
provided as an argument named “absolute_threshold” for a given residue. The function
should be named as follows, with the same arguments definition:'''
    
    

    protein_counter = 0
    target_protein_counter = 0

    for protein in FASTA_iterator(filename):

        protein_counter += 1  # count new protein
        target_residue_counter = 0 # reset target residue counter
        residue_counter = 0 # reset residue counter

        for amino_acid in protein[1]:

            residue_counter += 1       # count residue

            if amino_acid == residue:
                target_residue_counter += 1    # count target residue
        
        if target_residue_counter >= absolute_threshold and target_residue_counter/residue_counter >= relative_threshold:
            target_protein_counter += 1

    
    return(target_protein_counter/protein_counter)




#### print sequence summary ####

def print_sequence_summary(filename ,
                            output_filename,
                            first_n = 10,
                            last_m = 10):
    
    '''Given a protein FASTA file (filename), save on a output file named output_filename the
protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency
in the protein of all the aminoacids found in the protein (the aminoacids that do not appear
in the protein should not be shown). The fields must be separated by a tabulator, and one
protein by line.'''
    
    with open(output_filename , "w") as wf:
     
        for protein in FASTA_iterator(filename):  # iterate over proteins in fasta
            
            frequency_dict = {} # create residue frequency dictionary

            for amino_acid in protein[1]:  # iterate over residues of the sequence

                if amino_acid not in frequency_dict:

                    frequency_dict[amino_acid] = 1 

                else:
                    frequency_dict[amino_acid] += 1
            
            frequency_dict_str = str(frequency_dict)
            frequency_dict_str = frequency_dict_str[frequency_dict_str.index("{") + 1:frequency_dict_str.index("}")] + "\n"
            wf.write(protein[0] + "\t" + protein[1][:first_n] + "\t" + protein[1][-last_m:] + "\t" + frequency_dict_str)
    
    return





#### get max sequence length ####

def get_max_sequence_length_from_FASTA_file(fasta_filename):
    '''A function that, given a multiline FASTA file, returns the length of the sequence with the maximum length'''

    length_dict = {}

    for protein in FASTA_iterator(fasta_filename):

        length_dict[protein[0]] = len(protein[1])

    return(max(length_dict.values()))


#### get min sequence length ####

def get_min_sequence_length_from_FASTA_file(fasta_filename):
    '''A function that, given a multiline FASTA file, returns the length of the sequence with the minimum length'''
    length_dict = {}

    for protein in FASTA_iterator(fasta_filename):

        length_dict[protein[0]] = len(protein[1])
    
    
    return(min(length_dict.values()))



#### get longest sequences ####

def get_longest_sequences_from_FASTA_file(fasta_filename):

    '''A function that, given a FASTA file, returns a list of tuples (identifier, sequence) corresponding to the sequence(s) with maximum length. The list must be sorted by the identifier (case insensitive sorted).'''
    max_length = get_max_sequence_length_from_FASTA_file(fasta_filename)
    
    my_list = [protein for protein in FASTA_iterator(fasta_filename) if len(protein[1]) == max_length]
    my_list.sort(key = lambda protein: protein[0].lower())
    
    return(my_list) 






#### get shortest sequences ####

def get_shortest_sequences_from_FASTA_file(fasta_filename):

    '''A function that, given a FASTA file, returns a list of tuples (identifier, sequence) corresponding to the sequence(s) with minimum length. The list must be sorted by the identifier (case insensitive sorted).'''
    min_length = get_min_sequence_length_from_FASTA_file(fasta_filename)
    
    my_list = [protein for protein in FASTA_iterator(fasta_filename) if len(protein[1]) == min_length]
    my_list.sort(key = lambda protein: protein[0].lower())
    
    return(my_list) 




#### get molecular weights as dictionary ####

def get_molecular_weights(fasta_filename):

    '''A function that, given a protein FASTA file, returns a dictionary with the molecular weights of all the proteins in the file. The dictionary keys must be the protein identifiers and the associated values must be a float corresponding to the molecular weight.'''
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
    
    
    protein_dict = {}

    for protein in FASTA_iterator(fasta_filename):

        protein_dict[protein[0]] = 0 
        
        for residue in protein[1]:
            if amino_acid_weights.get(residue):

                protein_dict[protein[0]] += amino_acid_weights.get(residue)

            else:
                print(f"The letter {residue} is not considered an amino acid by this program")
        

    return(protein_dict) 





#### get sequence with lowest mw ####

def get_sequence_with_min_molecular_weight(fasta_filename):
    '''A function that, given a protein FASTA file, returns a tuple with (identifier, sequence) of the protein with the lowest molecular weight. If there are two or more proteins having the minimum molecular weight, just return the first one.'''
    
    weights_dict = get_molecular_weights(fasta_filename)

    min_weight = min(weights_dict.values())  # get minimal weight

    for item in weights_dict.items():

        if abs(item[1] - min_weight) < 0.00001:       # if weight of protein is = minimal
            
            for protein in FASTA_iterator(fasta_filename):  # look into the id,sequence tuples

                if protein[0] == item[0]:  # find the protein with the id of minimal weight
                    return_tuple = protein # return it
                    break                  # break loop if protein is found, to return only the first one


    return(return_tuple)




#### mean mol weight ####

def get_mean_molecular_weight(fasta_filename):

    '''A function that, given a protein FASTA file, returns the mean of the molecular weights of all proteins'''

    weights_dict = get_molecular_weights(fasta_filename)
    sum_of_weights = 0
    for weight in weights_dict.values():

        sum_of_weights += weight
    
    mean_weight = sum_of_weights/len(weights_dict)

    return(mean_weight)