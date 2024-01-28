# FASTA_iterator

def FASTA_iterator(fasta_filename):
    '''Generator function that reads a fasta file. In each iteration, the function
    returns a tuple with the following format: (identifier , sequence).'''
    with open(fasta_filename , 'r') as f:

        seq = ""

        for line in f:

            if line.startswith('>'):                           # if header
                
                if seq:
                    
                    yield (raw_header , seq)
                    seq = ""                                   # reset sequence string
            
                
                raw_header = line.strip()                      # save header of current protein
                raw_header = raw_header.replace(">" , "")


            else:                                              # if line is protein sequence

                seq += line.strip()                            # append lines of the sequence to "seq"

        
        yield(raw_header , seq)                                # include last sequence




#  compare_fasta_file_identifiers
        
def compare_fasta_file_identifiers(fasta_filenames_list):

    '''Given a list of FASTA files, create a function that returns a dictionary that contains
    the 4 following keys with the associated values:
    - "intersection": a set with the common indentifiers found in all the files
    - "union": a set with all the identifiers (unique) found in all the files.
    - "frequency": a dictionary with all the identifiers as keys and the number of files in which it appears as values. 
    - "specific": a dictionary with the name of the input files as keys and a set with the specific identifiers as values.'''


    # create necessary dictionaries

    file_ids = {}
    output_dict = {}
    frequency_dict = {}
    specific_dict = {}

    for file in fasta_filenames_list:            # iterate over files

        ids_set = set([])

        for protein in FASTA_iterator(file):    # iterate over proteins in the file
            
            # create set with all the protein identifiers
            ids_set.add(protein[0].upper())
            
            # create "frequency" dictionary: count the number of files in which ID is present

            if protein[0].upper() not in frequency_dict:
                frequency_dict[protein[0].upper()] = 1

            else:
                frequency_dict[protein[0].upper()] += 1


        # add IDs to dictionary with format: {"file_name": IDs in filename}
                
        file_ids[file] = ids_set                # add file to file_ids as key, and ids_set as value
        specific_dict[file] = set([])
    
    
    list_of_sets = list(file_ids.values())      # create list of sets from the dictionary
    
    

    # detect file-specific IDs

    file_counter = 0 # initialize file counter
    for file_name in specific_dict:
        
        list_rest_of_sets = list_of_sets[:file_counter] + list_of_sets[file_counter+1:]  # create list of all sets except one
        union_of_rest_of_sets = list_rest_of_sets[0].union(*list_rest_of_sets[1:])       # compute union of all sets  from the list


        specific_dict[file_name] = file_ids[file_name].difference(union_of_rest_of_sets) # compute difference of left out set from the rest
        file_counter +=1 
        

    # fill output_dict

    output_dict["intersection"] = list_of_sets[0].intersection(*list_of_sets[1:])  # save intersection of all sets
    output_dict["union"] = list_of_sets[0].union(*list_of_sets[1:])                # save union of sets
    output_dict["frequency"] = frequency_dict
    output_dict["specific"] = specific_dict
    
    
    return(output_dict)
    




