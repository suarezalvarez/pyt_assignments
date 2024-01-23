#### get_proteins_ratio_by_residue_threshold

def get_proteins_ratio_by_residue_threshold(filename,
residue,relative_threshold=0.03, absolute_threshold=10):
    '''Given a multi-line protein FASTA file (stored in a file with path 
    defined filename), returns a float corresponding to the ratio of proteins 
    in the fasta file having a relative frequency higher or equal than a given 
    threshold provided as an argument named “relative_threshold” and having an 
    absolute frequency of the same residue higher or equal than a given threshold 
    provided as an argument named “absolute_threshold” for a given residue. 
    The function should be named as follows, with the same arguments 
    definition.'''

    with open(filename,  'r') as filename:
        relative_freq = 0                  # initialize relative frequency 
        protein_counter = 0                # initialize protein counter
        target_protein_counter = 0         # initialize target protein counter
        protein_counter = 0                # new protein
        residue_counter = 0                # initialize residue counter
        target_residue_counter = 0         # initialize target residue counter
        
        for line in filename:
            if line.startswith(">") and residue_counter != 0:
                
                protein_counter += 1       # new protein
                relative_freq = target_residue_counter/residue_counter                                      # calculate relative frequency of the target residue

                if relative_freq >= relative_threshold and target_residue_counter >= absolute_threshold:    # if thresholds are reached, count protein and set target achieved as True (stop iterating over sequence)

                    target_protein_counter += 1        

                
                residue_counter = 0        # reset residue counter
                target_residue_counter = 0 # reset target residue counter
                
            elif not line.startswith(">"):                                                                                           # run only if the thresholds have not been reached
                
                for char in line.strip():                                                                   # iterate over the characters of the lines that belong to the sequence, not the header
                    residue_counter += 1                                                                    # count residues
                
                    if char == residue:
                        target_residue_counter += 1                                                         # count target residues
                
    
        protein_counter += 1                                                                        # take into account last protein
        relative_freq = target_residue_counter/residue_counter                                      # calculate relative frequency of the target residue
        if relative_freq >= relative_threshold and target_residue_counter >= absolute_threshold:    # if thresholds are reached, count protein and set target achieved as True (stop iterating over sequence)
                
            target_protein_counter += 1        

    return(target_protein_counter/protein_counter)



#### print_sequence_summary


def print_sequence_summary(filename , output_filename , first_n = 10 , last_m = 10):
    '''Given a protein FASTA file (filename), save on a output file named output_filename the protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency in the protein of all the aminoacids found in the protein (the aminoacids that do not appear in the protein should not be shown). The fields must be separated by a tabulator, and one protein by line.'''
    
    with open(filename , 'r') as read_file:
        with open(output_filename , 'w') as written_tab:
            

            seq = ""                 # initialice sequence
            raw_header = ""          # initialice header

            for line in read_file:                                                          # iterate over lines
                
                if line.startswith('>'):                                                    # if header
                    if seq:
                        start = seq[:first_n] + "\t"                                        # take first residues of previous protein
                        end = seq[-last_m:] + "\t"                                          # take last residues of previous protein
                        
                        count_dict = {}
                        for aa in seq:                                                      # count residues of previous protein
                            if aa in count_dict:
                                count_dict[aa] += 1
                            
                            else:
                                count_dict[aa] = 1
                        
                        
                        count = str(count_dict)
                        count = count[count.index("{")+1:count.index("}")] + "\n"           # remove curly brackets and convert to string
                        count = count.replace("\'" ,"")
                        written_tab.write(raw_header+start+end+count)                       # write to output file
                        
                        seq = ""                                                            # reset sequence string
                
                    
                    raw_header = line.strip()                                               # save header of current protein
                    raw_header = raw_header.replace(">" , "") + "\t"                        # remove ">" symbol
                
                else:                                                                       # if line is protein sequence

                    seq += line.strip()                                                     # append lines of the sequence to "seq"

            start = seq[:first_n] + "\t"                                        # take first residues of previous protein
            end = seq[-last_m:] + "\t"                                          # take last residues of previous protein
                        
            count_dict = {}
            for aa in seq:                                                      # count residues of last protein
                if aa in count_dict:
                    count_dict[aa] += 1
                            
                else:
                    count_dict[aa] = 1
                        
            count = str(count_dict)        
            count = count[count.index("{")+1:count.index("}")] + "\n"           # remove curly brackets and convert to string
            count = count.replace("\'" ,"")
            written_tab.write(raw_header+start+end+count)                       # write to output file

