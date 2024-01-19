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
                
            else:                                                                                           # run only if the thresholds have not been reached
                for char in line.strip():                                                                   # iterate over the characters of the lines that belong to the sequence, not the header
                    residue_counter += 1                                                                    # count residues
                    if char == residue:
                        target_residue_counter += 1                                                         # count target residues
                
    
        protein_counter += 1                                                                        # take into account last protein
        relative_freq = target_residue_counter/residue_counter                                      # calculate relative frequency of the target residue
        if relative_freq >= relative_threshold and target_residue_counter >= absolute_threshold:    # if thresholds are reached, count protein and set target achieved as True (stop iterating over sequence)
                
            target_protein_counter += 1        

    return(target_protein_counter/protein_counter)


'''Given a protein FASTA file (filename), save on a output file named output_filename the
protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency
in the protein of all the aminoacids found in the protein (the aminoacids that do not appear
in the protein should not be shown). The fields must be separated by a tabulator, and one
protein by line.'''

N = 4
M = 5

from collections import Counter


with open('example_fasta_file.fa' , 'r') as read_file:
    with open('tab_fasta_example.tab' , 'w') as written_tab:
        

        seq = ""
        raw_header = ""

        for line in read_file:
            if line.startswith('>'):
                if seq:
                    start = seq[:N] + "\t"
                    end = seq[-M:] + "\t"

                    count = str(Counter(seq))
                    count = count[count.index("{")+1:count.index("}")] + "\n"

                    written_tab.write(raw_header+start+end+count)

                    
                    
                
                
                raw_header = line.strip()
                raw_header = raw_header.replace(">" , "") + "\t"
                print("hello")
                print(raw_header)

                
            
            else:

                seq += line.strip()
                



seq = "AAAAAABBBBB"
Counter(seq)