#def get_proteins_ratio_by_residue_threshold(filename,
#residue,relative_threshold=0.03, absolute_threshold=10):
#    '''Given a multi-line protein FASTA file (stored in a file with path 
 #   defined filename), returns a float corresponding to the ratio of proteins 
  #  in the fasta file having a relative frequency higher or equal than a given 
   # threshold provided as an argument named “relative_threshold” and having an 
    #absolute frequency of the same residue higher or equal than a given threshold 
    #provided as an argument named “absolute_threshold” for a given residue. 
    #The function should be named as follows, with the same arguments 
    #definition.'''

   # with open(filename, 'r') as fd:
   #     for line in filename:
            

residue = "A"
relative_threshold = 0.001
absolute_threshold = 20



with open('/home/martin/master/2nd_trim/pyt/assignment2/example_fasta_file.fa', 'r') as filename:
    relative_freq = 0
    residue_counter = 0
    protein_counter = 0


    for line in filename:

        if line.startswith(">"):    # if line is header
            
            #if (relative_freq >= relative_threshold) and (residue_counter >= absolute_threshold):
                
            
            residue_counter = 0 # reset counter
            
            total_residue_counter = 0 # reset counter
        
        else:
            for char in line.strip():
                total_residue_counter += 1  # count all residues
                print(total_residue_counter)
                

                if char == residue:
                    residue_counter += 1  # count residue of interest
                    
                relative_freq = residue_counter/total_residue_counter

                

total_residue_counter


















    
            