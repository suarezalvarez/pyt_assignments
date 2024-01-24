def calculate_aminoacid_frequencies(fasta_filename, 
                                    subsequences_filename, 
                                    number_of_repetitions, 
                                    output_filename):
    
    '''Given a multi-line FASTA file (fasta_filename) and a "sub-sequences" file (subsequences_filename)
    (one sequence in each line), calculates the proportion of proteins in the FASTA file containing at least N times (number_of repetitions)
    each of the sub-sequences. Output is a tabular file sorted by proportion in descending order.'''



    ### create dictionary of subsequences ###

    subsequences = {} # empty dictionary to be filled with subsequences
    subsequence_counter = 0

    with open(subsequences_filename , 'r') as subsequences_filename:
        
        
        for subsequence in subsequences_filename:
            
            subsequence_counter += 1
            subsequences[subsequence.strip()] = 0    # add subsequences as keys





    ### count subsequences in each sequence ###

    with open(fasta_filename , 'r') as fasta_filename:
              
        protein_counter = 0             # initialize counter
        seq = ""                        # define seq



        for line in fasta_filename:                                                 # iterate over lines in fasta


            if line.startswith(">") and seq:                                        # if line is header
                protein_counter += 1

                for subsequence in subsequences: 

                    if seq.count(subsequence.strip()) >= number_of_repetitions:                    
                        subsequences[subsequence.strip()] += 1

                seq = ""




            elif not line.startswith(">"):

                seq += line.strip()                     # join lines into a single sequence string



        # take into account last protein
                
        protein_counter += 1

        for subsequence in subsequences:                
            if seq.count(subsequence.strip()) >= number_of_repetitions:                    
                subsequences[subsequence.strip()] += 1



    ### Formatting output ###
                
    sorted_subsequences = sorted(subsequences.items() ,                 # sort subsequences by proportion in descending order
                                key = lambda item: item[1] , 
                                reverse=True)


    with open(output_filename , 'w') as f:

        print("#Number of proteins:" + str(protein_counter).rjust(40-len("#Number of proteins:")),                      # header "number of proteins"
            file = f)
        print("#Number of subsequences" + str(subsequence_counter).rjust(40 - len("#Number of subsequences")),          # header "number of subsequences"
            file = f)
        print("#Subsequence proportions:",                                                                              # header "subsequence proportions"
            file = f)

        for index in range(len(sorted_subsequences)):                                                                   # print subsequence, absolute and relative proportions, right aligned in position 20 and 40
            print(sorted_subsequences[index][0] + 
                str(sorted_subsequences[index][1]).rjust(20-len(sorted_subsequences[index][0]))+
                str("%.4f" %(sorted_subsequences[index][1]/protein_counter)).rjust(20),
                file = f)
            
    
    return
