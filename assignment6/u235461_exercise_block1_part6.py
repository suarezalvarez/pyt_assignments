#!/usr/bin/env python
import math
import sys
import tempfile
import fileinput

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path = sys.stdin):
    with open(pdb_file_path , "r") as structure_pdb:

        output_dict = {}        # create dictionary with all chains, residues, atoms and coordinates



        for line in structure_pdb:      
            if line.startswith('ATOM'):                         # select lines that start with "ATOM"
            
                # save features of the line
                atom_name = line[13:17].strip()               # unique atom label
                residue = line[17:21].strip()                 # residue name
                chain = line[21:22].strip()                   # chain id
                res_num = line[23:27].strip()                 # residue number in the chain
                coords = {'X': float(line[31:39].strip()),    # coordinates of each atom
                            'Y': float(line[39:47].strip()),
                            'Z': float(line[47:55].strip())}

                # if new chain
                if chain not in output_dict:                

                    # id of previous chain is a key in output_dict, with chain_dict as value
                    chain_dict = {}
                    output_dict[chain] = chain_dict    




                if residue+res_num not in chain_dict:
                    residue_dict = {}
                # residues of the chain


                residue_dict[atom_name] = coords            # store atom coordinates in residue
                chain_dict[residue+res_num] = residue_dict          # store residues in chain
                output_dict[chain] = chain_dict  


    chain_distance = {}             # empty dictionary for the distances inside a chain
    for chain in output_dict:       # iterate over chains

        residues_distance = {}                          # empty dictionary for the distances in a residue

        residues_list = list(output_dict[chain].keys())   

        for i in range(len(residues_list)):

            residue1 = residues_list[i]

            for residue2 in residues_list[i+1:]:


                atoms_distance = []                 
                for atom1 in output_dict[chain][residue1]:

                    for atom2 in output_dict[chain][residue2]:


                        coordinates1 = output_dict[chain][residue1][atom1]
                        coordinates2 = output_dict[chain][residue2][atom2]

                        atoms_distance.append(math.sqrt((coordinates1['X']-coordinates2['X'])**2+
                                                (coordinates1['Y']-coordinates2['Y'])**2 +
                                                (coordinates1['Z']-coordinates2['Z'])**2))


                residues_distance[residue1+"_"+residue2] = min(atoms_distance)


        chain_distance[chain] = residues_distance




    for chain in chain_distance:
        sum_dist = 0
        count_pairs = 0

        for pair in chain_distance[chain]:
                count_pairs += 1
                sum_dist += chain_distance[chain][pair]

        chain_mean_distance = sum_dist/count_pairs

        chain_distance[chain] = float(str("%.4f" %chain_mean_distance))


    # formatting

    for item in chain_distance.items():

        output_string = str(item)
        output_string = output_string.replace("(" , "" , -1)
        output_string = output_string.replace(")" , "" , -1)
        output_string = output_string.replace("'" , "" , -1)
        output_string = output_string.replace("," , ":" , -1)
        print(output_string , file = sys.stdout)

    return 



if __name__ == '__main__':

    if len(sys.argv) > 1:
        calculate_pdb_chain_mean_minimum_distances(sys.argv[1])

    else:
        # Assuming pdb_input is a file-like object
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            for line in fileinput.input():
                temp_file.write(line)
            
            temp_file_path = temp_file.name
        calculate_pdb_chain_mean_minimum_distances(temp_file_path)