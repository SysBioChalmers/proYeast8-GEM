from Bio.PDB import *
import os    ##for directory
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

#function
def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    ''' an example
    residue_one = residue1
    residue_two = residue2'''
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord ##some errors: KeyError: 'CA'
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""

    ''' an example
    chain_one = chain
    chain_two = chain'''
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

#set the directory
os.chdir('/Users/luho/PycharmProjects/pdb/code')
print(os.getcwd()) #obtain the present directory


#get the paired distance
p = PDBParser()

#input the geneID or PDB ID
pdbID = 'YLL015W'
chainID = 'A'

#input the relative coordinated of PDB structure
start0 = 279
end0 = 1556
coordinate = list(range(start0,end0))
length0 = len(coordinate) + 1

# set directory for the input and output
infile = '../data/test_pdb/' + pdbID + '.pdb'
outfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/ResidueDistance_' + pdbID + '.txt'

structure = p.get_structure(pdbID, infile)
model = structure[0]
chain = model[chainID]  # should input the chain information for each structure
ss = calc_dist_matrix(chain, chain)

# how to add quality control for the distance before save the file for the downstream analysis
# First step we can compare the dimension with the relative length
# the dimension is 390 for chain A, while the relative coordinates is  25-417
# thus it can be found that the actual amino acid length (417-25)=393 is not equal to the dimension of matrix, thus need manual check
# In the manual check part,
# we can parse the pdb structure and obtain the residue information
# then we can compare the dimension of the distance matrix and the length of residue which has coordinates
# Finally for a strict comparison, the residue amino acids sequence should be compared the sequence from structure

dimension1 = list(ss.shape)
if dimension1[0] == length0:
    np.savetxt(outfile, ss, delimiter=',')
else:
    raise ValueError("wrong dimension of distance matrix")



# batch process
# the followed code were used to process all the pdb structure at one time
# it will return the pdb id which need further manual check
#input all the pdb information
infile0 = '../data/test_pdb/' + 'pdb_inf_distance.xlsx'
pdb_inf = pd.read_excel(infile0)
PDB_check = []
for i in range(0, len(pdb_inf['PDB_name'])):
    print(i)
    # get the paired distance
    p = PDBParser()

    # input the geneID or PDB ID
    pdbID = pdb_inf['PDB_name'][i]
    chainID = pdb_inf['chainID'][i]

    # input the relative coordinated of PDB structure
    start0 = pdb_inf['start'][i]
    end0 = pdb_inf['end'][i]
    coordinate = list(range(start0, end0))
    length0 = len(coordinate) + 1

    # set directory for the input and output
    infile = '../data/test_pdb/' + pdbID + '.pdb'
    outfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/ResidueDistance_' + pdbID + '.txt'

    structure = p.get_structure(pdbID, infile)
    model = structure[0]
    chain = model[chainID]  # should input the chain information for each structure
    ss = calc_dist_matrix(chain, chain)

    # how to add quality control for the distance before save the file for the downstream analysis
    # First step we can compare the dimension with the relative length
    # the dimension is 390 for chain A, while the relative coordinates is  25-417
    # thus it can be found that the actual amino acid length (417-25)=393 is not equal to the dimension of matrix, thus need manual check
    # In the manual check part,
    # we can parse the pdb structure and obtain the residue information
    # then we can compare the dimension of the distance matrix and the length of residue which has coordinates
    # Finally for a strict comparison, the residue amino acids sequence should be compared the sequence from structure

    dimension1 = list(ss.shape)
    if dimension1[0] == length0:
        np.savetxt(outfile, ss, delimiter=',')
    else:
        PDB_check.append(pdbID)
        continue