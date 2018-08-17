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



#pre-process the pdb meta information before calculate the residue distance
#set the directory
os.chdir('/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/PDB quality analysis/result')
#read meta data for one group of structure
pdb_sce = pd.read_csv('pdb_homo for PDB structure without experiment pdb.txt', sep="\t")
pdb_sce['coordinate_id0'] = pdb_sce['from'].apply(str) + '_' + pdb_sce['to'].apply(str)+\
                                        '_' + pdb_sce['template'] + \
                                        '_' + pdb_sce['coordinate_id'] + '.pdb'

#split the template to obtain the chainID for each homo pdb file
pdb_sce0 = pdb_sce.join(pdb_sce['template'].str.split('.', expand=True).rename(columns={0:'A', 1:'B', 2:'chainID'}))

#check whether the pdb_sce existed in the aimed file
pdbfile = '/Users/luho/Documents/pdb file/Swiss model information_2018_7_20/pdb_homo/'
homo_pdb_all = os.listdir(pdbfile)
np.setdiff1d(pdb_sce['coordinate_id0'], homo_pdb_all)





#one example to calculate the residue distance
#get the paired distance
p = PDBParser()
index = 0
#input the geneID or PDB ID
pdbID = pdb_sce['coordinate_id0'][index]
chainID = pdb_sce0 ['chainID'][index]

#input the relative coordinated of PDB structure
start0 = pdb_sce['from'][index]
end0 = pdb_sce['to'][index]
coordinate = list(range(start0,end0))
length0 = len(coordinate) + 1

# set directory for the input and output
infile = pdbfile + pdbID
outfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/' + pdbID + '.txt'

structure = p.get_structure(pdbID, infile)
model = structure[0]

#first obtain the chainID for the model
chainID0 = []
for chain in model:
    chainID0.append(chain.get_id())

if chainID in chainID0:
    chain = model[chainID]
    ss = calc_dist_matrix(chain, chain)
    dimension1 = list(ss.shape)

# how to add quality control for the distance before save the file for the downstream analysis
# First step we can compare the dimension with the relative length
# the dimension is 390 for chain A, while the relative coordinates is  25-417
# thus it can be found that the actual amino acid length (417-25)=393 is not equal to the dimension of matrix, thus need manual check
# In the manual check part,
# we can parse the pdb structure and obtain the residue information
# then we can compare the dimension of the distance matrix and the length of residue which has coordinates
# Finally for a strict comparison, the residue amino acids sequence should be compared the sequence from structure
    if dimension1[0] == length0:
        np.savetxt(outfile, ss, delimiter=',')
    else:
        raise ValueError("wrong dimension of distance matrix")
else:
    print ("Oops!  ChainID is not right")
    pass
len(ss)




# batch process
# the followed code were used to process all the pdb structure at one time
# it will return the pdb id which need further manual check
#input all the pdb information
pdb_inf = pdb_sce0
PDB_check = []
chain_error = []
for i in range(0, len(pdb_inf['coordinate_id0'])):
    print(i)
    # get the paired distance
    p = PDBParser()

    # input the geneID or PDB ID
    pdbID = pdb_inf['coordinate_id0'][i]
    chainID = pdb_inf['chainID'][i]

    # input the relative coordinated of PDB structure
    start0 = pdb_inf['from'][i]
    end0 = pdb_inf['to'][i]
    coordinate = list(range(start0, end0))
    length0 = len(coordinate) + 1

    # set directory for the input and output
    infile = pdbfile + pdbID
    outfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/' + pdbID + '.txt'

    structure = p.get_structure(pdbID, infile)
    model = structure[0]
    # first obtain the chainID for the model
    chainID0 = []
    for chain in model:
        chainID0.append(chain.get_id())

    if  chainID in chainID0:
        chain = model[chainID]
        ss = calc_dist_matrix(chain, chain)
        dimension1 = list(ss.shape)

    # how to add quality control for the distance before save the file for the downstream analysis
    # First step we can compare the dimension with the relative length
    # the dimension is 390 for chain A, while the relative coordinates is  25-417
    # thus it can be found that the actual amino acid length (417-25)=393 is not equal to the dimension of matrix, thus need manual check
    # In the manual check part,
    # we can parse the pdb structure and obtain the residue information
    # then we can compare the dimension of the distance matrix and the length of residue which has coordinates
    # Finally for a strict comparison, the residue amino acids sequence should be compared the sequence from structure
        if dimension1[0] == length0:
            np.savetxt(outfile, ss, delimiter=',')
        else:
            PDB_check.append(pdbID)
            continue
    else:
        print("Oops!  ChainID is not right")
        chain_error.append(pdbID)
        pass

'''chain error analysis
In the above analysis, we could find some chainID is wrong and distance can't be calculated,
Thus we need find why chainID error occur
After we find the right chainID, we can re-run the above code to obtain all the residue distance'''

len(chain_error)
chain_error = pd.Series(chain_error)
writer = pd.ExcelWriter('/Users/luho/PycharmProjects/pdb/result/error in chainID for manual check.xlsx')
chain_error.to_excel(writer,'Sheet1')
writer.save()

