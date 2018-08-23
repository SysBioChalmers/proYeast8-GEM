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

#function to preprocess residue from experiment PDB files, sometimes the residue sequence in the structure is different from the coordinates
#thus we should remove the redundant residue based on the given coordinates
def preprocessResidueEXP(chain0, start1, end1):
    '''
    :param chain0: aimed chain with the coordinates
    :param start1: start1 mean that in the blast analysis, start1 represent the first residue which could find in original protein
    :param end1: end1 mean that in the blast analysis, end1 represent the last residue which could find in original protein. start1 and
    end1 is not consistent with the residue order recorded in pdb file
    :return: a refined chain which could mapping onto the original protein
    '''
    chain_one = chain0
    len(chain_one)
    row0 = []
    residue_one0 = []
    for row, residue_one in enumerate(chain_one):
        row0.append(row)
        residue_one0.append(residue_one)

    #obtain the aimed target residue sequence
    target_residue = list(range(start1-1,end1))
    residue_one1 = [residue_one0[i] for i in target_residue]
    return residue_one1


#this function could help to check why keyError: 'CA' happened
def calc_residue_dist_at_two_pos(chain_check, s1, s2):
    chain_one = chain_check
    ss0 = []
    for row, residue_one in enumerate(chain_one):
        ss0.append(residue_one)
    return calc_residue_dist(ss0[s1],ss0[s2])

#pre-process the pdb meta information before calculate the residue distance
#set the directory
os.chdir('/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/PDB quality analysis/result')
#read meta data for one group of structure
pdb_sce = pd.read_csv('pdb_Ex refine for final residue distance calculation.txt', sep="\t")
pdb_sce['coordinate_id0'] = pdb_sce['template'] + '.pdb'

#split the template to obtain the chainID for each homo pdb file
pdb_sce0 = pdb_sce

#check whether the pdb_sce existed in the aimed file
pdbfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/PDB quality analysis/pdb_ex_right_format/'
exp_pdb_all = os.listdir(pdbfile)
np.setdiff1d(pdb_sce['coordinate_id0'], exp_pdb_all)


# batch process
# the followed code were used to process all the pdb structure at one time
# it will return the pdb id which need further manual check
#input all the pdb information
pdb_inf = pdb_sce0
PDB_check = []
chain_error = []
key_error= []
for i in range(0, len(pdb_inf['coordinate_id0'])):
    print(i)
    # get the paired distance
    p = PDBParser()
    # input the geneID or PDB ID
    pdbID0 = pdb_sce['template'][i]
    pdbID = pdb_inf['coordinate_id0'][i]
    chainID = pdb_inf['chain_new'][i]

    # input the relative coordinated of PDB structure
    # the start0, end0 are quite different from that used in the homology model
    start0 = pdb_inf['qstart2'][i]
    end0 = pdb_inf['qend2'][i]
    coordinate = list(range(start0, end0))
    length0 = len(coordinate) + 1

    # set directory for the input and output
    infile = pdbfile + pdbID
    outfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/' + pdbID0 + '@' + chainID + '.txt'

    structure = p.get_structure(pdbID, infile)
    model = structure[0]
    # first obtain the chainID for the model
    chainID0 = []
    for chain in model:
        chainID0.append(chain.get_id())
    if  chainID in chainID0:
        chain11 = model[chainID]
        chain_filter = preprocessResidueEXP(chain0=chain11, start1=start0, end1=end0)
        try:
            ss = calc_dist_matrix(chain_filter, chain_filter)
        except KeyError:
            key_error.append(i)
            pass
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

len(key_error)
key_error = pd.Series(key_error)
writer = pd.ExcelWriter('/Users/luho/PycharmProjects/pdb/result/error in residue sequence for experimental pdb file.xlsx')
key_error.to_excel(writer,'Sheet1')
writer.save()



