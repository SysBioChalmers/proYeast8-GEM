from Bio.PDB import *
import os    ##for directory
import numpy
import pandas as pd
#set the directory
os.chdir('/Users/luho/PycharmProjects/python learning/venv/project2_pdb/code')
print (os.getcwd()) #obtain the present directory

#function
def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    ''' an example
    residue_one = residue1
    residue_two = residue2'''
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord ##some errors: KeyError: 'CA'
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""

    ''' an example
    chain_one = chain
    chain_two = chain'''
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer


#get the paired distance
p = PDBParser()
gene = 'YBR115C'
infile = '../data/' + gene + '.pdb'
outfile = '/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/ResidueDistance_' + gene + '.xlsx'

structure = p.get_structure(gene, infile)
model = structure[0]
chain = model['A']
ss = calc_dist_matrix(chain,chain)
ResidueDistance = pd.DataFrame(ss)


#get sequence
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
ppb=PPBuilder()
seq = []
for pp in ppb.build_peptides(structure):
    seq.append(pp.get_sequence())
    print(pp.get_sequence())

ResidueDistance.columns = list(seq[0]) #sometimes the seq number is not equal to the distance number
ResidueDistance['row'] = list(seq[0])
ResidueDistance0 = ResidueDistance.set_index('row')

writer = pd.ExcelWriter(outfile)
ResidueDistance0.to_excel(writer,'Sheet1')
writer.save()











'''code to check why there exist error for protein structure YIL160C'''
lines = open('../data/YIL160C.pdb', 'r').readlines()
index1 = []
aa = []
for line in lines:
    column = line.split(' ')
    column0 = [x for i, x in enumerate(column) if x !='']
    print(column0)
    if 'CA' in column0[2]:
        index1.append(column0[5])
        aa.append(column0[3])

import numpy as np
[i for i, x in enumerate(index1) if x =='417']
s1 = index1[0:390]
s2 = [str(x) for x in range(25,418)]
x = np.array(s1)
y = np.array(s2)
np.setdiff1d(y,x)


s3 = index1[390:]
z = np.array(s3)
np.setdiff1d(z,y)

ResidueDistance.columns = list(seq[2]) #sometimes the seq number is not equal to the distance number
ResidueDistance['row'] = list(seq[2])





'''example to calculate the distance'''
#residues = structure.get_residues()
# Iterate over all atoms in a structure
#for atom in structure.get_atoms():
#    print(atom)

# Iterate over all residues in a model
#for residue in model.get_residues():
#    print(residue)


#residue = chain[100]
#atom = residue['CA']
#calculate the distance
#residue1 = chain[99]
#residue2 = chain[101]
#distance = residue1['CA'] - residue2['CA']

#use the function to get the distance between any two residues