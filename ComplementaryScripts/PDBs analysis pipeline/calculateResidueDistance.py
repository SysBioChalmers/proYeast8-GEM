'''Function to calculate the distance between two residues in a protein structure
Hongzhong Lu 2018-07-17'''

#package

from Bio.PDB import *
import os    ##for directory
import numpy

# function
def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    ''' an example
    residue_one = residue1
    residue_two = residue2'''
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
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


# main
# set the directory
os.chdir('/Users/luho/PycharmProjects/python learning/venv/project2_pdb/code')
print (os.getcwd()) #obtain the present directory

# input an example
p = PDBParser()
structure = p.get_structure('1n8p', '../data/1n8p.pdb')

for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom)


# Iterate over all atoms in a structure
for atom in structure.get_atoms():
    print(atom)

# Iterate over all residues in a model
for residue in model.get_residues():
    print(residue)

model = structure[0]
chain = model['A']
residue = chain[100]
atom = residue['CA']

#calculate the distance
residue1 = chain[99]
residue2 = chain[101]
distance = residue1['CA'] - residue2['CA']

#use the function to get the distance between any two residues
ss = calc_dist_matrix(chain,chain)
print(ss)




