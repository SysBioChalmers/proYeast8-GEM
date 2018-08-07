"""print the chain and amino acid sequence for each pdb structure
these sequence is just recorded in the pdb files
"""

import os
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser

#set the directory
os.chdir('/Users/luho/PycharmProjects/pdb/code')

#get the paired distance
p = PDBParser()
gene = 'YOR270C'
infile = '../data/' + gene + '.pdb'

handle = open(infile, "rU") # a pdb structure  as the input
chainInf ={}
chainLen =[]
for record in SeqIO.parse(handle, "pdb-seqres"):
    chainInf[record.id] = record.seq
    chainLen.append(len(record.seq))
    print(">" + record.id + "\n" + record.seq)