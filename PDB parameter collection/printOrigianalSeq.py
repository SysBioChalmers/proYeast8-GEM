"""print the chain and amino acid sequence for each pdb structure
these sequence is just recorded in the pdb web database
https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=3J9T&compressionType=uncompressed
actually the residues contained in the pdb files with coordinates can be different from the above residues sequence
"""

import os
from Bio import SeqIO #input biopython
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