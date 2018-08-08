from Bio.PDB import *
import os    ##for directory
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder,three_to_one

#set the directory
os.chdir('/Users/luho/PycharmProjects/pdb/code')

"""print the amino acid sequence for specific chain which has coordinated"""
def getChainInf(oneChain, pdbID):
    '''
    :param oneChain:
    :return:
    a list contains the chainID, the amino acid sequence and the coordinate
    '''
    s1 = []
    for aa in oneChain:
        s1.append(aa)

    seq1 = []
    order1 = []
    s0 = None

    for x in s1:
        t1 = str(x).split(' ')
        try:
            s0 = three_to_one(t1[1])
        except KeyError:
            pass
        seq1.append(s0)
        t2 = [x for i, x in enumerate(t1) if 'resseq' in x]
        t3 = ''.join(t2).replace('resseq=', '')
        order1.append(t3)

    seq10 = np.array(seq1)
    order10 = np.array(order1)
    index = seq10 != np.array(None)
    seq11 = seq10[index]
    order11 = order10[index]

    seq12 = ''.join(seq11)
    order12 = ','.join(order11)
    chainInf = [pdbID + '.' + oneChain.id, seq12, order12]
    return chainInf

"""save the fasta information for each chain ID"""
def saveFasta(id):
    '''
    :param id: pdb id
    :return:  save the fasta infomation for each chain of this PDB structure
    '''
    infile = '../data/pdb_ex_right_format/' + id + '.pdb'
    outfile = '../result/pdb_ex_seq/' + id
    fastaPDB = open(outfile, 'w')
    p = PDBParser()
    structure = p.get_structure(id, infile)
    chainID = []
    chainInf = []
    for model in structure :
        for chain in model:
            s = getChainInf(chain,id)
            print(s)
            s1 = '\n'.join(s[0:2])
            chainID.append(s)
            chainInf.append(s1)

    for item in chainInf:
        fastaPDB.write("%s\n" % item)


# example to display the residue sequence with coordinates
# single input
id = 'YLL015W'
infile = '../data/test_pdb/' + id + '.pdb'
outfile = '../result/' + id
fastaPDB = open(outfile, 'w')
p = PDBParser()
structure = p.get_structure(id, infile)
chainID = []
chainInf = []
for model in structure :
    for chain in model:
        s = getChainInf(chain,id)
        print(s)
        s1 = '\n'.join(s[0:2])
        chainID.append(s)
        chainInf.append(s1)

for item in chainInf:
    fastaPDB.write("%s\n" % item)
fastaPDB.close()
#change chainID into a dataframe for easy observation
#the chainSum could be used to check the calculated distance
chainSum = pd.DataFrame(chainID)
chainSum.columns = ['chainID','residue','order']

#example code to find difference for different chains
len(chainSum['residue'][0])
len(chainSum['residue'][1])
chainA = chainSum['order'][0]
chainB = chainSum['order'][1]
chainA0 = np.array(chainA.split(','))
chainB0 = np.array(chainB.split(','))
np.setdiff1d(chainB0,chainA0)



#batch process
#this process could produce the residue (amino acids) sequence for all the experimental structures
pdb_all = os.listdir('../data/pdb_ex_right_format')
pdb_all = [x.replace('.pdb','') for i,x in enumerate(pdb_all)]

for x in pdb_all:
    saveFasta(x)

#merge all files
filenames = os.listdir('../result/pdb_ex_seq')

with open('../result/pdb_ex_seq_summary', 'w') as outfile:
    for fname in filenames:
        fname0 = '../result/pdb_ex_seq/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())

#change the pdb seq file into fasta format
pdb_seq = open('../result/pdb_ex_seq_summary','r').readlines()
pdb_seq0 = open('../result/pdb_ex_seq_fasta','w')
for line in pdb_seq:
    if '.' in line:
        s1 = '>' + line
    else:
        s1 = line
    pdb_seq0.writelines(s1)
pdb_seq0.close()