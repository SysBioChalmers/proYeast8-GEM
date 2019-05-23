""""This script is mainly used to merge the hotspot analysis results
it contains the cluster for a lot of genes
9-4-2018, Hongzhong Lu"""

import os    ##for directory

#set the directory
os.chdir('/Users/luho/PycharmProjects/pdb/code')

#merge all files
filenames = os.listdir('../data/hotspot from pdb_homo for g4')
with open('../result/hotspot_homo_g4', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_homo for g4/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())


#merge all files
filenames = os.listdir('../data/hotspot from pdb_ex for PDETOH_high')
with open('../result/hotspot from pdb_ex for PDETOH_high', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_ex for PDETOH_high/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())



filenames = os.listdir('../data/hotspot from pdb_homo for PDETOH_high')
with open('../result/hotspot from pdb_homo for PDETOH_high', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_homo for PDETOH_high/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())



#merge all files
filenames = os.listdir('../data/hotspot from pdb_ex for PDETOH_medium')
with open('../result/hotspot from pdb_ex for PDETOH_medium', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_ex for PDETOH_medium/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())



filenames = os.listdir('../data/hotspot from pdb_homo for PDETOH_medium')
with open('../result/hotspot from pdb_homo for PDETOH_medium', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_homo for PDETOH_medium/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())



#merge all files
filenames = os.listdir('../data/hotspot from pdb_ex for PDETOH_low')
with open('../result/hotspot from pdb_ex for PDETOH_low', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_ex for PDETOH_low/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())



filenames = os.listdir('../data/hotspot from pdb_homo for PDETOH_low')
with open('../result/hotspot from pdb_homo for PDETOH_low', 'w') as outfile:
    for fname in filenames:
        fname0 = '../data/hotspot from pdb_homo for PDETOH_low/' + fname
        with open(fname0) as infile:
            outfile.write(infile.read())



