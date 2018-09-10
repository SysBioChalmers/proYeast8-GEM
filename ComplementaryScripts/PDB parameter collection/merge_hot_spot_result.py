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