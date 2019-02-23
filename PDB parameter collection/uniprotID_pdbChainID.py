"""this code is showing two methods to establish the relations between uniprotID and pdb chainID"""
import os    ##for directory
import pypdb
import pprint

#First method
os.chdir('/Users/luho/PycharmProjects/pdb/code')
all_info = pypdb.get_all_info('6fai')
pprint.pprint(all_info)

#second method
#download all the Structure-chainID-UniprotID mapping from http://www.ebi.ac.uk/pdbe/docs/sifts/