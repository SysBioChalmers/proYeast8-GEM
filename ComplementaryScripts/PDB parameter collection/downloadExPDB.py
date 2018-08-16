"""the code is mainly used to download the PDB files based on link provided by
swiss database
8-14-2018 Hongzhong Lu
"""

import os
from urllib import request
os.chdir('/Users/luho/PycharmProjects/pdb/code')

#comparison between swiss and pdb
#Two examples to download homology and experimental pdb files
infile ='../result/PDB_ex/'
id1 = '5b3a9bd4a03a1b3eeac03f41' #coordinate id provided by swiss model database
id2 = '5b20b1b602efd01abb492e00'
link1 = 'https://swissmodel.expasy.org/repository/uniprot/P07703.pdb?from=1&to=335&template=5fja.1.C&provider=swissmodel' #link
link2 = 'https://swissmodel.expasy.org/repository/uniprot/P07703.pdb?from=1&to=335&template=5w5y&provider=pdb'

testfile2 = request.urlretrieve(link1, infile+id1+'.pdb')
testfile3 = request.urlretrieve(link2, infile+id2+'.pdb')
