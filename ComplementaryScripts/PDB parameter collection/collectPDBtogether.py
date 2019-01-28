'''
Copy all swiss model files together.
2019-01-28 Hongzhong Lu
'''

import shutil
import os

#the dir contains the homology pdb files downloaded from swiss database
source = '/Users/luho/Documents/pdb file/Swiss model information_2018_7_20/SWISS-MODEL_Repository'

#the dir contains the homology pdb files we used to calculate the residue distance
dest1 = '/Users/luho/Documents/pdb file/Swiss model information_2018_7_20/pdb_homo'
#list the files in the directory
for dirName, subdirList, fileList in os.walk(source):
    print('Found directory: %s' % dirName)
    for fname in fileList:
        print('\t%s' % dirName + '/' + fname)
        if fname.endswith('DS_Store'):
            pass
        else:
            shutil.copy(dirName + '/' + fname, dest1)