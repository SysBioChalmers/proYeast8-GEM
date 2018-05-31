'''
Copy all swiss model files together.

Gang Li
2018-04-15
'''

import os

infile = '../Results/directories.text'
outdir = '../pdb_files/'
if not os.path.exists(outdir): os.mkdir(outdir)

for line in open(infile):
    if 'pdb' in line: os.system('mv ../{} {}'.format(line.strip(),outdir))
os.system('mv INDEX* ..')
os.system('mv README ..')
os.system('rm -r *')

os.system('mv ../INDEX* .')
os.system('mv ../README .')
os.system('mv ../pdb_files .')




## code from Hongzhong
import shutil
import os
import shutil
import os
source = '/Users/luho/PycharmProjects/python learning/venv/project2_pdb/data'
dest1 = '/Users/luho/PycharmProjects/python learning/venv/project2_pdb/pdb_files'
#list the files in the directory
for dirName, subdirList, fileList in os.walk(source):
    print('Found directory: %s' % dirName)
    for fname in fileList:
        print('\t%s' % dirName + '/' + fname)
        if fname.endswith('DS_Store'):
            pass
        else:
            shutil.copy(dirName + '/' + fname, dest1)












