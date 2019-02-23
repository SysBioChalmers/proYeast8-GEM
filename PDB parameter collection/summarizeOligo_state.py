"""the code is mainly used to parse the oligo state
of each experimental pdb structure for file obtained by GangLi.
The origial pdb oligo state information is downloaded from pdb database
https://www.rcsb.org/structure/3qht
Hongzhong Lu 13-8-2018
"""

import os
import numpy as np
import pandas as pd
os.chdir('/Users/luho/PycharmProjects/pdb/code')
infile = "../data/oligo_state_raw.txt"
oligo_state_PDB_ex = open(infile,'r+', encoding="utf-8").readlines()

#split the text based on the unique sign "//"
index0 = []
for i,x in enumerate(oligo_state_PDB_ex):
    if "//" in x:
        index0.append(i)
    else:
        continue

#summarize the paramter for each item
#this function is used to summarize the oligo state for each pdb file
def getOligoState(s1,s2,ss = oligo_state_PDB_ex):
    oligo = dict()
    i, j = s1,s2
    item1 = ss[i:j]
    len0 = len(item1)
    key = item1[1].strip('\n')
    value = item1[2:len0]
    value =[x.strip('\n') for i,x in enumerate(value)]
    oligo[key] = value
    return oligo

#loop to get the ogligo state information of each experiment file
pdb_oligo = {}
for i, j in zip(index0[0:(len(index0)-1)], index0[1:len(index0)]):
    print(i,j)
    s = getOligoState(i,j)
    pdb_oligo.update(s)

#post-process
#this function is used to remove the redundent information of desciption
def simpleOligo(key,dict0):
    #key = '1dpj'
    #dict0 = pdb_oligo
    s = dict0[key]
    s0 = []
    s1 = []
    for x in s:
        x1 = x.split(" -")
        if "Biological assembly" in x1[1]:
            x2 = x1[0] + '@' + 'BA' #BA means Biological assemble
            x2 = x2.replace('Global Stoichiometry: ','')
        else:
            x2 = x1[0]
            x2 = x2.replace('Global Stoichiometry: ', '')

        s0.append(x2)
        s1 = np.unique(np.array(s0)).tolist()
    return s1
#test
simpleOligo(key='3epk',dict0=pdb_oligo)


#loop of dict to simplify the description
pdb_oligo0 ={}
for key0 in pdb_oligo:
    s = simpleOligo(key=key0,dict0=pdb_oligo)
    pdb_oligo0[key0] = s


#this function is used to remove the duplicated oligo state of pdb structure
def removeDuplicatedOligo(key,dict0):
    #key = '4jpo'
    #dict0 = pdb_oligo0
    s = dict0[key]
    s2 = []
    s3 = []
    index2 = []
    index3 = []
    s4 = []
    index4 = []
    for i, x in enumerate(s):
        if '@' in x:
            x0 = x.split('@')
            s2.append(x0[0])
            index2.append(i)
        else:
            s3.append(x)
            index3.append(i)
    s2 = np.unique(s2)
    s3 = np.unique(s3)
    index4 = np.in1d(s3,s2)
    s2 = [x +'@BA' for i, x in enumerate(s2)]
    s4 = np.append(s2,s3[~index4])
    s40 = s4.tolist()
    s40 = ';'.join(s40)
    return s40

#test example
removeDuplicatedOligo('4jpo',dict0=pdb_oligo0)
removeDuplicatedOligo('1pjd',dict0=pdb_oligo0)


#loop to remove the duplicated oligo state
pdb_oligo1 ={}
for key0 in pdb_oligo0:
    s = removeDuplicatedOligo(key=key0,dict0=pdb_oligo0)
    pdb_oligo1[key0] = s
    print(key0, s)

#change the dict into the dataframe
pdb_oligo2 = pd.Series(pdb_oligo1).to_frame()

#save the results in xlsx format
writer = pd.ExcelWriter('../result/Oligo state of experimental structure.xlsx')
pdb_oligo2.to_excel(writer,'Sheet1')
writer.save()