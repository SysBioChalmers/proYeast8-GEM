"""summary of pdb information for models simulated using swiss web service
input information is the report.html of the simulation results of each protein
export the summary of pdb parameters.
12th,November, 2018
Hongzhong Lu
"""

# import libraries
import urllib.request
from urllib.error import HTTPError
from bs4 import BeautifulSoup
import re
import pandas as pd
import numpy as np
import os    ##for directory
print (os.getcwd()) #obtain the present directory

# set the directory
os.chdir('/Users/luho/PycharmProjects/pdb/code')

def summarizeProteinStucture(gene0):
   '''
   This function was used to obtain all the paramters about the simulated pdb files from swiss model database
   :param gene0:
   :return:
   '''
   # input protein structure model information
   # gene0 ='YBR294W_2018-03-25'
   url = "file:///Users/luho/Documents/pdb%20file/swiss_model_manual%20simulation/" + gene0 + "/report.html"
   response = urllib.request.urlopen(url)
   html = response.read()
   soup0 = BeautifulSoup(html,"lxml")


   """ title summary"""
   header = soup0.find_all('th')
   #len(header)
   ss = [None] * (len(header))
   for i in range(len(header)):
       ss[i] = soup0.find_all('th')[i].get_text()


   ss_r = remove_values_from_list(ss, 'Ligand')

   x1 = [i for i,x in enumerate(ss_r) if x == 'Added to Model']
   x2 = [x+1 for x in x1]

   ss_new = [v for i, v in enumerate(ss_r) if i not in x1+x2]


   """remove the last table title"""
   tt1 = [i for i,x in enumerate(ss_new) if x == 'Template']
   tt2 = [i for i,x in enumerate(ss_new) if x == 'Description']
   tt3 = list(range(tt1[len(tt1)-1]-1,tt2[len(tt2)-1],1))
   ss_new = [v for i, v in enumerate(ss_new) if i not in tt3]
   #tt4 = [i for i, v in enumerate(ss_new) if v =='QSQE']
   #tt4 = [i-1 for i in tt4]
   ss_new = [v for i, v in enumerate(ss_new) if v!='QSQE']


   """data part summary"""
   header1 = soup0.find_all('td')
    #len(header1)
   ss1 = [None] * (len(header1))
   for i in range(len(header1)):
      ss1[i] = soup0.find_all('td')[i].get_text()

   y1 = [i for i, s in enumerate(ss1) if ' - ' in s]
   y2 = [y+3 for y in y1]

   ss1_simple = ss1[0:(y2[len(y2)-1]+2)] ##remove the content after the last 'pdb' position


   z1 = [i for i, s in enumerate(ss1) if 'PDB' in s]
   z2 = z1[1:] #148 328 end
   y3 = y2[0:(len(y2)-1)] #39 183  start

   mm = len(z2)
   WW = list()
   for i in range(0,mm,1):
      WW = list(range(y3[i],(z2[i]-1),1))+ WW


   ##39:148  --- 183:328
   #t1 = list(range(39,(148-1),1))
   #t2 = list(range(183,(328-1),1))

   ss1_simple0 = [v for i, v in enumerate(ss1_simple) if i not in WW]


   """ further refine the data results"""
   n1 = [i for i,x in enumerate(ss1_simple0) if x == '']
   n2 = [i for i,x in enumerate(ss1_simple0) if x == 'Cβ']
   n3 = [x+1 for x in n2]
   n4 = [i for i,x in enumerate(ss1_simple0) if x == 'All Atom']
   n5 = [x+1 for x in n4]
   n6 = [i for i,x in enumerate(ss1_simple0) if x == 'Torsion']
   n7 = [x+1 for x in n6]
   n8 = [i for i,x in enumerate(ss1_simple0) if x == 'QMEAN']
   n9 = [x+1 for x in n8]
   n10 = [i for i,x in enumerate(ss1_simple0) if x == 'Solvation']
   n11 = [x+1 for x in n10]
   n12 = [i for i,x in enumerate(ss1_simple0) if x == 'QSQE']
   NN = n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11+n12

   ss1_simple01 = [v for i, v in enumerate(ss1_simple0) if i not in NN]

   n13 = [i for i,x in enumerate(ss1_simple01) if 'QMEAN' in x]
   ss1_simple01 = [v for i, v in enumerate(ss1_simple01) if i not in n13]

   #remove the first item and the last two item
   ss1_final = ss1_simple01[1:(len(ss1_simple01)-2)]

   """remove data under QSQE"""
   #ss1_final = [v for i, v in enumerate(ss1_final) if i not in tt4]
   n14 = [i for i,v in enumerate(ss1_final) if v=='HHblits']
   n15 = [i-1 for i in n14]
   n16 = []
   for i in n15:
      if hasNumbers(ss1_final[i]):
         n16= [n16,i]

   n17 = list(flatten(n16))
   ss1_final = [v for i,v in enumerate(ss1_final) if i not in n17]



   """get the data from the above summary"""
   """create dataframe"""
   my_list = ss_new
   df0 = pd.DataFrame(np.array(my_list))
   df0['parameter'] = [None] * len(df0)
   #df0.ix[:,0]

   pp1 = [i for i,x in enumerate(ss_new) if 'Model' not in x]


   for i, j in zip(pp1, ss1_final):
      df0['parameter'][i] = j


   """save the data into dict"""

   start0 = [i for i,x in enumerate(my_list) if 'Model' in x]

   def getData(start, ss, gene):
      df1 = [None] * 17
      for i, j in zip(range(0, 17, 1), ss):
         df1[i] = df0.iloc[j, 1]
      df1[0] = df0.iloc[start, 0]
      df1.insert(len(df1), gene)
      return df1

   # crear dict
   results = dict()
   for i in start0:
      results[str(i)+'_'+gene0] = getData(start=i,ss=range(i,i+17,1),gene=gene0)
   print(results)
   return results


def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]


### function to check whether the string contains number
def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)


def flatten(L):
    for item in L:
        try:
            yield from flatten(item)
        except TypeError:
            yield item

print (hasNumbers('234.12 a22')) # returns True
print (hasNumbers('0')) # returns True




#process the gene name information
proteinWithout3D = pd.read_excel('../data/protein_gene_without3D.xlsx')
proteinWithSeq = proteinWithout3D[proteinWithout3D['sequence'].notnull()]
proteinWithSeq['geneName2'] = proteinWithSeq['geneNames'].str.replace('-','')
geneName = proteinWithSeq['geneName2']
print(geneName)

geneWith3D = os.listdir('/Users/luho/Documents/pdb file/swiss_model_manual simulation')
geneWith3D = [v for i,v in enumerate(geneWith3D) if v !='.DS_Store']
geneWith3D_standard = [i.split('_', 1)[0] for i in geneWith3D]
geneWith3D_standard = [v for i, v in enumerate(geneWith3D_standard) if v !='.DS']


""" find protein without simulation"""
geneNeedCheck = [x for x in geneName if x not in geneWith3D_standard]

proteinStructure = {}
# batch process
for i in range(len(geneWith3D)):
   print(i)
   proteinStructure.update(summarizeProteinStucture(gene0=geneWith3D[i]))

newResults = pd.DataFrame.from_items(proteinStructure.items())
newResults = newResults.T

writer = pd.ExcelWriter('../data/proteinStructure2.xlsx')
newResults.to_excel(writer,'Sheet1')
writer.save()

