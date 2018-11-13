# -*- coding: utf-8 -*-
'''this code is to read latest yeastGEM, estabolish the gene-protein-reaction relation, then the protein information can be merged into this dataframe
based on the geneID mapping
12th, Nov, 2018
Hongzhong Lu
'''

# Import packages
import pandas as pd
import os    ##for directory
from cobra.io import read_sbml_model

os.chdir('/Users/luho/PycharmProjects/model/cobrapy/code')


# import self function
def getRXNgeneMapping(rxn0, gpr0):
    '''this function is used to split the GPR;
    input, for example rxn0=['r1','g2']
    gpr0=['a or c','a and b']
    output, each rxn related with each gene'''
    s1 = rxn0
    s2 = gpr0
    s2 = s2.str.replace('and','@')
    s2 = s2.str.replace('or','@')
    s2 = s2.str.replace('\\( ','')
    s2 = s2.str.replace('\\(\\( ','')
    s2 = s2.str.replace('\\(', '')
    s2 = s2.str.replace('\\(\\(', '')
    s2 = s2.str.replace(' \\)','')
    s2 = s2.str.replace(' \\)\\) ','')
    s2 = s2.str.replace('\\)', '')
    s2 = s2.str.replace('\\)\\) ', '')
    s3 = splitAndCombine(s2,s1,sep0="@")
    s3['V2'] = s3['V2'].str.strip()
    s3.columns = ['rxnID', 'gene']
    return s3

def correctSomeWrongFormat(model0):
  """
  This function is used to correct some wrong format when read yeastGEM model from cobratoolbox
  """
  # Correct metabolite ids:
  for met in model0.metabolites:
    met.id = met.id.replace('__91__', '_')
    met.id = met.id.replace('__93__', '')
  #for reaction in model0.reactions:
  #    reaction.gene_reaction_rule = reaction.gene_reaction_rule.replace('__45__','-')
  for gene in model0.genes:
      gene.id = gene.id.replace('__45__', '-')

  return model0



def saveExcel(infile, outfile):
    '''
    function to save the dataframe into xlsx format
    :param infile:
    :param outfile:
    :return:
    '''
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer,'Sheet1')
    writer.save()



# input the subsystem information
gem_dataframe = pd.read_excel('/Users/luho/PycharmProjects/model/model_correction/result/yeastGEM_with subsystem.xlsx')

# input yeast8 for every update from yeastGEM repo
GEM_nov = read_sbml_model('../data/yeastGEM_nov.xml')
GEM_nov= correctSomeWrongFormat(GEM_nov)
#produce the dataframe for the metabolites and the rxn

gem_rxn_nov = produceRxnList(GEM_nov)


#establish rxn-gene mapping
proYeast_DataFrame = getRXNgeneMapping(gem_rxn_nov['rxnID'], gem_rxn_nov['GPR'])
proYeast_DataFrame0 = pd.merge(proYeast_DataFrame, gem_rxn_nov, on='rxnID')
saveExcel(proYeast_DataFrame0, '../result/proYeast_DataFrame_from_yeastGEM_november.xlsx')












