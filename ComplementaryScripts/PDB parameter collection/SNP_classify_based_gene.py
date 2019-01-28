import os    ##for directory
import pandas as pd
import numpy as np

#function
def singleMapping (description, item1, item2, dataframe=True):
    """get the single description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result

gene_feature['cvfName'] = singleMapping(chrom['vcfName'],chrom['annotateName'],gene_feature['chromosome'])
gene_feature.cvfName.notnull()
gene_feature0 = gene_feature[gene_feature.cvfName.notnull()]
len(gene_feature0.cvfName)
gene_feature1 = gene_feature0[['cvfName', 'start', 'end', 'locus_tag']]


def getGeneName(sample, annotation = gene_feature1):
    #sample = SNP_v7[0]
    s1 = sample.split('\t')
    chrom0 = s1[1]
    pos = float(s1[2])
    #exist0 = annotation.query(('cvfName == @chrom0') and ('start <= @pos') and ('end >= @pos')) # the search is faster using this method
    exist0 = annotation[(annotation['cvfName'] == chrom0) & (annotation['start'] <= pos) & (annotation['end'] >= pos)]
    if len(exist0) > 0:
        ss = exist0['locus_tag']
        ss = ss.tolist()[0]
    else:
        ss ='INTEGENIC'
    s1.insert(1, ss)
    s2= '\t'.join(s1)
    return s2

def calculateSNPperGENE(s1):
    #s1 = 'YAL012W'
    s10 = open('../result/' + s1, 'r').readlines()
    snp_num = len(s10)
    return snp_num

# set the directory
# input dataset
os.chdir('/home/hongzhong/PycharmProjects/bigData/code')
gene_feature = pd.read_csv('../data/gene_feature0.txt', sep='\t')
gene_feature['chromosome'] = gene_feature['chromosome'].str.lower()

chrom = pd.read_csv('../data/chrom', header=None, sep='\t')
chrom.columns = ['annotateName','vcfName']
chrom['annotateName'] = chrom['annotateName'].str.lower()

geneGEM = pd.read_excel('../data/gene_list_yeastGEM.xlsx')


# get gene name based on the gene coordinates in the genome
# here in the gene_feature1, it mainly contains the coordinates information for these genes which are functional and could be translated into proteins
out10 = open('../result/SNP_v8', 'w')
with open('../result/SNP_v7', 'r') as SNP_v7:
    for line in SNP_v7:
        line0 = getGeneName(line)
        print(line0)
        out10.writelines(line0)
out10.close()


#classify the sample based on the gene name
def classifyGene(gene):
    #gene ='YAL063C'
    outfile = '../result/' + gene
    out11 = open(outfile, 'w+')
    num = 0
    with open('../result/SNP_v8', 'r') as  SNP_v8:
        for line in SNP_v8:
            line0 = line.split('\t')
            if line0[1] == gene:
                num += 1
                print(num, line)
                out11.writelines(line)
            else:
                continue
    out11.close()




geneGEM0 = geneGEM['geneNames'].tolist()
for x in geneGEM0:
    classifyGene(x)

#summarize the gene SNP occured in the intergenic
notGene = 'INTEGENIC'
classifyGene(notGene)

#check how many genes have been summarized
gene_with_SNP = os.listdir('../result')
gene_ALL = geneGEM['geneNames'].tolist()
gene_no_SNP = [x for x in gene_ALL if x not in gene_with_SNP]


#summarize the SNP number occured in each gene
geneGEM['SNP_NUM'] = [None]*len(geneGEM['geneNames'])
for i, x in geneGEM.iterrows():
    geneGEM['SNP_NUM'][i] = calculateSNPperGENE(geneGEM['geneNames'][i])


export_csv = geneGEM.to_csv('../result/geneGEM with SNP number.csv', index = None, header=True)



#########################################################################
# process all gene

allGene = gene_feature1['locus_tag'].tolist()
for x in allGene:
    classifyGene(x)
len(gene_feature1['locus_tag'])

#check which gene don't have snp data
genListWithSNP = pd.read_csv('../data/gene list with SNP', sep='\t')
genListWithSNP0 = genListWithSNP['x'].tolist()

# find gene with no SNP data
# these genes with SNP data but may be not functional
geneNeedCheck = set(allGene)-set(genListWithSNP0)