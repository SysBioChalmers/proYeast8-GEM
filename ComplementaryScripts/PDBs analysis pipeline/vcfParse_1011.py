import os    ##for directory
import pandas as pd

# set the directory
os.chdir('/Users/luho/PycharmProjects/python learning/venv/project2_pdb/data')
print (os.getcwd()) #obtain the present directory

#first step
#fiter the snp from vcf file
#find snp insertion deletion in vcf file
lines = open('1011Matrix.gvcf', 'r').readlines()
out = open('snp','w')
out1 = open('insertion','w')
out2 = open('deletion','w')
flag1 = 0
flag2 = 0
flag3 = 0
for line in lines:
    if line[0] == '#':
        continue
    column = line.split('\t')
    alt = column[4]
    ref = column[3]
    if len(ref) ==1 and alt == 'A':
        out.write(line)
        flag1 +=1
    elif len(ref) ==1 and alt == 'T':
        out.write(line)
        flag1 += 1
    elif len(ref) ==1 and alt == 'C':
        out.write(line)
        flag1 += 1
    elif len(ref) == 1 and alt == 'G':
        out.write(line)
        flag1 += 1
    elif len(ref) == 1 and ',' in alt: # this should be changed into: elif len(ref)==1 and ',' in alt
        out.write(line)
        flag1 += 1
    elif len(alt) > 1 and ',' not in alt:
        out1.write(line)
        flag2 += 1
    elif len(ref) >1:
        out2.write(line)
        flag3 += 1

print('SNPs is:',flag1)
print('Insertion is:', flag2)
print('Deletion is:', flag3)
out.close()
out1.close()
out2.close()





#second step
#input the SNP file
#output file will be the SNP under each sample
#filter sample based on the GT
SNPs = open('snp_strain_name', 'r').readlines()
header0 = SNPs[0].split('\t')
strain_name0 = header0[9:]

def getHomoSample(sample_all, common_row, strain_name=strain_name0):

    index_new = []
    for i, x in enumerate(sample_all):
        if '1/1' in x:
            index_new.append(i)
    sample_all1 = pd.Series(sample_all)
    strain_name1 = pd.Series(strain_name)

    parameter_new = sample_all1[index_new]
    sample_new = strain_name1[index_new]

    parameter = parameter_new.tolist()
    sample = sample_new.tolist()

    commom_row1 = '\t'.join(common_row)
    newRecord = []
    for i in range(len(sample)):
        record = sample[i] +'\t'+ parameter[i] + '\t' + commom_row1
        newRecord.append(record)
    return newRecord

#read and establish a new file to store the result
SNPs = open('snp', 'r').readlines()
out3 = open('SNP_target','w')

for snp in SNPs:
    snp0 = snp.split('\t')
    common_row2 = snp0[0:9]
    sample_all2 = snp0[9:]
    ss = getHomoSample(sample_all2, common_row2) # loop based on row and then should choose which column be used
    print(ss)
    if ss ==[]:
        continue
    else:
        for item in ss:
            out3.write(item)
            out3.write('\n')

out3.close()


#check the quality of the SNP_target
