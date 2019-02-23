import os    ##for directory
import pandas as pd

# set the directory
os.chdir('/home/hongzhong/PycharmProjects/bigData/code')
infile = '../data/1011Matrix.gvcf'
outfile = '../result/'

#first step
#fiter the snp from vcf file
#find snp insertion deletion in vcf file

out = open(outfile + 'snp', 'w')
out1 = open(outfile + 'insertion', 'w')
out2 = open(outfile + 'deletion', 'w')
flag1 = 0
flag2 = 0
flag3 = 0

with open(infile) as infile0:
  for line in infile0:
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
    elif len(ref) == 1 and ',' in alt:   # this should be changed into: elif len(ref)==1 and ',' in alt
        out.write(line)
        flag1 += 1
    elif len(alt) > 1 and ',' not in alt:
        out1.write(line)
        flag2 += 1
    elif len(ref) >1:
        out2.write(line)
        flag3 += 1

print('SNPs is:', flag1)
print('Insertion is:', flag2)
print('Deletion is:', flag3)
out.close()
out1.close()
out2.close()





#second step
#input the SNP file
#output file will be the SNP under each sample
#filter sample based on the GT
SNPs = open('../data/snp_strain_name', 'r').readlines()
header0 = SNPs[0].split('\t')
strain_name0 = header0[9:]

def getHomoSample(sample_all, common_row, strain_name=strain_name0):
    index_new = []
    for i, x in enumerate(sample_all):
        if '1/1' in x:
            index_new.append(i)
        elif '2/2' in x:
            index_new.append(i)
        elif '3/3' in x:
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

out3 = open('../result/SNP_target', 'w')

with open('../result/snp') as SNP0:
    for snp in SNP0:
        snp0 = snp.split('\t')
        common_row2 = snp0[0:9]
        sample_all2 = snp0[9:]
        ss = getHomoSample(sample_all2, common_row2) # loop based on row and then should choose which column be used
        print(ss)
        if ss == []:
            continue
        else:
            for item in ss:
                out3.write(item)
                out3.write('\n')

out3.close()



#refine the dataformat
out4 = open('../result/snp1', 'w')
out5 = open('../result/snp2', 'w')
with open('../result/SNP_target', 'r') as TEST:
    for line in TEST:
        print(line)
        if len(line) < 20:
            out5.write(line)
        elif line.startswith("\t"):
            out5.write(line)
        else:
            out4.write(line)

out4.close()
out5.close()




#process snp2
out6 = open('../result/snp20', 'w')
SNP2 = open('../result/snp2', 'r').readlines()
index0 = []
num = 0
for i, line in enumerate(SNP2):
    num += 1
    if len(line) < 10:
        index0.append(i)

for j in index0:
    print(j)
    newstring = SNP2[j:(j+3)]
    newline = newstring[0].strip('\n') + newstring[1].strip('\n') + newstring[2]
    print(newline)
    out6.write(newline)
out6.close()




#check the quality of the SNP_target
# check in genotype of '2/2, 3/3', whether the length of ALT ==1
out7 = open('../result/SNP_v3', 'w')
with open('../result/snp1', 'r') as lines1:
    for s0 in lines1:
        s1 = s0.split('\t')
        if '1/1' in s1[1]:
            out7.write(s0)
        if '2/2' in s1[1]:
            ALT2 = s1[6].split(',')
            if len(ALT2[1]) == 1:
                out7.write(s0)
        elif '3/3' in s1[1]:
            ALT3 = s1[6].split(',')
            if len(ALT3[2]) == 1:
                out7.write(s0)
        print(s0)
out7.close()



#check the quality of the SNP_target
# check in genotype of '2/2, 3/3', whether the length of ALT ==1
out8 = open('../result/SNP_v3_2', 'w')
lines1 = open('../result/snp20', 'r').readlines()
    for s0 in lines1:
        s1 = s0.split('\t')
        if '1/1' in s1[1]:
            out8.write(s0)
        if '2/2' in s1[1]:
            ALT2 = s1[6].split(',')
            if len(ALT2[1]) == 1:
                out8.write(s0)
        elif '3/3' in s1[1]:
            ALT3 = s1[6].split(',')
            if len(ALT3[2]) == 1:
                out8.write(s0)
        print(s0)
out8.close()



#merge the snp1 and snp20 to form the file SNP_target0
filenames = ['../result/SNP_v3_2', '../result/SNP_v3']
with open('../result/SNP_v4', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                print(line)
                outfile.write(line)


#add filter
#site information: QD < 2 , MQ < 40
#sample information: GQ < 30, DP < 5
#sample informaiton: GT:AD:DP:GQ:PGT:PID:PL
out9 = open('../result/SNP_v5','w')
def FilterResult(sample):
    #sample = SACE_YDO[1]
    DP = 0
    GQ = 0
    QD = 0
    MQ = 0
    s2 = sample.split('\t')
    formatInf = s2[1].split(':')
    #[0]GT, [1]AD, [2]DP, [3]GQ, [4]PGT, [5]PID, [6]PL
    try:
        DP =  float(formatInf[2])
    except ValueError:
        pass
    try:
        GQ =  float(formatInf[3])
    except ValueError:
        pass

    #site information
    s2[9] = s2[9].replace('GQ_STDDEV',';GQ_STDDEV')
    siteInf = s2[9].split(';')
    siteInf = [x for i, x in enumerate(siteInf) if x !='']
    siteInf0 = [x.split('=')for i,x in enumerate(siteInf)]
    siteInf1 = dict(siteInf0)

    if siteInf1.get('QD') is None:
        pass
    else:
        try:
            QD = float(siteInf1.get('QD'))
        except ValueError:
            pass
    if siteInf1.get('MQ') is None:
        pass
    else:
        try:
            MQ = float(siteInf1.get('MQ'))
        except ValueError:
            pass

    if QD >= 2  and MQ >= 40 and DP >= 5 and GQ >= 30:
        return True
    else:
        return False


num = 0
with open ('../result/SNP_v4', 'r') as SNP_v4:
    for line in SNP_v4:
        highQuality = FilterResult(line)
        if highQuality:
            num +=1
            out9.writelines(line)
        else:
            continue
        print(highQuality,num)


#further simplify the information in each line
out9 = open('../result/SNP_v6','w')

def simpleSampleInf(sample):
    s1 = sample.split('\t')
    formatInf = s1[1].split(':')
    genotype = formatInf[0]
    if '2/2' in genotype:
        s1[6] = s1[6].split(',')[1]
    elif '3/3' in genotype:
        s1[6] = s1[6].split(',')[2]
    elif '4/4' in genotype:
        s1[6] = s1[6].split(',')[3]
    else:
        s1[6] = s1[6].split(',')[0]

    s2 = s1[0] + '\t' + s1[2] + '\t' + s1[3] + '\t' + s1[5] + '\t' +s1[6] + '\n'
    return s2

out9 = open('../result/SNP_v6','w')
with open('../result/SNP_v5','r') as SNP_v5:
    for line in SNP_v5:
        newline = simpleSampleInf(line)
        print(newline)
        out9.writelines(newline)


#further check the snp quality and remove the snp with alt lenght >2
out10 = open('../result/SNP_v7','w')
SNP_v6 = open('../result/SNP_v6','r').readlines()
num = []
for i, line in enumerate(SNP_v6):
    ss = line.split('\t')
    if len(ss) == 5 and len(ss[4].strip('\n')) == 1:
        out10.writelines(line)
    else:
        num.append(i)
        continue
    print(line)