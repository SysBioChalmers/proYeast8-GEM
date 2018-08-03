'''
https://github.com/bbuchfink/diamond/blob/master/diamond_manual.pdf
give the column name information for the blast analysis
Gang Li  2018-05-01
'''


import xlrd

infile_map = '../Results/yeast_3D_information.xls'
balst_file = '../PDB/pdb_against_orf.tsv'



# 1. load blast score
score = dict()
for line in open(balst_file):
    cont = line.split('\t')
    score[(cont[1],cont[0].split('_')[0])] = line
# 1. get pdb list with gene list
book = xlrd.open_workbook(infile_map)
sh = book.sheet_by_index(0)

outfile = '../Results/map_exppdb_prot_id.tsv'
fhand = open(outfile,'w')
fhand.write('#qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalues\tbitscore\n')
for i in range(sh.nrows-1):
    if str(sh.cell(i+1,4).value) == 'FALSE': continue
    key = (str(sh.cell(i+1,1).value),str(sh.cell(i+1,3).value))
    if key in score.keys():
        fhand.write(score[key])
fhand.close()