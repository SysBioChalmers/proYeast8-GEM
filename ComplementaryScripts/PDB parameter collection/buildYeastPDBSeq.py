'''
build the fasta file for pdb structures form S. cerevisiae

Gang Li
2018-04-28
'''
import xlrd
from Bio import SeqIO

infile = 'Experiment and template PDB id_v2.xlsx'
book = xlrd.open_workbook(infile)
sh = book.sheet_by_index(0)

isYeast = dict()
for i in range(sh.nrows-1):
    isYeast[sh.cell(i+1,0).value] = True

outfile = '../PDB/yeast_pdb.fasta'
fhand = open(outfile,'w')
for rec in SeqIO.parse('../PDB/pdb_seqres.txt','fasta'):
    pdb = rec.id.split('_')[0]
    if not isYeast.get(pdb,False): continue
    fhand.write('>{}\n'.format(rec.id))
    fhand.write('{}\n'.format(rec.seq))
fhand.close()
