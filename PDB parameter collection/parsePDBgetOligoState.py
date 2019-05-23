'''
Parse pdb websites to get oligo state of given list of pdb ids
Output: text file

Gang Li
2018-04-28
'''
import urllib2
from bs4 import BeautifulSoup
import xlrd

infile = 'Experiment and template PDB id_v2.xlsx'
book = xlrd.open_workbook(infile)
sh = book.sheet_by_index(0)

outfile = '../Results/oligo_state_raw.txt'
all_done_pdb = list()
for item in open(outfile).read().split('//'):
    if len(item)<1: continue
    all_done_pdb.append(item.split('\n')[1])
print (all_done_pdb)


fhand = open(outfile.replace('raw','raw1'),'w')
for i in range(sh.nrows-1):
    pdb_id = sh.cell(i+1,0).value
    if pdb_id in all_done_pdb: continue
    url = 'https://www.rcsb.org/structure/'+pdb_id
    page = urllib2.urlopen(url)
    soup = BeautifulSoup(page,'html.parser')
    contents = soup.find_all('br')
    fhand.write('//\n{}\n'.format(pdb_id))
    print (pdb_id)
    for item in contents:
        cont = item.text.encode('utf-8')
        if 'Global Stoichiometry' in cont and '//' not in cont:
            fhand.write(cont+'\n')
            print (cont)
    print
fhand.close()
