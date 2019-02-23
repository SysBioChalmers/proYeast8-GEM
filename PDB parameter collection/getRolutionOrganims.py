'''
Get the source organism and resolution information for a given list of proteins

Gang Li
2018-04-28
'''

import xlrd
import pypdb

infile = 'Experiment and template PDB id_v2.xlsx'
book = xlrd.open_workbook(infile)
sh1 = book.sheet_by_index(0)
sh2 = book.sheet_by_index(2)

outfile = '../Results/resolution_source_organism_template.txt'
fhand = open(outfile,'w')

def get_pdb_list(shs):
    pdb_list = list()
    for sh in shs:
        for i in range(sh.nrows-1):
            pdb_id = sh.cell(i+1,0).value
            pdb_list.append(pdb_id)
    return pdb_list

pdb_list = get_pdb_list([sh1,sh2])
print(len(pdb_list))

# get source organism dict
orgs = dict()
for line in open('../PDB/source.idx'):
    line = line.replace('\n','')
    cont = line.split('\t')
    if len(cont)<2: continue
    orgs[cont[0].lower()] = cont[1].lower().title()

# get resolution dictionary
resols = dict()
for line in open('../PDB/resolu.idx'):
    if ';' not in line: continue
    cont = line.split(';')
    try:
        resols[cont[0].strip().lower()] = float(cont[1].strip())
    except: print(cont[1].strip())

# get deposition_date

'''
for line in open('../PDB/status_query.csv'):
    cont = line.split(',')
    dates[cont[1].replace("'",'').lower()] = cont[3].replace("'",'')
    print cont[1].replace("'",'').lower(),cont[3].replace("'",'')
'''

# save Results
fhand.write('#pdb_id\tresolution\tsource_org\tdeposition_date\n')
for i in range(len(pdb_list)):
    pdb_id = pdb_list[i]
    date = pypdb.describe_pdb(pdb_id).get('deposition_date','None')
    print(pdb_id,date)
    fhand.write('{}\t{}\t{}\t{}\n'.format(pdb_id,
                                          resols.get(pdb_id,'None'),
                                          orgs.get(pdb_id,'None'),
                                          date))
fhand.close()
