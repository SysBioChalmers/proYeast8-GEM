'''
    Get the ligands and resolution information from PDB.
    1. pdb ids are from swiss model.
    2. export a text file
    #pdb_id,resolution,ligands

Gang Li
2018-02-11
'''

import pypdb

swiss_model_file = '../SWISS-MODEL_Repository/INDEX'
outfile = '../PDB/ligandsAndresolution.txt'

fhand = open(outfile, 'w')
fhand.write('#pdb_id,resolution,ligands\n')
pdb_ids = dict()
for line in open(swiss_model_file):
    if line.startswith('#'): continue
    if line.startswith('UniProtKB_ac'): continue
    cont = line.split('\t')
    pdb_id = cont[8].split('.')[0]
    if pdb_ids.get(pdb_id,0)==0:
        pdb_ids[pdb_id] =1
        chain_id = cont[8].split('.')[-1]
        des = pypdb.describe_pdb(pdb_id)
        res = des.get('resolution',' ')

        ligands_info = pypdb.get_ligands(pdb_id)
        ligands = ''
        try:
            for item in ligands_info['ligandInfo']['ligand']:
                ligands += item.get('@chemicalID',' ')+';'
        except:None
        fhand.write('{0}\t{1}\t{2}\n'.format(pdb_id,res,ligands[:-1]))
        print(pdb_id,res,ligands[:-1])

fhand.close()
