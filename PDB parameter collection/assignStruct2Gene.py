'''# python 3
    Structure information of yeast genes. Structure info is from SWISS-Model.
    Export an excel file.

    sheet:
    (1) Protein structures
        seq_uniprot,seq_len,pdb_id,struct_is_experimental,from,to,struct_percent_seq_cov,struct_resolution,struct_chemicals,qmean,qmean_norm,url

    (2) Domain structures
        seq_uniprot,locus,domain,from,to,pdb_id,struct_is_experimental,struct_resolution,struct_chemicals,qmean,qmean_norm,url

Gang Li
2018-02-10
'''
import xlwt, xlrd


swiss_model_file = '../SWISS-MODEL_Repository/INDEX'
uniport_id_file = '../Uniprot/uniprot-saccharomyces+cerevisiae_detailed version (ID in difference database).xlsx'
struct_info_file = '../PDB/ligandsAndresolution.txt'
pfam_domain_file = '../Pfam/559292.tsv'

book = xlwt.Workbook()
protein_sheet = book.add_sheet('Protein')
domain_sheet = book.add_sheet('Domain')

# 1. write head lines
protein_heads = ['seq_uniprot',
                 'locus',
                 'seq_len',
                 'pdb_id',
                 'struct_is_experimental',
                 'from',
                 'to',
                 'struct_percent_seq_cov',
                 'struct_resolution',
                 'struct_chemicals',
                 'qmean',
                 'qmean_norm',
                 'url']
for i in range(len(protein_heads)): protein_sheet.write(0,i,protein_heads[i])


# 2. Load uniprot_id and locus id
wb = xlrd.open_workbook(uniport_id_file)
sh = wb.sheet_by_index(0)
uniprot2locus = dict()
for i in range(sh.nrows-1):
    i +=1
    uniprot2locus[sh.cell(i,0).value] = sh.cell(i,6).value

# 3. Load ligands and resolution of structures
ligands, resolutions = dict(), dict()
for line in open(struct_info_file):
    if line.startswith('#'): continue
    cont = line.split('\t')
    ligands[cont[0]] = cont[2]
    resolutions[cont[0]] = cont[1]


# 5. Write the protein structure information
k = 1
struct_info = dict()
# struct_info={uniprot_id:{(from,to):['pdb_id','struct_is_experimental','struct_resolution','struct_chemicals','qmean','qmean_norm','url']}}
for line in open(swiss_model_file):
    if line.startswith('#'): continue
    if line.startswith('UniProtKB_ac'): continue
    cont = line.split('\t')

    if cont[4] == 'PDB': is_exp = 'TRUE'
    else: is_exp = 'FALSE'

    pdb_id = cont[8].split('.')[0]

    uniprot_id = cont[0].split('-')[0]

    info_dict = struct_info.get(cont[0],{})
    info_dict[(int(cont[5]),int(cont[6]))] = [cont[8],is_exp,resolutions[pdb_id],ligands[pdb_id],cont[9],cont[10],cont[11]]
    struct_info[uniprot_id] = info_dict
    record_protein = [cont[0],                      #'seq_uniprot'
                      uniprot2locus[uniprot_id],    #'Locus'
                      cont[2],                      #'seq_len'
                      cont[8],                      #'pdb_id'
                      is_exp,                       #'struct_is_experimental'
                      cont[5],                      #'from'
                      cont[6],                      #'to'
                      cont[7],                      #'struct_percent_seq_cov'
                      resolutions[pdb_id],          #'struct_resolution',
                      ligands[pdb_id],              #'struct_chemicals',
                      cont[9],                      #'qmean'
                      cont[10],                     #'qmean_norm'
                      cont[11]                      #'url'
                      ]
    for i in range(len(record_protein)): protein_sheet.write(k,i,record_protein[i])
    k+=1

# 6. Write domain information
domain_heads_2 = [  'pdb_id',
                    'struct_is_experimental',
                    'struct_resolution',
                    'struct_chemicals',
                    'qmean',
                    'qmean_norm',
                    'url']
k = 1
for line in open(pfam_domain_file):
    if line.startswith('#') and 'seq id' not in line: continue
    if line.startswith('#'):
        cont = line.strip().replace('>','').split('<')
        col = 0
        for item in cont[1:]:
            domain_sheet.write(0,col,item.strip())
            col += 1
        for item in domain_heads_2:
            domain_sheet.write(0,col,item)
            col += 1
        continue

    # get structures mapped to this domains
    mapped_pdb = list()
    cont = line.strip().split('\t')
    if struct_info.get(cont[0],None) is not None:
        for pos, info in struct_info[cont[0]].items():
            if int(cont[1]) >= pos[0] and int(cont[2]) <= pos[1]:
                mapped_pdb.append(info)
    if len(mapped_pdb) == 0:
        for i in range(len(cont)):
            domain_sheet.write(k,i,cont[i])
        k += 1
    else:
        for pdb_info in mapped_pdb:
            col = 0
            for item in cont:
                domain_sheet.write(k,col,item)
                col += 1
            for item in pdb_info:
                domain_sheet.write(k,col,item)
                col += 1
            k+=1

book.save('../Results/yeast_3D_information.xls')
