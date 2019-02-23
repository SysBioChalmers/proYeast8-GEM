'''
Extract the model and template information from atom files.

Gang Li
2018-04-15
'''

import xlwt
import os

index_file = '../SWISS-MODEL_Repository/INDEX'
atom_dir = '../SWISS-MODEL_Repository/pdb_files/'

# 1. get coordinate id from INDEX file
coordinate2uniprot= dict()
for line in open(index_file):
    if line.startswith('UniProtKB_ac') or line.startswith('#'): continue
    cont = line.split('\t')
    coordinate2uniprot[cont[3]] = cont[0]

# 2. get coordinate id from swiss model file names
coordinate_ids = list()
for name in os.listdir(atom_dir):
    if name.endswith('pdb'): coordinate_ids.append(name.split('_')[-1].split('.')[0])
print len(coordinate2uniprot),len(coordinate_ids)

# 3. get model infomation from atom_files
info = dict()
model_attrs = list()
temp_attrs = list()
# info = {'coordinate_id':[{model_attr},{temp_attr}]}
for name in os.listdir(atom_dir):
    if name.endswith('pdb'):
        coor_id = name.split('_')[-1].split('.')[0]
        model_attr, temp_attr = dict(),dict()
        model, temp = False, False
        for line in open(os.path.join(atom_dir,name)):
            if line.startswith('ATOM'): break
            if 'REMARK   3 MODEL INFORMATION' in line:
                model = True
                continue
            if 'REMARK   3 MODEL LIGAND' in line:
                model = False
                continue
            if 'REMARK   3 TEMPLATE 1' in line:
                model = False
                temp = True
                continue
            if 'LIGND' in line: continue

            cont = line.strip().split()
            if len(cont) == 2: continue
            if model:
                model_attr[cont[2]] = cont[3]
                if cont[2] not in model_attrs: model_attrs.append(cont[2])

            if temp:
                temp_attr[cont[2]] = cont[3]
                if cont[2] not in temp_attrs: temp_attrs.append(cont[2])
        info[coor_id] = [model_attr,temp_attr]

# 4. Save tesults
# UniProtKB_ac,coordinate_id, model_attrs,temp_attrs

wb = xlwt.Workbook()
ws = wb.add_sheet('model_info')
ws.write(0,0,'UniProtKB_ac')
ws.write(0,1,'coordinate_id')
k = 2
for item in model_attrs:
    ws.write(0,k,item)
    k += 1
k += 1
for item in temp_attrs:
    ws.write(0,k,item)
    k += 1

row = 1
for key,data in info.items():
    ws.write(row,0,coordinate2uniprot[key])
    ws.write(row,1,key)
    col = 2
    for item in model_attrs:
        ws.write(row,col,data[0].get(item,''))
        col += 1
    col += 1

    for item in temp_attrs:
        ws.write(row,col,data[1].get(item,''))
        col += 1

    row += 1
wb.save('../Results/model_temp_info.xls')
