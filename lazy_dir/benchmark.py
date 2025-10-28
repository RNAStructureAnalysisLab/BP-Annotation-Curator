import pandas as pd

all_annotations_new = pd.read_csv("base_pairing_agreements.csv")[['All Annotations', 'Expected Contact Type']]
all_annotations_old = pd.read_csv("base_pairing_agreements_old.csv")['All Annotations']

def standardize_contact_types(row, mode):
    if mode == "old":
        row = row.split(',')[1:]
    else:
        row = row.split(',')
    for i in range(len(row)):
        contact_type = row[i]
        if contact_type.startswith('n') and contact_type != 'nbp':
            row[i] = contact_type[1] + contact_type[2:].upper()
        elif contact_type.startswith(('c', 't')):
            row[i] = contact_type[0] + contact_type[1:].upper()
    return row

all_annotations_new['All Annotations'] = all_annotations_new.apply(
    lambda entry: standardize_contact_types(entry['All Annotations'], "new"),
    axis=1
)
all_annotations_old = all_annotations_old.apply(
    lambda entry: standardize_contact_types(entry, "old")
)

rows = []
index_to_tool = {0: 'CL', 1: 'FR', 2: 'MC', 3: 'RV', 4: 'DSSR'}
for contact_type in [
        'cWW', 'cWH', 'cWS', 'cHW', 'cHH', 'cHS', 'cSW', 'cSH', 'cSS',
        'tWW', 'tWH', 'tWS', 'tHW', 'tHH', 'tHS', 'tSW', 'tSH', 'tSS', 'nbp'
]:
    row = {'expected contact type': contact_type, 'CL': 0, 'FR': 0, 'MC': 0, 'RV': 0, 'DSSR': 0}
    for i in range(len(all_annotations_new)):
        expected_contact_type = all_annotations_new.iloc[i]['Expected Contact Type']
        if expected_contact_type != contact_type:
            continue
        
        entry = all_annotations_new.iloc[i]['All Annotations']
        for j in range(len(entry)):
            description = entry[j]
            if expected_contact_type == description:
                row[index_to_tool[j]] += 1
                
            
            #if j == 0 and description == "nbp" and expected_contact_type != description:
                #input("huh?")
    rows.append(row)
pd.DataFrame(rows).to_csv('how_often_tool_agreed_with_expected_bp.csv')



all_annotations_new = all_annotations_new['All Annotations']



rows = []
index_to_tool = {0: 'CL', 1: 'FR', 2: 'MC', 3: 'RV', 4: 'DSSR'}
for contact_type in [
        'cWW', 'cWH', 'cWS', 'cHW', 'cHH', 'cHS', 'cSW', 'cSH', 'cSS',
        'tWW', 'tWH', 'tWS', 'tHW', 'tHH', 'tHS', 'tSW', 'tSH', 'tSS', 'nbp'
]:
    row = {'contact type': contact_type, 'total BP': 0, 'rows all nbp': 0, 'CL': 0, 'FR': 0, 'MC': 0, 'RV': 0, 'DSSR': 0}
    for entry in all_annotations_new:
        if contact_type == 'nbp':
            all_nbp = True
        row['total BP'] += 1
        for i in range(len(entry)):
            description = entry[i]
            if contact_type == description:
                row[index_to_tool[i]] += 1
            else:
                all_nbp = False
        if all_nbp:
            row['rows all nbp'] += 1
            
    rows.append(row)
pd.DataFrame(rows).to_csv('new.csv')
    
rows = []
index_to_tool = {0: 'CL', 1: 'FR', 2: 'MC', 3: 'RV', 4: 'DSSR'}
for contact_type in [
        'cWW', 'cWH', 'cWS', 'cHW', 'cHH', 'cHS', 'cSW', 'cSH', 'cSS',
        'tWW', 'tWH', 'tWS', 'tHW', 'tHH', 'tHS', 'tSW', 'tSH', 'tSS', 'nbp'
]:
    row = {'contact type': contact_type, 'total BP': 0, 'rows all nbp': 0, 'CL': 0, 'FR': 0, 'MC': 0, 'RV': 0, 'DSSR': 0}
    for entry in all_annotations_old:
        if contact_type == 'nbp':
            all_nbp = True
        row['total BP'] += 1
        for i in range(len(entry)):
            description = entry[i]
            if contact_type == description:
                row[index_to_tool[i]] += 1
            else:
                all_nbp = False
        if all_nbp:
            row['rows all nbp'] += 1
    rows.append(row)
pd.DataFrame(rows).to_csv('old.csv')


