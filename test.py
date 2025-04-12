import os

file_path = os.path.join('Data', 'Raw', 'RCSB', 'PDB_Files', '6MKN.pdb')

with open(file_path, 'r') as file:
    for line in file:
        print(line.strip())
        
print(file_path)