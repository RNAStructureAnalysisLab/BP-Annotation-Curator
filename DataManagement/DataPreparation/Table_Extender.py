# AUTHOR: Kristopher Church

import json
import csv
import os
import pandas as pd
from Bio.PDB import PDBParser

class Table_Extender:
    MOTIF_CLUSTER_DIRECTORY = os.path.join(
        'Data', 'Raw', 'R3DMA', 'cluster_tables_3.95'
    )
    JSON_ANNOTATION_DIRECTORY = os.path.join(
        'Data', 'Preprocessed', 'JSON_annotations'
    )
    EXTENDED_TABLES_DIRECTORY = os.path.join('extended_tables')
    USED_PDBS_FILE_PATH = os.path.join(
        'Data', 'Raw', 'RCSB', 'used_pdb_ids.txt'
    )
    PDB_DIRECTORY = os.path.join('Data', 'Raw', 'RCSB', 'PDB_Files')
    PARSER = PDBParser(QUIET=True)
    
    model_cache = {} # key = path to PDB file, value = model of PDB structure
    
    def run():
        for file_name in os.listdir(Table_Extender.MOTIF_CLUSTER_DIRECTORY):
            motif_cluster = pd.read_csv(
                os.path.join(Table_Extender.MOTIF_CLUSTER_DIRECTORY, file_name)
            )
            
            Table_Extender._prepare_dataframe(motif_cluster)
            Table_Extender._fill_table(motif_cluster)
            
            # TODO
            
    # INPUT: A pandas dataframe representing motif cluster table
    # ACTION: Modify dataframe in-place by removing rows where
    # PDBs that weren't used, where the residue_id contains '|', or the chain
    # name is 'ASIT' or 'VVV'
    # Initializes each contact type to be a 6-tuple. Each element in 
    # the 6-tuple a contact type reported from a specific tool in the order
    # of: R3DMA, CL, # TODO finish this
    def _prepare_dataframe(motif_cluster):
        used_pdb_ids = []
        with open(Table_Extender.USED_PDBS_FILE_PATH, 'r') as pdb_ids_file:
            for pdb_id in pdb_ids_file:
                used_pdb_ids.append(pdb_id.strip())
        
        rows_to_drop = []
        for i in range(len(motif_cluster)):
            # Check first if there is anything that invalidates the format
            row = motif_cluster.iloc[i]
            if (
                row['PDB'] not in used_pdb_ids or
                'VVV' in row['Chain(s)'] or
                'ASIT' in row['Chain(s)'] or
                any('|' in str(entry) for entry in row)
            ):
                rows_to_drop.append(i)
                continue
            
            # Convert the contact types into a 6-tuple format
            for column_name in motif_cluster.columns[::-1]:
                if '-' in column_name:
                    contact_type = motif_cluster.at[i, column_name]
                    if isinstance(contact_type, float): # blanks are floats
                        contact_type = 'nbp' # nbp: 'Not a Base Pair'
                    six_tuple = contact_type + ',nbp,nbp,nbp,nbp,nbp'
                    motif_cluster.at[i, column_name] = six_tuple
                else:
                    break
        
        # TODO make a csv showing which rows are dropped
        motif_cluster.drop(rows_to_drop, inplace=True)
     
    # TODO
    def _fill_table(motif_cluster):
        column_names = motif_cluster.columns
        segments = Table_Extender._get_segments(column_names)
        for i in range(len(motif_cluster)):
            row = motif_cluster.iloc[i]
            mapping = Table_Extender._map_chains_to_segments(row, segments)
            
            # Load the JSON files
            json_files = {}
            for tool in ['CL', 'FR', 'MC', 'MO', 'RV']:
                file_path = os.path.join(
                    Table_Extender.JSON_ANNOTATION_DIRECTORY, f'{tool}.json'
                )
                with open(file_path, 'r') as file:
                    json_files[tool] = json.load(file)
            
            # Create a dictionary mapping the residues to their column label
            residue_to_column = {}
            for chain, chain_segments in mapping.items():
                for segment in chain_segments:
                    for column_index in range(segment[0], segment[1] + 1, 2):
                        residue_to_column[
                            f'{chain}{row.iloc[column_index]}' + 
                            f'{row.iloc[column_index - 1]}'
                        ] = column_names[column_index]
            
            # Add residues found in the JSON files
            for chain, chain_segments in mapping.items():
                for segment in chain_segments:
                    for column_index in range(segment[0], segment[1] + 1, 2):
                        residue = (
                            f'{chain}{row.iloc[column_index]}' +
                            f'{row.iloc[column_index - 1]}'      
                        )
                        Table_Extender._add_to_table(
                            motif_cluster, residue_to_column, row['PDB'], 
                            residue, json_files, column_index
                        )
       
    # INPUT: A list with the column_names of a motif cluster table
    # OUTPUT: A list of integer pairs. Each integer is a column index. First
    # element in the pair is when a chain starts, and second is when that chain
    # ends. Each pair is a segment of the chain that is part of the motif
    def _get_segments(column_names):
        pair_list = []
        start_index = -1 # -1 represents 'uninitialized'
        for i in range(len(column_names)):
            name = column_names[i]
            if name.isdigit() and start_index == -1:
                start_index = i
            elif 'break' in name:
                pair_list.append((start_index, i - 1))
                start_index = -1
            elif '-' in name and start_index != -1:
                pair_list.append((start_index, i - 1))
                return pair_list
          
    # INPUT: Pandas series and a list of integer pairs
    # OUTPUT: A dictionary where the key is a chain, and the value is a list of
    # integer pairs. Each pair represents the range of columns in the dataframe
    # corresponding to that chain
    def _map_chains_to_segments(row, segments):
        chains = row['Chain(s)'].split('+')
        chain_count = len(chains)
        
        if chain_count == 1:
            return {chains[0]: segments}
        else:
            mapping = {}
            chain_index = 0
            segment_index = 0
            already_compared = set()
            while segment_index < len(segments):
                chain = chains[chain_index]
                segment = segments[segment_index]
                
                # Logic for handling cases where no mapping is found
                comparison = (chain, segment)
                if comparison in already_compared:
                    raise Exception('Chain-Segment mapping failed')
                else:
                    already_compared.add(comparison)
                
                if Table_Extender._in_pdb(row, chain, segment):
                    mapping[chain] = segment
                    segment_index += 1
                chain_index = (chain_index + 1) % chain_count
                
            return mapping
        # TODO: might break down in rare cases. Consider four segments
        # 0, 1, 2, 3 where the correct chain mapping is in the order 
        # A, B, B, A. If segments 0 and 2 by coincidence happen to be the same,
        # then the algorithm will have mapped it as A, B, A, A because it will 
        # have checked with A before checking with B so A claimed it first.
        # Regardless, how would I ever know the correct mapping here since it
        # can be anything in the form x1, B, x2, A where x1 and x2 are A or B?
        # !!!CONSIDER 8VFS in IL_00225.12.csv (the chains could map to either)
          
    def _in_pdb(row, chain, segment):
        pdb_path = os.path.join(
            Table_Extender.PDB_DIRECTORY, f"{row['PDB']}.pdb"
        )
        if pdb_path not in Table_Extender.model_cache:
            model = Table_Extender.PARSER.get_structure('PDB', pdb_path)[0]
            Table_Extender.model_cache[pdb_path] = model
        else:
            model = Table_Extender.model_cache[pdb_path]
        
        if chain not in model:
            return False
        
        # Check if every residue is in the PDB
        for j in range(segment[0], segment[1] + 1, 2):
            found = False
            residue_id = row.iloc[j]
            nucleotide = row.iloc[j - 1]
            
            for residue in model[chain]:
                residue_sequence_number = residue.id[1]
                if (
                    int(residue_sequence_number) == int(residue_id) and 
                    residue.get_resname().upper() == nucleotide
                ):
                    found = True
                    break
                
            if not found:
                return False
               
        return True
        # TODO: perhaps in the preprocessor create a JSON for the models to
        # allow O(1) look-up
        
    # TODO
    def _add_to_table(
        motif_cluster, residue_to_column, pdb, residue, json_files, 
        column_index
    ):
        for tool, file in json_files.items():
            info = file[pdb].get(residue)