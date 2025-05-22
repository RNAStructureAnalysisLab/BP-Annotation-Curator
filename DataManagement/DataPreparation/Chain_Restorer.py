# AUTHOR: Kristopher Church
# -----------------------------------------------------------------------------
# PURPOSE:
# Converts the dssr annotations into CSV files and stores them into
# 'Data/AnnotationTools/Annotations'. When converted to a CSV, the chain IDs
# found in the dssr file, are replaced with the original chain ID found in the
# corresponding PDBx file. Applies this chain ID replacement to non dssr
# annotations as well
# NOTE: does not transfer over entries that contain modified nucleosides
# because it seems that other tools don't include those either, so cases where
# a modified nucleoside exists can't be fairly compared across tools since the
# other tools will always fail.
# -----------------------------------------------------------------------------

import os
import re
import pandas as pd

class Chain_Restorer:
    ANNOTATION_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    DSSR_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'DSSR_Annotations'
    )
    PDBX_DIRECTORY = os.path.join('Data', 'Raw', 'RCSB', 'PDBx_Files')
    RESIDUE_TEMPLATE = r'([a-zA-Z]+|\d[a-zA-Z]+|\d)(-?\d+)'
    ORIGINAL_PDB_CHAINS = {} # initialized by restore()
    
    @staticmethod
    def restore():
        Chain_Restorer._load_original_chains(
            Chain_Restorer.ORIGINAL_PDB_CHAINS
        )
        Chain_Restorer._update_non_dssr_annotations()
        Chain_Restorer._add_dssr_annotations()
    
    # INPUT: Receives an empty dictionary 'chain_dictionary' meant to store a
    # PDB ID string as a key, and a nested dictionary as the value. The nested
    # dictionary has strings as the key and value where the key is current
    # chain ID found in the DSSR annotation, and the value is the original 
    # chain ID it must be replaced with.
    # ACTION: Modifies in place 'chain_dictionary' to have the aforementioned
    # contents
    @staticmethod
    def _load_original_chains(chain_dictionary):
        input_file_name = os.path.join(
            'Data', 'Raw', 'RCSB', 'original_pdbx_chains.txt'
        )
        
        with open(input_file_name, 'r') as file:
            for line in file:
                line = line.strip()
                pdb_id, nested_dictionary = line.split(':', 1)
                values = nested_dictionary.strip()[1:-1].split(',')
                nested_dictionary = {}
                
                for chain_id_pair in values:
                    current_id, original_id = chain_id_pair.split(':')
                    current_id = current_id.split("'")[1]
                    original_id = original_id.split("'")[1]
                    nested_dictionary[current_id] = original_id
                    
                chain_dictionary[pdb_id] = nested_dictionary
          
    # ACTION: Gets a list of files that were PDBx. Then visits
    # 'Data/Raw/AnnotationTools/Annotations' and replaces the chain IDs of PDB
    # files that used to be PDBx
    @staticmethod
    def _update_non_dssr_annotations():
        pdbx_file_names = os.listdir(Chain_Restorer.PDBX_DIRECTORY)
        for pdbx in pdbx_file_names:
            pdbx = os.path.splitext(pdbx)[0] # get only the PDBx ID
            if Chain_Restorer._was_remapped(pdbx):
                for tool in ['CL', 'MC', 'FR', 'MO']:
                    pdb_file_name = os.path.join(
                        Chain_Restorer.ANNOTATION_DIRECTORY, 
                        f"{pdbx}_{tool}.csv"
                    )
                    if not os.path.exists(pdb_file_name):
                        continue  # Skip if file doesn't exist
                    df = pd.read_csv(
                        pdb_file_name, dtype={'residue1': str, 'residue2': str}
                    )
                    for i in range(len(df)):
                        Chain_Restorer._convert_chain(i, df, 'residue1', pdbx)
                        Chain_Restorer._convert_chain(i, df, 'residue2', pdbx)
                    df.to_csv(pdb_file_name, index=False)
                        
    
    @staticmethod
    def _add_dssr_annotations():
        pdbx_ids = [
            pdbx.split('.')[0] for pdbx in os.listdir(
                Chain_Restorer.PDBX_DIRECTORY
            )
        ]
        header = r'List of (\d+) base pairs'
        dssr_residue_template = r'([AGCU])(-?\d+)' # disregard modified nucleosides
        dssr_file_names = os.listdir(Chain_Restorer.DSSR_DIRECTORY)
        for pdb in dssr_file_names:
            pdb_input_file_name = os.path.join(
                Chain_Restorer.DSSR_DIRECTORY, pdb
            )
            pdb = pdb.split('.')[0] # Reduce to only the PDB ID
            pdb_output_file_name = os.path.join(
                Chain_Restorer.ANNOTATION_DIRECTORY, f"{pdb}_DSSR.csv"
            )
            rows = []
            
            # Read file line by line. Store relevant base pairing information 
            # in 'rows'
            with open(pdb_input_file_name, 'r') as dssr_file:
                found_base_pairs = False
                base_pair_count = 0
                number_of_base_pairs = None
                for line in dssr_file:
                    if not found_base_pairs:
                        match = re.match(header, line)
                        if match:
                            number_of_base_pairs = int(match.group(1))
                            found_base_pairs = True
                    elif base_pair_count > number_of_base_pairs:
                        break
                    elif base_pair_count > 0:
                        parts = [
                            part.strip() for part in line.split(' ') if part
                        ]
                        residue1 = parts[1].split('.')
                        residue2 = parts[2].split('.')
                        match1 = re.match(dssr_residue_template, residue1[1])
                        match2 = re.match(dssr_residue_template, residue2[1])
                       
                        try:
                            if pdb in pdbx_ids:
                                residue1 = (
                                    Chain_Restorer.ORIGINAL_PDB_CHAINS[pdb][
                                        residue1[0]
                                    ] + match1.group(2)
                                )
                                residue2 = (
                                    Chain_Restorer.ORIGINAL_PDB_CHAINS[pdb][
                                        residue2[0]
                                    ] + match2.group(2)
                                )
                            else:
                                residue1 = residue1[0] + match1.group(2)
                                residue2 = residue2[0] + match2.group(2)
                        except AttributeError:
                            # Should be a modified nucleoside so disregard
                            base_pair_count += 1
                            continue
                        
                        n_type = parts[3].replace('-', '').replace('+', '')
                        contact_type = parts[-1]
                        rows.append({
                            "residue1": residue1, "residue2": residue2, 
                            "n_type": n_type, "description": contact_type
                        })
                        base_pair_count += 1
                    else:
                        base_pair_count += 1
                        
            # Write 'rows' as a csv file in 'ANNOTATION_DIRECTORY'
            df = pd.DataFrame(rows)
            df.to_csv(pdb_output_file_name, index=False)
    
    @staticmethod
    def _was_remapped(pdb_id):
        for key, value in Chain_Restorer.ORIGINAL_PDB_CHAINS[pdb_id].items():
            if key != value:
                return True
        return False
    
    @staticmethod
    def _convert_chain(i, df, column_name, pdb_id):
        residue = df.loc[i][column_name]
        groups = re.match(Chain_Restorer.RESIDUE_TEMPLATE, residue)
        try:
            df.at[i, column_name] = (
                f"{Chain_Restorer.ORIGINAL_PDB_CHAINS[pdb_id][groups.group(1)]}" + 
                f"{groups.group(2)}"
            )
        # This happens when it was already remapped from a previous run
        except KeyError:
            return