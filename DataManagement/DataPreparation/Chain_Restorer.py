# AUTHOR: Kristopher Church
# -----------------------------------------------------------------------------
# PURPOSE:
# TODO: update the purpose description, somewhat outdated
# Creates DSSR CSV files that are stored into 
# 'Data/AnnotationTools/Annotations'. These files are renamed to be consistent
# with the names of the other files, and the chain IDs in the originall DSSR
# file are remapped to the original chain IDs
# -----------------------------------------------------------------------------

import os

class Chain_Restorer:
    OUTPUT_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    INPUT_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'DSSR_Annotations'
    )
    ORIGINAL_PDB_CHAINS = {} # initialized by parse()
    
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
                
    @staticmethod
    def _update_non_dssr_annotations():
        pass
    
    @staticmethod
    def _add_dssr_annotations():
        pass