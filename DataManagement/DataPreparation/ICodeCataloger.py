#IMPORTANT
# run step 7. Then run this. Then run the rest of the steps. Then run ValidateCount

import os
import re
import json
from Bio.PDB import PDBParser, Model

class ICodeCataloger:
    PDB_DIRECTORY = os.path.join(
        '..', '..', 'Data', 'Raw', 'RCSB', 'PDB_Files'
    )
    JSON_DIRECTORY = os.path.join(
        '..', '..', 'Data', 'Preprocessed', 'JSON_Annotations'
    )
    RESIDUE_IDS_TXT = os.path.join(
        '.', 'residue_ids.txt'
    )
    
    @staticmethod
    def run() -> None:
        i_code_residue_pairs = set()
        residue_ids = set()
        parser = PDBParser(QUIET=True)
        for pdb in os.listdir(ICodeCataloger.PDB_DIRECTORY):
            pdb_model = ICodeCataloger._load_pdb(pdb, parser)
            ICodeCataloger._add_residues(residue_ids, pdb_model, pdb)
        ICodeCataloger._remove_residues_from_all(residue_ids, i_code_residue_pairs)
        ICodeCataloger._write_residues(i_code_residue_pairs)
    
    #==========================================================================
    # Helper Methods Below
    #==========================================================================

    @staticmethod
    def _load_pdb(pdb: str, parser: PDBParser) -> Model:
        pdb_id = os.path.splitext(pdb)[0]
        structure = parser.get_structure(
            pdb_id, os.path.join(
                ICodeCataloger.PDB_DIRECTORY, pdb
            )
        )
        return structure[0]
    
    @staticmethod
    def _add_residues(residue_ids: set[str], pdb_model: Model, pdb: str) -> None:
        nucleotides = ['A', 'G', 'U', 'C'] #only standard nucleotides
        pdb_id = os.path.splitext(pdb)[0]
        for chain in pdb_model:
            for residue in chain:
                _, sequence_number, insertion_code = residue.id
                residue_name = residue.get_resname()
                if residue_name in nucleotides and insertion_code != ' ':
                    residue_ids.add(f'{pdb_id}_{chain.id}_{sequence_number}')
                    
    @staticmethod
    def _remove_residues_from_all(residue_ids: set[str], i_code_residue_pairs: set[str]) -> None:
        for json_filename in os.listdir(ICodeCataloger.JSON_DIRECTORY):
            if 'rejected' not in json_filename:
                data = ICodeCataloger._get_json(json_filename)
                ICodeCataloger._remove_residues(data, residue_ids, i_code_residue_pairs)
                ICodeCataloger._dump_data(data, json_filename)
                
    @staticmethod
    def _get_json(json_filename: str) -> dict:
        path_to_json = os.path.join(
            ICodeCataloger.JSON_DIRECTORY, json_filename
        )
        with open(path_to_json, 'r') as json_file:
            return json.load(json_file)
        
    @staticmethod
    def _remove_residues(data: dict, residue_ids_to_remove: set[str], i_code_residue_pairs: set[str]) -> None:
        keys_to_remove = set()
        pattern = r'(\d*[A-Za-z]+|\d)(-?\d+)'
        for pdb_id, pdb_info in data.items():
            for first_residue, residue_info in pdb_info.items():
                first_res_matches = re.match(pattern, first_residue)
                first_chain = first_res_matches.group(1)
                first_residue_id = first_res_matches.group(2)
                for second_residue in residue_info.keys():
                    second_res_matches = re.match(pattern, second_residue)
                    second_chain = second_res_matches.group(1)
                    second_residue_id = second_res_matches.group(2)
                    if (
                        f'{pdb_id}_{first_chain}_{first_residue_id}' in residue_ids_to_remove or
                        f'{pdb_id}_{second_chain}_{second_residue_id}' in residue_ids_to_remove
                    ):
                        keys_to_remove.add(f'{pdb_id}_{first_residue}_{second_residue}')
                        i_code_residue_pairs.add(f'{pdb_id}_{first_residue}_{second_residue}') #TODO this only writes if there are still i_code residue paris in the json data, which after one use of this script there won't be any and will be a blank text file, fix this
           
        for key in keys_to_remove:
            parsed_key = key.split('_')
            del data[parsed_key[0]][parsed_key[1]][parsed_key[2]]
            
    @staticmethod
    def _dump_data(data: dict, json_filename: str) -> None:
         file_path = os.path.join(ICodeCataloger.JSON_DIRECTORY, json_filename)
         with open(file_path, 'w') as json_file:
             json.dump(data, json_file, indent=4)
             
    @staticmethod
    def _write_residues(residue_ids: set[str]) -> None:
        with open(ICodeCataloger.RESIDUE_IDS_TXT, 'w') as f:
            for residue_id in sorted(residue_ids):
                f.write(f'{residue_id}\n')
            
ICodeCataloger.run()