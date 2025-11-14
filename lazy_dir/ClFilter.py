import os
import json

class ClFilter():
    JSON_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'JSON_Annotations'
    )
    THRESHOLD = 0.5
    
    @staticmethod
    def run() -> None:
        data = ClFilter._load_json(os.path.join(ClFilter.JSON_DIRECTORY, 'CL.json'))
        entries_below_threshold = ClFilter._get_below_threshold(data)
        
        '''
        for json_file in os.listdir(ClFilter.JSON_DIRECTORY):
            if 'rejected' in json_file:
                continue
            json_path = os.path.join(ClFilter.JSON_DIRECTORY, json_file)
            data = ClFilter._load_json(json_path)
            ClFilter._remove_from_data(entries_below_threshold, data)
            ClFilter._write_json(json_path, data)
        '''
        
        ClFilter._remove_from_data(entries_below_threshold, data)
        ClFilter._write_json(os.path.join(ClFilter.JSON_DIRECTORY, 'CL.json'), data)
        
    @staticmethod
    def _load_json(path: str) -> dict:
        with open(path, 'r') as f:
            return json.load(f)
        
    @staticmethod
    def _get_below_threshold(data: dict) -> dict:
        entries_below_threshold = {}
        for pdb, pdb_info in data.items():
            for first_residue, first_residue_info in pdb_info.items():
                for second_residue, descriptions in first_residue_info.items():
                    for description in descriptions:
                        if description == 'nbp':
                            continue
                        _, weight = description.split(' ')
                        if float(weight) < ClFilter.THRESHOLD:
                            if pdb not in entries_below_threshold:
                                entries_below_threshold[pdb] = {}
                            if first_residue not in entries_below_threshold[pdb]:
                                entries_below_threshold[pdb][first_residue] = []
                            entries_below_threshold[pdb][first_residue].append(second_residue)
                            
        return entries_below_threshold
    
    @staticmethod
    def _remove_from_data(entries_below_threshold: dict, data: dict) -> None:
        for pdb, pdb_info in entries_below_threshold.items():
            if pdb not in data:
                continue
            for first_residue, second_residues in pdb_info.items():
                if first_residue not in data[pdb]:
                    continue
                for second_residue in second_residues:
                    if second_residue in data[pdb][first_residue]:
                        del data[pdb][first_residue][second_residue]
                    
    @staticmethod
    def _write_json(json_path: str, data: dict) -> None:
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=4)
        
ClFilter.run()