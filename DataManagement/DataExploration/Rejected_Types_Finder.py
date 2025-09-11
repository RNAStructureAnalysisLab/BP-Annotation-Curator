import os
import json

class Rejected_Types_Finder:
    REJECTED_JSON_DIRECTORY = os.path.join('Data', 'Preprocessed', 'JSON_Annotations')
    DATA_EXPLORATION_DIRECTORY = os.path.join('Data', 'ExplorationFindings')
    
    rejected_residue_interactions = set()
    
    # --- PUBLIC FUNCTIONS ---

    @staticmethod
    def run():
        for json_file_name in os.listdir(Rejected_Types_Finder.REJECTED_JSON_DIRECTORY):
            if "_rejected" not in json_file_name:
                continue
            json_file_path = os.path.join(Rejected_Types_Finder.REJECTED_JSON_DIRECTORY, json_file_name)
            with open(json_file_path, "r") as f:
                json_data = json.load(f)
                Rejected_Types_Finder._extract_strings(json_data)
                
            Rejected_Types_Finder._write_to_text_file()

    # --- PRIVATE FUNCTIONS ---

    @staticmethod
    def _extract_strings(obj):
        """
        Recursively traverse JSON-like data structures (dicts/lists/scalars).
        If a list of strings is found, add its contents to the set.
        """
        if isinstance(obj, dict):
            for value in obj.values():
                Rejected_Types_Finder._extract_strings(value)
        elif isinstance(obj, list):
            if all(isinstance(item, str) for item in obj):
                Rejected_Types_Finder.rejected_residue_interactions.update(obj)
            else:
                for item in obj:
                    Rejected_Types_Finder._extract_strings(item)
                    
    @staticmethod
    def _write_to_text_file():
        os.makedirs(Rejected_Types_Finder.DATA_EXPLORATION_DIRECTORY, exist_ok=True)
        output_path = os.path.join(
            Rejected_Types_Finder.DATA_EXPLORATION_DIRECTORY, "rejected_types.txt"
        )
        with open(output_path, "w") as f:
            for item in sorted(Rejected_Types_Finder.rejected_residue_interactions):
                f.write(f"{item}\n")