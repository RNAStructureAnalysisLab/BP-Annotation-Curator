import os
import json
import pandas as pd

# class converts the contact types found in a residue pair column into its
# absolute RC notation and stores as a json file for fast lookup

class RC_Converter:
    # directory configuration
    extended_tables_directory = os.path.join(
        'Data', 'Preprocessed', 'Extended_Tables'
    )
    output_directory = os.path.join(
        'Data', 'ExplorationFindings'    
    )
    
    # PUBLIC METHODS ----------------------------------------------------------
    
    @staticmethod
    def run():
        json_data = {}
        for cluster_table_name in os.listdir(
            RC_Converter.extended_tables_directory
        ):
            cluster_df = pd.read_csv(os.path.join(
                RC_Converter.extended_tables_directory, cluster_table_name
            ))
            RC_Converter._parse_for_json(
                cluster_table_name, cluster_df, json_data
            )
        output_file = os.path.join(RC_Converter.output_directory, 'absolute_rc_notation.json')
        with open(output_file, 'w') as f:
            json.dump(json_data, f, indent=4)
    
    # PRIVATE METHODS ---------------------------------------------------------
    
    @staticmethod
    def _parse_for_json(cluster_table_name, cluster_df, json_data):
        for residue_pair_column in cluster_df.columns[::-1]:
            if '-' not in residue_pair_column:
                break
            column_data = cluster_df[residue_pair_column].str.split(',')
            column_data = RC_Converter._standardize(column_data)
            rc, cc = RC_Converter._get_counts(column_data)
            for contact_type in [
                'cWW', 'cWH', 'cWS', 'cSW', 'cSH', 'cSS', 'cHW', 'cHH', 'cHS',
                'tWW', 'tWH', 'tWS', 'tSW', 'tSH', 'tSS', 'tHW', 'tHH', 'tHS',
                'nbp', 'REJECT'
            ]:
                absolute_rc_notation = f"r{rc.get(contact_type, 0)}c{cc.get(contact_type, 0)}"
                
                if cluster_table_name not in json_data:
                    json_data[cluster_table_name] = {}
                if residue_pair_column not in json_data[cluster_table_name]:
                    json_data[cluster_table_name][residue_pair_column] = {}
                
                json_data[cluster_table_name][residue_pair_column][contact_type] = absolute_rc_notation
        
    @staticmethod
    def _standardize(column_data: pd.Series) -> pd.Series:
        def transform_segment(segment: str) -> str:
            if segment.startswith("nc"):
                return "nc" + segment[2:].upper()
            elif segment.startswith("c"):
                return "c" + segment[1:].upper()
            return segment
    
        return column_data.apply(lambda lst: [transform_segment(s) for s in lst])
    
    @staticmethod
    def _get_counts(column_data: pd.Series):
        rf, cf = {}, {}
    
        for lst in column_data:
            if not lst:  # skip empty lists
                continue
    
            # --- Raw frequency count ---
            for elem in lst:
                rf[elem] = rf.get(elem, 0) + 1
    
            # --- Mode frequency count ---
            # Count frequencies inside the current list
            local_counts = {}
            for elem in lst:
                local_counts[elem] = local_counts.get(elem, 0) + 1
    
            # Find mode(s) in this list
            max_count = max(local_counts.values())
            modes = [k for k, v in local_counts.items() if v == max_count]
    
            # Increment cf for each mode
            for mode in modes:
                cf[mode] = cf.get(mode, 0) + 1
    
        return rf, cf