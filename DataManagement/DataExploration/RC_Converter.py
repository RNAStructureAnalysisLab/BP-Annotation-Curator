import os
import json
import pandas as pd

# class converts the contact types found in a residue pair column into its
# absolute RC notation and stores as a json file for fast lookup

# NOTE: the nbp subtypes have overlap. IE: a single nbp might be counted as both
# !cWW and !tSH.

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
        json_data_with_all_tools, json_data_without_R3DMA, json_data_without_FR = {}, {}, {}
        for cluster_table_name in os.listdir(
            RC_Converter.extended_tables_directory
        ):
            cluster_df = pd.read_csv(os.path.join(
                RC_Converter.extended_tables_directory, cluster_table_name
            ))
            RC_Converter._parse_for_json(
                cluster_table_name, cluster_df, json_data_with_all_tools, json_data_without_R3DMA, json_data_without_FR
            )
        with open(os.path.join(RC_Converter.output_directory, 'absolute_rc_notation.json'), 'w') as f:
            json.dump(json_data_with_all_tools, f, indent=4)
        with open(os.path.join(RC_Converter.output_directory, 'absolute_rc_notation_without_R3DMA.json'), 'w') as f:
            json.dump(json_data_without_R3DMA, f, indent=4)
        with open(os.path.join(RC_Converter.output_directory, 'absolute_rc_notation_without_FR.json'), 'w') as f:
            json.dump(json_data_without_FR, f, indent=4)
    
    # PRIVATE METHODS ---------------------------------------------------------
    
    @staticmethod
    def _parse_for_json(cluster_table_name, cluster_df, json_data_with_all_tools, json_data_without_R3DMA, json_data_without_FR):
        for residue_pair_column in cluster_df.columns[::-1]:
            if '-' not in residue_pair_column:
                break
            
            column_data = cluster_df[residue_pair_column].str.split(',')
            column_data_without_R3DMA = column_data.str[1:]
            column_data_without_FR = column_data.str[:2] + column_data.str[3:]
            
            column_data = RC_Converter._standardize(column_data)
            column_data_without_R3DMA = RC_Converter._standardize(column_data_without_R3DMA)
            column_data_without_FR = RC_Converter._standardize(column_data_without_FR)
            
            rc, cc = RC_Converter._get_counts(column_data)
            rc_R3DMA, cc_R3DMA = RC_Converter._get_counts(column_data_without_R3DMA)
            rc_FR, cc_FR = RC_Converter._get_counts(column_data_without_FR)
            
            for contact_type in [
                '!cW.', '!cS.', '!cH.', '!tW.', '!tS.', '!tH.',
                '!c.W', '!c.S', '!c.H', '!t.W', '!t.S', '!t.H',
                'cW.', 'cS.', 'cH.', 'tW.', 'tS.', 'tH.',
                'c.W', 'c.S', 'c.H', 't.W', 't.S', 't.H',
                #'INCOMPATIBLE'
                'REJECT'
            ]:
                absolute_rc_notation = f"r{rc.get(contact_type, 0)}c{cc.get(contact_type, 0)}"
                absolute_rc_notation_R3DMA = f"r{rc_R3DMA.get(contact_type, 0)}c{cc_R3DMA.get(contact_type, 0)}"
                absolute_rc_notation_FR = f"r{rc_FR.get(contact_type, 0)}c{cc_FR.get(contact_type, 0)}"
                
                if cluster_table_name not in json_data_with_all_tools:
                    json_data_with_all_tools[cluster_table_name] = {}
                if residue_pair_column not in json_data_with_all_tools[cluster_table_name]:
                    json_data_with_all_tools[cluster_table_name][residue_pair_column] = {}
                json_data_with_all_tools[cluster_table_name][residue_pair_column][contact_type] = absolute_rc_notation
                
                if cluster_table_name not in json_data_without_R3DMA:
                    json_data_without_R3DMA[cluster_table_name] = {}
                if residue_pair_column not in json_data_without_R3DMA[cluster_table_name]:
                    json_data_without_R3DMA[cluster_table_name][residue_pair_column] = {}
                json_data_without_R3DMA[cluster_table_name][residue_pair_column][contact_type] = absolute_rc_notation_R3DMA
                
                if cluster_table_name not in json_data_without_FR:
                    json_data_without_FR[cluster_table_name] = {}
                if residue_pair_column not in json_data_without_FR[cluster_table_name]:
                    json_data_without_FR[cluster_table_name][residue_pair_column] = {}
                json_data_without_FR[cluster_table_name][residue_pair_column][contact_type] = absolute_rc_notation_FR
        
    @staticmethod
    def _standardize(column_data: pd.Series) -> pd.Series:
        def transform_contact_type(contact_type: str) -> str:
            if contact_type.startswith('nc') or contact_type.startswith('nt'):
                return contact_type[1] + contact_type[2:].upper()
            elif contact_type.startswith('c') or contact_type.startswith('t'):
                return contact_type[0] + contact_type[1:].upper()
            else:
                return contact_type
        def transform_row(lst: list) -> str | list:
            #if any("REJECT" in c for c in lst):
            #    return "INCOMPATIBLE"
            return [transform_contact_type(c) for c in lst]
    
        return column_data.apply(transform_row)
        
    
    @staticmethod
    def _get_counts(column_data: pd.Series):
        rf, cf = {}, {}
    
        for lst in column_data:
            if not lst:  # skip empty lists
                continue
    
            # Count frequencies inside the current list
            nbp_local_count = 0
            contact_type_local_counts = {}
            for elem in lst:
                if elem == 'nbp':
                    nbp_local_count += 1
                if elem == "REJECT":
                    contact_type_local_counts[elem] = contact_type_local_counts.get(elem, 0) + 1
                else:
                    first, second = RC_Converter._extract_edge(elem) # IE: turns cWW to cW. and c.W
                    contact_type_local_counts[first] = contact_type_local_counts.get(first, 0) + 1
                    contact_type_local_counts[second] = contact_type_local_counts.get(second, 0) + 1
                    
            # --- Raw frequency count ---
            nbp_subtypes = []
            for contact_type, local_count in contact_type_local_counts.items():
                rf[contact_type] = rf.get(contact_type, 0) + local_count
                if nbp_local_count > 0:
                    nbp_subtype = f"!{contact_type}"
                    rf[nbp_subtype] = rf.get(nbp_subtype, 0) + nbp_local_count
                    nbp_subtypes.append(nbp_subtype)
            # For rows with only nbp, add all subtypes since denying all
            if len(nbp_subtypes) == 0 and nbp_local_count > 0:
                for subtype in ['!cW.', '!cS.', '!cH.', '!tW.', '!tS.', '!tH.',
                                '!c.W', '!c.S', '!c.H', '!t.W', '!t.S', '!t.H']:
                    nbp_subtypes.append(subtype)
                    rf[subtype] = rf.get(subtype, 0) + nbp_local_count
    
            # --- Mode frequency count ---
            # Find mode(s) in this list
            modes = []
            if len(contact_type_local_counts.values()) == 0:
                max_count = 0
            else:
                max_count = max(contact_type_local_counts.values())
            if nbp_local_count >= max_count:
                max_count = nbp_local_count
                modes = nbp_subtypes
            for k, v in contact_type_local_counts.items():
                if v == max_count:
                    modes.append(k)
    
            # Increment cf for each mode
            for mode in modes:
                cf[mode] = cf.get(mode, 0) + 1
    
        return rf, cf   

    @staticmethod
    def _extract_edge(elem: str) -> (str, str):
        return (elem[0] + elem[1] + '.', elem[0] + '.' + elem[2])