import os
import re
import json
import shutil
import pandas as pd
from collections import Counter

class Tool_Consensus_2:
    extended_tables_directory = os.path.join(
        'Data', 'Preprocessed', 'Extended_Tables'    
    )
    rc_json_path = os.path.join(
        'Data', 'ExplorationFindings', 'absolute_rc_notation.json'
    )
    consensus_directory = os.path.join(
        'Data', 'Consensus', 'Tool_Consensus'
    )
    mode_consensus_directory = os.path.join(
        consensus_directory, 'Mode_Based'
    )
    
    tie_instances = []
    
    # PUBLIC METHODS ----------------------------------------------------------
    
    @staticmethod
    def run():
        if os.path.exists(Tool_Consensus_2.mode_consensus_directory):
            shutil.rmtree(Tool_Consensus_2.mode_consensus_directory)
        os.makedirs(Tool_Consensus_2.mode_consensus_directory)
        
        rc_dictionary = Tool_Consensus_2._load_rc_notation()
        for table_name in os.listdir(Tool_Consensus_2.extended_tables_directory):
            cluster_df = pd.read_csv(os.path.join(
                Tool_Consensus_2.extended_tables_directory, table_name    
            ))
            consensus_cluster_df = Tool_Consensus_2._get_consensus(
                cluster_df, table_name, rc_dictionary
            )
            Tool_Consensus_2._output_df(consensus_cluster_df, table_name)
        df = pd.DataFrame(Tool_Consensus_2.tie_instances)
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'tie_instances.csv'
        ), index=False)
    
    # PRIVATE METHODS ---------------------------------------------------------
    @staticmethod
    def _load_rc_notation():
        with open(Tool_Consensus_2.rc_json_path, 'r') as f:
            return json.load(f)
    
    @staticmethod
    def _get_consensus(cluster_df, table_name, rc_dictionary):
        consensus_df = cluster_df.copy()
        for column_name in cluster_df.columns[::-1]:
            if '-' not in column_name:
                break
            consensus_df[column_name] = cluster_df[column_name].apply(
                lambda cell: Tool_Consensus_2._get_mode(cell, table_name, rc_dictionary, column_name)
            )
        return consensus_df
    
    @staticmethod
    def _get_mode(cell, table_name, rc_dictionary, column_name):
        items = [x for x in cell.split(",")]
        if not items:
            return None

        counts = Counter(items)
        mode, max_frequency = counts.most_common(1)[0]
        top_contact_types = [contact_type for contact_type, frequency in counts.items() if frequency == max_frequency]
        if len(top_contact_types) > 1:
            mode = Tool_Consensus_2._resolve_tie(top_contact_types, table_name, rc_dictionary, column_name)
        
        return mode
    
    @staticmethod
    def _output_df(consensus_cluster_df, table_name):
        consensus_cluster_df.to_csv(
            os.path.join(Tool_Consensus_2.mode_consensus_directory, table_name), index=False
        )
        
    @staticmethod
    def _resolve_tie(top_contact_types, table_name, rc_dictionary, column_name):
        mode = None
        max_count = 0
        for i in [True, False]: # just loop through twice
            checking_consensus_count = i
            for contact_type in top_contact_types:
                if contact_type == 'INCOMPATIBLE' or contact_type == 'REJECT':
                    continue
                rc = rc_dictionary[table_name][column_name][contact_type]
                counts = re.match(r"r(\d+)c(\d+)", rc)
                if checking_consensus_count:
                    consensus_count = int(counts.group(2))
                    if consensus_count == max_count:
                        mode = None
                        break
                    if  consensus_count > max_count:
                        max_count = consensus_count
                        mode = contact_type
                else:
                    report_count = int(counts.group(1))
                    if report_count == max_count:
                        mode = None
                        break
                    if report_count > max_count:
                        max_count = report_count
                        mode = contact_type
            if mode == None:
                max_count = 0
            else:
                Tool_Consensus_2.tie_instances.append({
                    'cluster': table_name, 'PDB': pdb, 'column': column_name, 
                    'tied_contact_types': top_contact_types
                })
                return mode
        
        input("need to add further logic for resolving ties")
        # TODO Further logic for resolving ties
        
        return None