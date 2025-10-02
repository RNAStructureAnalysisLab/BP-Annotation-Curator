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
    rc_R3DMA_json_path = os.path.join(
        'Data', 'ExplorationFindings', 'absolute_rc_notation_without_R3DMA.json'
    )
    rc_FR_json_path = os.path.join(
        'Data', 'ExplorationFindings', 'absolute_rc_notation_without_FR.json'
    )
    consensus_directory = os.path.join(
        'Data', 'Consensus', 'Tool_Consensus'
    )
    mode_consensus_directory = os.path.join(
        consensus_directory, 'Mode_Based'
    )
    
    # PUBLIC METHODS ----------------------------------------------------------
    
    @staticmethod
    def run():
        tie_instances = []
        unresolved_tie_instances = []
        
        tie_instances_without_R3DMA = []
        unresolved_tie_instances_without_R3DMA = []
        
        tie_instances_without_FR = []
        unresolved_tie_instances_without_FR = []
        
        if os.path.exists(Tool_Consensus_2.mode_consensus_directory):
            shutil.rmtree(Tool_Consensus_2.mode_consensus_directory)
        os.makedirs(Tool_Consensus_2.mode_consensus_directory)
        os.makedirs(os.path.join(Tool_Consensus_2.mode_consensus_directory, 'AllTools'))
        os.makedirs(os.path.join(Tool_Consensus_2.mode_consensus_directory, 'WithoutR3DMA'))
        os.makedirs(os.path.join(Tool_Consensus_2.mode_consensus_directory, 'WithoutFR'))
        
        #rc_dictionary = Tool_Consensus_2._load_rc_notation(Tool_Consensus_2.rc_json_path)
        rc_dictionary_R3DMA = Tool_Consensus_2._load_rc_notation(Tool_Consensus_2.rc_R3DMA_json_path)
        #rc_dictionary_FR = Tool_Consensus_2._load_rc_notation(Tool_Consensus_2.rc_FR_json_path)
        
        for table_name in os.listdir(Tool_Consensus_2.extended_tables_directory):
            cluster_df = pd.read_csv(os.path.join(
                Tool_Consensus_2.extended_tables_directory, table_name    
            ))
            
            #consensus_cluster_df = Tool_Consensus_2._get_consensus(
                #cluster_df, table_name, rc_dictionary, len(cluster_df), tie_instances, unresolved_tie_instances
            #)
            consensus_cluster_R3DMA_df = Tool_Consensus_2._get_consensus(
                cluster_df, table_name, rc_dictionary_R3DMA, len(cluster_df), tie_instances_without_R3DMA, unresolved_tie_instances_without_R3DMA
            )
            #consensus_cluster_FR_df = Tool_Consensus_2._get_consensus(
                #cluster_df, table_name, rc_dictionary_FR, len(cluster_df), tie_instances_without_FR, unresolved_tie_instances_without_FR
            #)
            
            #Tool_Consensus_2._output_df(consensus_cluster_df, consensus_cluster_R3DMA_df, consensus_cluster_FR_df, table_name)
            Tool_Consensus_2._output_df(consensus_cluster_R3DMA_df, table_name)
        df = pd.DataFrame(tie_instances)
        '''
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'tie_instances.csv'
        ), index=False)
        df = pd.DataFrame(unresolved_tie_instances)
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'unresolved_tie_instances.csv'
        ), index=False)
        '''
        
        df = pd.DataFrame(tie_instances_without_R3DMA)
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'tie_instances_without_R3DMA.csv'
        ), index=False)
        df = pd.DataFrame(unresolved_tie_instances_without_R3DMA)
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'unresolved_tie_instances_without_R3DMA.csv'
        ), index=False)
        
        '''
        df = pd.DataFrame(tie_instances_without_FR)
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'tie_instances_without_FR.csv'
        ), index=False)
        df = pd.DataFrame(unresolved_tie_instances_without_FR)
        df.to_csv(os.path.join(
            Tool_Consensus_2.consensus_directory, 'unresolved_tie_instances_without_FR.csv'
        ), index=False)
        '''
    
    # PRIVATE METHODS ---------------------------------------------------------
    @staticmethod
    def _load_rc_notation(path):
        with open(path, 'r') as f:
            return json.load(f)
    
    @staticmethod
    def _get_consensus(cluster_df, table_name, rc_dictionary, rows_in_df, tie_instances, unresolved_tie_instances):
        consensus_df = cluster_df.copy()
        for column_name in cluster_df.columns[::-1]:
            if '-' not in column_name:
                break
            consensus_df[column_name] = cluster_df.apply(
                lambda row: Tool_Consensus_2._get_mode(
                    row[column_name], table_name, rc_dictionary, column_name, row["PDB"], rows_in_df, tie_instances, unresolved_tie_instances
                ),
                axis=1
            )

        return consensus_df
    
    @staticmethod
    def _get_mode(cell, table_name, rc_dictionary, column_name, pdb, rows_in_df, tie_instances, unresolved_tie_instances):
        competing_contact_types = [Tool_Consensus_2._standardize(contact_type) for contact_type in cell.split(",")][1:]
        if not competing_contact_types:
            return None
        '''
        if 'REJECT' in competing_contact_types: # TODO why are there still rows with REJECT?
            #print("why is this here")
            return 'INCOMPATIBLE'
        '''

        counts = Counter(competing_contact_types) 
        mode, max_frequency = counts.most_common(1)[0]
        top_contact_types = [contact_type for contact_type, frequency in counts.items() if frequency == max_frequency]
        if len(top_contact_types) > 1:
            mode = Tool_Consensus_2._resolve_tie(top_contact_types, table_name, rc_dictionary, column_name, pdb, rows_in_df, tie_instances, unresolved_tie_instances)

        return mode
    
    @staticmethod
    def _output_df(consensus_cluster_R3DMA_df, table_name):
        '''
        consensus_cluster_df.to_csv(
            os.path.join(Tool_Consensus_2.mode_consensus_directory, 'AllTools', table_name), index=False
        )
        '''
        consensus_cluster_R3DMA_df.to_csv(
            os.path.join(Tool_Consensus_2.mode_consensus_directory, 'WithoutR3DMA', table_name), index=False
        )
        '''
        consensus_cluster_FR_df.to_csv(
            os.path.join(Tool_Consensus_2.mode_consensus_directory, 'WithoutFR', table_name), index=False
        )
        '''
        
    @staticmethod
    def _resolve_tie(top_contact_types, table_name, rc_dictionary, column_name, pdb, rows_in_df, tie_instances, unresolved_tie_instances):
        for combine_contact_types in [False, True]:
            mode = None
            max_count = None
            for checking_consensus_count in [True, False]:
                for contact_type in top_contact_types:
                    if 'INCOMPATIBLE' in contact_type or 'REJECT' in contact_type:
                        continue
                    rc, edge_based = Tool_Consensus_2._get_rc(combine_contact_types, rc_dictionary, table_name, column_name, contact_type)
                    counts = re.match(r"r(\d+)c(\d+)", rc)
                        
                    if checking_consensus_count:
                        consensus_count = int(counts.group(2))
                        if max_count is None or consensus_count > max_count:
                            max_count = consensus_count
                            mode = contact_type
                        elif consensus_count == max_count:
                            mode = None
                            break
                    else: # comparing report counts instead
                        report_count = int(counts.group(1))
                        if max_count is None or report_count > max_count:
                            max_count = report_count
                            mode = contact_type
                        elif report_count == max_count:
                            mode = None
                            break
                if mode == None:
                    max_count = None
                else:
                    if checking_consensus_count and not edge_based:
                        case = "contact type motif-wide consensus count"
                    elif checking_consensus_count and edge_based:
                        case = "edge-based motif-wide consensus count"
                    elif not edge_based:
                        case = "contact type report count"
                    else:
                        case = "edge-based report count"
                    tie_instances.append({
                        'cluster': table_name, 'PDB': pdb, 'column': column_name, 
                        'tied_contact_types': top_contact_types, 'consensus': mode,
                        'resolved_with': case, 'rows_in_cluster': rows_in_df
                    })
                    return mode
             
        # all cases failed, let's check whether the tie is only between nbp and something else
        non_nbp_count = 0
        non_nbp_contact_type = None
        for tied_contact_type in top_contact_types:
            if tied_contact_type != "nbp":
                non_nbp_count += 1
                non_nbp_contact_type = tied_contact_type
        if non_nbp_count < 2:
            tie_instances.append({
                'cluster': table_name, 'PDB': pdb, 'column': column_name,
                'tied_contact_types': top_contact_types, 'consensus': non_nbp_contact_type,
                'resolved_with': "prioritized non-nbp", 'rows_in_cluster': rows_in_df
            })
            return non_nbp_contact_type
                
                
        unresolved_tie_instances.append({
            'cluster': table_name, 'PDB': pdb, 'column': column_name,
            'tied_contact_types': top_contact_types, 'consensus': mode,
            'rows_in_cluster': rows_in_df
        })
        
        return "unresolved tie"
    
    @staticmethod
    def _get_rc(combine_contact_types, rc_dictionary, table_name, column_name, contact_type):
        if not combine_contact_types and 'nbp' not in contact_type:
            return (rc_dictionary[table_name][column_name][contact_type], False)
        if not combine_contact_types and 'nbp' in contact_type:
            return (Tool_Consensus_2._get_most_frequent_nbp_subtype(rc_dictionary, table_name, column_name), False)
        if 'nbp' in contact_type:
            return (Tool_Consensus_2._get_most_frequent_nbp_subtype(rc_dictionary, table_name, column_name), True)
        else:
            first_edge = contact_type[1]
            second_edge = contact_type[2]
            geometry = contact_type[0]
            counts_by_first_edge = [0,0] # initialized to r0c0
            counts_by_second_edge = [0,0] # initialized to r0c0
            for contact_type_key in rc_dictionary[table_name][column_name].keys():
                first_edge_match = re.match(fr"{geometry}{first_edge}.", contact_type_key)
                second_edge_match = re.match(fr"{geometry}.{second_edge}", contact_type_key)
                if first_edge_match is not None:
                    first_counts = re.match(r"r(\d+)c(\d+)", rc_dictionary[table_name][column_name][first_edge_match.group()])
                    counts_by_first_edge[0] += int(first_counts.group(1))
                    counts_by_first_edge[1] += int(first_counts.group(2))
                if second_edge_match is not None:
                    second_counts = re.match(r"r(\d+)c(\d+)", rc_dictionary[table_name][column_name][second_edge_match.group()])
                    counts_by_second_edge[0] += int(second_counts.group(1))
                    counts_by_second_edge[1] += int(second_counts.group(2))
            # For now only return the one with the highest consensus count, but using the second one later could help resolve ties
            if counts_by_first_edge[1] > counts_by_second_edge[1]: # TODO could be a slight bias towards second edge since it happens when equal or greater than
                return (f"r{counts_by_first_edge[0]}c{counts_by_first_edge[1]}", True)
            elif counts_by_first_edge[1] < counts_by_second_edge[1]:
                return (f"r{counts_by_second_edge[0]}c{counts_by_second_edge[1]}", True)
            if counts_by_first_edge[0] > counts_by_second_edge[0]:
                return (f"r{counts_by_first_edge[0]}c{counts_by_first_edge[1]}", True)
            elif counts_by_first_edge[0] < counts_by_second_edge[0]:
                return (f"r{counts_by_second_edge[0]}c{counts_by_second_edge[1]}", True)
            return (f"r{counts_by_second_edge[0]}c{counts_by_second_edge[1]}", True) # if equal
        
    @staticmethod
    def _standardize(contact_type):
        if contact_type.startswith('nc') or contact_type.startswith('nt'):
            return contact_type[1] + contact_type[2:].upper()
        elif contact_type.startswith('c') or contact_type.startswith('t'):
            return contact_type[0] + contact_type[1:].upper()
        else:
            return contact_type
        
    @staticmethod
    def _get_most_frequent_nbp_subtype(rc_dictionary, table_name, column_name):
        max_report_count = 0
        max_motif_wide_consensus_count = 0
        for contact_type_key, rc_counts in rc_dictionary[table_name][column_name].items():
            if '!' not in contact_type_key: #disregard non nbp subtypes
                continue
            rc_counts = re.match(r"r(\d+)c(\d+)", rc_counts)
            if int(rc_counts.group(1)) > max_report_count:
                max_report_count = int(rc_counts.group(1))
            if int(rc_counts.group(2)) > max_motif_wide_consensus_count:
                max_motif_wide_consensus_count = int(rc_counts.group(2))
        return f"r{max_report_count}c{max_motif_wide_consensus_count}"