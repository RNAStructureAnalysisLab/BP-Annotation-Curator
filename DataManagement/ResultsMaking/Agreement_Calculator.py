# AUTHOR: Kristopher Church

import os
import shutil
import pandas as pd
from collections import Counter

class Agreement_Calculator:
    CONSENSUS_DIRECTORY = os.path.join( # TODO use cluster consensus instead when finished
        'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based_Profiles'
    )
    EXTENDED_TABLE_DIRECTORY = os.path.join(
        'Data', 'Preprocessed', 'Extended_Tables'
    )
    OUTPUT_DIRECTORY = os.path.join('Data', 'ResultsAnalysis')
    TOOLS = ['R3DMA', 'CL', 'FR', 'MC', 'RV', 'DSSR']
    
    consensus = pd.DataFrame
    extended_table = pd.DataFrame
        
    # table with each base pair in the extended table, shows whether it was 
    # also reported by an annotation tool, and whether it's contact type 
    # matched with that annotation tool.
    intersections = []
    
    @staticmethod
    def run():
        # Reset the contents of OUTPUT_DIRECTORY folder to be blank
        if os.path.exists(Agreement_Calculator.OUTPUT_DIRECTORY):
            shutil.rmtree(Agreement_Calculator.OUTPUT_DIRECTORY)
        os.makedirs(Agreement_Calculator.OUTPUT_DIRECTORY)
        
        for file_name in os.listdir(
                Agreement_Calculator.CONSENSUS_DIRECTORY
        ):
            Agreement_Calculator.consensus = pd.read_csv(
                os.path.join(
                    Agreement_Calculator.CONSENSUS_DIRECTORY, file_name
                )
            )
            Agreement_Calculator.extended_table = pd.read_csv(
                os.path.join(
                    Agreement_Calculator.EXTENDED_TABLE_DIRECTORY, file_name
                )
            )
            
            Agreement_Calculator._calculate(os.path.splitext(file_name)[0])
        
        Agreement_Calculator._export_table()
            
    @staticmethod
    def _calculate(cluster_name):
        # Note that consensus and extended_table have same dimensions and labels
        for row_label in Agreement_Calculator.consensus.index:
            pdb = Agreement_Calculator.consensus.loc[row_label, 'PDB']
            chains = Agreement_Calculator.consensus.loc[row_label, 'Chain(s)']
            for column_name in reversed(
                    Agreement_Calculator.consensus.columns
            ):
                if '-' in column_name:
                    consensus, consensus_count_cluster, consensus_proportion, consensus_count_row = Agreement_Calculator._get_consensus_information(row_label, column_name)
                    residue_ids, nucleotides, contact_type, all_annotations = (
                        Agreement_Calculator._get_base_pair(
                            row_label, column_name, consensus
                        )
                    )
                    if contact_type == None:
                        continue
                    presence_list = (
                        Agreement_Calculator._tools_presence(
                            row_label, column_name, contact_type
                        )
                    )
                    
                    # Add to intersections
                    Agreement_Calculator.intersections.append({
                        'Cluster': cluster_name, 'PDB': pdb, 'Chain(s)': chains,
                        'Residue IDs': residue_ids, 'Nucleotides': nucleotides,
                        'Consensus Contact Type': contact_type, 
                        'Agree BP (R3DMA)': presence_list[0][0],
                        'Agree BP (CL)': presence_list[1][0],
                        'Agree BP (FR)': presence_list[2][0],
                        'Agree BP (MC)': presence_list[3][0],
                        'Agree BP (RV)': presence_list[4][0],
                        'Agree BP (DSSR)': presence_list[5][0],
                        'Agree BP Count': Agreement_Calculator._sum(
                            presence_list, 0
                        ),
                        'Matches (R3DMA)': presence_list[0][1],
                        'Matches (CL)': presence_list[1][1],
                        'Matches (FR)': presence_list[2][1],
                        'Matches (MC)': presence_list[3][1],
                        'Matches (RV)': presence_list[4][1],
                        'Matches (DSSR)': presence_list[5][1],
                        'Matches Count': Agreement_Calculator._sum(
                            presence_list, 1
                        ),
                        'Cluster Consensus': consensus,
                        'Cluster Consensus Count': consensus_count_cluster,
                        'Cluster Consensus Proportion': consensus_proportion,
                        'Cluster Consensus Count in that BP': consensus_count_row,
                        'All Annotations': all_annotations,
                        'Comments': None
                    })
                    
                else:
                    break
        
    @staticmethod
    def _export_table():
        df = pd.DataFrame(Agreement_Calculator.intersections)
        output_path = os.path.join(
            Agreement_Calculator.OUTPUT_DIRECTORY, 'base_pairing_agreements.csv'
        )
        df.to_csv(output_path, index=False)
    
    @staticmethod
    def _get_base_pair(row_label, column_name, consensus):
        '''
        contact_type = Agreement_Calculator.consensus.loc[
            row_label, column_name
        ].split(' ')[0]
        '''
        all_annotations = Agreement_Calculator.extended_table.loc[
            row_label, column_name
        ].split(' ')[0]
        if all_annotations == 'INCOMPATIBLE':
            return (None, None, None, None)
        contact_type = Agreement_Calculator._resolve_contact_type(
            all_annotations, consensus
        )
        column_name1, column_name2 = column_name.split('-')
        columns = list(Agreement_Calculator.consensus.columns)
        
        # Get residue ids
        residue1 = Agreement_Calculator.consensus.loc[
            row_label, column_name1
        ]
        residue2 = Agreement_Calculator.consensus.loc[
            row_label, column_name2
        ]
        residue_ids = (residue1, residue2)
        
        # Get nucleotides
        column_position = columns.index(column_name1)
        nucleotide_column_label = columns[column_position - 1]
        nucleotide1 = Agreement_Calculator.consensus.loc[
            row_label, nucleotide_column_label
        ]
        column_position = columns.index(column_name2)
        nucleotide_column_label = columns[column_position - 1]
        nucleotide2 = Agreement_Calculator.consensus.loc[
            row_label, nucleotide_column_label
        ]
        nucleotides = (nucleotide1, nucleotide2)
        
        return (residue_ids, nucleotides, contact_type, all_annotations)
        

    @staticmethod
    def _tools_presence(row_label, column_name, contact_type):
        presence_list = [] # List of boolean pairs
        contact_types = Agreement_Calculator.extended_table.loc[
            row_label, column_name
        ].split(',')
        for i in range(len(Agreement_Calculator.TOOLS)):
            presence_list.append(
                Agreement_Calculator._presence(
                    contact_types[i], contact_type
                )
            )
            
        return presence_list
    
    # OUTPUT: A boolean pair. First element describes whether the tool and
    # extended tanble agree that there is a base pairing interaction  or not 
    # for that pair of residues. Second element describes whether they agree on
    # what that base pairing interaction is if there is one. 
    # Recall 'nbp' stands for 'Not a Base Pair'
    @staticmethod
    def _presence(tool_contact_type, consensus_contact_type):
        # Instances where consensus believes there is NO base pairing
        # interaction
        if tool_contact_type == 'nbp' and consensus_contact_type == 'nbp': 
            return (True, True) # vacously they agree on bp interaction since there isn't one
        elif tool_contact_type == 'nbp':
            return (False, False)
        
        # Instances where consensus believes there IS a base pairing
        # interaction
        elif tool_contact_type == consensus_contact_type:
            return (True, True) # Both agree is bp, and on what
        elif consensus_contact_type != 'nbp':
            return (True, False) # Both agree is bp, but disagree on what
        else: # consensus_contact_type == 'nbp'
            return (False, False) # Both disagree on whether residues are bp
        # !!!
        # TODO would like to do deeper analysis on instances where consensus is
        # nbp, and instances where it is not. This will tell us whether a tool
        # might be over reporting or under reporting base pair interactions
        # compared to the othe tools
        
    @staticmethod
    def _sum(presence_list, index):
        boolean_sum = 0
        
        for pair in presence_list:
            first, second = pair
            if index == 0:
                boolean_sum += int(first)
            else:
                boolean_sum += int(second)
        
        return boolean_sum
    
    # TODO: this is fine, but redesign the invoking function so that it iterate
    # through the entire column at once
    @staticmethod
    def _get_consensus_information(row_label, column_name):
        counts = {} # counts of the contact type in the entire column
        column = Agreement_Calculator.extended_table[column_name]
        most_frequent_contact_type = '?'
        max_count = 0
        total_count = 0
        
        for contact_types in column:
            contact_types = contact_types.strip().split(',')
            for contact_type in contact_types:
                if contact_type == 'REJECT':
                    continue
                total_count += 1
                
                # add to total counts
                if contact_type not in counts:
                    counts[contact_type] = 1
                else:
                    counts[contact_type] = counts[contact_type] + 1
                
                current_count = counts[contact_type]
                if current_count > max_count:
                    most_frequent_contact_type = contact_type
                    max_count = current_count
                elif current_count == max_count:
                    most_frequent_contact_type = '?'
        
        row_values = column.loc[row_label]
        row_values = row_values.strip().split(',')
        row_count = 0
        for contact_type in row_values:
            if contact_type == most_frequent_contact_type:
                row_count += 1
        
        return (most_frequent_contact_type, max_count, max_count/total_count, row_count)
   
    #==========================================================================
    # RESOLVE CONTACT TYPE Section Below:
    #==========================================================================
    
    @staticmethod
    def _resolve_contact_type(all_annotations, cluster_consensus):
        all_annotations = all_annotations.split(',')
        Agreement_Calculator._ensure_lw(all_annotations)
        most_reported = Agreement_Calculator._most_reported_contact_types(all_annotations)
        if 'REJECT' in most_reported:
            most_reported.remove('REJECT')
        contact_type = Agreement_Calculator._tool_consensus(most_reported, cluster_consensus)
        if isinstance(contact_type, list):
            contact_type = contact_type[0]
        return contact_type
    
    @staticmethod
    def _ensure_lw(all_annotations):
        for i in range(len(all_annotations)):
            annotation = all_annotations[i]
            if 'n' in annotation and annotation != 'nbp':
                annotation = annotation[1] + annotation[2:].upper()
            elif annotation != 'nbp':
                annotation = annotation[0] + annotation[1:].upper()
            all_annotations[i] = annotation
    
    @staticmethod
    def _most_reported_contact_types(all_annotations):
        counts = Counter(all_annotations)
        max_count = max(counts.values())
        return [item for item, count in counts.items() if count == max_count]
    
    @staticmethod
    def _tool_consensus(most_reported, cluster_consensus):
        amount_tied = len(most_reported)
        if amount_tied == 2:
            contact_type = Agreement_Calculator._get_tie_from_two(most_reported, cluster_consensus)
        elif amount_tied > 2:
            contact_type = Agreement_Calculator._cluster_based_tie_breaker(most_reported, cluster_consensus)
        else: # amount_tied == 1
            contact_type = most_reported
        return contact_type
    
    @staticmethod
    def _get_tie_from_two(most_reported, cluster_consensus):
        first = most_reported[0]
        second = most_reported[1]
        result = 'nbp'
        if first != 'nbp' and second != 'nbp':
            result = Agreement_Calculator._cluster_based_tie_breaker(most_reported, cluster_consensus)
        else:
            if first != 'nbp':
                result = first
            if second != 'nbp':
                result = second
        
        return result
    
    @staticmethod
    def _cluster_based_tie_breaker(most_reported, cluster_consensus):
        if cluster_consensus in most_reported:
            return cluster_consensus
        else:
            return '?'
    
    #==========================================================================
                