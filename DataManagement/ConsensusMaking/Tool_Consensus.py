# AUTHOR: Kristopher Church

import os
import shutil
import pandas as pd
from collections import Counter

# currently uses a basic mode-based conensus algorithm.
# Currently will consider 'near' contact types as being normal. IE: ncSH will
# merely be considered as cSH.
class Tool_Consensus:
    TABLE_DIRECTORY = os.path.join('Data', 'Preprocessed', 'Extended_Tables')
    OUTPUT_DIRECTORY = os.path.join(
        'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based_Profiles'
    )
    AMBIGUITIES_FILE_PATH = os.path.join(
        'Data', 'Consensus', 'Tool_Consensus', 'ambiguities.txt'
    )
    
    table = pd.DataFrame
    ambiguities = []
    
    @staticmethod
    def run():
        # Reset the contents of Mode_Based_Profiles folder to be blank
        if os.path.exists(Tool_Consensus.OUTPUT_DIRECTORY):
            shutil.rmtree(Tool_Consensus.OUTPUT_DIRECTORY)
        os.makedirs(Tool_Consensus.OUTPUT_DIRECTORY)
        
        # Access each extended table
        for csv_file_name in os.listdir(Tool_Consensus.TABLE_DIRECTORY):
            Tool_Consensus.table = pd.read_csv(
                os.path.join(Tool_Consensus.TABLE_DIRECTORY, csv_file_name)
            )
            
            Tool_Consensus._find_consensus(csv_file_name)
            Tool_Consensus._export_table(csv_file_name)
            
        # Write a file with the ambiguities
        with open(Tool_Consensus.AMBIGUITIES_FILE_PATH, 'w') as file:
            for line in Tool_Consensus.ambiguities:
                file.write(line)
    
    @staticmethod         
    def _find_consensus(csv_file_name):
        for row_label, _ in Tool_Consensus.table.iterrows():
            
            for column_name in Tool_Consensus.table.columns[::-1]:
                if '-' in column_name:
                    entry = Tool_Consensus.table.loc[row_label, column_name]
                    entry = Tool_Consensus._condense(
                        entry, csv_file_name, row_label, column_name
                    )
                    Tool_Consensus.table.loc[row_label, column_name] = entry
                else:
                    break
    
    # TODO: still in progress, comments might need to be streamlined as they 
    # reflect different intents based on different moments during development
    @staticmethod
    def _condense(entry, csv_file_name, row_label, column_name):
        entry = entry.split(',')
        
        # First get the consensus if 'near' types are left as is
        # This is because we want to quantify the impact of treating 'near'
        # types as normal, specifically by seeing when it would have changed
        # the consensus contact type to a completely different one. IE: find
        # cases where originally the consensus might be cSS, but by considering
        # 'near' contact types as being normal, it makes the consensus cSH now.
        # Call these case 'ambiguities'
        counts = Counter(entry)
        if len(counts) == 1:
            contact_type, _ = counts.most_common(1)[0]
            weight = 1.00
        else:
            top1, top2 = counts.most_common(2)
            if top1[1] == top2[1]: # If there was a tie
                contact_type = '?'
                weight = 0.0
            else:
                contact_type, amount = top1
                weight = amount / len(entry)
        if contact_type.startswith('n') and contact_type != 'nbp':
            contact_type = contact_type[1:]
        original_consensus = (contact_type, weight)
        
        # Now actually determine the consensus we will be using for mode-based
        entry = [ # convert 'near' types first
            s[1:] if s.startswith('n') and s != 'nbp' else s for s in entry
        ]
        counts = Counter(entry)
        if len(counts) == 1:
            contact_Type, _ = counts.most_common(1)[0]
            weight = 1.00
        else:
            top1, top2 = counts.most_common(2)
            if top1[1] == top2[1]: #  If there was a tie
                contact_type = '?'
                weight = 0.0
            else:
                contact_type, amount = top1
                weight = amount / len(entry)
        
        # Document the ambiguities
        if contact_type == '?' and original_consensus[0] != '?':
            Tool_Consensus.ambiguities.append(
                f'CASE 1 -- {csv_file_name} at row {row_label} and column ' + 
                f'{column_name}: the ORIGINAL contact type ' + 
                f'{original_consensus[0]} is probably more accurate as the ' + 
                'tool consensus.\n'
            )
            return f'{original_consensus[0]} {original_consensus[1]}'
        elif contact_type != '?' and original_consensus[0] == '?':
            Tool_Consensus.ambiguities.append(
                f'CASE 2 -- {csv_file_name} at row {row_label} and column ' + 
                f'{column_name}: the NEW contact type {contact_type} ' +
                'is probably more accurate as the tool consensus.\n'
            )
        elif contact_type != '?': # and original_consensus[0] != '?' also
            Tool_Consensus.ambiguities.append(
                f'CASE 3 --{csv_file_name} at row {row_label} and column ' + 
                f'{column_name}: inspect these more carefully. Since \'near\' ' +
                'means that contact type is not as certain, it is important ' +
                'to know to what extent it adds certainty towards the final ' +
                'consensus. it may or may not warrant the NEW contact type.\n'
            )
        
        return f'{contact_type} {weight:.2f}'
    
    @staticmethod
    def _export_table(csv_file_name):
        file_path = os.path.join(Tool_Consensus.OUTPUT_DIRECTORY, csv_file_name)
        Tool_Consensus.table.to_csv(file_path, index=False)
        
        