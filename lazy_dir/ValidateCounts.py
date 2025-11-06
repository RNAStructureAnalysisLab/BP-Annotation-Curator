import os
import re
import json
import pandas as pd

class ValidateCount:
    RAW_DATA_DIRECTORY = os.path.join(
        '..', 'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    PREPROCESSED_DATA_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'JSON_Annotations'
    )
    R3DMA_DIRECTORY = os.path.join(
        '..', 'Data', 'Raw', 'R3DMA', 'cluster_tables_3.95'
    )
    EXTENDED_TABLE_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'Extended_Tables'    
    )
    CONSENSUS_DIRECTORY = os.path.join(
        '..', 'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based', 'WithoutR3DMA'    
    )
    I_CODE_RESIDUES_PATH = os.path.join(
        '..', 'DataManagement', 'DataPreparation', 'residue_ids.txt'
    )
    OUTPUT_DIRECTORY = '.'
    
    @staticmethod
    def run():
        
        residue_counts = ValidateCount._configure_counts()
        ValidateCount._raw(residue_counts)
        ValidateCount._r3dma(residue_counts)
        ValidateCount._print(residue_counts, "RAW BP COUNTS")
        
        residue_counts = ValidateCount._configure_counts()
        ValidateCount._preprocessed(residue_counts)
        ValidateCount._print(residue_counts, "PREPROCESSED BP COUNTS")
        
        i_code_residues = ValidateCount._load_i_code_residues()
        ValidateCount._extended_tables(i_code_residues)
        ValidateCount._consensus_table(i_code_residues)
        
        
        '''
        raw_unique_residues = {
            'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 'DSSR': set()
        }
        preprocessed_unique_residues = {
            'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 'DSSR': set()
        }
        only_in_raw = {
            'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 'DSSR': set()
        }
        only_in_preprocessed = {
            'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 'DSSR': set()
        }
        
        ValidateCount._raw(raw_unique_residues)
        ValidateCount._preprocessed(preprocessed_unique_residues)
        ValidateCount._get_differences(only_in_raw, only_in_preprocessed, raw_unique_residues, preprocessed_unique_residues)
        ValidateCount._analyze(only_in_raw, only_in_preprocessed, raw_unique_residues, preprocessed_unique_residues)
        ValidateCount._write_differences(only_in_raw, only_in_preprocessed)
        '''
        
    @staticmethod
    def _raw(raw_unique_residues):
        for pdb_tool_csv in os.listdir(ValidateCount.RAW_DATA_DIRECTORY):
            pdb, tool = os.path.splitext(pdb_tool_csv)[0].split('_')
            if tool == 'MO':
                continue
            
            try:
                df = pd.read_csv(os.path.join(ValidateCount.RAW_DATA_DIRECTORY, pdb_tool_csv), dtype={
                    'residue1': str, 'residue2': str})
            except pd.errors.EmptyDataError:
                #TODO justify why we have empty annotations
                continue
            for i in range(len(df)):
                residue1 = df.loc[i, 'residue1']
                residue2 = df.loc[i, 'residue2']
                try:
                    residue1, residue2 = ValidateCount._sort_residues(
                        residue1, residue2
                    )
                except TypeError: # no chain in the residues probably
                    continue
                raw_unique_residues[tool].add(f'{pdb}_{residue1}_{residue2}')
    
    @staticmethod
    def _preprocessed(preprocessed_unique_residues):
        for tool_json in os.listdir(ValidateCount.PREPROCESSED_DATA_DIRECTORY):
            if 'rejected' in tool_json:
                continue
            tool = os.path.splitext(tool_json)[0]
            file_path = os.path.join(
                ValidateCount.PREPROCESSED_DATA_DIRECTORY, tool_json
            )
            with open(file_path, 'r') as json_file:
                data = json.load(json_file)
                #sum_descriptions = 0
                for pdb, residue_pairs in data.items():
                    for first_residue, pair_information in residue_pairs.items():
                        first_residue = first_residue[:-1]
                        for second_residue, descriptions in pair_information.items():
                            #sum_descriptions += len(descriptions)
                            second_residue = second_residue[:-1]
                            first_residue_new, second_residue_new = ValidateCount._sort_residues(first_residue, second_residue) # The sorting method I used for preprocessed is a little different
                            preprocessed_unique_residues[tool].add(f'{pdb}_{first_residue_new}_{second_residue_new}')
                #print(sum_descriptions)
        
    @staticmethod
    def _sort_residues(residue1: str, residue2: str) -> None:
        try:
            pattern = r'(\d*[A-Za-z]+|\d)(-?\d+)'
            r1_matches = re.match(pattern, residue1)
            r2_matches = re.match(pattern, residue2)
            r1_chain = r1_matches.group(1)
            r1_res_id = r1_matches.group(2)
            r2_chain = r2_matches.group(1)
            r2_res_id = r2_matches.group(2)
        except AttributeError:
            input(f'{residue1} {residue2}')
                        
        if ValidateCount._chain_less_than(r1_chain, r2_chain):
            return residue1, residue2
        elif ValidateCount._chain_less_than(r2_chain, r1_chain):
            return residue2, residue1
        elif int(r1_res_id) < int(r2_res_id):
            return residue1, residue2
        else:
            return residue2, residue1
        
    @staticmethod
    def _chain_less_than(chain1: str, chain2: str) -> bool:
        c1_digit, c1_letter = ValidateCount._parse_chain(chain1)
        c2_digit, c2_letter = ValidateCount._parse_chain(chain2)
        
        if c1_digit != '' and c2_digit != '':
            if c1_digit == c2_digit:
                if c1_letter != '' and c2_letter != '':
                    return ValidateCount._letter_comparison(c1_letter, c2_letter)
                return c2_letter != ''
            return int(c1_digit) < int(c2_digit)
        if c1_digit != '':
            return False
        if c2_digit != '':
            return True
        return ValidateCount._letter_comparison(c1_letter, c2_letter)
    
    @staticmethod
    def _parse_chain(chain):
        pattern = r'(\d*)([A-Za-z]*)'
        matches = re.match(pattern, chain)
        
        if len(matches.groups()) == 2:
            return (matches.group(1), matches.group(2))
        if matches.group(1).isdigit():
            return (matches.group(1), '')
        
        return ('', matches.group(1))
    
    @staticmethod
    def _letter_comparison(chain1, chain2):
        # If chain was a series of letters
        i = 0
        while i < len(chain1) and i < len(chain2):
            if ord(chain1[i]) < ord(chain2[i]):
                return True
            elif ord(chain1[i]) > ord(chain2[i]):
                return False
            i += 1
        # If tie, consider the shorter one smaller
        return len(chain1) < len(chain2)
    
    @staticmethod
    def _get_differences(only_in_raw, only_in_preprocessed, raw_unique_residues, preprocessed_unique_residues):
        for tool in ['CL', 'FR', 'MC', 'RV', 'DSSR']:
            only_in_raw[tool] = raw_unique_residues[tool].difference(preprocessed_unique_residues[tool])
            only_in_preprocessed[tool] = preprocessed_unique_residues[tool].difference(raw_unique_residues[tool])
            
    @staticmethod
    def _analyze(only_in_raw, only_in_preprocessed, raw_unique_residues, preprocessed_unique_residues):
        for tool in ['CL', 'FR', 'MC', 'RV', 'DSSR']:
            print(f'raw {tool} {len(raw_unique_residues[tool])}')
            print(f'preprocessed {tool} {len(preprocessed_unique_residues[tool])}')
            print()
        '''
        print(len(raw_unique_residues['CL']))
        print(len(preprocessed_unique_residues['CL']))
        print(len(only_in_raw['CL']))
        print(len(only_in_preprocessed['CL']))
        '''
        
        '''
        print("==============================================================")
        print(only_in_raw['CL'])
        print('\n\n\n')
        print("==============================================================")
        print(only_in_preprocessed['CL'])
        '''
        
    @staticmethod
    def _write_differences(only_in_raw, only_in_preprocessed):
        for tool in ['FR']: #'FR', 'MC', 'RV', 'DSSR']:
            path = os.path.join(ValidateCount.OUTPUT_DIRECTORY, f'r-p_{tool}.txt')

            only_raw_sorted = sorted(only_in_raw[tool])
            only_pre_sorted = sorted(only_in_preprocessed[tool])

            with open(path, 'w', encoding='utf-8') as f:
                for item in only_raw_sorted:
                    f.write(item + '\n')
            
            path = os.path.join(ValidateCount.OUTPUT_DIRECTORY, f'p-r_{tool}.txt')
            with open(path, 'w', encoding='utf-8') as f:
                for item in only_pre_sorted:
                    f.write(item + '\n')
     
    ###########################################################################
    # VERSION 2 FUNCTIONS BELOW
    ###########################################################################               
     
    @staticmethod
    def _configure_counts() -> dict[set[str]]:
        return {
            'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 'DSSR': set(), 'R3DMA': 0
        }
    
    @staticmethod
    def _print(residue_counts: dict[set[str]], step: str) -> None:
        count_sum = 0
        print(f'{'='*5} {step} {'='*5}')
        for tool, item in residue_counts.items():
            if isinstance(item, int):
                print(f'{tool}: {item}')
                count_sum += item
            else:
                print(f'{tool}: {len(item)}')
                count_sum += len(item)
            
        print(f'TOTAL: {count_sum}\n\n\n')
        
    @staticmethod
    def _r3dma(residue_counts: dict[set[str]]) -> None:
        for cluster_table in os.listdir(ValidateCount.R3DMA_DIRECTORY):
            df = pd.read_csv(os.path.join(
                ValidateCount.R3DMA_DIRECTORY, cluster_table
            ))
            for column in df.columns[::-1]:
                if '-' not in column:
                    break
                residue_counts['R3DMA'] += len(df[column])
     
    @staticmethod 
    def _count_bp(bp_counts, contact_types, first_res_id, second_res_id, chain, pdb_id, i_code_residues):
        tool = {0: 'R3DMA', 1: 'CL', 2: 'FR', 3: 'MC', 4: 'RV', 5: 'DSSR'}
        pattern = r'(\d*)' # try not to grab i_code ar anything else by accident
        first_res_sequence = re.match(pattern, str(first_res_id))
        second_res_sequence = re.match(pattern, str(second_res_id))
        
        chain = chain.split('+')
        if len(chain) == 2:
            first_chain, second_chain = chain
        elif len(chain):
            first_chain, second_chain = (chain[0], chain[0])
        else:
            return # TODO is it oaky to not count entries where there are more than 2 chains?
            
            
        formatted_pair = f'{pdb_id}_{first_chain}{first_res_sequence.group(1)}_{second_chain}{second_res_sequence.group(1)}'

        if formatted_pair in i_code_residues:
            return

        contact_types = contact_types.split(',')
        for i in range(len(contact_types)):
            if contact_types[i] != 'REJECT':
                bp_counts[tool[i]] += 1
                
    @staticmethod 
    def _count_consensus_bp(bp_counts, contact_type, first_res_id, second_res_id, chain, pdb_id, i_code_residues) -> int:
        pattern = r'(\d*)' # try not to grab i_code ar anything else by accident
        first_res_sequence = re.match(pattern, str(first_res_id))
        second_res_sequence = re.match(pattern, str(second_res_id))
        
        chain = chain.split('+')
        if len(chain) == 2:
            first_chain, second_chain = chain
        elif len(chain):
            first_chain, second_chain = (chain[0], chain[0])
        else:
            return bp_counts # TODO is it oaky to not count entries where there are more than 2 chains?
            
            
        formatted_pair = f'{pdb_id}_{first_chain}{first_res_sequence.group(1)}_{second_chain}{second_res_sequence.group(1)}'
        reversed_formatted_pair = f'{pdb_id}_{second_chain}{second_res_sequence.group(1)}_{first_chain}{first_res_sequence.group(1)}'

        if formatted_pair in i_code_residues or reversed_formatted_pair in i_code_residues:
            return bp_counts
        
        if contact_type != 'REJECT':
            return bp_counts + 1
        return bp_counts
                
     
    @staticmethod
    def _extended_tables(i_code_residues: set[str]) -> None:
        bp_counts = {'CL': 0, 'FR': 0, 'MC': 0, 'RV': 0, 'DSSR': 0, 'R3DMA': 0}
        for extended_table in os.listdir(ValidateCount.EXTENDED_TABLE_DIRECTORY):
            df = pd.read_csv(os.path.join(ValidateCount.EXTENDED_TABLE_DIRECTORY, extended_table))
            ValidateCount._get_bp_counts_table(df, bp_counts, i_code_residues)

        ValidateCount._print_bp_counts(bp_counts, 'EXTENDED TABLES')
        
    @staticmethod
    def _get_bp_counts_table(df: pd.DataFrame, bp_counts: dict[int], i_code_residues: set[str]) -> None:
        for column in df.columns[::-1]:
            if '-' not in column:
                break
            first_res, second_res = column.split('-')
            for contact_types, first_res_id, second_res_id, chain, pdb_id in zip(df[column], df[first_res], df[second_res], df['Chain(s)'], df['PDB']):
                ValidateCount._count_bp(bp_counts, contact_types, first_res_id, second_res_id, chain, pdb_id, i_code_residues)
                
    @staticmethod
    def _get_bp_counts_consensus_table(df: pd.DataFrame, bp_counts: int, i_code_residues: set[str]) -> int:
        for column in df.columns[::-1]:
            if '-' not in column:
                break
            first_res, second_res = column.split('-')
            for contact_type, first_res_id, second_res_id, chain, pdb_id in zip(df[column], df[first_res], df[second_res], df['Chain(s)'], df['PDB']):
                bp_counts = ValidateCount._count_consensus_bp(bp_counts, contact_type, first_res_id, second_res_id, chain, pdb_id, i_code_residues)
        return bp_counts
        
    @staticmethod
    def _consensus_table(i_code_residues: set[str]) -> None:
        bp_counts = 0
        for consensus_table in os.listdir(ValidateCount.CONSENSUS_DIRECTORY):
            df = pd.read_csv(os.path.join(ValidateCount.CONSENSUS_DIRECTORY, consensus_table))
            bp_counts = ValidateCount._get_bp_counts_consensus_table(df, bp_counts, i_code_residues)
        
        print(f'{'='*5} CONSENSUS TABLES {'='*5}')
        print(f'TOTAL: {bp_counts}')
    
    @staticmethod
    def _print_bp_counts(bp_counts: dict[int], step: str) -> None:
        count_sum = 0
        print(f'{'='*5} {step} {'='*5}')
        for tool, count in bp_counts.items():
            count_sum += count
            print(f'{tool}: {count}')
            
        print(f'TOTAL: {count_sum}\n\n\n')
        
    @staticmethod
    def _load_i_code_residues() -> set[str]:
        i_code_residues = set()
        with open(ValidateCount.I_CODE_RESIDUES_PATH, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    i_code_residues.add(line)
        
        return i_code_residues
    
ValidateCount.run()