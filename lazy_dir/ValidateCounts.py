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
    OUTPUT_DIRECTORY = '.'
    
    @staticmethod
    def run():
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

    
ValidateCount.run()