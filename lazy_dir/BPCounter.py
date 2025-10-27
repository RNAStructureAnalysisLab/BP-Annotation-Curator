import os
import re
import json
import pandas as pd

#TODO
# in the raw data, there are some csv files that are empty. Why?
# in the raw data, there are CSVs like 4V9F_RV.csv where there seems to be no
    # chain.
# Each tool did not annotate the same number of files. 4V9F only exists for RV 
    # and DSSR (TURNS OUT 4V9F was originally PDBx)
    # The PDB file has chains '0' and '9' yet it didn't get annotated, could it
    # be because it didn't have headers?
        # Also I don't think I saw those chains in the PDBx so I want to double check

class BPCounter:
    ###########################################################################
    # Config
    RAW_DATA_DIRECTORY = os.path.join(
        '..', 'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    PREPROCESSED_DATA_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'JSON_Annotations'
    )
    EXTENDED_COUNTS_DATA = os.path.join(
        '..', 'Data', 'Preprocessed', 'Extended_Tables'    
    )
    
    all_annotations = {
        'CL': {}, 'FR': {}, 'MC': {}, 'RV': {}, 'DSSR': {}
    }
    tools = {1: 'CL', 2: 'FR', 3: 'MC', 4: 'RV', 5: 'DSSR'}
    
    @staticmethod
    def run() -> None:
        BPCounter._initial_counts()
        BPCounter._reset_counts()
        BPCounter._preprocessor_counts()
        BPCounter._reset_counts()
        BPCounter._extended_counts()
        BPCounter._reset_counts()
    
    # =========================================================================
    # Pipeline Methods Below
    # =========================================================================
    
    @staticmethod
    def _initial_counts() -> None:
        for file in os.listdir(BPCounter.RAW_DATA_DIRECTORY):
            pdb = file.split('_')[0]
            tool = file.split('_')[1].split('.')[0]
            if tool not in BPCounter.all_annotations.keys():
                continue
            file_path = os.path.join(BPCounter.RAW_DATA_DIRECTORY, file)
            try:
                df = pd.read_csv(file_path, dtype={
                    'residue1': str, 'residue2': str}
                )
            except pd.errors.EmptyDataError:
                #TODO justify why we have empty annotations
                continue
            for i in range(len(df)):
                residue1 = df.loc[i, 'residue1']
                residue2 = df.loc[i, 'residue2']
                try:
                    residue1, residue2 = BPCounter._sort_residues(
                        residue1, residue2
                    )
                except TypeError: # no chain in the residues probably
                    continue
                if pdb not in BPCounter.all_annotations[tool]:
                    BPCounter.all_annotations[tool][pdb] = set()
                BPCounter.all_annotations[tool][pdb].add(f'{residue1}{residue2}')
        
        BPCounter._print("Residue Pair Counts in Raw Data")
        
    @staticmethod
    def _preprocessor_counts() -> None:
        for file in os.listdir(BPCounter.PREPROCESSED_DATA_DIRECTORY):
            if 'rejected' in file:
                continue
            tool = file.split('.')[0]
            file_path = os.path.join(
                BPCounter.PREPROCESSED_DATA_DIRECTORY, file
            )
            with open(file_path, 'r') as json_file:
                data = json.load(json_file)
                for _, residue_pairs in data.items():
                    for second_residue in residue_pairs.values():
                        BPCounter.all_annotations[tool] += len(second_residue)
        
        BPCounter._print("Residue Pair Counts in Preprocessed Data")
        
    @staticmethod
    def _extended_counts() -> None:
        for file in os.listdir(BPCounter.EXTENDED_COUNTS_DATA):
            file_path = os.path.join(
                BPCounter.EXTENDED_COUNTS_DATA, file
            )
            extended_table_df = pd.read_csv(file_path)
            for column in extended_table_df.columns[::-1]:
                if '-' not in column:
                    break
                extended_table_df[column].apply(
                    lambda entry: BPCounter._count_extended_residues(entry)
                )
            
        BPCounter._print("Residue Pair Counts in Extended Tables")
    
    ###########################################################################
    # Auxiliary functions below
    ###########################################################################
      
    @staticmethod
    def _count_extended_residues(entry):
        entry = entry.split(',')
        for i in range(len(entry)):
            if i == 0: # skip R3DMA
                continue
            if entry[i] != 'REJECT':
                tool = BPCounter.tools[i]
                BPCounter.all_annotations[tool] += 1
          
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
                        
        if BPCounter._chain_less_than(r1_chain, r2_chain):
            return residue1, residue2
        elif BPCounter._chain_less_than(r2_chain, r1_chain):
            return residue2, residue1
        elif int(r1_res_id) < int(r2_res_id):
            return residue1, residue2
        else:
            return residue2, residue1
           
    #TODO verify the format because we have '1K' as an example of a chain
    @staticmethod
    def _chain_less_than(chain1: str, chain2: str) -> bool:
        c1_digit, c1_letter = BPCounter._parse_chain(chain1)
        c2_digit, c2_letter = BPCounter._parse_chain(chain2)
        
        if c1_digit != '' and c2_digit != '':
            if c1_digit == c2_digit:
                if c1_letter != '' and c2_letter != '':
                    return BPCounter._letter_comparison(c1_letter, c2_letter)
                return c2_letter != ''
            return int(c1_digit) < int(c2_digit)
        if c1_digit != '':
            return False
        if c2_digit != '':
            return True
        return BPCounter._letter_comparison(c1_letter, c2_letter)
    
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
    def _reset_counts() -> None:
        for tool in BPCounter.all_annotations.keys():
            BPCounter.all_annotations[tool] = 0
            
    @staticmethod
    def _print(headline: str) -> None:
        total_residue_pairs = 0
        print('='*60 + '\n' + headline + '\n' + '.'*60)
        for tool, item in BPCounter.all_annotations.items():
            if isinstance(item, dict):
                count = sum(len(pdb_set) for pdb_set in item.values())
                print(f'{tool}: {count}')
                total_residue_pairs += count
            elif isinstance(item, int):
                print(f'{tool} {item}')
                total_residue_pairs += item
        print(f'TOTAL: {total_residue_pairs}')
    
BPCounter.run()