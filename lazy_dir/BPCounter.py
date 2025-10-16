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
    
    all_annotations = {
        'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 'DSSR': set()
    }
    
    @staticmethod
    def run() -> None:
        BPCounter._initial_counts()
        BPCounter._reset_counts()
        BPCounter._preprocessor_counts()
        BPCounter._reset_counts()
    
    # =========================================================================
    # Pipeline Methods Below
    # =========================================================================
    
    @staticmethod
    def _initial_counts() -> None:
        for file in os.listdir(BPCounter.RAW_DATA_DIRECTORY):
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
                BPCounter.all_annotations[tool].add(f'{residue1}{residue2}')
        
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
                    BPCounter.all_annotations[tool] += len(residue_pairs.values())
        
        BPCounter._print("Residue Pair Counts in Preprocessed Data")
    
    ###########################################################################
    # Auxiliary functions below
    ###########################################################################
                
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
        #######################################################################
        # If chain was a single digit
        if chain1[0].isdigit() and chain2[0].isdigit():
            return int(chain1) < int(chain2)
        elif chain1[0].isdigit():
            return True
        elif chain2[0].isdigit():
            return False
        
        #######################################################################
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
            if isinstance(item, set):
                print(f'{tool}: {len(item)}')
                total_residue_pairs += len(item)
            elif isinstance(item, int):
                print(f'{tool} {item}')
                total_residue_pairs += item
        print(f'TOTAL: {total_residue_pairs}')
    
BPCounter.run()