import os
import copy
import pandas as pd

class ClFalsePositives():
    EXTENDED_TABLES_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'Extended_Tables'
    )
    CONSENSUS_TABLES_DIRECTORY = os.path.join(
        '..', 'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based', 'WithoutR3DMA'
    )
    
    @staticmethod
    def run() -> None:
        confusion_matrix_bp = {
            'CL': {"tp": 0, "fp": 0, "tn": 0, "fn": 0},
            'FR': {"tp": 0, "fp": 0, "tn": 0, "fn": 0},
            'MC': {"tp": 0, "fp": 0, "tn": 0, "fn": 0},
            'RV': {"tp": 0, "fp": 0, "tn": 0, "fn": 0},
            'DSSR': {"tp": 0, "fp": 0, "tn": 0, "fn": 0}
        }
        tp_bp = {
            'CL': {
                'cWW': 0, 'cWH': 0, 'cWS': 0, 'cHW': 0, 'cHH': 0, 'cHS': 0, 'cSW': 0, 'cSH': 0, 'cSS': 0, 
                'tWW': 0, 'tWH': 0, 'tWS': 0, 'tHW': 0, 'tHH': 0, 'tHS': 0, 'tSW': 0, 'tSH': 0, 'tSS': 0
            },
            'FR': {
                'cWW': 0, 'cWH': 0, 'cWS': 0, 'cHW': 0, 'cHH': 0, 'cHS': 0, 'cSW': 0, 'cSH': 0, 'cSS': 0, 
                'tWW': 0, 'tWH': 0, 'tWS': 0, 'tHW': 0, 'tHH': 0, 'tHS': 0, 'tSW': 0, 'tSH': 0, 'tSS': 0
            },
            'MC': {
                'cWW': 0, 'cWH': 0, 'cWS': 0, 'cHW': 0, 'cHH': 0, 'cHS': 0, 'cSW': 0, 'cSH': 0, 'cSS': 0, 
                'tWW': 0, 'tWH': 0, 'tWS': 0, 'tHW': 0, 'tHH': 0, 'tHS': 0, 'tSW': 0, 'tSH': 0, 'tSS': 0
            },
            'RV': {
                'cWW': 0, 'cWH': 0, 'cWS': 0, 'cHW': 0, 'cHH': 0, 'cHS': 0, 'cSW': 0, 'cSH': 0, 'cSS': 0, 
                'tWW': 0, 'tWH': 0, 'tWS': 0, 'tHW': 0, 'tHH': 0, 'tHS': 0, 'tSW': 0, 'tSH': 0, 'tSS': 0
            },
            'DSSR': {
                'cWW': 0, 'cWH': 0, 'cWS': 0, 'cHW': 0, 'cHH': 0, 'cHS': 0, 'cSW': 0, 'cSH': 0, 'cSS': 0, 
                'tWW': 0, 'tWH': 0, 'tWS': 0, 'tHW': 0, 'tHH': 0, 'tHS': 0, 'tSW': 0, 'tSH': 0, 'tSS': 0
            }
        }
        confusion_matrix_ct = copy.deepcopy(confusion_matrix_bp)
        fp_bp = copy.deepcopy(tp_bp)
        tp_ct = copy.deepcopy(tp_bp)
        fp_ct = copy.deepcopy(tp_bp)
        
        for consensus_file in os.listdir(ClFalsePositives.CONSENSUS_TABLES_DIRECTORY):
            consensus = pd.read_csv(os.path.join(ClFalsePositives.CONSENSUS_TABLES_DIRECTORY, consensus_file))
            extended = pd.read_csv(os.path.join(ClFalsePositives.EXTENDED_TABLES_DIRECTORY, consensus_file))
            assert len(consensus) == len(extended)
            ClFalsePositives._compare_tables(consensus, extended, confusion_matrix_bp, confusion_matrix_ct, tp_bp, fp_bp, tp_ct, fp_ct)
        
        ClFalsePositives._store_counts(confusion_matrix_bp, 'confusion_matrix_bp.csv')
        ClFalsePositives._store_counts(confusion_matrix_ct, 'confusion_matrix_ct.csv')
        ClFalsePositives._store_counts(tp_bp, 'true_positives_bp.csv')
        ClFalsePositives._store_counts(fp_bp, 'false_positives_bp.csv')
        ClFalsePositives._store_counts(tp_ct, 'true_positives_ct.csv')
        ClFalsePositives._store_counts(fp_ct, 'false_positives_ct.csv')
    
    @staticmethod
    def _compare_tables(consensus: pd.DataFrame, extended: pd.DataFrame, confusion_matrix_bp: dict, confusion_matrix_ct: dict, tp_bp: dict, fp_bp: dict, tp_ct: dict, fp_ct: dict) -> None:
        for column in consensus.columns[::-1]:
            if '-' not in column:
                break
            for c_entry, e_entry in zip(consensus[column], extended[column]):
                ClFalsePositives._count_base_pairing(c_entry, e_entry, confusion_matrix_bp, tp_bp, fp_bp)
                ClFalsePositives._count_contact_types(c_entry, e_entry, confusion_matrix_ct, tp_ct, fp_ct)
    
    @staticmethod
    def _count_base_pairing(c_entry: str, e_entry: str, confusion_matrix: dict, tp: dict, fp: dict) -> None:
        if c_entry == 'REJECT':
            return
        index_to_tool = {0: 'CL', 1: 'FR', 2: 'MC', 3: 'RV', 4: 'DSSR'}
        
        e_entry = e_entry.split(',')[1:]
        for i in range(len(e_entry)):
            contact_type = e_entry[i]
            if contact_type == 'REJECT':
                continue
            
            is_base_pairing = c_entry != 'nbp'
            tool = index_to_tool[i]
            if contact_type == 'nbp' and is_base_pairing: # false negative
                confusion_matrix[tool]["fn"] += 1
            elif contact_type == 'nbp' and not is_base_pairing: # true negative
                confusion_matrix[tool]["tn"] += 1
            elif contact_type != 'nbp' and is_base_pairing: # true positive
                confusion_matrix[tool]["tp"] += 1
                tp[tool][contact_type] += 1
            else: # false positive
                confusion_matrix[tool]["fp"] += 1
                fp[tool][contact_type] += 1
                
    @staticmethod
    def _count_contact_types(c_entry: str, e_entry: str, confusion_matrix: dict, tp: dict, fp: dict) -> None:
        if c_entry == 'REJECT':
            return
        index_to_tool = {0: 'CL', 1: 'FR', 2: 'MC', 3: 'RV', 4: 'DSSR'}
        
        e_entry = e_entry.split(',')[1:]
        for i in range(len(e_entry)):
            contact_type = e_entry[i]
            if contact_type == 'REJECT':
                continue
            
            is_base_pairing = c_entry != 'nbp'
            tool = index_to_tool[i]
            if contact_type == c_entry and is_base_pairing: # true positive
                confusion_matrix[tool]["tp"] += 1
                tp[tool][contact_type] += 1
            elif contact_type == c_entry and not is_base_pairing: # true negative
                confusion_matrix[tool]["tn"] += 1
            elif contact_type != c_entry and is_base_pairing: # false negative
                confusion_matrix[tool]["fn"] += 1
            else: # false positive
                confusion_matrix[tool]["fp"] += 1
                fp[tool][contact_type] += 1
                
    @staticmethod
    def _store_counts(confusion_matrix: dict, file_name: str) -> None:
        df = pd.DataFrame(confusion_matrix)
        df.to_csv(file_name)
    
ClFalsePositives.run()