import os
import copy
import pandas as pd

class CTCounter():
    EXTENDED_TABLES_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'Extended_Tables'
    )
    CONSENSUS_TABLES_DIRECTORY = os.path.join(
        '..', 'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based', 'WithoutR3DMA'
    )
    
    @staticmethod
    def run() -> None:
        TOOLS = ['CL', 'FR', 'MC', 'RV', 'DSSR']
        CONTACT_TYPES_ROWS = [
            'cWW','cWH','cWS','cHW','cHH','cHS','cSW','cSH','cSS',
            'tWW','tWH','tWS','tHW','tHH','tHS','tSW','tSH','tSS','nbp', 'REJECT', '?'
        ]
        CONTACT_TYPES_COLS = [
            'cW.', 'c.W', 'cS.', 'c.S', 'cH.', 'c.H',
            'tW.', 't.W', 'tS.', 't.S', 'tH.', 't.H', 'nbp', 'REJECT', '?'
        ]
        tables_ct = {
            tool: pd.DataFrame(0, index=CONTACT_TYPES_ROWS, columns=CONTACT_TYPES_COLS)
            for tool in TOOLS
        }
        tables_bp = copy.deepcopy(tables_ct)
        
        for consensus_file in os.listdir(CTCounter.CONSENSUS_TABLES_DIRECTORY):
            consensus = pd.read_csv(os.path.join(CTCounter.CONSENSUS_TABLES_DIRECTORY, consensus_file))
            extended = pd.read_csv(os.path.join(CTCounter.EXTENDED_TABLES_DIRECTORY, consensus_file))
            assert len(consensus) == len(extended)
            CTCounter._compare_tables(consensus, extended, tables_bp, tables_ct)
        
        CTCounter._write_tables(tables_ct)
            
    @staticmethod
    def _compare_tables(consensus: pd.DataFrame, extended: pd.DataFrame, tables_bp: dict[pd.DataFrame], tables_ct: dict[pd.DataFrame]) -> None:
        for column in consensus.columns[::-1]:
            if '-' not in column:
                break
            for c_entry, e_entry in zip(consensus[column], extended[column]):
                CTCounter._update_count(tables_bp, tables_ct, c_entry, e_entry)
    
    @staticmethod
    def _update_count(tables_bp: dict[pd.DataFrame], tables_ct: dict[pd.DataFrame], c_entry: str, e_entry: str) -> None:
        if c_entry == 'REJECT':
            return
        
        e_entry = e_entry.split(',')[1:]
        index_to_tool = {0: 'CL', 1: 'FR', 2: 'MC', 3: 'RV', 4: 'DSSR'}
        for i in range(len(e_entry)):
            if c_entry == 'unresolved tie':
                c_entry = "?"
            tool = index_to_tool[i]
            tool_contact_type = e_entry[i]
            
            if 'nbp' != tool_contact_type and 'REJECT' != tool_contact_type and '.' not in tool_contact_type:
                first = tool_contact_type[0] + tool_contact_type[1] + '.'
                second = tool_contact_type[0] + '.' + tool_contact_type[2]
                tables_ct[tool].loc[c_entry, first] += 1
                tables_ct[tool].loc[c_entry, second] += 1
                    
            else:
                tables_ct[tool].loc[c_entry, tool_contact_type] += 1
            
    @staticmethod
    def _write_tables(tables_ct: dict[pd.DataFrame]) -> None:
        for tool, df in tables_ct.items():
            path = os.path.join('lazy_data', f'{tool}_confusion_matrix.csv')
            df.to_csv(path)
    
CTCounter.run()
