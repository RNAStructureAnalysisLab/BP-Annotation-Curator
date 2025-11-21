import os
import pandas as pd

class TableMaker():
    R3DMA_DIR = os.path.join(
        '..', 'Data', 'Raw', 'R3DMA', 'cluster_tables_3.95'
    )
    EXTENDED_TABLES_DIR = os.path.join(
        '..', 'Data', 'Preprocessed', 'Extended_Tables'
    )
    COMBINED_TABLE_PATH = os.path.join(
        '..', 'Data', 'ResultsAnalysis', 'BasePairingAgreementsSummary', 'base_pairing_agreements.csv'
    )
    OUTPUT_DIR = os.path.join(
        '..', 'Data', 'ResultsAnalysis', 'BenchmarkDatasetInfo'
    )
    
    @staticmethod
    def run() -> None:
        #r3dma_df = TableMaker._generate_r3dma_table()
        #extended_df = TableMaker._generate_extended_table()
        benchmark_dataset, benchmark_statistics = TableMaker._generate_benchmark_info()
        #TableMaker._output({'r3dma.csv': r3dma_df, 'extended.csv': extended_df})
        TableMaker._output({'benchmark.csv': benchmark_dataset, 'benchmark_statistics.csv': benchmark_statistics})
        
    @staticmethod
    def _generate_r3dma_table() -> pd.DataFrame:
        loop_types = ["HL", "IL", "Junction", "Total"]
        columns = ["Cluster Count", "Motif Instance Count", "Basepair Count"]
        df = pd.DataFrame(0, index=loop_types, columns=columns)
        for cluster_table in os.listdir(TableMaker.R3DMA_DIR):
            cluster_type = cluster_table.split('_')[0]
            if cluster_type[0] == 'J':
                cluster_type = "Junction"
            df.loc[cluster_type, 'Cluster Count'] += 1
            df.loc['Total', 'Cluster Count'] += 1
                
            path = os.path.join(
                TableMaker.R3DMA_DIR, cluster_table
            )
            cluster_df = pd.read_csv(path)
            
            df.loc[cluster_type, 'Motif Instance Count'] += len(cluster_df)
            df.loc['Total', 'Motif Instance Count'] += len(cluster_df)
            
            for col in cluster_df.columns[::-1]:
                if '-' not in col:
                    break
                for entry in cluster_df[col]:
                    if entry != "":
                        df.loc[cluster_type, 'Basepair Count'] += 1
                        df.loc['Total', 'Basepair Count'] += 1
            
        return df.reset_index(names="Loop Type")  
    
    @staticmethod
    def _generate_extended_table() -> pd.DataFrame:
        loop_types = ["HL", "IL", "Junction", "Total"]
        columns = ["Cluster Count", "Motif Instance Count", "Residue Pair Count"]
        df = pd.DataFrame(0, index=loop_types, columns=columns)
        for extended_table in os.listdir(TableMaker.EXTENDED_TABLES_DIR):
            cluster_type = extended_table.split('_')[0]
            if cluster_type[0] == 'J':
                cluster_type = "Junction"
            df.loc[cluster_type, 'Cluster Count'] += 1
            df.loc['Total', 'Cluster Count'] += 1
            path = os.path.join(
                TableMaker.EXTENDED_TABLES_DIR, extended_table
            )
            extended_df = pd.read_csv(path)
            df.loc[cluster_type, 'Motif Instance Count'] += len(extended_df)
            df.loc['Total', 'Motif Instance Count'] += len(extended_df)
            
            for col in extended_df.columns[::-1]:
                if '-' not in col:
                    break
                num_residue_pairs = len(extended_df[col])
                df.loc[cluster_type, 'Residue Pair Count'] += num_residue_pairs
                df.loc['Total', 'Residue Pair Count'] += num_residue_pairs
            
        return df.reset_index(names="Loop Type")  
    
    @staticmethod
    def _generate_benchmark_info() -> (pd.DataFrame, pd.DataFrame):
        combined_df = pd.read_csv(TableMaker.COMBINED_TABLE_PATH)
        df = combined_df[['PDB', 'Chain(s)', 'Residue IDs', 'Nucleotides', 'Expected Contact Type']]
        df = df[(df['Expected Contact Type'] != "unresolved tie") & (df['Expected Contact Type'] != "REJECT")]
        df['Expected Contact Type'] = df['Expected Contact Type'].replace("nbp", "None")
        
        contact_types = [
            "cWW", "cWH", 'cWS', 'cHW', 'cHH', 'cHS', 'cSW', 'cSH', 'cSS',
            "tWW", 'tWH', 'tWS', 'tHW', 'tHH', 'tHS', 'tSW', 'tSH', 'tSS', 'None'
        ]
        columns=["Count", "All Agree", "4 Agree", "3 Agree", "2 Agree", "1 Agree"]
        statistics_df = pd.DataFrame(0, index=contact_types, columns=columns)
        for agree_counts, consensus in zip(combined_df['Expected Matches Count'], combined_df['Expected Contact Type']):
            if consensus == "unresolved tie" or consensus == "REJECT":
                continue
            if consensus == "nbp":
                consensus = "None"
            statistics_df.loc[consensus, 'Count'] += 1
            if agree_counts == 5:
                agree_counts = "All"
            statistics_df.loc[consensus, f'{agree_counts} Agree'] += 1
        
        return (df, statistics_df.reset_index(names="Contact Type"))
    
    @staticmethod
    def _output(dataframes: dict[str, pd.DataFrame]) -> None:
        os.makedirs(TableMaker.OUTPUT_DIR, exist_ok=True)
        for name, df in dataframes.items():
            path = os.path.join(TableMaker.OUTPUT_DIR, name)
            df.to_csv(path, index=False)
    
TableMaker.run()