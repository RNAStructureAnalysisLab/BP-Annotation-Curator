import os
import pandas as pd

class Tool_Benchmarker:
    DATASET_PATHS = [
        os.path.join(
            'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based', 'AllTools'
        ),
        os.path.join(
            'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based', 'WithoutR3DMA'
        ),
        os.path.join(
            'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based', 'WithoutFR'
        )
    ]
    
    EXTENDED_TABLES_DIRECTORY = os.path.join(
        'Data', 'Preprocessed', 'Extended_Tables'
    )
    
    OUTPUT_PATH = os.path.join(
        'Data', 'ResultsAnalysis', 'tool_counts.csv'
    )
    
# --- PUBLIC METHODS ---
    
    def run():
        benchmarked_rows = []
        index_to_tool = {0: "R3DMA", 1: "CL", 2: "FR", 3: "MC", 4: "RV", 5: "DSSR"}
        for path in Tool_Benchmarker.DATASET_PATHS:
            benchmark_row = {'Dataset': os.path.basename(path), 'Total BP': 0, 'R3DMA': 0, 'CL': 0, 'FR': 0, 'MC': 0, 'RV': 0, 'DSSR': 0}
            for table_name in os.listdir(path):
                consensus_df = pd.read_csv(os.path.join(path, table_name))
                cluster_df = pd.read_csv(os.path.join(Tool_Benchmarker.EXTENDED_TABLES_DIRECTORY, table_name)) # could optomize by caching this cluster rather than fetching it upon each loop iteration
                
                for column_name in consensus_df.columns[::-1]:
                    if '-' not in column_name:
                        break
                    for row_idx in range(len(consensus_df)):
                        consensus_entry = consensus_df.at[row_idx, column_name]
                        cluster_entry = cluster_df.at[row_idx, column_name]
                        Tool_Benchmarker._update_benchmark_row(
                            cluster_entry, consensus_entry, benchmark_row, index_to_tool
                        )
                        
            benchmarked_rows.append(benchmark_row)
            
        result_df = pd.DataFrame(benchmarked_rows)
        result_df.to_csv(Tool_Benchmarker.OUTPUT_PATH)

# --- PRIVATE METHODS ---

    def _update_benchmark_row(cluster_entry, consensus_entry, benchmark_row, index_to_tool):
        benchmark_row['Total BP'] += 1
        cluster_entry = cluster_entry.split(',')
        for i in range(len(cluster_entry)):
            if cluster_entry[i] == consensus_entry:
                benchmark_row[index_to_tool[i]] += 1



















'''
We want to create a row called "Total Residue Pairs", and then a row for each
tool counting how many residue pairs that tool agreed with for the consensus.
One column is "all tools", the other is "without R3DMA", and the last is
"without FR"
'''