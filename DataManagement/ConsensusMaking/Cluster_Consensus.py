# AUTHOR: Kristopher Church

import os
import pandas as pd

class Cluster_Consensus:
    INPUT_DIRECTORY = os.path.join(
        'Data', 'Consensus', 'Tool_Consensus', 'Mode_Based_Profiles'
    )
    
    profile = pd.DataFrame
    
    @staticmethod
    def run():
        # Access each tool consensus profile
        for csv_file_name in os.listdir(Cluster_Consensus.INPUT_DIRECTORY):
            csv_file_path = os.path.join(
                Cluster_Consensus.INPUT_DIRECTORY, csv_file_name
            )
            Cluster_Consensus.profile = pd.read_csv(csv_file_path)
            
            Cluster_Consensus._clusterwise_consensus()
            Cluster_Consensus._export()
            
    @staticmethod
    def _clusterwise_consensus():
        
        
    @staticmethod
    def _export():
        pass        