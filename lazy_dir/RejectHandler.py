import os
import re
import json
import pandas as pd

class RejectHandler:
    INPUT_PATH = os.path.join(
        '..', 'Data', 'ResultsAnalysis', 'BasePairingAgreementsSummary', 'base_pairing_agreements.csv'
    )
    JSON_DIRECTORY = os.path.join(
        '..', 'Data', 'Preprocessed', 'JSON_Annotations'  
    )
    OUTPUT_PATH = os.path.join('base_pairing_agreements.csv')
    REJECTED_CASES = os.path.join('rejects.csv')
    TOOLS = ['CL', 'DSSR', 'FR', 'MC', 'RV']
    
    @staticmethod
    def run():
        rejects_df = pd.DataFrame(columns=['Cluster', 'Tool', 'PDB', 'Chain(s)', 'Residue IDs', 'Description', 'Complication'])
        json_files = RejectHandler._get_json()
        df = pd.read_csv(RejectHandler.INPUT_PATH)
        df.apply(
            lambda row: RejectHandler._fill_reject_entry(
                row['All Annotations'], row['Cluster'], row['PDB'],
                row['Chain(s)'], row['Residue IDs'], row['Nucleotides'], json_files,
                rejects_df
            ),
            axis=1
        )
        df.to_csv(RejectHandler.OUTPUT_PATH)
        rejects_df.to_csv(RejectHandler.REJECTED_CASES)
    
    # =========================================================================
    # helper methods below
    # =========================================================================
    
    @staticmethod
    def _get_json():
        json_files = {}
        for json_file in os.listdir(RejectHandler.JSON_DIRECTORY):
            if '_' in json_file:
                continue
            path = os.path.join(RejectHandler.JSON_DIRECTORY, json_file)
            tool = json_file.split('.')[0]
            with open(path, 'r') as f:
                json_files[tool] = json.load(f)
        
        return json_files
    
    @staticmethod
    def _fill_reject_entry(
            all_annotations, cluster, pdb, chain, residue_id, nucleotides, json_files,
            rejects_df
    ):
        all_annotations = all_annotations.split(',')
        residue_id, nucleotides = RejectHandler._parse_data(residue_id, nucleotides)
        for i in range(len(all_annotations)):
            if all_annotations[i] == 'REJECT':
                json = json_files[RejectHandler.TOOLS[i]]
                try:
                    if len(chain) == 1:
                        contact_type = json[pdb][f'{chain}{residue_id[0]}{nucleotides[0]}'][f'{chain}{residue_id[1]}{nucleotides[1]}']
                    else:
                        contact_type = json[pdb][f'{chain[0]}{residue_id[0]}{nucleotides[0]}'][f'{chain[2]}{residue_id[1]}{nucleotides[1]}']
                    if len(contact_type) > 1: 
                        raise AttributeError
                    contact_type = contact_type[0].split(" ")[0]
                    complication = None
                except Exception as e:
                    if isinstance(e, AttributeError):
                        complication = "more than one possible contact type"
                        contact_type = "REJECT"
                    elif isinstance(e, KeyError): #This shouldn't happen at all
                        complication = "KeyError, couldn't find residue in json"
                        contact_type = 'nbp'
                    else:
                        raise e
                all_annotations[i] = contact_type
                new_row = {'Cluster': cluster, 'Tool': RejectHandler.TOOLS[i], 'PDB': pdb, 'Chain(s)': chain, 'Residue IDs': residue_id, 'Description': contact_type, 'Complication': complication}
                rejects_df.loc[len(rejects_df)] = new_row
        all_annotations = ','.join(all_annotations)
        return all_annotations
    
    @staticmethod
    def _parse_data(residue_id, nucleotides):
        residue_id = re.findall(r'\d+', residue_id)
        nucleotides = re.findall(r'[A-Za-z]', nucleotides)
        return (residue_id, nucleotides)
    
RejectHandler.run()