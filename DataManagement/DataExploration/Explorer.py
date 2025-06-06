# AUTHOR: Kristopher Church

import os
import json
import pandas as pd
from pandas.errors import EmptyDataError

class Explorer:
    ANNOTATION_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    R3DMA_DIRECTORY = os.path.join(
        'Data', 'Raw', 'R3DMA', 'cluster_tables_3.95'
    )
    OUTPUT_DIRECTORY = os.path.join('Data', 'ExplorationFindings')
    
    descriptions = {
        'R3DMA': set(), 'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 
        'DSSR': set()
    }
    
    @staticmethod
    def document_descriptions():
        file_names = os.listdir(Explorer.ANNOTATION_DIRECTORY)
        
        # Get nucleosides from the annotation tools
        for file_name in file_names:
            tool_name = os.path.splitext(file_name)[0].split('_')[1]
            if tool_name == "MO":
                continue
            try:
                annotation = pd.read_csv(
                    os.path.join(
                        Explorer.ANNOTATION_DIRECTORY, file_name
                    )
                )
                Explorer._add_descriptions(annotation, tool_name)
            # All base pairs have modified nucleoside or has no base pairs
            except EmptyDataError: 
                continue
        
        # Get descriptions (contact types) from R3DMA
        file_names = os.listdir(Explorer.R3DMA_DIRECTORY)
        for file_name in file_names:
            cluster_table = pd.read_csv(
                os.path.join(
                    Explorer.R3DMA_DIRECTORY, file_name
                )
            )
            Explorer._add_descriptions(cluster_table, 'R3DMA')
        
        # Write the descriptions as a json
        descriptions_serializable = {
            k: sorted(list(v)) for k, v in Explorer.descriptions.items()
        }
        output_path = os.path.join(Explorer.OUTPUT_DIRECTORY, "descriptions.json")
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(descriptions_serializable, f, indent=2, ensure_ascii=False)
        
    @staticmethod
    def _add_descriptions(data, tool_name):
        if tool_name != 'R3DMA':
            Explorer.descriptions[tool_name].update(
                data['description'].dropna()
            )
        else:
            reversed_columns = data.columns[::-1]
            for column_name in reversed_columns:
                if '-' not in column_name:
                    break
                Explorer.descriptions[tool_name].update(
                    data[column_name].dropna()
                )
                