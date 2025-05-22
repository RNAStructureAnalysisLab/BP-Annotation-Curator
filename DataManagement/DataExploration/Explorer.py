# AUTHOR: Kristopher Church

import os
import pandas as pd

class Explorer:
    ANNOTATION_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    R3DMA_DIRECTORY = os.path.join(
        'Data', 'Raw', 'R3DMA', 'cluster_tables_3.95'
    )
    OUTPUT_DIRECTORY = os.path.join('Data', 'ExplorationFindings')
    
    nucleoside_types = {
        'R3DMA': set(), 'CL': set(), 'FR': set(), 'MC': set(), 'RV': set(), 
        'DSSR': set()
    }
    
    @staticmethod
    def explore():
        annotation_file_names = os.listdir(Explorer.ANNOTATION_DIRECTORY)
        
        # Get nucleosides from the annotation tools
        for file_name in annotation_file_names:
            tool_name = os.path.splitext(file_name)[0].split('_')[1]
            annotation = pd.read_csv(file_name)
            Explorer._document_nucleosides(annotation, tool_name)
        
        # Get nucleosides from R3DMA
        
        # TODO: write the nucleoside types
        
    @staticmethod
    def _document_nucleosides(annotation, tool_name):
        n_types = annotation['n_type']
        # TODO