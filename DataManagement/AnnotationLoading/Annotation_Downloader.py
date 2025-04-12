# AUTHOR: Sameer Dingore and Kristopher Church

import os
import shutil
from DataManagement.AnnotationLoading.AnnotationUtilities.Threader import Threader
from DataManagement.AnnotationLoading.AnnotationUtilities.Job_Manager import Job_Manager

class Annotation_Downloader:
    INPUT_FILENAME = os.path.join('Data', 'Raw', 'RCSB', 'used_pdb_ids.txt')
    RETRY_INPUT_FILENAME = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'failed_pdb_ids'
    )
    OUTPUT_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    
    @staticmethod
    def download():
        # Reset the contents of the AnnotationTools folder to be blank
        if os.path.exists(Annotation_Downloader.OUTPUT_DIRECTORY):
            shutil.rmtree(Annotation_Downloader.OUTPUT_DIRECTORY)
        os.makedirs(Annotation_Downloader.OUTPUT_DIRECTORY)
        
        Threader.create_threads(
            Annotation_Downloader.INPUT_FILENAME, Job_Manager.create_job
        )
        
    @staticmethod 
    def retry_download():
        Threader.create_threads(
            Annotation_Downloader.RETRY_INPUT_FILENAME, Job_Manager.create_job
        )
        
    @staticmethod 
    def if_failed_pdb_ids():
        failed_pdb_ids = Threader.get_failed_pdb_ids()
        with open(Threader.RETRY_INPUT_FILENAME, 'w') as file:
            for pdb_id in failed_pdb_ids:
                file.write(f'{pdb_id}\n')
        
        if failed_pdb_ids: # If the list is not empty
            return True
        return False
        
        
        