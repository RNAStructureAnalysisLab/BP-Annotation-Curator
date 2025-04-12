# AUTHOR: Sameer Dingore

import concurrent.futures
import os
import time
from DataManagement.AnnotationLoading.AnnotationUtilities.Browser import Browser
from DataManagement.AnnotationLoading.AnnotationUtilities.Tools import Logger

class Threader:
    MAX_WORKERS = 10
    SLEEP_TIME = 5
    FAILED_ATTEMPTS_FILE = os.path.join(
        'DataManagement', 'AnnotationLoading', 'failed_attempts.txt'
    )
    INPUT_DIRECTORY = os.path.join('Data', 'Raw', 'RCSB', 'PDB_Files')
    
    handler_results = {}
    failed_pdb_ids = []
    
    # INPUT:
    # ACTION:
    # OUTPUT:
    def create_threads(input_file_name, processing_function):
        successfully_processed = set()
        failed_attempts_file = open(Threader.FAILED_ATTEMPTS_FILE, 'w')
        
        try:
            with open(input_file_name, 'r') as file:
                # Create a list of PDBs that need to be processed with a job
                pdb_ids = file.readlines()
                pdb_ids = [pdb_id.strip() for pdb_id in pdb_ids]
                
                # Create jobs for each PDB ID. Store results of successful
                # jobs and remove them from list of active jobs
                with concurrent.futures.ThreadPoolExecutor(
                        max_workers = Threader.MAX_WORKERS
                ) as executor:
                    active_processing = [] # pairs of (pdb_id, future)
                    while pdb_ids or active_processing:
                        for pdb_id, future in active_processing:
                            if future.done():
                                active_processing.remove((pdb_id, future))
                                pdb_id_result = future.result()
                                if pdb_id_result:
                                    successfully_processed.add(pdb_id_result)
                                else:
                                    failed_attempts_file.write(
                                        f'''
                                        Failed to process PDB ID: 
                                        {pdb_id_result}\n
                                        '''
                                    )
                                    Threader.failed_pdb_ids.append(pdb_id)
                                    
                        # Create jobs for any idle workers
                        available_slots = (
                            Threader.MAX_WORKERS - len(active_processing)
                        )
                        if available_slots > 0 and pdb_ids:
                            new_pdb_ids = pdb_ids[:available_slots]
                            pdb_ids = pdb_ids[available_slots:]
                            
                            # Begin processing the additional PDB IDs
                            for pdb_id in new_pdb_ids:
                                future = executor.submit(
                                    Threader.process_pdb_file, pdb_id, 
                                    processing_function
                                )
                                active_processing.append((pdb_id, future))
                                
                        time.sleep(Threader.SLEEP_TIME)
        
        except Exception as e:
            Logger.log(f'An error occurred: {e}')
            failed_attempts_file.write(f'Failed with error: {str(e)}\n')
            
        failed_attempts_file.close()
        
        
    # INPUT:
    # ACTION:
    # OUTPUT:
    def process_pdb_file(pdb_id, processing_function):
        pdb_file_path = os.path.join(Threader.INPUT_DIRECTORY, f'{pdb_id}.pdb')
        
        with concurrent.futures.ThreadPoolExecutor() as executor:
            if os.path.exists(pdb_file_path):
                job_id = Browser.get_job_id(pdb_file_path)
                Logger.log("Job is created with job id: " + job_id)
                Threader.handler_results[job_id] = executor.submit(
                    processing_function, [job_id, pdb_id]
                )
                for job_id, future in Threader.handler_results.items():
                    if future.result():
                        Logger.log(f'Threader: Job {job_id} is finished.')
                        return True
                    Logger.log(f'Threader: Job {job_id} is still processing.')
                    return False
            return False
        
    # ACTION: temporarily copies the list of failed PDB IDs and resets it to be
    # blank.
    # OUTPUT: returns the temporary copy of failed_PDB IDs
    def get_failed_pdb_ids():
        temp = Threader.failed_pdb_ids
        Threader.failed_pdb_ids = []
        
        return temp