# AUTHOR: Sameer Dingore

import concurrent.futures
import time
from DataManagement.AnnotationLoading.AnnotationUtilities.Tools import Logger, Redirecter, Downloader

class Job_Manager:
    
    @staticmethod
    def create_job(ids):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(Job_Manager._worker, ids)
            result = future.result() # Blocking. Waiting for thread to finish
            
        return result
       
    # INPUT: a 2-tuple of the job_id and the pdb_id
    @staticmethod
    def _worker(ids):
        job_id, pdb_id = ids
        while True:
            if Job_Manager._check_job_status(job_id): # If job has finished
                Job_Manager._finish_job(job_id, pdb_id)
                return True
            else:
                Logger.log(
                    f'Job {job_id} is still processing. Sleeping for 10 minutes'
                )
                time.sleep(600)
    
    # INPUT:
    # OUTPUT: returns 'True' if the job has finished, or 'False' if still being
    # processed
    @staticmethod
    def _check_job_status(job_id):
        url = 'https://genesilico.pl/clarna/alg/cl_job/' + job_id + '/status/'
        return Redirecter.check_redirect(url)
    
    # INPUT:
    # OUTPUT:
    @staticmethod
    def _finish_job(job_id, pdb_id):
        Downloader.download_clarna_results(job_id, pdb_id)
        Logger.log(f'Job {job_id} is finished')