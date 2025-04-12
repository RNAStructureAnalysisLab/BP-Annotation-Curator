# AUTHOR: Sameer Dingore

import requests
from DataManagement.AnnotationLoading.AnnotationUtilities.Tools import Logger

class Browser:
    def get_job_id(file_path):
        url = 'https://genesilico.pl/clarna/alg/cl_job_start/'
        session = requests.Session()
        response = session.get(url)
        
        if response.status_code != 200:
            Logger.log('Failed to access the website.')
            exit()
            
        csrf_token = response.cookies.get('csrftoken')
        data = {
            "csrfmiddlewaretoken": csrf_token,
            "pdb_as_text": "",
            "email": "",
            "use_rmsd": "1",
            "use_dist": "1",
            "show_fuzzy": "1",
        }
        
        def extract_job_id(url):
            start_str = "https://genesilico.pl/clarna/alg/cl_job/"
            end_str = "/status/"
            start_index = url.find(start_str)
            end_index = url.find(end_str)

            if start_index != -1 and end_index != -1:
                extracted_string = url[start_index + len(start_str) : end_index]
                return extracted_string
            else:
                return "error"
            
        with open(file_path, 'rb') as pdb_file:
            files = {'pdb_file': ('your_file.pdb', pdb_file)}
            response = session.post(url, data=data, files=files)
        
        if response.status_code == 200:
            Logger.log('Form submitted successfully for ' + file_path)
            redirect_url = response.url
            job_id = extract_job_id(redirect_url)
            return job_id
        else:
            Logger.log('Form submission failed for ' + file_path)
            return
        session.close()
            