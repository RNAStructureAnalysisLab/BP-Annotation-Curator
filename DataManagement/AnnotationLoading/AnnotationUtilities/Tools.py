# AUTHOR: Sameer Dingore

import datetime
import os

class Logger:
    def log(input_string):
        current_datetime = datetime.datetime.now()
        formatted_datetime = current_datetime.strftime("[%Y-%m-%d %H:%M:%S]")
        log_entry = f"{formatted_datetime} : {input_string}"

        # Print the string
        print(log_entry)

        # Append the string to 'log.txt' file
        file_path = os.path.join(
            'DataManagement', 'AnnotationLoading', 'log.txt'
        )
        with open(file_path, 'a') as file:
            file.write(log_entry + '\n')
            
# =============================================================================

import requests
        
class Redirecter:
    def check_redirect(url):
        try:
            response = requests.head(url, allow_redirects=True)
            return response.url != url
        except requests.exceptions.RequestException:
            return False

# =============================================================================

import os
import requests
from urllib.parse import urlparse

class Downloader:
    @staticmethod 
    def download_clarna_results(job_id, pdb_file):
        base_url = "https://genesilico.pl/clarna/alg/cl_job"
        classifiers = ["CL", "RV", "MC", "FR", "MO"]
        output_folder = os.path.join(
            "Data", "Raw", "AnnotationTools", "Annotations"
        )
        
        for classifier in classifiers:
            url = f"{base_url}/{job_id}/result/export/?classifier={classifier}&format=csv&imperfect=1&contact_types=bp,stacking,base_phosphate,base_ribose,other&score_tolerance=0.55"
            response = requests.get(url)
            if response.status_code == 200:
                content_disposition = response.headers.get("Content-Disposition")
                if content_disposition:
                    filename = content_disposition.split("filename=")[1]
                else:
                    parsed_url = urlparse(url)
                    filename = os.path.basename(parsed_url.path)
                extension = os.path.splitext(filename)[-1]
                new_filename = f"{pdb_file}_{classifier}{extension}"
                file_path = os.path.join(output_folder, new_filename)

                with open(file_path, "wb") as file:
                    file.write(response.content)
                    
                Logger.log(f"File for {classifier} downloaded and saved as {file_path}")
            else:
                Logger.log(
                    f"Failed to download file for {classifier}. HTTP status code: {response.status_code}"
                )