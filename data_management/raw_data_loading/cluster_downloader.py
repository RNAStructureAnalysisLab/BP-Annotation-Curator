import os
import csv
import sys
import time
import shutil
import requests
from bs4 import BeautifulSoup
from project_management import job as Job
from project_management import configuration as config

class ClusterDownloader:
    
    @classmethod
    def run(cls, r3dma_version: str) -> None:
        jobs = cls._create_jobs()
        filtered_jobs = cls._evaluate_jobs(jobs)
        cls._run_jobs(jobs)
        
    #==========================================================================
    # HELPER CLASS METHODS BELOW
    #==========================================================================
    
    @classmethod
    def _create_jobs(cls) -> set[Job]:
        jobs = set()
        
        # Get components to make the jobs
        inputs = cls._get_input_links()
        outputs = cls._get_output_links(inputs)
        
        # Make jobs using the components
        for input_paths, output_paths in zip(inputs, outputs):
            jobs.add(Job(input_paths, output_paths))
        
        return jobs
    
    @staticmethod
    def _evalutate_jobs(jobs: set[Job]) -> set[Job]:
        filtered_jobs = set()
        
        return filtered_jobs
    
    @classmethod
    def _run_jobs(cls, jobs: set[Job]) -> None:
        pass
    
    @staticmethod
    def _get_input_links(cls) -> set[str]:
        inputs = set()
        
        return inputs
    
    #==========================================================================
    # HELPER STATIC METHODS
    #==========================================================================
    
    @staticmethod
    def _get_links(url: str) -> set[str]:
        attempts = 0
        
        while attempts < config.max_attempts_1:
            session = requests.Session()
            response = session.get(url)
            
            if response.status_code != 200:
                attempts += 1
            else:
                page_html = BeautifulSoup(response.text, 'html.parser')
                return set([
                    anchor['href'] for anchor in page_html.find_all('a', href=True) 
                    if anchor['href'].startswith( # links for specific clusters
                        'https://rna.bgsu.edu/rna3dhub/motif/view/'
                    ) or (
                        anchor['href'].startswith( # links for motif types
                        'https://rna.bgsu.edu/rna3dhub/motifs/release/'
                        ) and 'btn' in anchor.get('class', [])
                    )
                ])
            
        config.failed_links_1.add(url)
        return set()
    
    @staticmethod
    def _get_output_links(input_links: set[str]) -> set[str]:
        outputs = set()
        
        return outputs
        
        