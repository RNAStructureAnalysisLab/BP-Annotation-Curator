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
    def run(cls) -> None:
        jobs = cls._create_jobs()
        filtered_jobs = cls._evaluate_jobs(jobs)
        cls._run_jobs(jobs)
        
    #==========================================================================
    # HELPER METHODS BELOW
    #==========================================================================
    
    @staticmethod
    def _create_jobs() -> list[Job]:
        jobs = []
        
        return jobs
    
    @staticmethod
    def _evalutate_jobs(jobs: list[Job]) -> list[Job]:
        filtered_jobs = []
        
        return filtered_jobs
    
    @classmethod
    def _run_jobs(cls, jobs: list[Job]) -> None:
        pass