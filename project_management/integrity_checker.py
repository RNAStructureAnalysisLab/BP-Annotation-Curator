# AUTHOR: Kristopher Church

import os
import hashlib
from typing import Optional
from project_management import configuration as config

'''
PURPOSE: Determine whether outputs of each pipeline step exist and ar valid
USE: IntegrityChecker.check_step(step_name, manifest_data) -> bool
'''   
class IntegrityChecker:
    
    @staticmethod
    def check_step(step_name: str, manifest_data: Optional[dict]) -> bool:
        checkers = {
            "cluster_downloader": IntegrityChecker._check_cluster_downloader,
            "pdb_downlaoder": IntegrityChecker._check_pdb_downloader,
        }
        
        if step_name not in checkers:
            return 