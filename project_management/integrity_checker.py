# AUTHOR: Kristopher Church

import os
import hashlib
from typing import Optional
from project_management.configuration import Configuration

'''
PURPOSE: Determine whether outputs of each pipeline step exist and ar valid
USE: IntegrityChecker.check_step(manifest_data) -> bool
'''   
class IntegrityChecker:
    FILE_HASHES_CACHE = {} # key is file name, value is its sha256 hash
    
    '''
    
    '''
    @staticmethod
    def get_corrupted_files(manifest_data: Optional[dict], config: Configuration) -> set[str]:
        if manifest_data is None:
            return set()
        
        return IntegrityChecker._compare_hashes(manifest_data, config)
    
    #==========================================================================
    # PRIVATE METHODS
    #==========================================================================
    
    @staticmethod
    def _compare_hashes(manifest_data: dict, config: Configuration) -> set[str]:
        corrupted_filepaths = set()
        
        for directory, directory_files in manifest_data.items():
            directory_path = config.get(directory)
            actual_filenames = set(os.listdir(directory_path))
            expected_filenames = set(directory_files.keys())
            all_filenames = actual_filenames | expected_filenames
            
            for filename in all_filenames:
                filepath = os.path.join(directory_path, filename)
                expected_hash = directory_files.get(filename)
                if expected_hash is None:
                    # Delete this file since not in manifest
                    os.remove(filepath)
                    continue
                
                actual_hash = IntegrityChecker._compute_hash(filepath)
                if expected_hash != actual_hash:
                    corrupted_filepaths.add(filepath)
        
        return corrupted_filepaths
            
    @classmethod
    def _compute_hash(cls, filepath: str) -> Optional[str]:
        if filepath in cls.FILE_HASHES_CACHE:
            return cls.FILE_HASHES_CACHE[filepath]
        
        if not os.path.isfile(filepath):
            return None
        
        sha256 = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        cls.FILE_HASHES_CACHE[filepath] = sha256.hexdigest()
        
        return sha256.hexdigest()