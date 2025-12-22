# AUTHOR: Kristopher Church

import os
import json
from datetime import datetime, timezone

'''
PURPOSE:
USE:
'''
class ManifestManager:
    ALLOWED_STATUSES = {"missing", "running", "failed", "success"}
    
    '''
    INPUT:
    ACTION:
    OUTPUT:
    '''
    @staticmethod
    def load(manifest_path: str) -> dict | None:
        if not os.path.exists(manifest_path):
            return None
        
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest_data = json.load(f)
            
        ManifestManager._validate(manifest_data)
        return manifest_data
        
    
    '''
    INPUT:
    ACTION:
    OUTPUT:
    '''
    @staticmethod
    def save(manifest_path: str, manifest_data: dict) -> None: 
        ManifestManager._validate(manifest_data)
        
        tmp_path = manifest_path + ".tmp"
        with tmp_path.open("w", encoding="utf-8") as f:
            json.dump(manifest_data, f, indent=2)
            
        # Ensures atomic writing of json data (safe on Windows/Linux)
        os.replace(tmp_path, manifest_path)
        
    '''
    '''
    @staticmethod 
    def update_manifest(manifest_data : dict, step_name: str, status: str, output: dict = None) -> dict:
        if status not in ManifestManager.ALLOWED_STATUSES:
            raise ValueError(f"Invalid status: {status!r}")
        
        if manifest_data is None:
            manifest_data = {
                "step": step_name,
                "status": None,
                "started_at": None,
                "completed_at": None,
                "output": None
            }
            
        if manifest_data.get("step") != step_name:
            raise ValueError(
                f"Manifest file name mismatch: {manifest_data.get('step')} != {step_name}"
            )
        
        manifest_data["status"] = status
        
        if status == "running":
            manifest_data["started_at"] = datetime.now(timezone.utc).isoformat()
        elif status == "success":
            manifest_data["completed_at"] = datetime.now(timezone.utc).isoformat()
            manifest_data["output"] = output
            
        return manifest_data
        
    #==========================================================================
    # PRIVATE HELPER METHODS BELOW
    #==========================================================================
    '''
    INPUT:
    ACTION:
    OUTPUT:
    '''
    @staticmethod
    def _validate(manifest_data: dict) -> None:
        status = manifest_data.get("status")
        if status is not None and status not in ManifestManager.ALLOWED_STATUSES:
            raise ValueError(f"Invalid manifest status: {status!r}")