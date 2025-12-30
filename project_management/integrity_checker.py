# AUTHOR: Kristopher Church

import os
import hashlib
from typing import Optional
from project_management import configuration as config

'''
PURPOSE: Determine whether outputs of each pipeline step exist and ar valid
USE: IntegrityChecker.check_step(manifest_data) -> bool
'''   
class IntegrityChecker:
    
    @staticmethod
    def check_step(manifest_data: Optional[dict], prev_manifest_data: Optional[dict]=None) -> bool:
        if manifest_data is None:
            return False
        
        step_name = manifest_data.get("step")
        
        if step_name not in IntegrityChecker.CHECKERS:
            return False
        
        status = manifest_data.get("status")
        if status in ("failed", "running"):
            return False
        
        
        return IntegrityChecker.CHECKERS[step_name](manifest_data, prev_manifest_data)
    
    #==========================================================================
    # PRIVATE CHECKER METHODS
    #==========================================================================
    
    @staticmethod
    def _check_cluster_downloader(manifest_data: dict, prev_manifest_data: Optional[dict]) -> bool:
        # Check if no data from RNA 3D Motif Atlas
        if not os.path.isdir(config.CLUSTERS_DIR):
            return False
        motif_cluster_filenames = [f for f in os.listdir(config.CLUSTERS_DIR) if f.endswith(".csv")]
        num_clusters = len(motif_cluster_filenames)
        if num_clusters == 0: 
            return False
        
        # Check if PDB info not extracted from RNA 3D Motif Atlas data
        if not os.path.isfile(config.USED_PDB_IDS_TXT):
            return False
        with open(config.USED_PDB_IDS_TXT, "r") as f:
            num_pdbs = sum(1 for line in f if line.strip())
        if num_pdbs == 0:
            return False
        
        # Check if counts of outputs match those reported in the manifest
        counts = manifest_data.get("counts")
        if counts is not None:
            expected_clusters = counts.get("clusters", 0)
            expected_pdbs = counts.get("pdb_ids", 0)
            if expected_clusters != num_clusters or expected_pdbs != num_pdbs:
                return False
        
        return True
    
    @staticmethod
    def _check_pdb_downloader(manifest_data: dict, prev_manifest_data: Optional[dict]) -> bool:
        # Ensure chronological consistency with previous step
        if prev_manifest_data is not None and prev_manifest_data.get("step") == "cluster_downloader":
            prev_completed = prev_manifest_data.get("completed_at")
            curr_started = manifest_data.get("started_at")
            if prev_completed and curr_started and curr_started < prev_completed:
                return False
        else:
            raise ValueError("pdb_downloader requires prev_manifest_data from cluster_downloader")
        
        # Check if no data from RCSB
        if not os.path.isdir(config.RCSB_DIR):
            return False
        pdb_filenames = [f for f in os.listdir(config.PDB_DIR) if f.endswith(".pdb")]
        pdbx_filenames = [f for f in os.listdir(config.PDBX_DIR) if f.endswith(".cif")]
        num_pdbs = len(pdb_filenames)
        num_pdbxs = len(pdbx_filenames)
        if num_pdbs == 0 and num_pdbxs == 0:
            return False
        
        # Check if counts of inputs match those reported in the manifest
        '''
        NOTE: we already know num_pdbs is nonzero and consistent with the reported output
            of the previous manifest, since _check_cluster_downloader() should always be
            executed before this method.
        '''
        with open(config.USED_PDB_IDS_TXT, "r") as f:
            num_pdb_ids = sum(1 for line in f if line.strip())
        if num_pdb_ids != num_pdbs + num_pdbxs:
            return False
        
        # Check if counts of outputs match those reported in the manifest
        counts = manifest_data.get("counts")
        if counts is not None:
            expected_pdbs = counts.get("pdb_ids", 0)
            expected_pdbxs = counts.get("pdbx_ids", 0)
            if expected_pdbs != num_pdbs or expected_pdbxs != num_pdbxs:
                return False
            
        return True

IntegrityChecker.CHECKERS = {
    "cluster_downloader": IntegrityChecker._check_cluster_downloader,
    "pdb_downloader": IntegrityChecker._check_pdb_downloader,
}