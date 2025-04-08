# AUTHOR: Kristopher Church

# Looks through INPUT_FILE which is a text file representing the list of all
# PDB IDs found in the hairpin loop, internal loop, and three-way junction
# motif clusters on RNA 3D Motif Atlas. For each PDB ID in this file, the
# script visits RCSB to download either a corresponding PDB or PDBx file.

# =============================================================================

import os
import requests
import time
import shutil

class PDB_Downloader:
    CRAWL_DELAY = 1 # No delay specified at https://www.rcsb.org/robots.txt
    INPUT_DIRECTORY = 'Data/Raw/RCSB'
    INPUT_FILE = f'{INPUT_DIRECTORY}/pdb_ids.txt'
    PDB_OUTPUT_DIRECTORY = 'Data/Raw/RCSB/PDB_Files'
    PDBX_OUTPUT_DIRECTORY = 'Data/Raw/RCSB/PDBx_Files'
    RCSB_URL = 'https://files.rcsb.org/download'
    SESSION = requests.Session()
    
    # INPUT: a text file path with the PDB IDs to download, the ones used in
    # RNA 3d Motif Atlas.
    # ACTION: Calls a helper function 'download_pdb' to download each
    # individual PDB. If that download was unsuccessful, stores that PDB ID in
    # a file called 'failed.txt'
    @staticmethod
    def download():
        # Remake directories upon each run, ensures empty directories
        if os.path.exists(PDB_Downloader.PDB_OUTPUT_DIRECTORY):
            shutil.rmtree(PDB_Downloader.PDB_OUTPUT_DIRECTORY)
        if os.path.exists(PDB_Downloader.PDBX_OUTPUT_DIRECTORY):
            shutil.rmtree(PDB_Downloader.PDBX_OUTPUT_DIRECTORY)
        os.makedirs(PDB_Downloader.PDB_OUTPUT_DIRECTORY)
        os.makedirs(PDB_Downloader.PDBx_OUTPUT_DIRECTORY)
        
        # Remake text file of PDB IDs that didn't get downloaded on each run
        failed_file_path = os.path.join(
            PDB_Downloader.INPUT_DIRECTORY, 'failed_ids.txt'
        )
        with open(failed_file_path, 'w') as file:
            pass
        
        # Attempt downloads for all PDB IDs in the text file
        with open(PDB_Downloader.INPUT_FILE, 'r') as pdb_list:
            with open(failed_file_path, 'a') as fails:
                for line in pdb_list:
                    pdb_id = line.strip()
                    # If download failed
                    if not PDB_Downloader.download_pdb(pdb_id):
                        fails.write(f'{pdb_id}\n')
                        
                    # Delay between requests
                    time.sleep(PDB_Downloader.CRAWL_DELAY)
          
        # TODO move this to the PDBx script
        # Create a text file with all PDB IDs that were downloaded as a PDB
        # PDB file, but sorted by decreasing size (will be used by
        # Annotation_Downloader)
        pdb_files = sorted(
            os.listdir(PDB_Downloader.PDB_OUTPUT_DIRECTORY),
            key=lambda f: os.path.getsize(
                os.path.join(PDB_Downloader.PDB_OUTPUT_DIRECTORY, f)
            )
        )
        # Write the sorted list
        with open(
                os.path.join(
                    PDB_Downloader.INPUT_DIRECTORY, 'used_pdb_ids.txt'
                ), 'w'
        ) as sorted_file:
            for pdb_file in pdb_files:
                sorted_file.write(f'{pdb_file.split('.')[0]}\n')
                
    # INPUT: a string representing the ID of a single PDB
    # ACTION: Stores the downloaded PDB file into a local directory, or stores
    # the downloaded PDBx file into another local directory
    # OUTPUT: 'True' if the download was successful, 'False' otherwise
    @staticmethod
    def download_pdb(pdb_id):
        # Attempts to download the file for a PDB ID. Up to 3 attempts
        for attempt in [1, 2, 3]:
            try:
                # Attempt to download the PDB ID as a PDB file
                url = os.path.join(PDB_Downloader.RCSB_URL, f'{pdb_id}.pdb')
                response = PDB_Downloader.SESSION.get(url)
                if response.status_code == 200: # Successfully downloaded
                    with open(
                            os.path.join(
                                PDB_Downloader.PDB_OUTPUT_DIRECTORY, 
                                f'{pdb_id}.pdb'
                            ),
                            'wb'
                    ) as pdb_file:
                        pdb_file.write(response.content)
                    return True
            
                # Attempt to download the PDB ID as a PDBx file instead
                url = os.path.join(PDB_Downloader.RCSB_URL, f'{pdb_id}.cif')
                response = PDB_Downloader.SESSION.get(url)
                if response.status_code == 200: # Successfully downloaded
                    with open(
                            os.path.join(
                                PDB_Downloader.PDBX_OUTPUT_DIRECTORY, 
                                f'{pdb_id}.cif'
                            ),
                            'wb'
                    ) as pdbx_file:
                        pdbx_file.write(response.content)
                    return True
                
                raise Exception('Could not access .pdb or .cif')
            
            except Exception:
                if attempt == 3: # If is final attempt
                    return False
                else:
                    time.sleep(PDB_Downloader.CRAWL_DELAY)