from DataManagement.RawDataLoading.Cluster_Downloader import Cluster_Downloader
from DataManagement.RawDataLoading.PDB_Downloader import PDB_Downloader

# =============================================================================
# Prompts user  for which step of the pipeline to begin at
print(
    '''
    --- BASE PAIRING PROJECT PIPELINE ---
    STEP 1: Retrieve motif cluster table data for each cluster in R3DMA
    STEP 2: Download PDB or PDBx files from RCSB for each PDB_ID in R3DMA
    STEP 3: Convert any downloaded PDBx files to the PDB file format
    STEP 4: Push all PDB files to ClaRNA for annotation
    STEP 5: (MANUAL STEP, Dr. Islam) Push all PDB files to DSSR for annotation
    STEP 6: Preprocess the annotation data
    STEP 7: Extend the motif cluster tables originally from R3DMA
    *more to come...
    -------------------------------------
    
    '''
)
while True:
    starting_step_number = input(
        'Enter a number representing which pipeline step to begin at: '
    ).strip()
    if starting_step_number.isdigit():
        starting_step_number = int(starting_step_number)
        break
    print('TypeError -- Please enter a number: ')
# =============================================================================

# =============================================================================
# STEP 1: Retrieve motif cluster table data for each cluster in R3DMA

if 1 >= starting_step_number:
    print('Beginning STEP 1:\n')
    Cluster_Downloader.download()
    
    # If some were not downloaded retry downloading those:
    while Cluster_Downloader.failed_links:
        response = input(
            'Some motif clusters were not downloaded. Can be found in ' + 
            'Data/Raw/R3DMA/failed_links.txt. Retry these? (y/n) '
        ).strip()
        
        if response == 'Y' or response == 'y':
            Cluster_Downloader.retry_download()
        elif response == 'N' or response == 'n':
            break
        else:
            print('Incorrect input, please enter a single letter.')
# =============================================================================

# =============================================================================
# STEP 2: Download PDB or PDBx files from RCSB for each PDB_ID in R3DMA

if 2 >= starting_step_number:
    print('Beginning STEP 2:\n')
    PDB_Downloader.download()
    
    # If some PDBs were not downloaded as PDB files or PDBx files
    while PDB_Downloader.failed_ids:
        response = input(
            'Some PDB IDs were not downloaded as either PDB or PDBx files ' +
            'Can be found in Data/Raw/RCSB/failed_ids.txt. Retry these? (y/n) ' 
        ).strip()
        
        if response == 'Y' or response == 'y':
            PDB_Downloader.retry_download()
        elif response == 'N' or response == 'n':
            break
        else:
            print('Incorrect input, please enter a single letter.')
# =============================================================================
# STEP 3: Convert any downloaded PDBx files to the PDB file format

if 3 >= starting_step_number:
    print('Beginning STEP 3:\n')
