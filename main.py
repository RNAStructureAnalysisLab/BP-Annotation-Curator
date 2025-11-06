# AUTHOR: Kristopher Chruch
# USE: This file represents the entire pipeline of the base pairing project.
# TODO finish description of use

# TODO: Within the R3DMA directory, create a text file 'current_version" which
# has the version, instead of it being in the name of the cluster directory
# =============================================================================

import sys
from DataManagement.RawDataLoading.Cluster_Downloader import Cluster_Downloader
from DataManagement.RawDataLoading.PDB_Downloader import PDB_Downloader
from DataManagement.DataPreparation.PDB_Maker import PDB_Maker
from DataManagement.AnnotationLoading.Annotation_Downloader import Annotation_Downloader
from DataManagement.DataPreparation.Chain_Restorer import Chain_Restorer
from DataManagement.DataExploration.Explorer import Explorer
from DataManagement.DataPreparation.Preprocessor import Preprocessor
from DataManagement.DataPreparation.Table_Extender import Table_Extender
from DataManagement.ConsensusMaking.Tool_Consensus import Tool_Consensus
from DataManagement.ResultsMaking.Agreement_Calculator import Agreement_Calculator
from DataManagement.ResultsMaking.Agreement_Analyzer import Agreement_Analyzer
from DataManagement.DataExploration.Catalogue import Catalogue
from DataManagement.DataExploration.RC_Converter import RC_Converter
from DataManagement.ConsensusMaking.Tool_Consensus_2 import Tool_Consensus_2
from DataManagement.DataExploration.Rejected_Types_Finder import Rejected_Types_Finder
from DataManagement.ResultsMaking.Tool_Benchmarker import Tool_Benchmarker
from DataManagement.ResultsMaking.Agreement_Calculator_Cluster_Consensus import Agreement_Calculator_Cluster_Consensus

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
    STEP 6: Perform data exploration on the annotations
    STEP 7: Preprocess the annotation data
    STEP 8: Extend the motif cluster tables originally from R3DMA
    STEP 9: Create tool consensus tables for each motif cluster from STEP 8
    STEP 10: Find the base pair intersections between the extended tables and
        the tools
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




            
# =============================================================================
# STEP 3: Convert any downloaded PDBx files to the PDB file format

if 3 >= starting_step_number:
    print('Beginning STEP 3:\n')
    PDB_Maker.convert_all()
    
# =============================================================================




    
# =============================================================================
# STEP 4: Push all PDB files to ClaRNA for annotation

if 4 >= starting_step_number:
    print('Beginning STEP 4:\n')
    Annotation_Downloader.download()
    
    while Annotation_Downloader.if_failed_pdb_ids():
        response = input(
            'Some PDB files were not annotated. Can be found in ' +
            'Data/Raw/AnnotationTools/failed_pdb_ids.txt. Retry these? (y/n) '
        ).strip()
        
        if response == 'Y' or response == 'y':
            Annotation_Downloader.retry_download()
        elif response == 'N' or response == 'n':
            break
        else:
            print('Incorrect input, please enter a single letter.')
            
# =============================================================================





# =============================================================================
# STEP 5: (MANUAL STEP, Dr. Islam) Push all PDB files to DSSR for annotation

if 5 >= starting_step_number:
    print("Beginning STEP 5:\n")
    response = input(
        "Ensure that DSSR annotations have been manually added into " +
        "'Data/Raw/AnnotationTools/DSSR_Annotations'. Otherwise the " + 
        "following steps in the pipeline will not make use of this data. " + 
        "Proceed (y/n): "
    )
    
    while True:
        if response == "n" or response == "N":
            print(
                "Rerun starting from STEP 5 once you have the data and " +
                "enter 'y' or 'Y' to proceed."
            )
            sys.exit()
        elif response == 'y' or response == "Y":
            Chain_Restorer.restore()
            break
        else:
            print("Incorrect input, please enter a single letter.")

# =============================================================================





# =============================================================================
# STEP 6: Perform data exploration on the annotations

# TODO: in progress, continue adding features as needed
if 6 >= starting_step_number:
    print('Beginning STEP 6:\n')
    Explorer.document_descriptions()
    
# =============================================================================
  




# =============================================================================
# STEP 7: Preprocess the annotation data

if 7 >= starting_step_number:
    print('Beginning STEP 7:\n')
    Preprocessor.convert_all()

# =============================================================================





# =============================================================================
# STEP 8: Extend the motif cluster tables originally from R3DMA

if 8 >= starting_step_number:
    print('Beginning STEP 8:\n')
    Table_Extender.run()

# =============================================================================




    
# =============================================================================
# STEP 9: Create tool consensus tables for each motif cluster from STEP 8
# TODO: incorporate other tool consensus paradigms

if 9 >= starting_step_number:
    print('Beginning STEP 9:\n')
    #Tool_Consensus.run()
    #Consensus_V1.run()

# =============================================================================
# STEP 10: Find the base pair intersections between the extended table and 
# tools annotations

if 10 >= starting_step_number:
    print('Beginning STEP 10:\n')
    #Agreement_Calculator.run()
if 11 >= starting_step_number: # should be part of step 10, temporary
    print('Beginning STEP 11:\n')
    #Agreement_Analyzer.run()
if 12 >= starting_step_number:
    print('Beginning STEP 12:\n')
    RC_Converter.run()
    Tool_Consensus_2.run()
    Rejected_Types_Finder.run()
if 13 >= starting_step_number:
    print('Beginning STEP 13:\n')
    Tool_Benchmarker.run()
    Agreement_Calculator_Cluster_Consensus.run()
    
# TODO: verify that pdb 8CRE is getting processed when we run the entire
# project again from the beginning. Originally is a PDBx file.
# For some reson it didn't end up in my list of files to convert to PDB and as
# such am missing a DSSR file for it. Some others missing too. Could be a 
# problem with one of the files, or a mistake when I was manually trying to 
# insert things (Same for 5J7L)
# NOTE: 4L47 and 4V9F are good ones to verify step 5
# TODO: 1F7U  and 2NZ4 dssr are strange. For example it reports the nucleoside 
# pseudouridine (P). P is not featured in the other annotations. a residue that
# is P seems to be completely omitted. Does R3DMA also omit these from the cluster
# or does it represent P as a 'U' since they are similar? If it represents it
# as a U, then this would explain the occurrence of "nbp" being reported for 
# certain residues, where it was a modified nucleoside

#TODO: IL_85033 has some unusual contact types reported from R3DMA
# NOTE: HL_28252.5 from version 3.95 is a good example of a cluster with
# nonstandard nucleosides
    