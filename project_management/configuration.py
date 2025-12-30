import os

PROJECT_ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#==============================================================================
# Paths for project_management (IE: python scripts)
#==============================================================================
PROJECT_MANAGEMENT_DIR = os.path.join(PROJECT_ROOT_DIR, "project_management")
CONFIGURATION_PY = os.path.join(PROJECT_MANAGEMENT_DIR, "configuration.py")

#==============================================================================
# Paths for state (IE: json)
#==============================================================================
STATE_DIR = os.path.join(PROJECT_ROOT_DIR, "state")

#==============================================================================
# Paths for data_management (IE: python scripts)
#==============================================================================
DATA_MANAGEMENT_DIR = os.path.join(PROJECT_ROOT_DIR, "data_management")
RAW_DATA_LOADING_DIR = os.path.join(DATA_MANAGEMENT_DIR, "raw_data_loading")
ANNOTATED_DATA_LOADING_DIR = os.path.join(RAW_DATA_LOADING_DIR, "annotated_data_loading")

#==============================================================================
# Paths for data (IE: PDBs, CSVs, etc)
#==============================================================================
DATA_DIR = os.path.join(PROJECT_ROOT_DIR, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
R3DMA_DIR = os.path.join(RAW_DIR, "r3dma")
CLUSTERS_DIR = os.path.join(R3DMA_DIR, "clusters")
USED_PDB_IDS_TXT = os.path.join(R3DMA_DIR, "used_pdb_ids.txt")
RCSB_DIR = os.path.join(RAW_DIR, "rcsb")
PDB_DIR = os.path.join(RCSB_DIR, "pdb")
PDBX_DIR = os.path.join(RCSB_DIR, "pdbx")
CONVERTED_PDBX_DIR = os.path.join(RCSB_DIR, "converted_pdbx")

#==============================================================================
# Paths for external data (IE: website links)
#==============================================================================
R3DMA_HOMEPAGE_URL = "https://rna.bgsu/rna3dhub/motifs"