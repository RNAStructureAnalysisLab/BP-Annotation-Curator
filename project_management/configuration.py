import os

class Configuration():
    def __init__(self, version: str):
        project_root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.version = version
        
        #======================================================================
        # Paths for state (manifests)
        #======================================================================
        state_dir = os.path.join(project_root_dir, "state", version)
        self.cluster_downloader_json = os.path.join(state_dir, "cluster_downloader.json")
        self.pdb_downloader_json = os.path.join(state_dir, "pdb_downloader.json")
        
        #======================================================================
        # Paths for data
        #======================================================================
        data_dir = os.path.join(project_root_dir, "data")
        raw_dir = os.path.join(data_dir, "raw")
        
        # R3DMA data (version-dependent)
        r3dma_dir = os.path.join(raw_dir, "r3dma", version)
        self.clusters_dir = os.path.join(r3dma_dir, "clusters")
        self.used_pdb_ids_txt = os.path.join(r3dma_dir, "used_pdb_ids.txt")
        
        # RCSB data (version-dependent)
        rcsb_dir = os.path.join(raw_dir, "rcsb", version)
        self.pdb_dir = os.path.join(rcsb_dir, "pdb")
        self.pdbx_dir = os.path.join(rcsb_dir, "pdbx")
        self.converted_pdbx_dir = os.path.join(rcsb_dir, "converted_pdbx")
        
        #======================================================================
        # Paths to external data
        #======================================================================
        self.r3dma_homepage_url = "https://rna.bgsu.edu/rna3dhub/motifs"
        
        #======================================================================
        # Variables for cluster_downloader.py (step 1)
        #======================================================================
        self.crawl_delay_1 = 10
        self.max_attempts_1 = 3
        self.failed_links_1 = set()