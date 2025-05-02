# AUTHORS: Jackie Luc and Kristopher Church

import os
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, Select

# Looks in within INPUT_DIRECTORY and converts all PDBx files into PDB. Where
# they are then stored into PDB_DIRECTORY
class PDB_Maker(Select): # Inherits from "Select" found in Bio.PDB
    INPUT_DIRECTORY = os.path.join("Data", "Raw", "RCSB", "PDBx_Files")
    PDB_DIRECTORY = os.path.join("Data", "Raw", "RCSB", "PDB_Files")
    RCSB_DIRECTORY = os.path.join("Data", "Raw", "RCSB")
    CLUSTER_TABLE_DIRECTORY = os.path.join(
        "Data", "Raw", "R3DMA", "cluster_tables_3.95"
    ) # TODO there is an issue when the version is updated
    RNA_NUCLEOTIDES = {"A", "U", "C", "G"}
    VALID_PDB_CHAIN_IDS = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"
    )
    PARSER = MMCIFParser(QUIET=True)
    
    # Dictionary that maps a PDB ID to another dictionary. This other
    # dictionary maps a single character chain to the original chain in the PDB
    # NOTE: PDBx files can have multi-character chain names, this must be
    # made single character for the conversion. The mappings are saved so
    # that the original chain IDs can be restored later in the pipeline
    to_original_chain = {}
    # Dictionary with the original chain ID as the key, and the new chain ID as
    # the value
    to_new_chain = {}
    chains_per_pdb = {} # PDB IDs as keys and values are sets of chain IDs
    
    # ACTION: Iterates through each PDBx file in INPUT_DIRECTORY, and converts
    # to a PDB file. It then finally stores it within PDB_DIRECTORY
    @staticmethod
    def convert_all():
        PDB_Maker.load_chains_per_pdb()
        PDB_Maker.load_chain_mappings()
        successfully_converted = [] # TODO fully implement why these variables are here
        unsuccessfully_converted = []
        
        # Loop for getting each PDBx file name
        for file_name in os.listdir(PDB_Maker.INPUT_DIRECTORY):
            pdb_id = os.path.splitext(file_name)[0]
            input_file_path = os.path.join(
                PDB_Maker.INPUT_DIRECTORY, file_name
            )
            output_file_path = os.path.join(
                PDB_Maker.PDB_DIRECTORY, pdb_id + ".pdb"
            )
            
            try:
                structure = PDB_Maker.PARSER.get_structure(pdb_id, input_file_path)
                PDB_Maker.to_new_chain = {
                    v:k for k, v in PDB_Maker.to_original_chain[structure.id].items()
                }
                
                # Remove unwanted chains, for reducing file size to allow PDB
                for model in structure:
                    for chain in list(model):
                        if chain.id not in PDB_Maker.to_new_chain:
                            model.detach_child(chain.id)
                            
                # Update chain IDs
                for model in structure:
                    for chain in model:
                        if chain.id in PDB_Maker.to_new_chain:
                            chain.id = PDB_Maker.to_new_chain[chain.id]
                            
                # Insert the converted PDBx files into the PDB directory
                io = PDBIO()
                io.set_structure(structure)
                io.save(output_file_path, PDB_Maker())
                successfully_converted.append(pdb_id)
            except Exception as e:
                unsuccessfully_converted.append((pdb_id, str(e)))
                print(str(e))
            
     
    @staticmethod 
    def load_chains_per_pdb():
        for cluster_csv in os.listdir(PDB_Maker.CLUSTER_TABLE_DIRECTORY):
            csv_path = os.path.join(
                PDB_Maker.CLUSTER_TABLE_DIRECTORY, cluster_csv
            )
            cluster = pd.read_csv(csv_path)
            
            for _, row in cluster.iterrows():
                pdb = row["PDB"]
                chains = str(row["Chain(s)"])
                split_chains = chains.split("+")
                for chain in split_chains:
                    if pdb in PDB_Maker.chains_per_pdb:
                        PDB_Maker.chains_per_pdb[pdb].add(chain)
                    else:
                        PDB_Maker.chains_per_pdb[pdb] = {chain}
     
    @staticmethod 
    def load_chain_mappings():
        for pdb_id, chains in PDB_Maker.chains_per_pdb.items():
            count = 0
            temp = {}
            for chain in chains:
                chain = str(chain)
                if len(chain) == 1: # If the chain is a single letter
                    temp[chain] = chain # Doesn't require remapping
                elif len(chain) > 1:
                    while PDB_Maker.VALID_PDB_CHAIN_IDS[count] in [
                            key for key, _ in temp.items()
                    ]: # TODO is this handling overflow if count is too big?
                        count += 1
                    temp[PDB_Maker.VALID_PDB_CHAIN_IDS[count]] = chain
            # Raise exception if the same pdb_id is found again
            if pdb_id in PDB_Maker.to_original_chain.keys():
                raise ValueError(f"Duplicate pdb_id encountered: {pdb_id}")
            PDB_Maker.to_original_chain[pdb_id] = temp
        
        # Store to_original_chain for later use in the pipeline
        output_path = os.path.join(
            PDB_Maker.RCSB_DIRECTORY, "original_pdbx_chains.txt"
        )
        with open(output_path, 'w') as file:
            for pdb_id, chain_mappings in PDB_Maker.to_original_chain.items():
                file.write(f"{pdb_id}: {chain_mappings}\n")
                
    # -------------------------------------------------------------------------
    # NOTE: Below are methods inherited from Select that are being overloaded
    # -------------------------------------------------------------------------
    
    # INPUT: A biopython residue object, and a dictionary with the original
    # chain ID as the key, and the new chain ID as the value
    # ACTION: Modifies the chain ID of the residue object to a new one
    # OUTPUT: a boolean that is True when the nucleotide type of the residue
    # is valid
    def accept_residue(self, residue):
        structure_id, model_id, chain_id, residue_id = residue.get_full_id()
        
        if chain_id in PDB_Maker.to_original_chain[structure_id].values():
            # Modify original chain
            residue.parent.id = PDB_Maker.to_new_chain[chain_id]
        return residue.get_resname() in PDB_Maker.RNA_NUCLEOTIDES