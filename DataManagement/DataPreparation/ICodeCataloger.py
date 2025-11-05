import os
from Bio.PDB import PDBParser, Model

class ICodeCataloger:
    PDB_DIRECTORY = os.path.join(
        '..', '..', 'Data', 'Raw', 'RCSB', 'PDB_Files'
    )
    
    @staticmethod
    def run() -> None:
        residue_ids = set()
        parser = PDBParser(QUIET=True)
        for pdb in os.listdir(ICodeCataloger.PDB_DIRECTORY):
            pdb_model = ICodeCataloger._load_pdb(pdb, parser)
            ICodeCataloger._add_residues(residue_ids, pdb_model, pdb)
        print(residue_ids)
    
    #==========================================================================
    # Helper Methods Below
    #==========================================================================

    @staticmethod
    def _load_pdb(pdb: str, parser: PDBParser) -> Model:
        pdb_id = os.path.splitext(pdb)[0]
        structure = parser.get_structure(
            pdb_id, os.path.join(
                ICodeCataloger.PDB_DIRECTORY, pdb
            )
        )
        return structure[0]
    
    @staticmethod
    def _add_residues(residue_ids: set[str], pdb_model: Model, pdb: str) -> None:
        nucleotides = ['A', 'G', 'U', 'C'] #only standard nucleotides
        pdb_id = os.path.splitext(pdb)[0]
        for chain in pdb_model:
            for residue in chain:
                _, sequence_number, insertion_code = residue.id
                residue_name = residue.get_resname()
                if residue_name in nucleotides and insertion_code != ' ':
                    residue_ids.add(f'{pdb_id}_{chain.id}_{sequence_number}')
    
    # load a pdb, go through all residue IDs and store the ones containing
    # an I-code into a set().
    # the information stored in the set should be {pdb}{chain}{residue_id_without_i_code}
    # so remove the icode when stored in this set.
    # Once we have this set, we can remove certain entries from the JSON files
    
ICodeCataloger.run()