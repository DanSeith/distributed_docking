from Bio import PDB
from Bio.PDB import *
import os
from subprocess import run

def prepare_7cn1():
    """
    Prepare the 7CN1 hERG structure for docking
    """
    # Download structure
    pdbl = PDB.PDBList()
    pdbl.retrieve_pdb_file('7CN1', pdir='./proteins', file_format='pdb')
    
    # Initialize parser
    parser = PDB.PDBParser()
    structure = parser.get_structure('7CN1', './proteins/pdb7cn1.ent')
    
    # Clean structure
    # Remove water and non-essential ions (keeping K+ in selectivity filter)
    for model in structure:
        for chain in model:
            for residue in list(chain):
                if residue.id[0] != ' ':  # Remove heteroatoms
                    if residue.resname != 'K':  # Keep potassium ions
                        chain.detach_child(residue.id)
    
    # Save cleaned structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save('./proteins/7cn1_cleaned.pdb')
    
    # Convert to PDBQT using AutoDock Tools
    prepare_receptor_script = """
        from ADT.MolKit import Read
        from ADT.AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
        
        receptor = Read('./proteins/7cn1_cleaned.pdb')[0]
        receptor_prep = AD4ReceptorPreparation(receptor)
        receptor_prep.cleanup()
    """
    
    with open('prepare_receptor.py', 'w') as f:
        f.write(prepare_receptor_script)
    
    # Run preparation script
    run(['python', 'prepare_receptor.py'])

if __name__ == "__main__":
    prepare_7cn1()