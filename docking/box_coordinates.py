from Bio import PDB
import numpy as np

def get_box_coordinates(pdb_file, residues, box_size):
    """
    Calculate the docking box coordinates based on specified residues.
    
    Args:
    pdb_file (str): Path to the cleaned PDB file.
    residues (list): List of residue names and their positions (e.g., ['TYR652', 'PHE656']).
    box_size (float): Size of the box in Angstroms.
    
    Returns:
    dict: Center coordinates and box dimensions.
    """
    # Initialize parser
    parser = PDB.PDBParser()
    structure = parser.get_structure('7CN1', pdb_file)
    
    # Collect coordinates of specified residues
    coords = []
    for model in structure:
        for chain in model:
            for residue in residues:
                res_name = residue[:3]  # Get the residue name (e.g., 'TYR')
                res_id = int(residue[3:])  # Get the residue number (e.g., 652)
                try:
                    res = chain[(' ', res_id, ' ')]  # Using the position ID tuple (chain ID, residue number, insertion code)
                    if res.get_resname() == res_name:  # Verify it's the correct residue type
                        # Get coordinates of the alpha carbon (CA) atom
                        ca_vector = res['CA'].get_vector()
                        coords.append(np.array([ca_vector[0], ca_vector[1], ca_vector[2]]))
                except KeyError:
                    print(f"Residue {residue} not found in the structure.")
    
    # Calculate the center of the box
    if coords:
        center = np.mean(coords, axis=0)
    else:
        raise ValueError("No valid coordinates found for the specified residues.")
    
    # Define box dimensions
    box_coordinates = {
        'center_x': center[0],
        'center_y': center[1],
        'center_z': center[2],
        'size_x': box_size,
        'size_y': box_size,
        'size_z': box_size
    }
    
    return box_coordinates

if __name__ == "__main__":
    residues = ['TYR652', 'PHE656', 'THR623', 'SER624', 'VAL625', 'GLY648', 'PHE557']
    box_size = 20.0  # Define the size of the box
    box_coords = get_box_coordinates('./proteins/7cn1_cleaned.pdb', residues, box_size)
    print("Docking Box Coordinates:", box_coords)