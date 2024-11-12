import os
import subprocess
import tempfile
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from docking.box_coordinates import get_box_coordinates

def perform_docking(smiles, prepared_protein_path):
    """
    Perform docking of a compound against a prepared protein using AutoDock Vina.
    Args:
    smiles (str): SMILES representation of the ligand.
    prepared_protein_path (str): File path to the prepared protein (PDBQT format).
    Returns:
    dict: Docking results including binding affinity and pose details.
    """
    residues = ['TYR652', 'PHE656', 'THR623', 'SER624', 'VAL625', 'GLY648', 'PHE557']
    box_size = 20.0  # Define the size of the box
    
    try:
        ligand_pdbqt = convert_smiles_to_pdbqt(smiles)
    except Exception as e:
        logging.error(f"Docking failed: {str(e)}")
        raise ValueError(f"Docking failed: {str(e)}")
    if not ligand_pdbqt:
        raise ValueError("Failed to convert SMILES to PDBQT.")
    
    # Get box coordinates with error handling
    try:
        box_coords = get_box_coordinates('./proteins/7cn1_cleaned.pdb', residues, box_size)
    except Exception as e:
        raise RuntimeError(f"Failed to get box coordinates: {e}")
        
    with tempfile.TemporaryDirectory() as tmpdirname:
        dock_output_path = os.path.join(tmpdirname, "docking_result.pdbqt")
        log_output_path = os.path.join(tmpdirname, "vina.log")
        
        # Define the Vina configuration with box coordinates
        vina_config = {
            'receptor': prepared_protein_path,
            'ligand': ligand_pdbqt,
            'out': dock_output_path,
            'log': log_output_path,
            'config': {
                'center_x': box_coords['center_x'],
                'center_y': box_coords['center_y'],
                'center_z': box_coords['center_z'],
                'size_x': box_coords['size_x'],
                'size_y': box_coords['size_y'],
                'size_z': box_coords['size_z'],
                'cpu': 4,
                'exhaustiveness': 8
            }
        }
        # Run Vina docking
        success = run_vina(vina_config=vina_config)
        if not success:
            raise RuntimeError("Vina docking failed.")
        # Parse the results
        binding_affinity, poses = parse_vina_output(log_output_path, dock_output_path)
        return {
            'binding_affinity': binding_affinity,
            'poses': poses
        }

def convert_smiles_to_pdbqt(smiles: str) -> str:
    try:
        # Add logging to debug the process
        logging.info(f"Converting SMILES: {smiles}")
        
        # Check if required environment is set up
        temp_dir = os.environ.get('TMPDIR', '/tmp')
        logging.info(f"Using temporary directory: {temp_dir}")
        
        # Ensure temp directory is writable
        if not os.access(temp_dir, os.W_OK):
            raise PermissionError(f"Temporary directory {temp_dir} is not writable")
        
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        # Generate 3D coordinates
        mol = Chem.AddHs(mol)  # Add hydrogens
        AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D conformation
        AllChem.MMFFOptimizeMolecule(mol)  # Energy minimization
        
        # Create temporary PDB file
        temp_pdb = os.path.join(temp_dir, "temp_ligand.pdb")
        temp_pdbqt = os.path.join(temp_dir, "temp_ligand.pdbqt")
        
        # Write molecule to PDB file
        Chem.MolToPDBFile(mol, temp_pdb)
        
        # Convert PDB to PDBQT using Open Babel
        cmd = ["obabel", temp_pdb, "-O", temp_pdbqt, "-xh"]
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Read the PDBQT content
        with open(temp_pdbqt, 'r') as f:
            pdbqt_content = f.read()
        
        logging.info(f"Conversion successful. Output file: {temp_pdbqt}")
        return pdbqt_content
        
    except Exception as e:
        logging.error(f"SMILES conversion failed: {str(e)}")
        raise ValueError(f"Failed to convert SMILES to PDBQT: {str(e)}")

def run_vina(vina_config):
    """
    Execute AutoDock Vina with the provided configuration.
    Args:
    vina_config (dict): Configuration parameters for Vina.
    Returns:
    bool: True if docking is successful, False otherwise.
    """
    try:
        cmd = [
            "vina",
            "--receptor", vina_config['receptor'],
            "--ligand", vina_config['ligand'],
            "--out", vina_config['out'],
            "--log", vina_config['log'],
            "--centroid",
            f"{vina_config['config']['center_x']},{vina_config['config']['center_y']},{vina_config['config']['center_z']}",
            "--size",
            f"{vina_config['config']['size_x']},{vina_config['config']['size_y']},{vina_config['config']['size_z']}",
            "--cpu", str(vina_config['config']['cpu']),
            "--exhaustiveness", str(vina_config['config']['exhaustiveness'])
        ]
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Vina execution failed: {e}")
        return False

def parse_vina_output(log_path, dock_output_path):
    """
    Parse Vina log and output files to extract docking results.
    Args:
    log_path (str): Path to Vina log file.
    dock_output_path (str): Path to Vina output PDBQT file.
    Returns:
    tuple: Binding affinity (float) and list of poses.
    """
    binding_affinity = None
    poses = []
    # Parse log file for binding affinity
    try:
        with open(log_path, 'r') as log_file:
            for line in log_file:
                if line.startswith("Affinity:"):
                    parts = line.strip().split()
                    binding_affinity = float(parts[1])
                    break
    except Exception as e:
        print(f"Error parsing Vina log file: {e}")
    # Parse docked poses from PDBQT file
    try:
        with open(dock_output_path, 'r') as dock_file:
            current_pose = []
            for line in dock_file:
                if line.startswith("MODEL"):
                    current_pose = []
                elif line.startswith("ENDMDL"):
                    poses.append('\n'.join(current_pose))
                else:
                    current_pose.append(line.strip())
    except Exception as e:
        print(f"Error parsing Vina output file: {e}")
    return binding_affinity, poses
