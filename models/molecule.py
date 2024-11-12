from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional

class MoleculeProcessor:
    @staticmethod
    def smiles_to_3d(smiles: str) -> Optional[Chem.Mol]:
        """Convert SMILES to 3D structure using RDKit"""
        try:
            # Create RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D conformation
            success = AllChem.EmbedMolecule(mol, randomSeed=42)
            if success == -1:
                return None
                
            # Optimize the structure
            AllChem.MMFFOptimizeMolecule(mol)
            
            return mol
            
        except Exception:
            return None 