from rdkit import Chem
from monomer_library import MonomerLibrary
from fragment_processor import TerminalNormalizer


class MonomerMatcher:
    """
    Matches molecular fragments to monomers in the library.
    
    Uses exact matching strategies:
    - Exact SMILES matching against library
    - Normalization attempts for terminal variations
    """
    
    def __init__(self, monomer_library: MonomerLibrary):
        self.monomer_library = monomer_library

    def find_exact_match(self, fragment: Chem.Mol):
        """
        Find exact match for a fragment.
        
        Args:
            fragment: RDKit molecule object representing a fragment
        
        Returns:
            MonomerData object if match found, None otherwise
        """
        try:
            frag_smiles = Chem.MolToSmiles(fragment, canonical=True)
            if not frag_smiles:
                return None

            # Try exact SMILES match
            exact_match = self.monomer_library.find_monomer_by_smiles(frag_smiles)
            if exact_match:
                return exact_match

            # Try with normalization
            normalized_match = self._find_with_normalization(fragment, frag_smiles)
            if normalized_match:
                return normalized_match

            return None

        except Exception:
            return None

    def _find_with_normalization(self, fragment: Chem.Mol, frag_smiles: str):
        """Find monomer by trying different normalization strategies"""
        try:
            normalizer = TerminalNormalizer()

            normalized1 = normalizer.normalize_for_library(fragment, is_c_terminal=True)
            if normalized1:
                norm_smiles1 = Chem.MolToSmiles(normalized1, canonical=True)
                match1 = self.monomer_library.find_monomer_by_smiles(norm_smiles1)
                if match1:
                    return match1

            normalized2 = normalizer.normalize_for_library(fragment, is_c_terminal=False)
            if normalized2:
                norm_smiles2 = Chem.MolToSmiles(normalized2, canonical=True)
                match2 = self.monomer_library.find_monomer_by_smiles(norm_smiles2)
                if match2:
                    return match2

            return None

        except Exception:
            return None


