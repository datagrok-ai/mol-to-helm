from rdkit import Chem
import os
from monomer_library import MonomerLibrary, MonomerData
from target_processor import PrecisionFragmentProcessor, PrecisionMonomerMatcher, HELMGenerator

# Global variable for caching the library
_MONOMER_LIBRARY = None


def _load_monomer_library():
    global _MONOMER_LIBRARY
    if _MONOMER_LIBRARY is None:
        # Define path to library relative to current directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        library_path = os.path.join(project_root, "libraries", "HELMCoreLibrary.json")

        if not os.path.exists(library_path):
            return None

        print("Loading monomer library...")
        _MONOMER_LIBRARY = MonomerLibrary()
        _MONOMER_LIBRARY.load_from_helm_json(library_path)

        if not _MONOMER_LIBRARY.monomers:
            return None
        
        print(f"Monomer library loaded: {len(_MONOMER_LIBRARY.monomers)} monomers")

    return _MONOMER_LIBRARY


def preload_library():
    """
    Preload the monomer library to ensure it's loaded once at the start.
    Returns True if successful, False otherwise.
    """
    library = _load_monomer_library()
    return library is not None


def mol_to_helm(mol: Chem.Mol, is_cyclic: bool = False) -> str:
    """
    Main function to convert RDKit molecule to HELM notation.

    Args:
        mol: RDKit molecule
        is_cyclic: Cyclic structure flag

    Returns:
        HELM notation as a string
    """
    library = _load_monomer_library()
    if not library:
        return ""

    processor = PrecisionFragmentProcessor(library)
    matcher = PrecisionMonomerMatcher(library)
    helm_generator = HELMGenerator()

    fragments = processor.process_molecule(mol, is_cyclic=is_cyclic)

    matched_monomers = []
    unknown_count = 0

    for fragment in fragments:
        monomer = matcher.find_best_match(fragment)
        if monomer:
            matched_monomers.append(monomer)
        else:
            unknown_count += 1
            mock_monomer = MonomerData()
            mock_monomer.symbol = f"X{unknown_count}"
            mock_monomer.name = f"Unknown_{unknown_count}"
            matched_monomers.append(mock_monomer)

    if matched_monomers:
        helm_notation = helm_generator.generate_helm_notation(matched_monomers, is_cyclic=is_cyclic)
        return helm_notation
    else:
        return ""