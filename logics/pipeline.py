from rdkit import Chem
import os
from monomer_library import MonomerLibrary, MonomerData
from target_processor import PrecisionFragmentProcessor, PrecisionMonomerMatcher, HELMGenerator

# Global variables for caching
_MONOMER_LIBRARY = None
_PROCESSOR = None
_MATCHER = None
_HELM_GENERATOR = None


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


def _get_processors():
    """
    Get or create singleton instances of processors.
    Returns tuple: (processor, matcher, helm_generator)
    """
    global _PROCESSOR, _MATCHER, _HELM_GENERATOR
    
    if _PROCESSOR is None or _MATCHER is None or _HELM_GENERATOR is None:
        library = _load_monomer_library()
        if not library:
            return None, None, None
        
        _PROCESSOR = PrecisionFragmentProcessor(library)
        _MATCHER = PrecisionMonomerMatcher(library)
        _HELM_GENERATOR = HELMGenerator()
    
    return _PROCESSOR, _MATCHER, _HELM_GENERATOR


def preload_library():
    """
    Preload the monomer library and initialize processors once at the start.
    Returns True if successful, False otherwise.
    """
    library = _load_monomer_library()
    if library is None:
        return False
    
    # Initialize processors
    processor, matcher, generator = _get_processors()
    return processor is not None


def mol_to_helm(mol: Chem.Mol, is_cyclic: bool = False) -> str:
    """
    Main function to convert RDKit molecule to HELM notation.

    Args:
        mol: RDKit molecule
        is_cyclic: Cyclic structure flag

    Returns:
        HELM notation as a string
    """
    # Use shared processor instances
    processor, matcher, helm_generator = _get_processors()
    if not processor:
        return ""

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


def convert_molecules_batch(molfiles: list, is_cyclic: bool = False) -> list:
    """
    Convert a batch of molecules from molfile format to HELM notation.
    
    Args:
        molfiles: List of molfile strings
        is_cyclic: Whether the peptides are cyclic (default: False)
    
    Returns:
        List of tuples: (success: bool, helm_notation: str)
        success is True if molecule was successfully converted, False otherwise
    """
    # Ensure library and processors are initialized before batch processing
    global _PROCESSOR
    if _PROCESSOR is None:
        print("Initializing monomer library and processors...")
        if not preload_library():
            print("ERROR: Failed to load monomer library")
            return [(False, "Library initialization failed") for _ in molfiles]
        print()
    
    # Use shared processor instances
    processor, matcher, helm_generator = _get_processors()
    if not processor:
        return [(False, "") for _ in molfiles]
    
    results = []
    
    for molfile in molfiles:
        mol = Chem.MolFromMolBlock(molfile)
        if not mol:
            results.append((False, ""))
            continue
        
        try:
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
                results.append((True, helm_notation))
            else:
                results.append((False, ""))
        except Exception as e:
            results.append((False, f"Error: {str(e)}"))
    
    return results