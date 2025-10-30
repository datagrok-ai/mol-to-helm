from rdkit import Chem
import os
from monomer_library import MonomerLibrary, MonomerData
from fragment_processor import FragmentProcessor
from monomer_matcher import MonomerMatcher
from helm_generator import HELMGenerator

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
        
        _PROCESSOR = FragmentProcessor(library)
        _MATCHER = MonomerMatcher(library)
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


def convert_molecules_batch(molfiles: list, library_path: str = None) -> list:
    """
    Convert a batch of molecules from molfile format to HELM notation.
    
    Args:
        molfiles: List of molfile strings
        library_path: Optional path to monomer library JSON file (as string).
                     If None, uses default cached library from HELMCoreLibrary.json
    
    Returns:
        List of tuples: (success: bool, helm_notation: str)
        success is True if molecule was successfully converted, False otherwise
    """
    # Determine which library to use
    if library_path is None:
        # Use cached global library
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
    else:
        # Load custom library from provided path (no caching)
        if not os.path.exists(library_path):
            print(f"ERROR: Library file not found: {library_path}")
            return [(False, "Library file not found") for _ in molfiles]
        
        print(f"Loading custom library from: {library_path}")
        library = MonomerLibrary()
        library.load_from_helm_json(library_path)
        
        if not library.monomers:
            print("ERROR: No monomers loaded from custom library")
            return [(False, "Library loading failed") for _ in molfiles]
        
        print(f"Custom library loaded: {len(library.monomers)} monomers")
        
        # Create processor instances for this library
        processor = FragmentProcessor(library)
        matcher = MonomerMatcher(library)
        helm_generator = HELMGenerator()
    
    results = []
    
    for i in range(len(molfiles)):
        molfile = molfiles[i]
        mol = Chem.MolFromMolBlock(molfile)
        if not mol:
            results.append((False, ""))
            continue
        
        try:
            # Process molecule into fragment graph
            graph = processor.process_molecule(mol)
            
            # Match each fragment to a monomer using graph topology
            unknown_count = 0
            for node_id, node in graph.nodes.items():
                # Count connections for this node
                neighbors = graph.get_neighbors(node_id)
                num_connections = len(neighbors)
                
                # Find matching monomer
                monomer = matcher.find_exact_match(node.mol, num_connections)
                if monomer:
                    node.monomer = monomer
                else:
                    unknown_count += 1
                    mock_monomer = MonomerData()
                    mock_monomer.symbol = f"X{unknown_count}"
                    mock_monomer.name = f"Unknown_{unknown_count}"
                    node.monomer = mock_monomer
            
            # Try to recover unmatched fragments by merging with neighbors
            max_recovery_attempts = 3  # Prevent infinite loops
            for attempt in range(max_recovery_attempts):
                had_changes = processor.recover_unmatched_fragments(graph, matcher)
                if not had_changes:
                    break
            
            if len(graph.nodes) > 0:
                helm_notation = helm_generator.generate_helm_from_graph(graph)
                results.append((True, helm_notation))
            else:
                results.append((False, ""))
        except Exception as e:
            results.append((False, f"Error: {str(e)}"))
    
    return results