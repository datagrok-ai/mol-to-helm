from rdkit import Chem
from monomer_library import MonomerLibrary
from fragment_graph import FragmentGraph, FragmentNode


class MonomerMatcher:
    """
    Matches molecular fragments to monomers using graph-aware R-group analysis.
    
    Revolutionary approach:
    - No hardcoded mappings
    - No complex normalization
    - Direct string comparison of canonical SMILES
    - Graph topology determines which R-groups are capped
    """
    
    def __init__(self, monomer_library: MonomerLibrary):
        self.monomer_library = monomer_library

    def find_exact_match(self, fragment: Chem.Mol, num_connections: int = 0):
        """
        Find exact match for a fragment based on graph topology.
        
        Args:
            fragment: RDKit molecule object representing a fragment
            num_connections: Number of connections this fragment has in the graph
        
        Returns:
            MonomerData object if match found, None otherwise
        """
        try:
            # Get canonical SMILES of the fragment
            frag_smiles = Chem.MolToSmiles(fragment, canonical=True)
            if not frag_smiles:
                return None

            # Use the library's new graph-aware matching
            match = self.monomer_library.find_monomer_by_fragment_smiles(
                frag_smiles, num_connections
            )
            
            return match

        except Exception:
            return None
    
    def match_graph(self, graph: FragmentGraph):
        """
        Match all fragments in a graph to monomers.
        
        Args:
            graph: FragmentGraph with unmatched nodes
        
        Returns:
            Number of successfully matched nodes
        """
        matched_count = 0
        
        for node_id, node in graph.nodes.items():
            # Count connections for this node
            neighbors = graph.get_neighbors(node_id)
            num_connections = len(neighbors)
            
            # Find matching monomer
            monomer = self.find_exact_match(node.mol, num_connections)
            
            if monomer:
                node.monomer = monomer
                matched_count += 1
        
        return matched_count
