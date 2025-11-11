from rdkit import Chem
from typing import Optional, List, Dict, Tuple
from enum import Enum


class LinkageType(Enum):
    """Types of linkages between fragments"""
    PEPTIDE = "peptide"
    DISULFIDE = "disulfide"
    ESTER = "ester"
    ETHER = "ether"
    THIOETHER = "thioether"
    UNKNOWN = "unknown"


class FragmentNode:
    """Represents a single molecular fragment (amino acid/monomer)"""
    
    def __init__(self, fragment_id: int, mol: Chem.Mol):
        self.id = fragment_id
        self.mol = mol  # RDKit molecule object
        self.smiles = Chem.MolToSmiles(mol, canonical=True) if mol else ""
        self.monomer = None  # Will be filled by matcher - MonomerData object
        self.is_c_terminal = False
        self.is_n_terminal = False
    
    def __repr__(self):
        monomer_name = self.monomer.symbol if self.monomer else "Unknown"
        return f"FragmentNode(id={self.id}, monomer={monomer_name}, smiles={self.smiles[:20]}...)"


class FragmentLink:
    """Represents a connection between two fragments"""
    
    def __init__(
        self, 
        from_node_id: int, 
        to_node_id: int,
        linkage_type: LinkageType,
        from_atom_idx: Optional[int] = None,
        to_atom_idx: Optional[int] = None
    ):
        self.from_node_id = from_node_id
        self.to_node_id = to_node_id
        self.linkage_type = linkage_type
        self.from_atom_idx = from_atom_idx  # Atom index in from_node's molecule
        self.to_atom_idx = to_atom_idx      # Atom index in to_node's molecule
    
    def __repr__(self):
        return f"FragmentLink({self.from_node_id} --{self.linkage_type.value}--> {self.to_node_id})"


class FragmentGraph:
    """
    Graph structure representing a molecule as fragments and their connections.
    
    Supports:
    - Linear peptides (chain of peptide bonds)
    - Cyclic peptides (peptide bond from last to first)
    - Disulfide bridges (additional S-S links)
    - Branched structures (multiple links per fragment)
    """
    
    def __init__(self):
        self.nodes: Dict[int, FragmentNode] = {}  # node_id -> FragmentNode
        self.links: List[FragmentLink] = []
    
    def add_node(self, node: FragmentNode) -> int:
        """Add a fragment node to the graph"""
        self.nodes[node.id] = node
        return node.id
    
    def add_link(self, link: FragmentLink):
        """Add a linkage between two nodes"""
        if link.from_node_id not in self.nodes or link.to_node_id not in self.nodes:
            raise ValueError(f"Cannot add link: nodes {link.from_node_id} or {link.to_node_id} not in graph")
        self.links.append(link)
    
    def get_node(self, node_id: int) -> Optional[FragmentNode]:
        """Get a node by ID"""
        return self.nodes.get(node_id)
    
    def get_neighbors(self, node_id: int) -> List[Tuple[int, LinkageType]]:
        """Get all neighbors of a node with their linkage types"""
        neighbors = []
        for link in self.links:
            if link.from_node_id == node_id:
                neighbors.append((link.to_node_id, link.linkage_type))
            elif link.to_node_id == node_id:
                neighbors.append((link.from_node_id, link.linkage_type))
        return neighbors
    
    def get_ordered_nodes(self) -> List[FragmentNode]:
        """
        Get nodes in sequential order (for linear/cyclic peptides).
        For branched structures, returns a depth-first traversal.
        """
        if not self.nodes:
            return []
        
        # Find starting node (N-terminal for peptides)
        start_node_id = None
        for node_id, node in self.nodes.items():
            if node.is_n_terminal:
                start_node_id = node_id
                break
        
        # If no N-terminal found, use first node
        if start_node_id is None:
            start_node_id = min(self.nodes.keys())
        
        # Traverse the graph
        ordered = []
        visited = set()
        self._traverse_from_node(start_node_id, visited, ordered)
        
        return ordered
    
    def _traverse_from_node(self, node_id: int, visited: set, ordered: list):
        """Helper for depth-first traversal"""
        if node_id in visited:
            return
        
        visited.add(node_id)
        ordered.append(self.nodes[node_id])
        
        # Get peptide bond neighbors first (to maintain chain order)
        peptide_neighbors = []
        other_neighbors = []
        
        for link in self.links:
            if link.from_node_id == node_id and link.to_node_id not in visited:
                if link.linkage_type == LinkageType.PEPTIDE:
                    peptide_neighbors.append(link.to_node_id)
                else:
                    other_neighbors.append(link.to_node_id)
        
        # Visit peptide bonds first, then others
        for neighbor_id in peptide_neighbors + other_neighbors:
            self._traverse_from_node(neighbor_id, visited, ordered)
    
    def get_fragment_sequence(self) -> List[str]:
        """Get sequence of monomer symbols (for matched fragments)"""
        ordered_nodes = self.get_ordered_nodes()
        return [
            node.monomer.symbol if node.monomer else f"X{node.id}"
            for node in ordered_nodes
        ]
    
    def is_cyclic(self) -> bool:
        """
        Detect if the peptide is cyclic.
        A cyclic peptide has a peptide bond connecting the last residue back to near the beginning.
        Handles cases where N-terminal caps (like 'ac' from Lys_Ac) create an extra fragment at position 0.
        """
        if len(self.nodes) < 3:
            return False
        
        # Get ordered nodes
        ordered = self.get_ordered_nodes()
        if len(ordered) < 3:
            return False
        
        # Get the last node ID
        last_id = ordered[-1].id
        
        # For a cyclic peptide, the last residue should connect back to one of the first few residues
        # (usually first, but could be second if there's an N-terminal cap like 'ac')
        # Check if last node has a peptide bond to any of the first 3 nodes
        first_few_ids = [ordered[i].id for i in range(min(3, len(ordered)))]
        
        for link in self.links:
            if link.linkage_type == LinkageType.PEPTIDE:
                # Check if link connects last node to one of the first few nodes
                if (link.from_node_id == last_id and link.to_node_id in first_few_ids) or \
                   (link.to_node_id == last_id and link.from_node_id in first_few_ids):
                    return True
        
        return False
    
    def find_all_cycles(self) -> List[List[int]]:
        """
        Find all cycles in the graph using DFS.
        Returns list of cycles, where each cycle is a list of node IDs.
        """
        cycles = []
        visited = set()
        rec_stack = set()
        parent = {}
        
        def dfs(node_id: int, path: List[int]):
            visited.add(node_id)
            rec_stack.add(node_id)
            path.append(node_id)
            
            # Get peptide bond neighbors
            neighbors = [n_id for n_id, link_type in self.get_neighbors(node_id) 
                        if link_type == LinkageType.PEPTIDE]
            
            for neighbor_id in neighbors:
                if neighbor_id not in visited:
                    parent[neighbor_id] = node_id
                    dfs(neighbor_id, path[:])
                elif neighbor_id in rec_stack and neighbor_id != parent.get(node_id):
                    # Found a cycle - extract it from path
                    cycle_start_idx = path.index(neighbor_id)
                    cycle = path[cycle_start_idx:] + [neighbor_id]
                    # Normalize cycle (start from smallest ID)
                    min_idx = cycle.index(min(cycle[:-1]))  # Don't include duplicate last element
                    normalized = cycle[min_idx:-1] + cycle[:min_idx]
                    if normalized not in cycles:
                        cycles.append(normalized)
            
            rec_stack.remove(node_id)
        
        # Try starting DFS from each unvisited node
        for node_id in self.nodes.keys():
            if node_id not in visited:
                parent[node_id] = None
                dfs(node_id, [])
        
        return cycles
    
    def get_connected_components(self) -> List[List[int]]:
        """
        Find all connected components in the graph.
        Returns list of components, where each component is a list of node IDs.
        """
        visited = set()
        components = []
        
        def dfs_component(node_id: int, component: List[int]):
            visited.add(node_id)
            component.append(node_id)
            neighbors = self.get_neighbors(node_id)
            for neighbor_id, _ in neighbors:
                if neighbor_id not in visited:
                    dfs_component(neighbor_id, component)
        
        for node_id in self.nodes.keys():
            if node_id not in visited:
                component = []
                dfs_component(node_id, component)
                components.append(sorted(component))
        
        return components
    
    def __len__(self):
        return len(self.nodes)
    
    def __repr__(self):
        return f"FragmentGraph(nodes={len(self.nodes)}, links={len(self.links)})"
    
    def to_dict(self) -> dict:
        """Convert graph to dictionary for serialization"""
        return {
            "nodes": [
                {
                    "id": node.id,
                    "smiles": node.smiles,
                    "monomer": node.monomer.symbol if node.monomer else None,
                    "is_n_terminal": node.is_n_terminal,
                    "is_c_terminal": node.is_c_terminal
                }
                for node in self.nodes.values()
            ],
            "links": [
                {
                    "from": link.from_node_id,
                    "to": link.to_node_id,
                    "type": link.linkage_type.value,
                    "from_atom": link.from_atom_idx,
                    "to_atom": link.to_atom_idx
                }
                for link in self.links
            ]
        }

