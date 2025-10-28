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
        self.is_cyclic = False
    
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
    
    def __len__(self):
        return len(self.nodes)
    
    def __repr__(self):
        return f"FragmentGraph(nodes={len(self.nodes)}, links={len(self.links)}, cyclic={self.is_cyclic})"
    
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
            ],
            "is_cyclic": self.is_cyclic
        }

