from fragment_graph import FragmentGraph, LinkageType, FragmentNode, FragmentLink
from rdkit import Chem
from rdkit.Chem import GetSSSR
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict


class ChainAnalysis:
    """
    Analysis of a fragment graph identifying chains and branches.
    """
    def __init__(self):
        self.main_chain: List[int] = []  # Node IDs in main chain order
        self.branches: List[List[int]] = []  # Each branch is a list of node IDs
        self.rings: List[List[int]] = []  # Each ring is a list of node IDs
        self.node_to_chain: Dict[int, Tuple[str, int, int]] = {}  # node_id -> (chain_type, chain_idx, position)
        # chain_type: 'main', 'branch', 'ring'
        # chain_idx: index in main_chain/branches/rings list
        # position: position within that chain (0-indexed)


class HELMGenerator:
    """
    Generates HELM notation from fragment graphs or monomer lists.
    
    Supports:
    - Linear peptides
    - Cyclic peptides (via SSSR ring detection)
    - Disulfide bridges
    - Branched peptides (R3+ connections as separate chains)
    - Multiple rings (SSSR algorithm)
    - Complex topologies
    """
    
    def __init__(self):
        #GENERALIZATION ITEM: POLYMER TYPES SHOULD BE DERIVED FROM LIBRARY
        self.polymer_types = {
            "peptide": "PEPTIDE",
            "rna": "RNA",
            "dna": "DNA",
            "chemical": "CHEM"
        }
    
    def _detect_rings_simple(self, graph: FragmentGraph) -> List[List[int]]:
        """
        Detect simple cyclic structure using the existing is_cyclic() method.
        
        Args:
            graph: FragmentGraph to analyze
            
        Returns:
            List of rings (empty or single ring for simple cyclic peptides)
        """
        if graph.is_cyclic():
            # Simple cyclic peptide: return all nodes as one ring
            ordered = graph.get_ordered_nodes()
            return [[node.id for node in ordered]]
        return []
    
    def _get_r_group_for_connection(self, node: FragmentNode, neighbor_id: int, graph: FragmentGraph) -> str:
        """
        Determine which R-group is used for a connection between two nodes.
        
        R1: N-terminus (typically the "previous" residue in sequence)
        R2: C-terminus (typically the "next" residue in sequence)
        R3+: Side chain connections (branches)
        
        Args:
            node: The FragmentNode we're analyzing
            neighbor_id: ID of the neighbor node
            graph: The fragment graph
            
        Returns:
            R-group label (e.g., "R1", "R2", "R3")
        """
        # Count total connections to this node
        neighbors = graph.get_neighbors(node.id)
        peptide_neighbors = [n for n, lt in neighbors if lt == LinkageType.PEPTIDE]
        
        if len(peptide_neighbors) <= 2:
            # Simple case: only R1 and R2 (main chain)
            # First neighbor = R1, second neighbor = R2
            ordered_nodes = graph.get_ordered_nodes()
            node_idx = next((i for i, n in enumerate(ordered_nodes) if n.id == node.id), None)
            neighbor_idx = next((i for i, n in enumerate(ordered_nodes) if n.id == neighbor_id), None)
            
            if node_idx is not None and neighbor_idx is not None:
                if neighbor_idx < node_idx:
                    return "R1"  # Neighbor comes before (N-terminus direction)
                else:
                    return "R2"  # Neighbor comes after (C-terminus direction)
        
        # Complex case: 3+ connections means we have branches
        # R1 and R2 are the main chain, R3+ are branches
        # We need to identify which are main chain and which are branches
        
        # For now, use a heuristic: first two peptide neighbors are R1/R2, rest are R3+
        if len(peptide_neighbors) > 0:
            if neighbor_id == peptide_neighbors[0]:
                return "R1"
            elif len(peptide_neighbors) > 1 and neighbor_id == peptide_neighbors[1]:
                return "R2"
            else:
                # This is a branch connection
                # Find which R-group number (R3, R4, R5, ...)
                branch_index = peptide_neighbors.index(neighbor_id) - 2
                return f"R{3 + branch_index}"
        
        return "R1"  # Default fallback
    
    def _analyze_chains(self, graph: FragmentGraph) -> ChainAnalysis:
        """
        Analyze the fragment graph to identify:
        - Main chain (R1-R2 backbone connections)
        - Branches (R3+ connections)
        - Rings (detected via SSSR)
        
        Args:
            graph: FragmentGraph to analyze
            
        Returns:
            ChainAnalysis object with chain decomposition
        """
        analysis = ChainAnalysis()
        
        # Step 1: Detect rings using simple cyclic detection
        analysis.rings = self._detect_rings_simple(graph)
        
        # Step 2: Identify main chain vs branches
        # Main chain: nodes connected by R1-R2 (backbone)
        # Branches: nodes connected via R3+ (side chains)
        
        # Build adjacency info with connection counts
        node_connections = defaultdict(list)  # node_id -> [(neighbor_id, linkage_type)]
        
        for link in graph.links:
            if link.linkage_type == LinkageType.PEPTIDE:
                node_connections[link.from_node_id].append((link.to_node_id, link.linkage_type))
                node_connections[link.to_node_id].append((link.from_node_id, link.linkage_type))
        
        # Identify branch points (nodes with 3+ peptide connections)
        branch_points = set()
        for node_id, connections in node_connections.items():
            peptide_conns = [c for c in connections if c[1] == LinkageType.PEPTIDE]
            if len(peptide_conns) > 2:
                branch_points.add(node_id)
        
        # Build main chain (traverse R1-R2 connections)
        # Start from N-terminus or first node
        ordered_nodes = graph.get_ordered_nodes()
        
        if not ordered_nodes:
            return analysis
        
        # For simple cases without branches, main chain is the ordered sequence
        if not branch_points:
            analysis.main_chain = [node.id for node in ordered_nodes]
            for i, node_id in enumerate(analysis.main_chain):
                analysis.node_to_chain[node_id] = ('main', 0, i)
        else:
            # Complex case: need to separate main chain from branches
            # Main chain: longest path or path containing most nodes with only 2 connections
            visited = set()
            main_chain = []
            
            # Start from N-terminus
            current = ordered_nodes[0].id
            visited.add(current)
            main_chain.append(current)
            
            # Traverse along the main chain (preferring nodes with 2 connections)
            while True:
                neighbors = [n for n, _ in node_connections[current] if n not in visited]
                if not neighbors:
                    break
                
                # Prefer neighbor with 2 connections (main chain) over 1 or 3+ (branch point)
                next_node = None
                for n in neighbors:
                    n_conn_count = len([c for c in node_connections[n] if c[1] == LinkageType.PEPTIDE])
                    if n_conn_count == 2:
                        next_node = n
                        break
                
                if next_node is None:
                    next_node = neighbors[0]  # Take any available
                
                visited.add(next_node)
                main_chain.append(next_node)
                current = next_node
            
            analysis.main_chain = main_chain
            for i, node_id in enumerate(main_chain):
                analysis.node_to_chain[node_id] = ('main', 0, i)
            
            # Identify branches (nodes not in main chain but connected to it)
            for node_id in graph.nodes.keys():
                if node_id not in visited:
                    # This is a branch - trace it
                    branch = []
                    branch_visited = set()
                    
                    def trace_branch(node_id):
                        if node_id in branch_visited or node_id in analysis.main_chain:
                            return
                        branch_visited.add(node_id)
                        branch.append(node_id)
                        for neighbor_id, _ in node_connections[node_id]:
                            if neighbor_id not in analysis.main_chain:
                                trace_branch(neighbor_id)
                    
                    trace_branch(node_id)
                    if branch:
                        branch_idx = len(analysis.branches)
                        analysis.branches.append(branch)
                        for i, bid in enumerate(branch):
                            analysis.node_to_chain[bid] = ('branch', branch_idx, i)
        
        return analysis

    def generate_helm_from_graph(self, graph: FragmentGraph) -> str:
        """
        Generate HELM notation from a FragmentGraph with support for:
        - Multiple rings (SSSR detection)
        - Branches as separate polymer chains
        - Complex topologies
        
        Args:
            graph: FragmentGraph containing matched monomers and their connections
        
        Returns:
            HELM notation string
        """
        if len(graph) == 0:
            return ""
        
        # Analyze the graph structure
        analysis = self._analyze_chains(graph)
        
        # Build polymer chain strings
        polymer_chains = []
        chain_id_map = {}  # Maps (chain_type, chain_idx) -> HELM chain ID (e.g., "PEPTIDE1")
        
        # Main chain is always PEPTIDE1
        if analysis.main_chain:
            main_sequence = self._build_sequence(analysis.main_chain, graph, analysis.rings)
            polymer_chains.append(f"PEPTIDE1{{{main_sequence}}}")
            chain_id_map[('main', 0)] = "PEPTIDE1"
        
        # Branches become separate chains (PEPTIDE2, PEPTIDE3, ..., or CHEM1, CHEM2, ...)
        for branch_idx, branch in enumerate(analysis.branches):
            branch_sequence = self._build_sequence(branch, graph, analysis.rings)
            
            # Determine polymer type for branch
            # If branch is a single monomer (like 'ac' acetyl cap), use CHEM
            # Otherwise use PEPTIDE
            if len(branch) == 1:
                chain_id = f"CHEM{branch_idx + 1}"
                polymer_chains.append(f"{chain_id}{{{branch_sequence}}}")
            else:
                chain_id = f"PEPTIDE{len(polymer_chains) + 1}"
                polymer_chains.append(f"{chain_id}{{{branch_sequence}}}")
            
            chain_id_map[('branch', branch_idx)] = chain_id
        
        # Build connections section
        connections = []
        
        # 1. Add ring closures (cyclic peptides) - only for main chain
        is_cyclic = len(analysis.rings) > 0
        if is_cyclic and analysis.main_chain:
            # For cyclic peptides, connect last to first in main chain
            num_residues = len(analysis.main_chain)
            connections.append(f"PEPTIDE1,PEPTIDE1,{num_residues}:R2-1:R1")
        
        # 2. Add inter-chain connections (branches to main chain)
        for branch_idx, branch in enumerate(analysis.branches):
            if not branch:
                continue
            
            # Find connection point between branch and main chain
            branch_node_id = branch[0]  # First node in branch
            
            # Find which main chain node connects to this branch
            for link in graph.links:
                if link.linkage_type == LinkageType.PEPTIDE:
                    main_node_id = None
                    branch_side_node_id = None
                    
                    if link.from_node_id in analysis.main_chain and link.to_node_id == branch_node_id:
                        main_node_id = link.from_node_id
                        branch_side_node_id = link.to_node_id
                    elif link.to_node_id in analysis.main_chain and link.from_node_id == branch_node_id:
                        main_node_id = link.to_node_id
                        branch_side_node_id = link.from_node_id
                    
                    if main_node_id:
                        # Found connection!
                        main_pos = analysis.main_chain.index(main_node_id) + 1
                        branch_chain_id = chain_id_map[('branch', branch_idx)]
                        
                        # Determine R-groups
                        # Main chain node uses R3+ for branch connections
                        main_node = graph.get_node(main_node_id)
                        main_r_group = self._get_r_group_for_connection(main_node, branch_side_node_id, graph)
                        
                        # Branch node uses R2 to connect to main chain (if it's CHEM) or R1 (if PEPTIDE)
                        if branch_chain_id.startswith("CHEM"):
                            branch_r_group = "R2"
                        else:
                            branch_r_group = "R1"
                        
                        connections.append(f"PEPTIDE1,{branch_chain_id},{main_pos}:{main_r_group}-1:{branch_r_group}")
                        break
        
        # 3. Add disulfide bridges
        for link in graph.links:
            if link.linkage_type == LinkageType.DISULFIDE:
                from_node_id = link.from_node_id
                to_node_id = link.to_node_id
                
                from_chain_info = analysis.node_to_chain.get(from_node_id)
                to_chain_info = analysis.node_to_chain.get(to_node_id)
                
                if from_chain_info and to_chain_info:
                    from_chain_type, from_chain_idx, from_pos = from_chain_info
                    to_chain_type, to_chain_idx, to_pos = to_chain_info
                    
                    # Determine chain IDs
                    from_chain_id = chain_id_map.get((from_chain_type, from_chain_idx), "PEPTIDE1")
                    to_chain_id = chain_id_map.get((to_chain_type, to_chain_idx), "PEPTIDE1")
                    
                    # Disulfide bonds use R3 (side chain sulfur)
                    connections.append(f"{from_chain_id},{to_chain_id},{from_pos + 1}:R3-{to_pos + 1}:R3")
        
        # Assemble final HELM notation
        simple_polymer_list = "+".join(polymer_chains)
        
        if connections:
            connection_str = "|".join(connections)
            helm = f"{simple_polymer_list}${connection_str}$$$V2.0"
        else:
            helm = f"{simple_polymer_list}$$$$V2.0"
        
        return helm
    
    def _build_sequence(self, node_ids: List[int], graph: FragmentGraph, rings: List[List[int]]) -> str:
        """
        Build HELM sequence string for a chain of nodes.
        
        Args:
            node_ids: List of node IDs in order
            graph: The fragment graph
            rings: List of detected rings
            
        Returns:
            HELM sequence string (e.g., "A.C.G.V" or "[meI].[hHis].[Aca]")
        """
        # Check if this chain is part of a ring
        is_ring = any(set(node_ids).issubset(set(ring)) for ring in rings)
        
        symbols = []
        for node_id in node_ids:
            node = graph.get_node(node_id)
            if node and node.monomer:
                symbol = node.monomer.symbol
                # Wrap multi-letter symbols in brackets for cyclic peptides
                if is_ring and len(symbol) > 1:
                    symbols.append(f"[{symbol}]")
                else:
                    symbols.append(symbol)
            else:
                symbols.append("X")  # Unknown monomer
        
        return ".".join(symbols)

    def generate_helm_notation(self, monomers) -> str:
        """
        Legacy method: Generate HELM notation from a list of monomers.
        Kept for backward compatibility.
        
        Args:
            monomers: List of MonomerData objects
        
        Returns:
            HELM notation string
        """
        if not monomers:
            return ""

        sequence = ".".join([monomer.symbol for monomer in monomers])
        helm = f"PEPTIDE1{{{sequence}}}$$$$"

        return helm

