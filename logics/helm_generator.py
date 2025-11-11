from fragment_graph import FragmentGraph, LinkageType


class HELMGenerator:
    """
    Generates HELM notation from fragment graphs or monomer lists.
    
    Supports:
    - Linear peptides
    - Cyclic peptides
    - Multi-chain structures (BILN peptides)
    - Disulfide bridges
    - Custom linkages
    """
    
    def __init__(self):
        #GENERALIZATION ITEM: POLYMER TYPES SHOULD BE DERIVED FROM LIBRARY
        self.polymer_types = {
            "peptide": "PEPTIDE",
            "rna": "RNA",
            "dna": "DNA",
            "chemical": "CHEM"
        }

    def generate_helm_from_graph(self, graph: FragmentGraph) -> str:
        """
        Generate HELM notation from a FragmentGraph.
        
        Supports multi-chain structures (BILN peptides):
        - Detects all cycles (rings) using SSSR-like algorithm
        - Each cycle becomes a separate PEPTIDE chain
        - R1-R2 connections define backbone within each chain
        - R3 connections link chains together
        
        Args:
            graph: FragmentGraph containing matched monomers and their connections
        
        Returns:
            HELM notation string
        """
        if len(graph) == 0:
            return ""
        
        # Find all cycles in the graph (each cycle will be a separate PEPTIDE chain)
        cycles = graph.find_all_cycles()
        
        # Decision: Use multi-chain HELM only if:
        # 1. Multiple cycles exist (BILN-style structure), OR
        # 2. There are standalone nodes not in any cycle (attached fragments)
        
        if not cycles:
            # No cycles - simple linear peptide
            return self._generate_simple_helm(graph)
        
        if len(cycles) == 1:
            # Single cycle - check if there are standalone nodes attached
            nodes_in_cycle = set(cycles[0])
            standalone_nodes = [nid for nid in graph.nodes.keys() if nid not in nodes_in_cycle]
            
            if not standalone_nodes:
                # Simple single-cycle peptide - use simple generator
                return self._generate_simple_helm(graph)
        
        # Multi-chain structure detected (multiple cycles or cycle with attachments)
        return self._generate_multi_chain_helm(graph, cycles)
    
    def _generate_simple_helm(self, graph: FragmentGraph) -> str:
        """
        Generate HELM for simple linear or single-cycle peptides.
        This is the original implementation for backward compatibility.
        """
        # Get ordered sequence of monomers (backbone)
        ordered_nodes_raw = graph.get_ordered_nodes()
        
        # Check if cyclic
        is_cyclic = graph.is_cyclic()
        
        # Filter backbone: nodes that are part of R1-R2 chain are backbone
        # Nodes connected only via R3 (side chain) are branches
        backbone_nodes = []
        for i, node in enumerate(ordered_nodes_raw):
            is_branch = False
            
            if i == 0 and len(ordered_nodes_raw) > 1 and node.monomer:
                # Check if this first node lacks R1 (N-terminus)
                # If it has no R1, it's a cap that should be a branch
                has_r1 = 'R1' in node.monomer.r_groups
                
                if not has_r1:
                    # This is an N-terminal cap (like 'ac') at position 1
                    # It should be a branch, not part of the main backbone
                    is_branch = True
            
            if not is_branch:
                backbone_nodes.append(node)
        
        ordered_nodes = backbone_nodes
        sequence_symbols = [node.monomer.symbol if node.monomer else "X" for node in ordered_nodes]
        
        # Detect branch nodes (nodes not in backbone)
        ordered_node_ids = {node.id for node in ordered_nodes}
        branch_nodes = [(node_id, node) for node_id, node in graph.nodes.items() 
                       if node_id not in ordered_node_ids]
        
        # Generate sequence notation
        if is_cyclic:
            # Cyclic: wrap multi-letter monomers in brackets, single-letter ones stay as-is
            formatted_symbols = [f"[{symbol}]" if len(symbol) > 1 else symbol for symbol in sequence_symbols]
            sequence = ".".join(formatted_symbols)
        else:
            # Linear: no brackets
            sequence = ".".join(sequence_symbols)
        
        # Collect non-sequential connections (disulfide bridges, cyclic bonds, etc.)
        connections = []
        
        if is_cyclic:
            # Find the actual cyclic peptide bond (last residue connects back to beginning)
            last_id = ordered_nodes[-1].id
            first_few_ids = [ordered_nodes[i].id for i in range(min(3, len(ordered_nodes)))]
            
            for link in graph.links:
                if link.linkage_type == LinkageType.PEPTIDE:
                    # Check if this is the cyclic bond (last to one of first few)
                    is_cyclic_bond = False
                    from_id, to_id = None, None
                    
                    if link.from_node_id == last_id and link.to_node_id in first_few_ids:
                        from_id, to_id = link.from_node_id, link.to_node_id
                        is_cyclic_bond = True
                    elif link.to_node_id == last_id and link.from_node_id in first_few_ids:
                        from_id, to_id = link.to_node_id, link.from_node_id
                        is_cyclic_bond = True
                    
                    if is_cyclic_bond:
                        # Find positions (1-indexed)
                        from_pos = next((i + 1 for i, n in enumerate(ordered_nodes) if n.id == from_id), None)
                        to_pos = next((i + 1 for i, n in enumerate(ordered_nodes) if n.id == to_id), None)
                        
                        if from_pos and to_pos:
                            connections.append(f"PEPTIDE1,PEPTIDE1,{from_pos}:R2-{to_pos}:R1")
                            break
        
        # Add disulfide bridges
        for link in graph.links:
            if link.linkage_type == LinkageType.DISULFIDE:
                # Get positions in ordered sequence (1-indexed)
                from_pos = None
                to_pos = None
                for i, node in enumerate(ordered_nodes):
                    if node.id == link.from_node_id:
                        from_pos = i + 1
                    if node.id == link.to_node_id:
                        to_pos = i + 1
                
                if from_pos and to_pos:
                    # Format: PEPTIDE1,PEPTIDE1,from_pos:R3-to_pos:R3
                    connections.append(f"PEPTIDE1,PEPTIDE1,{from_pos}:R3-{to_pos}:R3")
        
        # Handle branch nodes (side chain modifications)
        # Create separate PEPTIDE chains for each branch
        branch_chains = []
        if branch_nodes:
            for branch_idx, (branch_node_id, branch_node) in enumerate(branch_nodes, start=2):
                branch_chain_name = f"PEPTIDE{branch_idx}"
                branch_symbol = branch_node.monomer.symbol if branch_node.monomer else f"X{branch_node_id}"
                
                # Format branch chain (single monomer)
                # In cyclic peptides, always use brackets for consistency with reference HELM
                if is_cyclic:
                    branch_chains.append(f"{branch_chain_name}{{[{branch_symbol}]}}")
                else:
                    branch_chains.append(f"{branch_chain_name}{{{branch_symbol}}}")
                
                # Find which backbone node this branch connects to
                for link in graph.links:
                    backbone_node_id = None
                    if link.from_node_id == branch_node_id and link.to_node_id in ordered_node_ids:
                        backbone_node_id = link.to_node_id
                    elif link.to_node_id == branch_node_id and link.from_node_id in ordered_node_ids:
                        backbone_node_id = link.from_node_id
                    
                    if backbone_node_id is not None:
                        # Find position of backbone node (1-indexed)
                        backbone_pos = next((i + 1 for i, n in enumerate(ordered_nodes) if n.id == backbone_node_id), None)
                        if backbone_pos:
                            # Determine which R-group the branch uses
                            branch_r_group = "R1"
                            if branch_node.monomer:
                                if 'R1' in branch_node.monomer.r_groups:
                                    branch_r_group = "R1"
                                elif 'R2' in branch_node.monomer.r_groups:
                                    branch_r_group = "R2"
                            
                            # Connection: backbone position R3 (side chain) to branch position 1 R-group
                            connections.append(f"PEPTIDE1,{branch_chain_name},{backbone_pos}:R3-1:{branch_r_group}")
                            break
        
        # Generate final HELM notation
        all_chains = [f"PEPTIDE1{{{sequence}}}"] + branch_chains
        helm_chains = "|".join(all_chains)
        
        if connections:
            connection_str = "|".join(connections)
            helm = f"{helm_chains}${connection_str}$$$V2.0"
        else:
            helm = f"{helm_chains}$$$$V2.0"
        
        return helm
    
    def _generate_multi_chain_helm(self, graph: FragmentGraph, cycles: list) -> str:
        """
        Generate HELM for multi-chain structures (BILN peptides).
        
        Strategy:
        1. Each cycle becomes a separate PEPTIDE chain
        2. Nodes not in cycles become additional chains
        3. R3 connections between chains are added as cross-links
        """
        # Identify which nodes belong to which cycles
        nodes_in_cycles = set()
        for cycle in cycles:
            nodes_in_cycles.update(cycle)
        
        # Find standalone nodes (not in any cycle)
        standalone_nodes = [nid for nid in graph.nodes.keys() if nid not in nodes_in_cycles]
        
        # Build PEPTIDE chains
        chains = []
        chain_node_map = {}  # Maps node_id -> (chain_idx, position_in_chain)
        
        # Add cyclic chains
        for cycle_idx, cycle in enumerate(cycles, start=1):
            chain_name = f"PEPTIDE{cycle_idx}"
            # Create sequence from cycle nodes
            sequence_symbols = []
            for pos, node_id in enumerate(cycle):
                node = graph.nodes[node_id]
                symbol = node.monomer.symbol if node.monomer else f"X{node_id}"
                sequence_symbols.append(symbol)
                chain_node_map[node_id] = (cycle_idx, pos + 1)  # 1-indexed position
            
            # Format with brackets for multi-letter symbols
            formatted = [f"[{s}]" if len(s) > 1 else s for s in sequence_symbols]
            sequence = ".".join(formatted)
            chains.append(f"{chain_name}{{{sequence}}}")
        
        # Add standalone chains (linear fragments not in cycles)
        next_chain_idx = len(cycles) + 1
        for node_id in standalone_nodes:
            chain_name = f"PEPTIDE{next_chain_idx}"
            node = graph.nodes[node_id]
            symbol = node.monomer.symbol if node.monomer else f"X{node_id}"
            chains.append(f"{chain_name}{{{symbol}}}")
            chain_node_map[node_id] = (next_chain_idx, 1)
            next_chain_idx += 1
        
        # Build connections
        connections = []
        
        # Add cyclic connections (R1-R2 within each cycle)
        for cycle_idx, cycle in enumerate(cycles, start=1):
            if len(cycle) >= 3:
                # Connect last to first
                chain_name = f"PEPTIDE{cycle_idx}"
                last_pos = len(cycle)
                connections.append(f"{chain_name},{chain_name},{last_pos}:R2-1:R1")
        
        # Add inter-chain connections (R3 links) and disulfide bridges
        processed_links = set()
        for link in graph.links:
            link_key = tuple(sorted([link.from_node_id, link.to_node_id]))
            if link_key in processed_links:
                continue
            
            from_chain_info = chain_node_map.get(link.from_node_id)
            to_chain_info = chain_node_map.get(link.to_node_id)
            
            if not from_chain_info or not to_chain_info:
                continue
            
            from_chain, from_pos = from_chain_info
            to_chain, to_pos = to_chain_info
            
            # Skip intra-cycle backbone peptide bonds (already handled by R1-R2 connection)
            if from_chain == to_chain and link.linkage_type == LinkageType.PEPTIDE:
                # Check if this is a sequential bond within the cycle
                cycle = cycles[from_chain - 1] if from_chain <= len(cycles) else []
                # Sequential bonds: adjacent positions or last-to-first
                if abs(from_pos - to_pos) == 1 or (from_pos == 1 and to_pos == len(cycle)) or (to_pos == 1 and from_pos == len(cycle)):
                    processed_links.add(link_key)
                    continue
            
            # Add cross-chain connections or intra-chain disulfide bridges
            if link.linkage_type == LinkageType.DISULFIDE:
                # Disulfide uses R3 (side chain cysteine)
                r_group = "R3"
            elif link.linkage_type == LinkageType.PEPTIDE:
                # Cross-chain peptide bond (side chain R3 connection)
                r_group = "R3"
            else:
                r_group = "R3"
            
            from_chain_name = f"PEPTIDE{from_chain}"
            to_chain_name = f"PEPTIDE{to_chain}"
            connections.append(f"{from_chain_name},{to_chain_name},{from_pos}:{r_group}-{to_pos}:{r_group}")
            processed_links.add(link_key)
        
        # Generate final HELM
        helm_chains = "|".join(chains)
        if connections:
            connection_str = "|".join(connections)
            helm = f"{helm_chains}${connection_str}$$$V2.0"
        else:
            helm = f"{helm_chains}$$$$V2.0"
        
        return helm

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
