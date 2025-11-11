from fragment_graph import FragmentGraph, LinkageType


class HELMGenerator:
    """
    Generates HELM notation from fragment graphs or monomer lists.
    
    Supports:
    - Linear peptides
    - Cyclic peptides
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
        
        Args:
            graph: FragmentGraph containing matched monomers and their connections
        
        Returns:
            HELM notation string
        """
        if len(graph) == 0:
            return ""
        
        # Get ordered sequence of monomers (backbone)
        ordered_nodes_raw = graph.get_ordered_nodes()
        
        # Check if cyclic
        is_cyclic = graph.is_cyclic()
        
        # Filter backbone: nodes that are part of R1-R2 chain are backbone
        # Nodes connected only via R3 (side chain) are branches
        # 
        # Logic: A node at position 1 is a branch if:
        # - It has no R1 (N-terminus) - meaning it's a cap like 'ac' that only has R2
        # - It only has 1 peptide connection (to the real backbone)
        # 
        # Example: [ac].K in cyclic peptide
        # - 'ac' has only R2, no R1 → it's a cap
        # - 'ac' connects to K's R3 (side chain), not K's R1 (backbone)
        # - So 'ac' should be PEPTIDE2, not part of PEPTIDE1
        
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
            # This handles cases where N-terminal caps (like 'ac') are at position 1
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
                
                # Format branch chain (single monomer, so no dots needed)
                if is_cyclic and len(branch_symbol) > 1:
                    branch_chains.append(f"{branch_chain_name}{{[{branch_symbol}]}}")
                else:
                    branch_chains.append(f"{branch_chain_name}{{{branch_symbol}}}")
                
                # Find which backbone node this branch connects to
                # Look for links connecting this branch to the main backbone
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
                            # If branch has R1, connect to R1; if only R2, connect to R2
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

