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
        
        # Get ordered sequence of monomers
        ordered_nodes = graph.get_ordered_nodes()
        sequence_symbols = [node.monomer.symbol if node.monomer else "X" for node in ordered_nodes]
        
        # Check if cyclic
        is_cyclic = graph.is_cyclic()
        
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
        
        # Generate final HELM notation
        if connections:
            connection_str = "|".join(connections)
            helm = f"PEPTIDE1{{{sequence}}}${connection_str}$$$V2.0"
        else:
            helm = f"PEPTIDE1{{{sequence}}}$$$$"
        
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

