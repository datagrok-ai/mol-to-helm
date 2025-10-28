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
        
        if graph.is_cyclic:
            # Cyclic peptide notation
            sequence = "].[".join(sequence_symbols)
            n = len(sequence_symbols)
            helm = f"PEPTIDE1{{[{sequence}]}}$PEPTIDE1,PEPTIDE1,{n}:R2-1:R1$$$V2.0"
        else:
            # Linear peptide notation
            sequence = ".".join(sequence_symbols)
            
            # Check for disulfide bridges or other non-peptide bonds
            has_special_bonds = any(
                link.linkage_type != LinkageType.PEPTIDE 
                for link in graph.links
            )
            
            if has_special_bonds:
                # Add connection notation for disulfide bridges
                connections = []
                for link in graph.links:
                    if link.linkage_type == LinkageType.DISULFIDE:
                        # Format: PEPTIDE1,PEPTIDE1,from_idx:R3-to_idx:R3
                        connections.append(
                            f"PEPTIDE1,PEPTIDE1,{link.from_node_id + 1}:R3-{link.to_node_id + 1}:R3"
                        )
                
                if connections:
                    connection_str = "|".join(connections)
                    helm = f"PEPTIDE1{{{sequence}}}${connection_str}$$$V2.0"
                else:
                    helm = f"PEPTIDE1{{{sequence}}}$$$$"
            else:
                helm = f"PEPTIDE1{{{sequence}}}$$$$"
        
        return helm

    def generate_helm_notation(self, monomers, is_cyclic: bool = False) -> str:
        """
        Legacy method: Generate HELM notation from a list of monomers.
        Kept for backward compatibility.
        
        Args:
            monomers: List of MonomerData objects
            is_cyclic: Whether the peptide is cyclic
        
        Returns:
            HELM notation string
        """
        if not monomers:
            return ""

        if is_cyclic:
            sequence = "].[".join([monomer.symbol for monomer in monomers])
            n = len(monomers)
            helm = f"PEPTIDE1{{[{sequence}]}}$PEPTIDE1,PEPTIDE1,{n}:R2-1:R1$$$V2.0"
        else:
            sequence = ".".join([monomer.symbol for monomer in monomers])
            helm = f"PEPTIDE1{{{sequence}}}$$$$"

        return helm

