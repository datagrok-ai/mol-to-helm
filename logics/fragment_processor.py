from rdkit import Chem
from fragment_graph import FragmentGraph, FragmentNode, FragmentLink, LinkageType


class BondDetector:
    #GENERALIZATION ITEM: BOND PATTERNS SHOULD BE DERIVED FROM LIBRARY
    def __init__(self):
        # True peptide bond: C and N both in backbone (each bonded to carbons)
        # Alpha carbons can be sp3 (X4) or sp2 (X3) for dehydroamino acids
        self.peptide_bond = Chem.MolFromSmarts('[C;X3,X4]-[C;X3](=[O;X1])-[N;X3]-[C;X3,X4]')
        # True disulfide bond: S-S where each S is bonded to carbon (cysteine residues)
        self.disulfide_bond = Chem.MolFromSmarts('[C;X4]-[S;X2]-[S;X2]-[C;X4]')
        # Primary amine at N-terminus (can be NH2 or NH3+), alpha-C can be sp3 or sp2
        self.primary_amine = Chem.MolFromSmarts('[N;H2,H3;X3,X4]-[C;X3,X4]')

    def find_cleavable_bonds(self, mol: Chem.Mol):
        """
        Find all cleavable bonds in the molecule.
        
        Returns:
            List of tuples: (atom1_idx, atom2_idx, LinkageType)
        """
        try:
            all_bonds = []

            # Find peptide bonds
            peptide_bonds = self._find_peptide_bonds(mol)
            all_bonds.extend([(bond[0], bond[1], LinkageType.PEPTIDE) for bond in peptide_bonds])
            
            # Find disulfide bonds
            disulfide_bonds = self._find_disulfide_bonds(mol)
            all_bonds.extend([(bond[0], bond[1], LinkageType.DISULFIDE) for bond in disulfide_bonds])

            # Order peptide bonds from N to C (keep disulfide bonds unordered)
            peptide_only = [(b[0], b[1]) for b in all_bonds if b[2] == LinkageType.PEPTIDE]
            ordered_peptide = self._order_bonds_from_n_to_c(mol, peptide_only)
            
            # Rebuild with types
            ordered_bonds = [(b[0], b[1], LinkageType.PEPTIDE) for b in ordered_peptide]
            ordered_bonds.extend([b for b in all_bonds if b[2] != LinkageType.PEPTIDE])

            return ordered_bonds

        except Exception:
            return []

    def _find_peptide_bonds(self, mol: Chem.Mol):
        bonds = []
        try:
            matches = mol.GetSubstructMatches(self.peptide_bond)
            for match in matches:
                if len(match) >= 5:
                    # Pattern: [C;X3,X4]-[C;X3](=[O;X1])-[N;X3]-[C;X3,X4]
                    # match[0]=alpha-C (sp2 or sp3), match[1]=carbonyl-C, match[2]=O, match[3]=N, match[4]=next-alpha-C (sp2 or sp3)
                    c_atom = match[1]  # Carbonyl carbon
                    n_atom = match[3]  # Nitrogen
                    bonds.append((c_atom, n_atom))
        except Exception:
            pass
        return bonds
    
    def _find_disulfide_bonds(self, mol: Chem.Mol):
        """Find disulfide bonds (S-S linkages)"""
        bonds = []
        try:
            matches = mol.GetSubstructMatches(self.disulfide_bond)
            for match in matches:
                if len(match) >= 4:
                    # Pattern: [C;X4]-[S;X2]-[S;X2]-[C;X4]
                    # match[0]=C, match[1]=S, match[2]=S, match[3]=C
                    s1_atom = match[1]  # First sulfur
                    s2_atom = match[2]  # Second sulfur
                    bonds.append((s1_atom, s2_atom))
        except Exception:
            pass
        return bonds

    def _order_bonds_from_n_to_c(self, mol: Chem.Mol, bonds):
        if not bonds:
            return bonds

        n_terminal = self._find_n_terminal(mol)
        if n_terminal is None:
            return bonds

        ordered = []
        visited = set()
        current = n_terminal

        while current is not None and len(ordered) < len(bonds):
            next_bond = None
            for bond in bonds:
                if bond not in visited and bond[1] == current:
                    next_bond = bond
                    break

            if next_bond is None:
                break

            ordered.append(next_bond)
            visited.add(next_bond)
            current = next_bond[0]

        for bond in bonds:
            if bond not in visited:
                ordered.append(bond)

        return ordered

    def _find_n_terminal(self, mol: Chem.Mol):
        try:
            matches = mol.GetSubstructMatches(self.primary_amine)
            if matches:
                return matches[0][0]

            max_h = -1
            n_term = None
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 7:
                    h_count = atom.GetTotalNumHs()
                    if h_count > max_h:
                        max_h = h_count
                        n_term = atom.GetIdx()
            return n_term

        except Exception:
            return None


class FragmentProcessor:
    def __init__(self, monomer_library):
        self.monomer_library = monomer_library
        self.bond_detector = BondDetector()

    def process_molecule(self, mol: Chem.Mol) -> FragmentGraph:
        """
        Process a molecule into a fragment graph.
        
        Args:
            mol: RDKit molecule object
        
        Returns:
            FragmentGraph object containing fragments and their connections
        """
        graph = FragmentGraph()
        # Store original molecule for fragment recovery
        graph.original_mol = mol
        
        try:
            bonds_to_cleave = self.bond_detector.find_cleavable_bonds(mol)

            if not bonds_to_cleave:
                # Single fragment (no cleavable bonds)
                node = FragmentNode(0, mol)
                node.is_n_terminal = True
                node.is_c_terminal = True
                graph.add_node(node)
                return graph

            # Extract bond info for fragmentation
            bond_indices = []
            bond_info = []  # (bond_idx, atom1, atom2, linkage_type)
            seen_bonds = set()  # Track which bonds we've already added
            
            for atom1, atom2, linkage_type in bonds_to_cleave:
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond:
                    bond_idx = bond.GetIdx()
                    if bond_idx not in seen_bonds:
                        bond_indices.append(bond_idx)
                        bond_info.append((bond_idx, atom1, atom2, linkage_type))
                        seen_bonds.add(bond_idx)
                    # Skip duplicate bonds silently
                # Skip invalid bonds silently

            if not bond_indices:
                # No valid bonds found
                node = FragmentNode(0, mol)
                node.is_n_terminal = True
                node.is_c_terminal = True
                graph.add_node(node)
                return graph

            # Fragment the molecule
            fragmented_mol = Chem.FragmentOnBonds(mol, bond_indices, addDummies=True)
            
            # Get fragments AND their atom mappings separately
            fragments_tuple = Chem.GetMolFrags(
                fragmented_mol, 
                asMols=True, 
                sanitizeFrags=True
            )
            fragments = list(fragments_tuple)
            
            # Store bond cleavage info for recovery - we'll use this to selectively re-fragment
            graph.cleaved_bond_indices = bond_indices
            graph.bond_info = bond_info
            print(f"DEBUG: Created {len(fragments)} fragments, cleaved {len(bond_indices)} bonds")

            # Create nodes for each fragment
            fragment_nodes = []
            for i, frag in enumerate(fragments):
                clean_frag = self._clean_fragment(frag)
                if clean_frag and clean_frag.GetNumAtoms() >= 3:
                    is_c_terminal = (i == len(fragments) - 1)
                    is_n_terminal = (i == 0)
                    # No normalization! Use fragment as-is
                    node = FragmentNode(i, clean_frag)
                    node.is_c_terminal = is_c_terminal
                    node.is_n_terminal = is_n_terminal
                    graph.add_node(node)
                    fragment_nodes.append((i, node))

            # Create links between fragments based on cleaved bonds
            # For sequential peptide bonds
            peptide_links = [b for b in bond_info if b[3] == LinkageType.PEPTIDE]
            for i in range(len(fragment_nodes) - 1):
                from_id, _ = fragment_nodes[i]
                to_id, _ = fragment_nodes[i + 1]
                link = FragmentLink(from_id, to_id, LinkageType.PEPTIDE)
                graph.add_link(link)
            
            # Add disulfide bridges (if any)
            # TODO: Track which fragments contain the S atoms for proper linking
            disulfide_links = [b for b in bond_info if b[3] == LinkageType.DISULFIDE]
            # For now, disulfide bonds require more complex atom tracking
            # This is a placeholder for future enhancement

            return graph

        except Exception as e:
            # Fallback: single node with original molecule
            node = FragmentNode(0, mol)
            node.is_n_terminal = True
            node.is_c_terminal = True
            graph.add_node(node)
            return graph

    def _clean_fragment(self, mol: Chem.Mol):
        try:
            mol_copy = Chem.Mol(mol)
            atoms_to_remove = []

            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atoms_to_remove.append(atom.GetIdx())

            atoms_to_remove.sort(reverse=True)
            if atoms_to_remove:
                emol = Chem.EditableMol(mol_copy)
                for atom_idx in atoms_to_remove:
                    emol.RemoveAtom(atom_idx)
                return emol.GetMol()

            return mol_copy

        except Exception:
            return None

    def _reconstruct_fragment(self, node_ids: list, graph: FragmentGraph) -> Chem.Mol:
        """
        Reconstruct a molecule by combining multiple fragment nodes.
        Re-fragments the original molecule, excluding bonds between the nodes to merge.
        """
        if not node_ids or not hasattr(graph, 'original_mol') or not hasattr(graph, 'cleaved_bond_indices'):
            return None
        
        try:
            # Sort node IDs to ensure consistent ordering
            sorted_nodes = sorted(node_ids)
            
            # Identify which bonds to exclude (bonds between consecutive merged nodes)
            bonds_to_exclude = set()
            for i in range(len(sorted_nodes) - 1):
                # We want to keep the bond between node i and node i+1
                # This bond would be at position sorted_nodes[i] in the cleaved_bond_indices
                if sorted_nodes[i] + 1 == sorted_nodes[i + 1]:
                    # Consecutive nodes - exclude the bond between them
                    if sorted_nodes[i] < len(graph.cleaved_bond_indices):
                        bonds_to_exclude.add(sorted_nodes[i])
            
            # Create new bond list excluding the bonds we want to keep
            new_bond_indices = [
                bond_idx for i, bond_idx in enumerate(graph.cleaved_bond_indices)
                if i not in bonds_to_exclude
            ]
            
            print(f"DEBUG reconstruct: Original had {len(graph.cleaved_bond_indices)} cleaved bonds, "
                  f"excluding {len(bonds_to_exclude)} bonds, new list has {len(new_bond_indices)} bonds")
            
            # Re-fragment with the modified bond list
            if not new_bond_indices:
                # No bonds to cleave - return whole molecule
                return graph.original_mol
            
            fragmented_mol = Chem.FragmentOnBonds(graph.original_mol, new_bond_indices, addDummies=True)
            fragments_tuple = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
            fragments = list(fragments_tuple)
            
            # Find which fragment corresponds to our merged nodes
            # The merged nodes should be at the position of the first node ID in sorted order
            target_idx = sorted_nodes[0]
            
            # Account for excluded bonds shifting fragment indices
            adjusted_idx = target_idx - sum(1 for excluded_idx in bonds_to_exclude if excluded_idx < target_idx)
            
            print(f"DEBUG reconstruct: Got {len(fragments)} fragments after re-fragmentation, "
                  f"target_idx={target_idx}, adjusted_idx={adjusted_idx}")
            
            if adjusted_idx < len(fragments):
                clean_frag = self._clean_fragment(fragments[adjusted_idx])
                return clean_frag if clean_frag else fragments[adjusted_idx]
            
            return None
        
        except Exception as e:
            print(f"DEBUG reconstruct: Exception: {e}")
            return None

    def _merge_nodes_in_graph(self, graph: FragmentGraph, nodes_to_merge: list, 
                              new_node: FragmentNode) -> None:
        """
        Remove old nodes, add new merged node, update all links.
        Preserves terminal flags from edge nodes.
        """
        if not nodes_to_merge:
            return
        
        # Sort node IDs to identify edge nodes
        sorted_nodes = sorted(nodes_to_merge)
        leftmost = sorted_nodes[0]
        rightmost = sorted_nodes[-1]
        
        # Preserve terminal flags
        if leftmost in graph.nodes:
            new_node.is_n_terminal = graph.nodes[leftmost].is_n_terminal
        if rightmost in graph.nodes:
            new_node.is_c_terminal = graph.nodes[rightmost].is_c_terminal
        
        # Update links: replace references to merged nodes
        updated_links = []
        nodes_to_merge_set = set(nodes_to_merge)
        
        for link in graph.links:
            from_in = link.from_node_id in nodes_to_merge_set
            to_in = link.to_node_id in nodes_to_merge_set
            
            # Skip internal links between merged nodes
            if from_in and to_in:
                continue
            
            # Update link if one end is being merged
            new_from = new_node.id if from_in else link.from_node_id
            new_to = new_node.id if to_in else link.to_node_id
            
            updated_links.append(FragmentLink(new_from, new_to, link.linkage_type))
        
        # Remove old nodes
        for node_id in nodes_to_merge:
            if node_id in graph.nodes:
                del graph.nodes[node_id]
        
        # Add new node and update links
        graph.add_node(new_node)
        graph.links = updated_links

    def recover_unmatched_fragments(self, graph: FragmentGraph, matcher) -> bool:
        """
        Try to recover unmatched fragments by merging with neighbors.
        Returns True if any merges were successful.
        """
        # Identify unmatched nodes
        unmatched_nodes = []
        for node_id, node in graph.nodes.items():
            if node.monomer and node.monomer.symbol.startswith("X"):
                unmatched_nodes.append(node_id)
        
        if not unmatched_nodes:
            return False
        
        print(f"DEBUG: Found {len(unmatched_nodes)} unmatched nodes: {unmatched_nodes}")
        
        had_changes = False
        
        # Try to recover each unmatched node
        for node_id in unmatched_nodes:
            # Check if node still exists (might have been merged already)
            if node_id not in graph.nodes:
                continue
            
            # Get neighbors
            neighbors = graph.get_neighbors(node_id)
            neighbor_ids = [n[0] for n in neighbors]
            
            if not neighbor_ids:
                continue
            
            # Separate left and right neighbors (assuming sequential order)
            left_neighbors = [n for n in neighbor_ids if n < node_id]
            right_neighbors = [n for n in neighbor_ids if n > node_id]
            
            # Try merge combinations: left only, right only, both
            merge_attempts = []
            
            if left_neighbors:
                merge_attempts.append([left_neighbors[0], node_id])
            if right_neighbors:
                merge_attempts.append([node_id, right_neighbors[0]])
            if left_neighbors and right_neighbors:
                merge_attempts.append([left_neighbors[0], node_id, right_neighbors[0]])
            
            # Try each merge combination
            for nodes_to_merge in merge_attempts:
                print(f"DEBUG: Trying to merge nodes {nodes_to_merge}")
                
                # Reconstruct combined molecule
                combined_mol = self._reconstruct_fragment(nodes_to_merge, graph)
                if not combined_mol:
                    print(f"DEBUG: Failed to reconstruct molecule for {nodes_to_merge}")
                    continue
                
                print(f"DEBUG: Reconstructed mol with {combined_mol.GetNumAtoms()} atoms")
                
                # Count expected connections for this merged fragment
                # Get all unique neighbors of the merged set
                all_neighbors = set()
                for nid in nodes_to_merge:
                    if nid in graph.nodes:
                        node_neighbors = graph.get_neighbors(nid)
                        for neighbor_id, _ in node_neighbors:
                            if neighbor_id not in nodes_to_merge:
                                all_neighbors.add(neighbor_id)
                
                num_connections = len(all_neighbors)
                print(f"DEBUG: Expecting {num_connections} connections")
                
                # Try to match the combined fragment
                monomer = matcher.find_exact_match(combined_mol, num_connections)
                
                if monomer:
                    print(f"DEBUG: SUCCESS! Matched to {monomer.symbol}")
                    # Success! Create new merged node
                    new_node_id = min(nodes_to_merge)  # Use lowest ID
                    new_node = FragmentNode(new_node_id, combined_mol)
                    new_node.monomer = monomer
                    
                    # Merge nodes in graph
                    self._merge_nodes_in_graph(graph, nodes_to_merge, new_node)
                    
                    had_changes = True
                    break  # Stop trying other combinations for this node
                else:
                    print(f"DEBUG: No match found for merge {nodes_to_merge}")
        
        return had_changes

