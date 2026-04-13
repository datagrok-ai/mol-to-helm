from rdkit import Chem
from fragment_graph import FragmentGraph, FragmentNode, FragmentLink, LinkageType


class BondDetector:
    #GENERALIZATION ITEM: BOND PATTERNS SHOULD BE DERIVED FROM LIBRARY
    def __init__(self):
        # True peptide bond: C and N both in backbone (each bonded to carbons)
        # First carbon can be aliphatic or aromatic (for amino acids like NMe2Abz)
        # Carbonyl carbon is sp2 (X3)
        # Exclude if carbonyl is in a small ring (r5 or r6) to avoid cleaving lactams like Pyr
        # !r5 = not in 5-membered ring, !r6 = not in 6-membered ring
        # This preserves lactams but allows large macrocycles and proline (C=O outside ring)
        # Nitrogen can be X2 (proline, imino) or X3 (standard amino, N-methyl)
        # N-C bond can be single (-) or double (=) for imine bonds in dehydro amino acids
        # Alpha carbon after N can be sp3 (X4) or sp2 (X3) for dehydroamino acids, or aromatic (#6 includes both)
        self.peptide_bond = Chem.MolFromSmarts('[#6]-[C;X3;!r5;!r6](=[O;X1])-[N;X2,X3]~[#6;X3,X4]')
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

            # Filter out internal amide bonds in CHEM linkers like FC01.
            # FC01 pattern: C(=O)-N-ArRing-N-C(=O) — two amide bonds connect to the
            # same aromatic ring via the alpha-C position (match[4]).
            # Real aromatic amino acids (3Abz) only have ONE such bond per ring.
            skip_indices = set()
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()

            # Map: ring_frozenset -> list of match indices where match[4] is aromatic on that ring
            # Only consider small rings (5-6 atoms) — large macrocycles should not be filtered
            ring_to_matches = {}
            for i, match in enumerate(matches):
                if len(match) < 5:
                    continue
                alpha_c_atom = mol.GetAtomWithIdx(match[4])
                if alpha_c_atom.GetIsAromatic():
                    for ring in rings:
                        if match[4] in ring and len(ring) <= 6:
                            ring_key = frozenset(ring)
                            if ring_key not in ring_to_matches:
                                ring_to_matches[ring_key] = []
                            ring_to_matches[ring_key].append(i)
                            break

            # If 2+ matches share an aromatic ring at their alpha-C position, skip them
            for ring_key, match_indices in ring_to_matches.items():
                if len(match_indices) >= 2:
                    skip_indices.update(match_indices)

            for i, match in enumerate(matches):
                if len(match) >= 5 and i not in skip_indices:
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


    def _find_staple_sidechain_bonds(self, mol, existing_bonds):
        """
        Find bonds to cleave in staple macrocycles.

        Detects large macrocycles (>10 atoms) that contain both quaternary
        alpha-methyl carbons (staple monomers like R8, S5) and non-peptide
        linker bonds. Cleaves at C=C double bonds within these macrocycles,
        which are the signature of olefin metathesis (RCMtrans/RCMcis) linkers.
        Also cleaves at C-S thioether bonds for FC01-type linkers.
        """
        ring_info = mol.GetRingInfo()
        large_rings = [set(ring) for ring in ring_info.AtomRings() if len(ring) > 10]
        if not large_rings:
            return []

        # Quick check: must have quaternary alpha-methyl C in a large ring
        quat_alpha = Chem.MolFromSmarts('[N;X2,X3]-[C;X4;H0](-[C;X3](=[O;X1]))-[CH3]')
        if not quat_alpha or not mol.HasSubstructMatch(quat_alpha):
            return []

        existing_atom_pairs = set()
        for a1, a2, _ in existing_bonds:
            existing_atom_pairs.add((min(a1, a2), max(a1, a2)))

        additional_bonds = []
        seen = set()

        for ring in large_rings:
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 not in ring or a2 not in ring:
                    continue

                # For C=C double bonds (RCMtrans/RCMcis): cleave the single bonds
                # one hop away on each side. RCMtrans structure: ...R3chain-CH2-C=C-CH2-R3chain...
                # Cleaving at CH2-R3chain boundaries keeps the correct R3 chain length.
                if (bond.GetBondTypeAsDouble() >= 2 and
                        mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and
                        mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                    # For each C=C atom, find its other ring neighbor and cleave that bond
                    for cc_atom_idx in (a1, a2):
                        other_cc = a2 if cc_atom_idx == a1 else a1
                        cc_atom = mol.GetAtomWithIdx(cc_atom_idx)
                        for nbr in cc_atom.GetNeighbors():
                            nbr_idx = nbr.GetIdx()
                            if nbr_idx == other_cc or nbr_idx not in ring:
                                continue
                            # nbr is the CH2 between C=C and R3 chain
                            # Find nbr's other ring neighbor (the R3 chain atom)
                            for nbr2 in nbr.GetNeighbors():
                                nbr2_idx = nbr2.GetIdx()
                                if nbr2_idx == cc_atom_idx or nbr2_idx not in ring:
                                    continue
                                pair = (min(nbr_idx, nbr2_idx), max(nbr_idx, nbr2_idx))
                                if pair not in existing_atom_pairs and pair not in seen:
                                    seen.add(pair)
                                    additional_bonds.append((nbr_idx, nbr2_idx, LinkageType.UNKNOWN))

                # Cleave C-S thioether bonds (FC01-type linker)
                atom1 = mol.GetAtomWithIdx(a1)
                atom2 = mol.GetAtomWithIdx(a2)
                if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 16) or
                        (atom1.GetAtomicNum() == 16 and atom2.GetAtomicNum() == 6)):
                    pair = (min(a1, a2), max(a1, a2))
                    if pair not in existing_atom_pairs and pair not in seen:
                        seen.add(pair)
                        additional_bonds.append((a1, a2, LinkageType.UNKNOWN))

        return additional_bonds

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

            # Detect R3 side-chain bonds for staple monomers (R8, S5, etc.)
            r3_bonds = self._find_staple_sidechain_bonds(mol, bonds_to_cleave)
            if r3_bonds:
                bonds_to_cleave.extend(r3_bonds)

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
            
            # Get fragments as molecules
            fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
            
            # Get atom mappings separately (which original atoms are in which fragment)
            atom_mappings = Chem.GetMolFrags(fragmented_mol, asMols=False, fragsMolAtomMapping=True)
            
            # Store bond cleavage info for recovery - we'll use this to selectively re-fragment
            graph.cleaved_bond_indices = bond_indices
            graph.bond_info = bond_info
            graph.atom_mappings = atom_mappings

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

            # Create links between fragments based on the actual cleaved bonds
            # Build mapping: original atom index → (fragment_idx, new_atom_idx_in_fragment)
            atom_to_fragment_and_idx = {}
            for frag_idx, original_atom_indices in enumerate(atom_mappings):
                for new_idx_in_frag, original_atom_idx in enumerate(original_atom_indices):
                    atom_to_fragment_and_idx[original_atom_idx] = (frag_idx, new_idx_in_frag)
            
            
            # For each cleaved bond, determine which fragments it connects
            link_count = 0
            for bond_idx, atom1_orig, atom2_orig, linkage_type in bond_info:
                # Find which fragments contain these atoms and their new indices
                frag1_info = atom_to_fragment_and_idx.get(atom1_orig)
                frag2_info = atom_to_fragment_and_idx.get(atom2_orig)
                
                if frag1_info is None or frag2_info is None:
                    print(f"DEBUG: Skipping bond atoms {atom1_orig}-{atom2_orig}: not found in fragments")
                    continue
                
                frag1, atom1_new = frag1_info
                frag2, atom2_new = frag2_info
                    
                # Create link even if both atoms are in same fragment (internal bond like in Phe_4Sdihydroorotamido)
                # This creates a "self-link" that will be used during recovery to reconstruct the monomer
                link = FragmentLink(frag1, frag2, linkage_type, 
                                   from_atom_idx=atom1_new, to_atom_idx=atom2_new)
                graph.add_link(link)
                link_count += 1
                
                if frag1 == frag2:
                    print(f"DEBUG: Link {link_count}: {linkage_type.value.upper()} SELF-LINK frag{frag1} "
                          f"orig_atoms({atom1_orig}<->{atom2_orig}) frag_atoms({atom1_new}<->{atom2_new})")
                else:
                    print(f"DEBUG: Link {link_count}: {linkage_type.value.upper()} frag{frag1}<->frag{frag2} "
                          f"orig_atoms({atom1_orig}<->{atom2_orig}) frag_atoms({atom1_new}<->{atom2_new})")
            
            print(f"DEBUG: Created {link_count} links total")

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

    def _reconstruct_fragment_with_links(self, node_ids: list, graph: FragmentGraph, 
                                         links_to_exclude: list) -> Chem.Mol:
        """
        Reconstruct a molecule by combining multiple fragment nodes, using link information.
        
        Args:
            node_ids: List of node IDs to merge
            graph: The fragment graph
            links_to_exclude: List of FragmentLink objects connecting the nodes to merge
        
        Returns:
            Combined RDKit molecule, or None if reconstruction fails
        """
        if not node_ids or not hasattr(graph, 'original_mol'):
            return None
        
        if not hasattr(graph, 'cleaved_bond_indices') or not hasattr(graph, 'bond_info'):
            return None
        
        try:
            # Find which bond indices correspond to the links we want to exclude
            bonds_to_exclude_indices = []
            
            for link in links_to_exclude:
                # Find the bond_info entry that matches this link's original atoms
                # We need to find which bond connected these fragments
                for bond_list_idx, (bond_idx, atom1, atom2, linkage_type) in enumerate(graph.bond_info):
                    # Check if this bond connects the fragments in this link
                    if hasattr(graph, 'atom_mappings'):
                        # Find which fragments contain these atoms
                        frag1 = None
                        frag2 = None
                        for frag_idx, atom_indices in enumerate(graph.atom_mappings):
                            if atom1 in atom_indices:
                                frag1 = frag_idx
                            if atom2 in atom_indices:
                                frag2 = frag_idx
                        
                        # If this bond connects the two fragments in the link, exclude it
                        if (frag1 == link.from_node_id and frag2 == link.to_node_id) or \
                           (frag1 == link.to_node_id and frag2 == link.from_node_id):
                            bonds_to_exclude_indices.append(bond_list_idx)
                            print(f"DEBUG: Excluding {linkage_type.value} bond at index {bond_list_idx} (atoms {atom1}<->{atom2})")
                            break
            
            # Create new bond list excluding the bonds we want to keep
            new_bond_indices = [
                bond_idx for i, bond_idx in enumerate(graph.cleaved_bond_indices)
                if i not in bonds_to_exclude_indices
            ]
            
            print(f"DEBUG reconstruct: Original had {len(graph.cleaved_bond_indices)} cleaved bonds, "
                  f"excluding {len(bonds_to_exclude_indices)} bonds, new list has {len(new_bond_indices)} bonds")
            
            # Re-fragment with the modified bond list
            if not new_bond_indices:
                # No bonds to cleave - return whole molecule
                return graph.original_mol
            
            fragmented_mol = Chem.FragmentOnBonds(graph.original_mol, new_bond_indices, addDummies=True)
            fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
            new_atom_mappings = Chem.GetMolFrags(fragmented_mol, asMols=False, fragsMolAtomMapping=True)
            
            # Find which new fragment contains atoms from our target nodes
            # Look for the fragment that contains atoms from the first node we want to merge
            sorted_nodes = sorted(node_ids)
            first_node_atoms = set(graph.atom_mappings[sorted_nodes[0]])
            
            target_fragment_idx = None
            for new_frag_idx, new_atoms in enumerate(new_atom_mappings):
                # Check if this new fragment contains any atoms from our first target node
                if first_node_atoms & set(new_atoms):
                    target_fragment_idx = new_frag_idx
                    break
            
            print(f"DEBUG reconstruct: Got {len(fragments)} fragments after re-fragmentation, "
                  f"target_fragment_idx={target_fragment_idx}")
            
            if target_fragment_idx is not None and target_fragment_idx < len(fragments):
                clean_frag = self._clean_fragment(fragments[target_fragment_idx])
                return clean_frag if clean_frag else fragments[target_fragment_idx]
            
            return None
        
        except Exception as e:
            print(f"DEBUG reconstruct: Exception: {e}")
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
        Try to recover unmatched fragments by merging with neighbors based on graph links.
        Returns True if any merges were successful.
        """
        # Identify unmatched nodes
        unmatched_nodes = []
        for node_id, node in graph.nodes.items():
            if node.monomer and node.monomer.symbol.startswith("X"):
                unmatched_nodes.append(node_id)
        
        if not unmatched_nodes:
            return False
        
        had_changes = False
        
        # Try to recover each unmatched node
        for node_id in unmatched_nodes:
            # Check if node still exists (might have been merged already)
            if node_id not in graph.nodes:
                continue
            
            # Get neighbors from graph links (returns list of (neighbor_id, linkage_type))
            neighbors = graph.get_neighbors(node_id)
            
            if not neighbors:
                continue
            
            # Try merging with each individual neighbor first
            for neighbor_id, linkage_type in neighbors:
                if neighbor_id not in graph.nodes:
                    continue
                    
                nodes_to_merge = sorted([node_id, neighbor_id])
                
                # Find the links between nodes we're merging
                links_to_exclude = []
                for link in graph.links:
                    from_in = link.from_node_id in nodes_to_merge
                    to_in = link.to_node_id in nodes_to_merge
                    if from_in and to_in:
                        links_to_exclude.append(link)
                
                # Reconstruct combined molecule
                combined_mol = self._reconstruct_fragment_with_links(nodes_to_merge, graph, links_to_exclude)
                if not combined_mol:
                    print(f"DEBUG: Failed to reconstruct molecule for {nodes_to_merge}")
                    continue
                
                print(f"DEBUG: Reconstructed mol with {combined_mol.GetNumAtoms()} atoms")
                
                # Count expected connections for this merged fragment
                # Get all unique neighbors of the merged set (excluding internal connections)
                all_neighbors = set()
                for nid in nodes_to_merge:
                    if nid in graph.nodes:
                        node_neighbors = graph.get_neighbors(nid)
                        for neighbor_id, _ in node_neighbors:
                            if neighbor_id not in nodes_to_merge:
                                all_neighbors.add(neighbor_id)
                
                num_connections = len(all_neighbors)
                
                # Try to match the combined fragment (exact match only)
                monomer = matcher.find_exact_match(combined_mol, num_connections)

                if monomer:
                    # Success! Create new merged node
                    new_node_id = min(nodes_to_merge)
                    new_node = FragmentNode(new_node_id, combined_mol)
                    new_node.monomer = monomer
                    
                    # Merge nodes in graph
                    self._merge_nodes_in_graph(graph, nodes_to_merge, new_node)
                    
                    had_changes = True
                    break  # Stop trying other neighbors for this node
            
            if had_changes:
                break  # Restart from beginning after a successful merge
        
        return had_changes
    
    def recover_unmatched_with_stereo_agnostic(self, graph: FragmentGraph, matcher) -> int:
        """
        Separate recovery procedure: Try to match remaining unmatched fragments 
        using stereochemistry-agnostic comparison.
        
        This handles poor quality input data where stereochemistry is not assigned.
        Only called after regular recovery attempts have finished.
        
        Args:
            graph: FragmentGraph with some unmatched nodes
            matcher: MonomerMatcher instance
        
        Returns:
            Number of fragments that were successfully matched
        """
        from rdkit import Chem
        
        # Find all unmatched nodes (nodes with mock/unknown monomers)
        unmatched_nodes = []
        for node_id, node in graph.nodes.items():
            if node.monomer and (node.monomer.symbol.startswith('X') or 
                                 node.monomer.name.startswith('Unknown')):
                unmatched_nodes.append(node_id)
        
        if not unmatched_nodes:
            return 0
        
        print(f"DEBUG: Attempting stereo-agnostic recovery for {len(unmatched_nodes)} unmatched nodes")
        
        matched_count = 0
        
        for node_id in unmatched_nodes:
            if node_id not in graph.nodes:
                continue
            
            node = graph.nodes[node_id]
            
            # Get fragment SMILES
            fragment_smiles = Chem.MolToSmiles(node.mol, canonical=True)
            
            # Count connections
            neighbors = graph.get_neighbors(node_id)
            num_connections = len(neighbors)
            
            # Try stereo-agnostic matching
            monomer = matcher.monomer_library.find_monomer_by_fragment_smiles_no_stereo(
                fragment_smiles, num_connections
            )
            
            if monomer:
                print(f"DEBUG: Stereo-agnostic match for node {node_id}: {monomer.symbol}")
                node.monomer = monomer
                matched_count += 1
            else:
                print(f"DEBUG: No stereo-agnostic match for node {node_id}")
        
        return matched_count

    def recover_unmatched_by_merging_stereo_agnostic(self, graph: FragmentGraph, matcher) -> bool:
        """
        Final recovery pass: merge pairs of BOTH-unmatched neighbor fragments and
        try stereo-agnostic matching on the combined result.

        This handles monomers like Phe_4Sdihydroorotamido that have internal amide
        bonds which get incorrectly cleaved, producing two unmatched fragments.

        Only merges when BOTH fragments in a pair are unmatched — never touches
        already-matched nodes to avoid regressions.

        Returns True if any merges were successful.
        """
        def _is_unmatched(node):
            return (node.monomer and
                    (node.monomer.symbol.startswith('X') or
                     node.monomer.name.startswith('Unknown')))

        unmatched_ids = [nid for nid, node in graph.nodes.items() if _is_unmatched(node)]
        if not unmatched_ids:
            return False

        had_changes = False

        for node_id in unmatched_ids:
            if node_id not in graph.nodes:
                continue
            if not _is_unmatched(graph.nodes[node_id]):
                continue

            neighbors = graph.get_neighbors(node_id)
            for neighbor_id, linkage_type in neighbors:
                if neighbor_id not in graph.nodes:
                    continue
                # Only merge with another unmatched neighbor
                if not _is_unmatched(graph.nodes[neighbor_id]):
                    continue

                nodes_to_merge = sorted([node_id, neighbor_id])

                # Find internal links between the merge candidates
                links_to_exclude = []
                for link in graph.links:
                    if (link.from_node_id in nodes_to_merge and
                            link.to_node_id in nodes_to_merge):
                        links_to_exclude.append(link)

                combined_mol = self._reconstruct_fragment_with_links(
                    nodes_to_merge, graph, links_to_exclude)
                if not combined_mol:
                    continue

                # Count external connections for merged fragment
                all_neighbors = set()
                for nid in nodes_to_merge:
                    if nid in graph.nodes:
                        for nbr_id, _ in graph.get_neighbors(nid):
                            if nbr_id not in nodes_to_merge:
                                all_neighbors.add(nbr_id)
                num_connections = len(all_neighbors)

                # Try exact match first, then stereo-agnostic
                monomer = matcher.find_exact_match(combined_mol, num_connections)
                if not monomer:
                    combined_smiles = Chem.MolToSmiles(combined_mol, canonical=True)
                    monomer = matcher.monomer_library.find_monomer_by_fragment_smiles_no_stereo(
                        combined_smiles, num_connections)

                if monomer:
                    new_node_id = min(nodes_to_merge)
                    new_node = FragmentNode(new_node_id, combined_mol)
                    new_node.monomer = monomer
                    self._merge_nodes_in_graph(graph, nodes_to_merge, new_node)
                    had_changes = True
                    break  # Restart from outer loop

            if had_changes:
                break

        return had_changes

