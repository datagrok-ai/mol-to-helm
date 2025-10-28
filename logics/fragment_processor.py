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
            
            for atom1, atom2, linkage_type in bonds_to_cleave:
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond:
                    bond_indices.append(bond.GetIdx())
                    bond_info.append((bond.GetIdx(), atom1, atom2, linkage_type))

            if not bond_indices:
                # No valid bonds found
                node = FragmentNode(0, mol)
                node.is_n_terminal = True
                node.is_c_terminal = True
                graph.add_node(node)
                return graph

            # Fragment the molecule
            fragmented_mol = Chem.FragmentOnBonds(mol, bond_indices, addDummies=True)
            
            # Get fragments with atom mapping
            fragments_tuple = Chem.GetMolFrags(
                fragmented_mol, 
                asMols=True, 
                sanitizeFrags=True,
                frags=None,
                fragsMolAtomMapping=None
            )
            fragments = list(fragments_tuple)

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

