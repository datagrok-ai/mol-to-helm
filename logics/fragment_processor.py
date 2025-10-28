from rdkit import Chem


class TerminalNormalizer:
    def __init__(self):
        self.carboxyl_pattern = Chem.MolFromSmarts('[C;X3](=[O;X1])[O;H1]')
        self.amide_pattern = Chem.MolFromSmarts('[C;X3](=[O;X1])[N;X3]')
        self.primary_amine = Chem.MolFromSmarts('[N;X3;H2]')

    def normalize_for_library(self, mol: Chem.Mol, is_c_terminal: bool = False) -> Chem.Mol:
        try:
            if mol is None:
                return mol

            mol_copy = Chem.Mol(mol)

            if is_c_terminal:
                mol_copy = self._convert_carboxyl_to_amide_like(mol_copy)

            mol_copy = self._normalize_primary_amine(mol_copy)

            smiles = Chem.MolToSmiles(mol_copy, canonical=True)
            return Chem.MolFromSmiles(smiles)

        except Exception:
            return mol

    def _convert_carboxyl_to_amide_like(self, mol: Chem.Mol) -> Chem.Mol:
        try:
            matches = mol.GetSubstructMatches(self.carboxyl_pattern)
            if not matches:
                return mol

            for match in matches:
                if len(match) >= 2:
                    c_atom_idx = match[0]
                    o_atom_idx = match[1]

                    c_atom = mol.GetAtomWithIdx(c_atom_idx)
                    o_atom = mol.GetAtomWithIdx(o_atom_idx)

                    if o_atom.GetTotalNumHs() == 1:
                        emol = Chem.EditableMol(mol)
                        emol.RemoveAtom(o_atom_idx)
                        new_mol = emol.GetMol()
                        return new_mol

            return mol
        except Exception:
            return mol

    def _normalize_primary_amine(self, mol: Chem.Mol) -> Chem.Mol:
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
            return Chem.MolFromSmiles(smiles)
        except Exception:
            return mol


class BondDetector:
    def __init__(self):
        self.peptide_bond = Chem.MolFromSmarts('[C;X3](=[O;X1])-[N;X3]')
        self.disulfide_bond = Chem.MolFromSmarts('[S;X2]-[S;X2]')

    def find_cleavable_bonds(self, mol: Chem.Mol):
        try:
            all_bonds = []

            peptide_bonds = self._find_peptide_bonds(mol)
            all_bonds.extend(peptide_bonds)

            ordered_bonds = self._order_bonds_from_n_to_c(mol, all_bonds)

            return ordered_bonds

        except Exception:
            return []

    def _find_peptide_bonds(self, mol: Chem.Mol):
        bonds = []
        try:
            matches = mol.GetSubstructMatches(self.peptide_bond)
            for match in matches:
                if len(match) >= 3:
                    c_atom = match[0]
                    n_atom = match[2]
                    bonds.append((c_atom, n_atom))
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
            primary_amine = Chem.MolFromSmarts('[N;H2;X3][C;X4]')
            matches = mol.GetSubstructMatches(primary_amine)
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
        self.normalizer = TerminalNormalizer()

    def process_molecule(self, mol: Chem.Mol, is_cyclic: bool = False):
        try:
            bonds_to_cleave = self.bond_detector.find_cleavable_bonds(mol)

            if not bonds_to_cleave:
                return [mol]

            bond_indices = []
            for atom1, atom2 in bonds_to_cleave:
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond:
                    bond_indices.append(bond.GetIdx())

            if not bond_indices:
                return [mol]

            fragmented_mol = Chem.FragmentOnBonds(mol, bond_indices, addDummies=True)
            fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)

            cleaned_fragments = []
            for i, frag in enumerate(fragments):
                clean_frag = self._clean_fragment(frag)
                if clean_frag and clean_frag.GetNumAtoms() >= 3:
                    is_c_terminal = (i == len(fragments) - 1)
                    normalized_frag = self.normalizer.normalize_for_library(clean_frag, is_c_terminal)
                    if normalized_frag:
                        cleaned_fragments.append(normalized_frag)

            return cleaned_fragments

        except Exception:
            return [mol]

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

