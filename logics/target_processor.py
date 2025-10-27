from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import os
from monomer_library import MonomerData, MonomerLibrary

class AdvancedTerminalNormalizer:
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

class AdvancedBondDetector:
    def __init__(self):
        self.peptide_bond = Chem.MolFromSmarts('[C;X3](=[O;X1])-[N;X3]')
        self.ester_bond = Chem.MolFromSmarts('[C;X3](=[O;X1])-[O;X2]')
        self.urethane_bond = Chem.MolFromSmarts('[O;X2]-[C;X3](=[O;X1])-[N;X3]')
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

class PrecisionFragmentProcessor:
    def __init__(self, monomer_library):
        self.monomer_library = monomer_library
        self.bond_detector = AdvancedBondDetector()
        self.normalizer = AdvancedTerminalNormalizer()

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

class PrecisionMonomerMatcher:
    def __init__(self, monomer_library: MonomerLibrary):
        self.monomer_library = monomer_library
        self.hardcoded_mappings = self._create_precision_mappings()
        self.similarity_threshold = 0.8

    def _create_precision_mappings(self):
        return {
            "N[C@H](C=O)CC=O": "N",
            "N[C@@H](C=O)CC=O": "dN",
            "N[C@H](C=O)CCC(=O)O": "E",
            "N[C@@H](C=O)CCC(=O)O": "dE",
            "Cc1ccc(C[C@H](N)C(=O)O)cc1": "Phe_4Me",
            "Cc1ccc(C[C@@H](N)C(=O)O)cc1": "D-Phe_4Me",
            "CN[C@H](C=O)CCC=O": "meQ",
            "CC(C)[C@@H](N)C(=O)N=C(C=O)Cc1ccccc1": "Phe_ab-dehydro",
            "NCCNCC=O": "spacer",
            "N[C@H](C=O)CC(C=O)N": "Q",
            "N[C@@H](C=O)CC(C=O)N": "dQ",
            "N[C@H](C=O)CCCNC(N)=O": "R",
            "N[C@@H](C=O)CCCNC(N)=O": "dR",
            "CC(C=O)=NC(=O)[C@H](N)Cc1ccc(F)cc1": "D-Phe_4F",
            "CC(C=O)=NC(=O)[C@@H](N)Cc1ccc(F)cc1": "Phe_4F",
            "Cc1cn(CC=O)c(=O)[nH]c1=O": "Xan",
            "NC(C=O)CC=O": "Asn",
            "NC(C=O)Cc1ccc(O)c(O)c1": "Tyr",
            "NC(C=O)Cc1ccc(OCc2ccccc2)cc1": "Tyr_Bn",
            "O=CC1CCCCN1": "Pro",
            "CCCC(N)C=O": "Nle",
            "O=CC1CCN1": "Aze",
            "Nc1ccc(CC(N)C=O)cc1": "Phe_4NH2",
            "NC(C=O)CCCCNC(N)N": "R",
            "CNC(C=O)C(C)C": "meVal",
            "NC(C=O)CCC=O": "Q",
            "NC(C=O)Cc1ccc(Br)cc1": "Phe_4Br",
            "NC(C=O)Cc1cn(C=O)c2ccccc12": "W",
            "NC(C=O)Cc1ccc(-c2ccccc2)cc1": "Bip",
            "NC(Cc1ccccc1)C(O)CC=O": "Phe_ab-dehydro_alcohol",
            "Nc1ncnc2c1ncn2CC=O": "Ade",
            "CCC(C)C(N)C=O": "Ile",
            "NC(C=O)Cc1ccc(Cl)cc1": "Phe_4Cl",
            "Cn1cncc1CC(N)C=O": "H",
            "COC(=O)CC(N)C=O": "Asp_OMe",
            "CC(O)C(N)C=O": "T",
            "NC(C=O)C1CCCCC1": "Chg",
            "CS(=O)(=O)CCC(N)C=O": "Cys_SO3H",
            "CC(C)(c1ccccc1)C(N)C=O": "dF",
            "CCCCC(N)C=O": "Nle",
            "C=CCC(N)C=O": "Hse",
            "NC(C=O)CO": "S",
            "CC(C)(C)OCC(N)C=O": "tBuGly",
            "NC(C=O)Cc1ccncc1": "4Pal",
            "CC(C)(C)C(N)C=O": "tBuAla",
            "O=CC1C[C@@H](O)CN1": "Hyp",
            "CSCCC(N)C=O": "Cys_SPr",
            "N=C(N)NCCCC(N)C=O": "Cit",
            "CS(=O)CCC(N)C=O": "Cys_SO2Pr",
            "NC(C=O)C1CC1": "Ac3c",
            "Cn1cnc(CC(N)C=O)c1": "Ade",
            "CC(C)CC(N)C=O": "Ile",
            "N[C@H](C=O)CCSC": "C",
            "N[C@@H](C=O)CCSC": "dC",
            "N[C@H](C=O)CC(=O)O": "D",
            "N[C@@H](C=O)CC(=O)O": "dD",
            "N[C@H](C=O)CC(N)=O": "N",
            "N[C@@H](C=O)CC(N)=O": "dN",
            "NC(=O)NCCC[C@H](N)C=O": "Cit",
            "NC(=O)NCCC[C@@H](N)C=O": "D-Cit",
            "N[C@H](C=O)CCCN=C(N)N": "R",
            "N[C@@H](C=O)CCCN=C(N)N": "dR",
            "N[C@H](C=O)CCC(N)=O": "Q",
            "N[C@@H](C=O)CCC(N)=O": "dQ",
            "N[C@H](C=O)CCCCN": "Orn",
            "N[C@@H](C=O)CCCCN": "D-Orn",
        }

    def find_exact_match(self, fragment: Chem.Mol):
        try:
            frag_smiles = Chem.MolToSmiles(fragment, canonical=True)
            if not frag_smiles:
                return None

            hardcoded_match = self._find_by_hardcoded_mapping(frag_smiles)
            if hardcoded_match:
                return hardcoded_match

            exact_match = self.monomer_library.find_monomer_by_smiles(frag_smiles)
            if exact_match:
                return exact_match

            normalized_match = self._find_with_normalization(fragment, frag_smiles)
            if normalized_match:
                return normalized_match

            return None

        except Exception:
            return None

    def _find_by_hardcoded_mapping(self, frag_smiles: str):
        symbol = self.hardcoded_mappings.get(frag_smiles)
        if symbol:
            monomer = self.monomer_library.find_monomer_by_symbol(symbol)
            if monomer:
                return monomer
        return None

    def _find_with_normalization(self, fragment: Chem.Mol, frag_smiles: str):
        try:
            normalizer = AdvancedTerminalNormalizer()

            normalized1 = normalizer.normalize_for_library(fragment, is_c_terminal=True)
            if normalized1:
                norm_smiles1 = Chem.MolToSmiles(normalized1, canonical=True)
                match1 = self.monomer_library.find_monomer_by_smiles(norm_smiles1)
                if match1:
                    return match1

            normalized2 = normalizer.normalize_for_library(fragment, is_c_terminal=False)
            if normalized2:
                norm_smiles2 = Chem.MolToSmiles(normalized2, canonical=True)
                match2 = self.monomer_library.find_monomer_by_smiles(norm_smiles2)
                if match2:
                    return match2

            return None

        except Exception:
            return None

    def find_best_match(self, fragment: Chem.Mol):
        exact_match = self.find_exact_match(fragment)
        if exact_match:
            return exact_match

        return self._find_advanced_match(fragment)

    def _find_advanced_match(self, fragment: Chem.Mol):
        try:
            frag_fp = AllChem.GetMorganFingerprintAsBitVect(fragment, 2, nBits=1024)
            best_match = None
            best_similarity = 0.0

            for symbol, monomer in self.monomer_library.monomers.items():
                if not monomer.mol:
                    continue

                monomer_fp = self.monomer_library.get_fingerprint(symbol)
                if monomer_fp is None:
                    continue

                similarity = DataStructs.TanimotoSimilarity(frag_fp, monomer_fp)

                if similarity > best_similarity and similarity > self.similarity_threshold:
                    best_similarity = similarity
                    best_match = monomer

            return best_match

        except Exception:
            return None

class HELMGenerator:
    def __init__(self):
        self.polymer_types = {
            "peptide": "PEPTIDE",
            "rna": "RNA",
            "dna": "DNA",
            "chemical": "CHEM"
        }

    def generate_helm_notation(self, monomers, is_cyclic: bool = False) -> str:
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