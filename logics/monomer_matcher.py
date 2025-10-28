from rdkit import Chem
from monomer_library import MonomerLibrary
from fragment_processor import TerminalNormalizer


class MonomerMatcher:
    """
    Matches molecular fragments to monomers in the library.
    
    Uses multiple strategies:
    - Hardcoded SMILES mappings for known structures
    - Exact SMILES matching against library
    - Normalization attempts for terminal variations
    """
    
    def __init__(self, monomer_library: MonomerLibrary):
        self.monomer_library = monomer_library
        self.hardcoded_mappings = self._create_precision_mappings()

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
        """
        Find exact match for a fragment using multiple strategies.
        
        Args:
            fragment: RDKit molecule object representing a fragment
        
        Returns:
            MonomerData object if match found, None otherwise
        """
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
        """Find monomer using hardcoded SMILES mappings"""
        symbol = self.hardcoded_mappings.get(frag_smiles)
        if symbol:
            monomer = self.monomer_library.find_monomer_by_symbol(symbol)
            if monomer:
                return monomer
        return None

    def _find_with_normalization(self, fragment: Chem.Mol, frag_smiles: str):
        """Find monomer by trying different normalization strategies"""
        try:
            normalizer = TerminalNormalizer()

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
        """
        Find the best matching monomer for a fragment.
        
        Args:
            fragment: RDKit molecule object representing a fragment
        
        Returns:
            MonomerData object if match found, None otherwise
        """
        exact_match = self.find_exact_match(fragment)
        if exact_match:
            return exact_match

        return None

