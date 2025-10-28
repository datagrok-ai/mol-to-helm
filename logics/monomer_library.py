from rdkit import Chem
from rdkit import RDLogger
from collections import defaultdict
from itertools import combinations
import json
import os

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

class MonomerData:
    def __init__(self):
        self.symbol = ""
        self.name = ""
        self.mol = None
        self.smiles = ""  # Original SMILES with R-groups
        self.r_groups = {}  # R-group label -> cap SMILES
        self.r_group_count = 0
        self.capped_smiles = {}  # frozenset of capped R-group labels -> canonical SMILES

    def __repr__(self):
        return f"Monomer({self.symbol}: {self.name}, R-groups: {self.r_group_count})"
    
    def generate_capped_smiles(self):
        """
        Generate all possible canonical SMILES with different R-group cappings.
        
        For a monomer with n R-groups, this generates 2^n - 1 versions:
        - Each R-group either capped (replaced with cap) or uncapped (removed)
        - We store all combinations except all-uncapped (meaningless)
        """
        if not self.mol or not self.r_groups:
            return
        
        r_group_labels = list(self.r_groups.keys())
        
        # Generate all non-empty combinations of R-groups to cap
        for r in range(1, len(r_group_labels) + 1):
            for capped_combo in combinations(r_group_labels, r):
                capped_set = frozenset(capped_combo)
                smiles = self._get_smiles_with_caps(capped_set)
                if smiles:
                    self.capped_smiles[capped_set] = smiles
    
    def _get_smiles_with_caps(self, caps_to_apply: frozenset) -> str:
        """
        Generate canonical SMILES with specific R-groups removed (uncapped).
        
        Args:
            caps_to_apply: Set of R-group labels to KEEP/CAP (others are removed)
        
        Returns:
            Canonical SMILES string with non-capped R-groups removed
        """
        try:
            # Simply remove dummy atoms that are NOT in caps_to_apply
            # This creates "open" ends where bonds were cleaved
            mol_copy = Chem.Mol(self.mol)
            
            # Find and remove dummy atoms for R-groups NOT in caps_to_apply
            atoms_to_remove = []
            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Dummy atom (R-group)
                    # Get isotope number which corresponds to R-group number
                    isotope = atom.GetIsotope()
                    if isotope > 0:
                        r_label = f"R{isotope}"
                        # If this R-group is NOT being capped, remove it
                        if r_label not in caps_to_apply:
                            atoms_to_remove.append(atom.GetIdx())
            
            # Remove dummy atoms in reverse order to preserve indices
            if atoms_to_remove:
                atoms_to_remove.sort(reverse=True)
                emol = Chem.EditableMol(mol_copy)
                for atom_idx in atoms_to_remove:
                    emol.RemoveAtom(atom_idx)
                mol_copy = emol.GetMol()
            
            if mol_copy:
                # Generate canonical SMILES
                return Chem.MolToSmiles(mol_copy, canonical=True)
            return ""
        except Exception as e:
            return ""


class MonomerLibrary:
    def __init__(self):
        self.monomers = {}
        self.smiles_to_monomer = {}
        self.name_to_monomer = {}
        self.symbol_to_monomer = {}

    def load_from_helm_json(self, json_path: str) -> None:
        if not os.path.exists(json_path):
            return

        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except Exception:
            return

        successful = 0
        for monomer_dict in data:
            try:
                monomer = self._parse_monomer(monomer_dict)
                if monomer and monomer.mol is not None:
                    self.monomers[monomer.symbol] = monomer
                    self.symbol_to_monomer[monomer.symbol] = monomer

                    clean_name = monomer.name.lower().replace(" ", "").replace("-", "").replace("_", "")
                    self.name_to_monomer[clean_name] = monomer

                    successful += 1
            except Exception:
                continue

    def _parse_monomer(self, monomer_dict: dict):
        monomer = MonomerData()
        monomer.symbol = monomer_dict.get('symbol', '')
        monomer.name = monomer_dict.get('name', '')

        if not monomer.symbol:
            return None

        smiles = monomer_dict.get('smiles', '')
        molfile = monomer_dict.get('molfile', '')

        if smiles:
            try:
                monomer.mol = Chem.MolFromSmiles(smiles)
                monomer.smiles = smiles
            except Exception:
                monomer.mol = None

        if monomer.mol is None and molfile:
            try:
                monomer.mol = Chem.MolFromMolBlock(molfile)
                if monomer.mol:
                    monomer.smiles = Chem.MolToSmiles(monomer.mol)
            except Exception:
                monomer.mol = None

        if monomer.mol is None:
            return None
        
        # Parse R-groups
        rgroups_list = monomer_dict.get('rgroups', [])
        for rgroup in rgroups_list:
            label = rgroup.get('label', '')
            cap_smiles = rgroup.get('capGroupSMILES', '')
            if label and cap_smiles:
                monomer.r_groups[label] = cap_smiles
        
        monomer.r_group_count = len(monomer.r_groups)
        
        # Generate all capped SMILES variations
        monomer.generate_capped_smiles()

        return monomer

    def find_monomer_by_fragment_smiles(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES against capped variations.
        
        Args:
            fragment_smiles: Canonical SMILES of the fragment  
            num_connections: Number of connections this fragment has in the graph
        
        Returns:
            MonomerData if match found, None otherwise
            
        Logic:
            - Fragment with N connections had N R-groups removed during fragmentation
            - capped_set contains R-groups that were KEPT (not removed)
            - Number of R-groups removed = total_r_groups - len(capped_set)
            - This should equal num_connections
        """
        # Search through all monomers
        for symbol, monomer in self.monomers.items():
            # Check if monomer has a capped SMILES variation matching this fragment
            for capped_set, capped_smiles in monomer.capped_smiles.items():
                # Calculate number of R-groups removed
                num_kept = len(capped_set)
                num_removed = monomer.r_group_count - num_kept
                
                # Match if: same SMILES and same number of removed R-groups
                if capped_smiles == fragment_smiles and num_removed == num_connections:
                    return monomer
        
        return None

    def find_monomer_by_symbol(self, symbol: str):
        return self.symbol_to_monomer.get(symbol)