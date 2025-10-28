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
        self.capped_smiles_cache = {}  # Cache: frozenset of removed R-groups -> canonical SMILES

    def __repr__(self):
        return f"Monomer({self.symbol}: {self.name}, R-groups: {self.r_group_count})"
    
    def get_capped_smiles_for_removed_rgroups(self, removed_rgroups: frozenset) -> str:
        """
        Get canonical SMILES with specific R-groups removed (lazy generation with caching).
        
        Args:
            removed_rgroups: frozenset of R-group labels that were removed (e.g., {'R1', 'R2'})
        
        Returns:
            Canonical SMILES with those R-groups removed, or empty string on error
            
        Example:
            For monomer with R1, R2:
            - get_capped_smiles_for_removed_rgroups({'R1'}) → SMILES with R1 removed, R2 kept
            - get_capped_smiles_for_removed_rgroups({'R2'}) → SMILES with R2 removed, R1 kept
            - get_capped_smiles_for_removed_rgroups({'R1', 'R2'}) → SMILES with both removed
        """
        # Check cache first
        if removed_rgroups in self.capped_smiles_cache:
            return self.capped_smiles_cache[removed_rgroups]
        
        # Generate on demand
        smiles = self._get_smiles_with_rgroups_removed(removed_rgroups)
        
        # Cache for future use
        self.capped_smiles_cache[removed_rgroups] = smiles
        
        return smiles
    
    def _get_smiles_with_rgroups_removed(self, removed_rgroups: frozenset) -> str:
        """
        Generate canonical SMILES with specific R-groups removed.
        
        Args:
            removed_rgroups: Set of R-group labels to REMOVE (e.g., {'R1', 'R2'})
        
        Returns:
            Canonical SMILES string with those R-groups removed
        """
        try:
            mol_copy = Chem.Mol(self.mol)
            
            # Find and remove dummy atoms for specified R-groups
            atoms_to_remove = []
            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Dummy atom (R-group)
                    isotope = atom.GetIsotope()
                    if isotope > 0:
                        r_label = f"R{isotope}"
                        # Remove if this R-group is in the removed set
                        if r_label in removed_rgroups:
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

        return monomer

    def find_monomer_by_fragment_smiles(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES with on-demand R-group removal.
        
        Args:
            fragment_smiles: Canonical SMILES of the fragment  
            num_connections: Number of connections this fragment has in the graph
        
        Returns:
            MonomerData if match found, None otherwise
            
        Logic:
            - Fragment with N connections → N R-groups were removed during fragmentation
            - For monomer with M R-groups, try all C(M,N) combinations of which N R-groups were removed
            - Generate SMILES for each combination on-demand (with caching)
            
        Example:
            Fragment has 1 connection, monomer has R1, R2:
            - Try removing R1 → check if SMILES matches
            - Try removing R2 → check if SMILES matches
        """
        # Search through all monomers
        for symbol, monomer in self.monomers.items():
            # Skip if monomer doesn't have enough R-groups
            if monomer.r_group_count < num_connections:
                continue
            
            # Generate all combinations of num_connections R-groups that could have been removed
            r_group_labels = list(monomer.r_groups.keys())
            
            # For each combination of R-groups that could have been removed
            for removed_combo in combinations(r_group_labels, num_connections):
                removed_set = frozenset(removed_combo)
                
                # Generate SMILES with these R-groups removed (lazy, cached)
                candidate_smiles = monomer.get_capped_smiles_for_removed_rgroups(removed_set)
                
                # Check if it matches the fragment
                if candidate_smiles == fragment_smiles:
                    return monomer
        
        return None

    def find_monomer_by_symbol(self, symbol: str):
        return self.symbol_to_monomer.get(symbol)