from rdkit import Chem
from rdkit import RDLogger
from collections import defaultdict
from itertools import combinations
import json
import os
import re

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

def remove_stereochemistry_from_smiles(smiles: str) -> str:
    """
    Remove stereochemistry markers from SMILES string.
    Converts [C@@H], [C@H] to C, etc.
    
    This is used for matching when input molecules don't have stereochemistry defined.
    """
    if not smiles:
        return smiles
    
    # Remove @ symbols (stereochemistry markers)
    # Pattern: [@]+ inside brackets
    smiles_no_stereo = re.sub(r'(@+)', '', smiles)
    
    # Also remove H when it's explicit in brackets like [C@@H] -> [C] -> C
    # But we need to be careful not to remove H from [H] or CH3
    # After removing @, we might have [CH] which should become C
    smiles_no_stereo = re.sub(r'\[([A-Z][a-z]?)H\]', r'\1', smiles_no_stereo)
    # Handle [C] -> C (single atoms in brackets with no other info)
    smiles_no_stereo = re.sub(r'\[([A-Z][a-z]?)\]', r'\1', smiles_no_stereo)
    
    return smiles_no_stereo

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
        Generate canonical SMILES with specific R-groups removed and others capped.
        
        Args:
            removed_rgroups: Set of R-group labels where bonds were broken (e.g., {'R1', 'R2'})
        
        Returns:
            Canonical SMILES string matching fragment structure
            
        Logic:
            - R-groups in removed_rgroups: Remove dummy atom (bond was broken)
            - R-groups NOT in removed_rgroups: Cap according to library (e.g., OH, H)
            - Final SMILES has NO [*:X] markers to match fragment SMILES
        """
        try:
            mol_copy = Chem.Mol(self.mol)
            
            # Identify which R-groups to cap vs remove
            kept_rgroups = set(self.r_groups.keys()) - removed_rgroups
            
            # Process each R-group
            # IMPORTANT: SMILES [*:1] uses atom map numbers, not isotopes!
            dummy_atoms_to_process = []
            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Dummy atom (R-group)
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        r_label = f"R{map_num}"
                        if r_label in removed_rgroups:
                            # Just remove this dummy atom
                            dummy_atoms_to_process.append((atom.GetIdx(), 'remove', r_label))
                        elif r_label in kept_rgroups:
                            # Need to cap this R-group
                            cap_smiles = self.r_groups.get(r_label, '')
                            dummy_atoms_to_process.append((atom.GetIdx(), 'cap', cap_smiles))
            
            # Apply caps to kept R-groups, remove others
            # Process in two passes: first cap, then remove
            # Cap R-groups: Replace [*:X] with the cap group (e.g., H or OH)
            for atom_idx, action, data in sorted(dummy_atoms_to_process, reverse=True):
                if action == 'cap':
                    cap_smiles = data
                    # For R1 cap '[*:1][H]', we just remove [*:1] (implicit H added)
                    # For R2 cap 'O[*:2]', we need to add O when removing [*:2]
                    # Simplified: check if cap has O
                    if 'O' in cap_smiles and '[*:' in cap_smiles:
                        # R2-like cap: need to add OH group
                        # Get the neighbor atom of the dummy
                        atom = mol_copy.GetAtomWithIdx(atom_idx)
                        neighbors = atom.GetNeighbors()
                        if neighbors:
                            neighbor = neighbors[0]
                            # Add OH to the neighbor before removing dummy
                            emol = Chem.EditableMol(mol_copy)
                            new_o_idx = emol.AddAtom(Chem.Atom(8))  # Oxygen
                            emol.AddBond(neighbor.GetIdx(), new_o_idx, Chem.BondType.SINGLE)
                            emol.RemoveAtom(atom_idx)
                            mol_copy = emol.GetMol()
                    else:
                        # R1-like cap: just remove dummy (implicit H)
                        emol = Chem.EditableMol(mol_copy)
                        emol.RemoveAtom(atom_idx)
                        mol_copy = emol.GetMol()
                elif action == 'remove':
                    # Just remove the dummy atom
                    emol = Chem.EditableMol(mol_copy)
                    emol.RemoveAtom(atom_idx)
                    mol_copy = emol.GetMol()
            
            if mol_copy:
                # Sanitize to add implicit hydrogens where needed
                Chem.SanitizeMol(mol_copy)
                # Generate canonical SMILES without any R-group markers
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
        # IMPORTANT: Only load PEPTIDE monomers (amino acids)
        # The library contains RNA, CHEM, etc. with overlapping symbols (A, C, G, T, U)
        polymer_type = monomer_dict.get('polymerType', '')
        if polymer_type != 'PEPTIDE':
            return None
        
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
                
                # Check if it matches the fragment (exact match only)
                if candidate_smiles == fragment_smiles:
                    return monomer
        
        return None
    
    def find_monomer_by_fragment_smiles_no_stereo(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES WITHOUT stereochemistry.
        Used only in recovery for handling poor quality input data.
        
        Args:
            fragment_smiles: Canonical SMILES of the fragment  
            num_connections: Number of connections this fragment has in the graph
        
        Returns:
            MonomerData if match found, None otherwise
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
                
                # Stereochemistry-agnostic comparison for poor quality data
                candidate_no_stereo = remove_stereochemistry_from_smiles(candidate_smiles)
                fragment_no_stereo = remove_stereochemistry_from_smiles(fragment_smiles)
                
                if candidate_no_stereo == fragment_no_stereo:
                    return monomer
        
        return None

    def find_monomer_by_symbol(self, symbol: str):
        return self.symbol_to_monomer.get(symbol)