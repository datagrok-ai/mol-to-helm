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
    Only strips brackets from SMILES organic subset atoms (B,C,N,O,P,S,F,Cl,Br,I).
    Atoms like Se, Te, etc. must keep their brackets to remain valid SMILES.
    """
    if not smiles:
        return smiles

    # SMILES organic subset: atoms that can appear without brackets
    organic_subset = {'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}

    # Remove @ symbols (stereochemistry markers)
    smiles_no_stereo = re.sub(r'(@+)', '', smiles)

    # Remove explicit H and brackets only for organic subset atoms
    # [C@@H] -> [CH] -> C, but [SeH] must stay as [SeH]
    def _simplify_bracket(match):
        atom = match.group(1)  # e.g. 'C', 'Se', 'N'
        has_h = match.group(2)  # 'H' or ''
        if atom in organic_subset:
            return atom  # Strip brackets (and H) for organic subset
        elif has_h:
            return f'[{atom}H]'  # Keep brackets and H for non-organic atoms
        else:
            return f'[{atom}]'  # Keep brackets for non-organic atoms

    smiles_no_stereo = re.sub(r'\[([A-Z][a-z]?)(H?)\]', _simplify_bracket, smiles_no_stereo)

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


def _canonicalize_no_stereo(smiles: str) -> str:
    """
    Remove stereochemistry and re-canonicalize through RDKit.
    This ensures consistent canonical SMILES regardless of how the molecule was constructed.
    String-only stereo removal can produce non-canonical SMILES.
    """
    no_stereo = remove_stereochemistry_from_smiles(smiles)
    mol = Chem.MolFromSmiles(no_stereo)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return no_stereo  # Fallback to string version if parse fails


class MonomerLibrary:
    def __init__(self):
        self.monomers = {}
        self.smiles_to_monomer = {}
        self.name_to_monomer = {}
        self.symbol_to_monomer = {}
        # Hash indices for O(1) matching (built after loading)
        self._smiles_index = {}         # canonical_smiles -> MonomerData
        self._smiles_no_stereo_index = {}  # stereo-free_smiles -> MonomerData

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

        # Build hash indices for O(1) matching
        self._build_smiles_indices()

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

    def _build_smiles_indices(self):
        """
        Pre-compute all possible capped SMILES for every monomer and build
        hash indices for O(1) lookup. Called once after loading all monomers.

        For each monomer with M R-groups, generates capped SMILES for all
        possible R-group removal combinations (up to 2^M - 1 entries, typically 1-7).

        Deduplicates monomers with identical SMILES+R-groups to avoid redundant
        capping computations (important for large libraries with variants).
        """
        self._smiles_index = {}
        self._smiles_no_stereo_index = {}

        # Dedup: group monomers by (smiles, r_group_keys) to avoid recomputing
        # identical capped forms for monomers with the same structure
        seen_structures = {}  # (smiles, r_group_frozenset) -> list of capped entries

        for symbol, monomer in self.monomers.items():
            if monomer.r_group_count == 0:
                continue

            r_group_labels = list(monomer.r_groups.keys())
            struct_key = (monomer.smiles, frozenset(monomer.r_groups.items()))

            if struct_key in seen_structures:
                # Reuse cached capped SMILES from an identical monomer
                for capped_smiles, n_removed in seen_structures[struct_key]:
                    key = (capped_smiles, n_removed)
                    if key not in self._smiles_index:
                        self._smiles_index[key] = monomer
                    ns_canonical = _canonicalize_no_stereo(capped_smiles)
                    if ns_canonical:
                        ns_key = (ns_canonical, n_removed)
                        if ns_key not in self._smiles_no_stereo_index:
                            self._smiles_no_stereo_index[ns_key] = monomer
                continue

            # First time seeing this structure — compute capped SMILES
            cached_entries = []

            for n_removed in range(1, monomer.r_group_count + 1):
                for removed_combo in combinations(r_group_labels, n_removed):
                    removed_set = frozenset(removed_combo)
                    capped_smiles = monomer.get_capped_smiles_for_removed_rgroups(removed_set)

                    if not capped_smiles:
                        continue

                    cached_entries.append((capped_smiles, n_removed))

                    key = (capped_smiles, n_removed)
                    if key not in self._smiles_index:
                        self._smiles_index[key] = monomer

                    ns_canonical = _canonicalize_no_stereo(capped_smiles)
                    if ns_canonical:
                        ns_key = (ns_canonical, n_removed)
                        if ns_key not in self._smiles_no_stereo_index:
                            self._smiles_no_stereo_index[ns_key] = monomer

            seen_structures[struct_key] = cached_entries

    def find_monomer_by_fragment_smiles(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES. O(1) hash lookup.
        """
        return self._smiles_index.get((fragment_smiles, num_connections))

    def find_monomer_by_fragment_smiles_no_stereo(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES WITHOUT stereochemistry.
        Used in recovery for handling poor quality input data. O(1) hash lookup.
        """
        ns_canonical = _canonicalize_no_stereo(fragment_smiles)
        if not ns_canonical:
            return None

        return self._smiles_no_stereo_index.get((ns_canonical, num_connections))

    def find_monomer_by_symbol(self, symbol: str):
        return self.symbol_to_monomer.get(symbol)