from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
from collections import defaultdict
import json
import os

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

class MonomerData:
    def __init__(self):
        self.symbol = ""
        self.name = ""
        self.mol = None
        self.smiles = ""
        self.r_groups = {}
        self.smiles_canonical = ""

    def __repr__(self):
        return f"Monomer({self.symbol}: {self.name})"

class MonomerLibrary:
    def __init__(self):
        self.monomers = {}
        self.smiles_to_monomer = {}
        self.name_to_monomer = {}
        self.symbol_to_monomer = {}
        self.fingerprint_cache = {}

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

                    if monomer.smiles_canonical:
                        self.smiles_to_monomer[monomer.smiles_canonical] = monomer
                        self.fingerprint_cache[monomer.symbol] = AllChem.GetMorganFingerprintAsBitVect(
                            monomer.mol, 2, nBits=1024
                        )

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
                monomer.smiles_canonical = self._create_canonical_smiles(monomer.mol)
            except Exception:
                monomer.mol = None

        if monomer.mol is None and molfile:
            try:
                monomer.mol = Chem.MolFromMolBlock(molfile)
                if monomer.mol:
                    monomer.smiles = Chem.MolToSmiles(monomer.mol)
                    monomer.smiles_canonical = self._create_canonical_smiles(monomer.mol)
            except Exception:
                monomer.mol = None

        if monomer.mol is None:
            return None

        return monomer

    def _create_canonical_smiles(self, mol: Chem.Mol) -> str:
        try:
            mol_copy = Chem.Mol(mol)
            atoms_to_remove = []
            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atoms_to_remove.append(atom.GetIdx())

            atoms_to_remove.sort(reverse=True)
            emol = Chem.EditableMol(mol_copy)
            for atom_idx in atoms_to_remove:
                emol.RemoveAtom(atom_idx)

            clean_mol = emol.GetMol()
            return Chem.MolToSmiles(clean_mol, canonical=True)
        except Exception:
            return ""

    def find_monomer_by_smiles(self, smiles: str):
        try:
            query_mol = Chem.MolFromSmiles(smiles)
            if not query_mol:
                return None

            query_smiles = Chem.MolToSmiles(query_mol, canonical=True)
            return self.smiles_to_monomer.get(query_smiles)
        except Exception:
            return None

    def find_monomer_by_symbol(self, symbol: str):
        return self.symbol_to_monomer.get(symbol)

    def get_fingerprint(self, symbol: str):
        return self.fingerprint_cache.get(symbol)