#!/usr/bin/env python3
"""Debug hArg matching issue in test 26"""

import pandas as pd
from rdkit import Chem
from monomer_library import MonomerLibrary, remove_stereochemistry_from_smiles
from fragment_processor import FragmentProcessor

# Load library
library = MonomerLibrary()
library.load_from_helm_json('../libraries/HELMCoreLibrary.json')

# Check if hArg is in library
if 'hArg' in library.monomers:
    monomer = library.monomers['hArg']
    print("hArg found in library!")
    print(f"  Symbol: {monomer.symbol}")
    print(f"  Name: {monomer.name}")
    print(f"  R-groups: {list(monomer.r_groups.keys())}")
    print(f"  Original SMILES: {monomer.smiles}")
    
    # Get capped SMILES (both R1 and R2 removed for 2 connections)
    removed = frozenset(['R1', 'R2'])
    lib_smiles = monomer.get_capped_smiles_for_removed_rgroups(removed)
    print(f"  Library SMILES (R1,R2 removed): {lib_smiles}")
    
    # Remove stereochemistry
    lib_smiles_no_stereo = remove_stereochemistry_from_smiles(lib_smiles)
    print(f"  Library SMILES (no stereo): {lib_smiles_no_stereo}")
else:
    print("ERROR: hArg not found in library!")

# Load test 26 and process it
df = pd.read_csv('../test-sets/HELM_CYCLIC.csv')
test_26 = df.iloc[25]  
molfile = test_26['molfile(sequence)']
mol = Chem.MolFromMolBlock(molfile)

print("\nProcessing test 26 molecule...")
processor = FragmentProcessor(library)
graph = processor.process_molecule(mol)

print(f"Total fragments: {len(graph.nodes)}")

# Find fragments that don't match (likely the hArg positions)
print("\nLooking for unmatched fragments that might be hArg:")
for node_id, node in graph.nodes.items():
    if not node.monomer or node.monomer.symbol.startswith('X'):
        neighbors = graph.get_neighbors(node_id)
        num_connections = len(neighbors)
        frag_smiles = Chem.MolToSmiles(node.mol, canonical=True)
        frag_smiles_no_stereo = remove_stereochemistry_from_smiles(frag_smiles)
        
        print(f"\nNode {node_id} (unmatched):")
        print(f"  Fragment SMILES: {frag_smiles}")
        print(f"  Fragment no stereo: {frag_smiles_no_stereo}")
        print(f"  Connections: {num_connections}")
        print(f"  Atoms: {node.mol.GetNumAtoms()}")
        
        # Check if it matches hArg when compared without stereo
        if 'hArg' in library.monomers:
            hArg_monomer = library.monomers['hArg']
            removed = frozenset(['R1', 'R2'])
            hArg_lib_smiles = hArg_monomer.get_capped_smiles_for_removed_rgroups(removed)
            hArg_lib_no_stereo = remove_stereochemistry_from_smiles(hArg_lib_smiles)
            
            if num_connections == 2:
                print(f"  hArg lib (no stereo): {hArg_lib_no_stereo}")
                print(f"  Match? {frag_smiles_no_stereo == hArg_lib_no_stereo}")

