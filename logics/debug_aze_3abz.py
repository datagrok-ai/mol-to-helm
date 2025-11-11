#!/usr/bin/env python3
"""Debug why Aze and 3Abz peptide bond is not detected"""

import pandas as pd
from rdkit import Chem
from fragment_processor import BondDetector

# Load test 30 from PROBLEMS.csv
df = pd.read_csv('../test-sets/PROBLEMS.csv')
test_30 = df.iloc[29]  # 0-indexed
molfile = test_30['molfile(sequence)']
mol = Chem.MolFromMolBlock(molfile)

print(f"Test 30 (PROBLEMS.csv) has {mol.GetNumAtoms()} atoms")
print(f"Expected sequence: {test_30['sequence']}")
print()

# Check current peptide bond pattern
detector = BondDetector()
print("Current peptide bond SMARTS:")
print(Chem.MolToSmarts(detector.peptide_bond))
print()

# Find peptide bonds
peptide_bonds = mol.GetSubstructMatches(detector.peptide_bond)
print(f"Detected {len(peptide_bonds)} peptide bonds")
print()

# Check Aze structure from library
print("=== Aze structure ===")
print("SMILES: O=C([C@@H]1CCN1[*:1])[*:2]")
print("Aze is a 4-membered azetidine ring")
print("Carbonyl is sp2 (X3), NOT in a ring (attached to ring)")
print()

# Check 3Abz structure from library  
print("=== 3Abz structure ===")
print("SMILES: O=C(c1cccc(N[*:1])c1)[*:2]")
print("3-Aminobenzoic acid: benzoic acid with NH2 at position 3")
print("The N is attached to position 3 of benzene ring (aromatic carbon)")
print()

# Try to identify the Aze-3Abz bond region in the molecule
# Look for azetidine ring pattern
azetidine_pattern = Chem.MolFromSmarts('C1CCN1')
azetidine_matches = mol.GetSubstructMatches(azetidine_pattern)
print(f"Found {len(azetidine_matches)} azetidine rings at atoms: {azetidine_matches}")

# Look for aminobenzoic acid pattern
aminobenzoic_pattern = Chem.MolFromSmarts('c1cccc(N)c1C(=O)')
aminobenzoic_matches = mol.GetSubstructMatches(aminobenzoic_pattern)
print(f"Found {len(aminobenzoic_matches)} aminobenzoic acid patterns at atoms: {aminobenzoic_matches}")
print()

# If we found both, check the bond between them
if azetidine_matches and aminobenzoic_matches:
    # Find the carbonyl carbon connected to azetidine
    for aze_atoms in azetidine_matches:
        aze_ring_carbons = [a for a in aze_atoms if mol.GetAtomWithIdx(a).GetSymbol() == 'C']
        # Check neighbors for carbonyl
        for c_idx in aze_ring_carbons:
            atom = mol.GetAtomWithIdx(c_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetTotalDegree() == 3:  # sp2 carbon
                    for n2 in neighbor.GetNeighbors():
                        if n2.GetSymbol() == 'O' and n2.GetTotalDegree() == 1:  # carbonyl O
                            carbonyl_c = neighbor.GetIdx()
                            print(f"Found Aze carbonyl carbon at atom {carbonyl_c}")
                            
                            # Check what it's bonded to
                            for bond_neighbor in neighbor.GetNeighbors():
                                if bond_neighbor.GetSymbol() == 'N':
                                    n_idx = bond_neighbor.GetIdx()
                                    print(f"  Carbonyl bonded to nitrogen at atom {n_idx}")
                                    
                                    # Check what N is bonded to
                                    for n_neighbor in bond_neighbor.GetNeighbors():
                                        if n_neighbor.GetIdx() != carbonyl_c:
                                            print(f"    N is bonded to atom {n_neighbor.GetIdx()} ({n_neighbor.GetSymbol()})")
                                            print(f"      Is aromatic? {n_neighbor.GetIsAromatic()}")
                                            print(f"      Hybridization: {n_neighbor.GetHybridization()}")
                                            print(f"      Total degree: {n_neighbor.GetTotalDegree()}")

# Try modified SMARTS patterns
print("\n=== Testing modified SMARTS patterns ===")

# Original pattern (should work)
pattern1 = Chem.MolFromSmarts('[#6]-[C;X3;!r5;!r6](=[O;X1])-[N;X2,X3]~[C;X3,X4]')
matches1 = mol.GetSubstructMatches(pattern1)
print(f"Original pattern: {len(matches1)} matches")

# Without ring constraints on carbonyl
pattern2 = Chem.MolFromSmarts('[#6]-[C;X3](=[O;X1])-[N;X2,X3]~[C;X3,X4]')
matches2 = mol.GetSubstructMatches(pattern2)
print(f"Without ring constraints: {len(matches2)} matches")

# With aromatic constraint on alpha-C
pattern3 = Chem.MolFromSmarts('[#6]-[C;X3;!r5;!r6](=[O;X1])-[N;X2,X3]~[c;X3]')
matches3 = mol.GetSubstructMatches(pattern3)
print(f"With aromatic alpha-C (lowercase c): {len(matches3)} matches")

# Check if the bond matches the simpler amide pattern
simple_amide = Chem.MolFromSmarts('[C](=[O])-[N]')
amide_matches = mol.GetSubstructMatches(simple_amide)
print(f"Simple amide C(=O)-N pattern: {len(amide_matches)} matches")

