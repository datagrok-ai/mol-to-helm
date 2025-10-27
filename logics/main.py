import pandas as pd
from rdkit import Chem
import os
from pipeline import mol_to_helm


def test_linear_peptides():
    print("=== TESTING LINEAR PEPTIDES ===")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    table_path = os.path.join(project_root, 'test-sets', 'HELM_LINEAR.csv')

    if not os.path.exists(table_path):
        print("File not found")
        return

    table = pd.read_csv(table_path)
    
    passed = 0
    failed_tests = []
    total = len(table)

    for index in range(len(table)):
        molfile = table['molfile(HELM)'][index]
        true_helm = table['HELM'][index]

        mol = Chem.MolFromMolBlock(molfile)
        if mol:
            predicted_helm = mol_to_helm(mol, is_cyclic=False)

            print(f"{index + 1:>4}. Refered:   {true_helm}")
            print(f"      Generated: {predicted_helm}")
            
            if true_helm == predicted_helm:
                passed += 1
            else:
                failed_tests.append(index + 1)
    
    print("\n" + "=" * 60)
    print(f"LINEAR PEPTIDES SUMMARY:")
    print(f"Passed: {passed}/{total} ({passed/total*100:.1f}%)")
    if failed_tests:
        print(f"Failed tests: {failed_tests}")
    else:
        print("All tests passed! ✓")
    print("=" * 60)



def test_biln_peptides():
    print("\n=== TESTING BILN PEPTIDES ===")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    table_path = os.path.join(project_root, 'test-sets', 'BILN_W_HELM.csv')

    if not os.path.exists(table_path):
        print("File not found")
        return

    table = pd.read_csv(table_path)
    
    passed = 0
    failed_tests = []
    total = len(table)

    for index in range(len(table)):
        molfile = table['molfile(BILN)'][index]
        true_helm = table['helm(BILN)'][index]
        biln_notation = table['BILN'][index]

        mol = Chem.MolFromMolBlock(molfile)
        if mol:
            predicted_helm = mol_to_helm(mol, is_cyclic=True)  # BILN typically represents cyclic peptides

            print(f"{index + 1:>4}. BILN:      {biln_notation}")
            print(f"      Refered:   {true_helm}")
            print(f"      Generated: {predicted_helm}")
            
            if true_helm == predicted_helm:
                passed += 1
            else:
                failed_tests.append(index + 1)
    
    print("\n" + "=" * 60)
    print(f"BILN PEPTIDES SUMMARY:")
    print(f"Passed: {passed}/{total} ({passed/total*100:.1f}%)")
    if failed_tests:
        print(f"Failed tests: {failed_tests}")
    else:
        print("All tests passed! ✓")
    print("=" * 60)



def test_cyclic_peptides():
    print("\n=== TESTING CYCLIC PEPTIDES ===")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    table_path = os.path.join(project_root, 'test-sets', 'HELM_cyclic.csv')

    if not os.path.exists(table_path):
        print("File not found")
        return

    table = pd.read_csv(table_path)
    
    passed = 0
    failed_tests = []
    total = len(table)

    for index in range(len(table)):
        molfile = table['molfile(sequence)'][index]
        true_helm = table['sequence'][index]

        mol = Chem.MolFromMolBlock(molfile)
        if mol:
            predicted_helm = mol_to_helm(mol, is_cyclic=True)

            print(f"{index + 1:>4}. Refered:   {true_helm}")
            print(f"      Generated: {predicted_helm}")
            
            if true_helm == predicted_helm:
                passed += 1
            else:
                failed_tests.append(index + 1)
    
    print("\n" + "=" * 60)
    print(f"CYCLIC PEPTIDES SUMMARY:")
    print(f"Passed: {passed}/{total} ({passed/total*100:.1f}%)")
    if failed_tests:
        print(f"Failed tests: {failed_tests}")
    else:
        print("All tests passed! ✓")
    print("=" * 60)




if __name__ == "__main__":
    test_linear_peptides()
    #test_cyclic_peptides()
    #test_biln_peptides()