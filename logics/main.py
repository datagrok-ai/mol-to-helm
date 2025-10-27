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

    for index in range(len(table)):
        molfile = table['molfile(HELM)'][index]
        true_helm = table['HELM'][index]

        mol = Chem.MolFromMolBlock(molfile)
        if mol:
            predicted_helm = mol_to_helm(mol, is_cyclic=False)

            print(f"{index + 1}. Refered:  {true_helm}")
            print(f"   Generated: {predicted_helm}")



def test_cyclic_peptides():
    print("\n=== TESTING CYCLIC PEPTIDES ===")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    table_path = os.path.join(project_root, 'test-sets', 'HELM_cyclic.csv')

    if not os.path.exists(table_path):
        print("File not found")
        return

    table = pd.read_csv(table_path)

    for index in range(len(table)):
        molfile = table['molfile(sequence)'][index]
        true_helm = table['sequence'][index]

        mol = Chem.MolFromMolBlock(molfile)
        if mol:
            predicted_helm = mol_to_helm(mol, is_cyclic=True)

            print(f"{index + 1}. Refered:  {true_helm}")
            print(f"   Generated: {predicted_helm}")




if __name__ == "__main__":
    test_linear_peptides()
    #test_cyclic_peptides()