import pandas as pd
import os
import time
from pipeline import convert_molecules_batch


def test_peptides(filename, test_name, molfile_column, helm_column, extra_column=None):
    """
    Generic function to test peptide conversion from molecule to HELM notation.
    
    Args:
        filename: Name of the CSV file in test-sets folder
        test_name: Display name for the test (e.g., "LINEAR PEPTIDES")
        molfile_column: Name of the column containing molfile data
        helm_column: Name of the column containing reference HELM notation
        extra_column: Optional column name to display additional info (e.g., BILN notation)
    """
    print(f"\n=== TESTING {test_name.upper()} ===")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    table_path = os.path.join(project_root, 'test-sets', filename)

    if not os.path.exists(table_path):
        print("File not found")
        return

    table = pd.read_csv(table_path)
    
    # Extract molfiles and reference HELM notations
    molfiles = table[molfile_column].tolist()
    reference_helms = table[helm_column].tolist()
    
    # Convert all molecules in batch
    print(f"Converting {len(molfiles)} molecules...")
    start_time = time.time()
    results = convert_molecules_batch(molfiles)
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # Process and display results
    passed = 0
    failed_tests = []
    total = len(table)

    for index, ((success, predicted_helm), true_helm) in enumerate(zip(results, reference_helms)):
        if success:
            # Display extra column if specified (e.g., BILN notation)
            if extra_column and extra_column in table.columns:
                extra_data = table[extra_column][index]
                print(f"{index + 1:>4}. {extra_column}:      {extra_data}")
            
            print(f"{index + 1:>4}. Refered:   {true_helm}")
            print(f"      Generated: {predicted_helm}")
            
            if true_helm == predicted_helm:
                passed += 1
            else:
                failed_tests.append(index + 1)
        else:
            # Molecule conversion failed
            if extra_column and extra_column in table.columns:
                extra_data = table[extra_column][index]
                print(f"{index + 1:>4}. {extra_column}:      {extra_data}")
            
            print(f"{index + 1:>4}. Refered:   {true_helm}")
            print(f"      Generated: FAILED - {predicted_helm if predicted_helm else 'Could not parse molecule'}")
            failed_tests.append(index + 1)
    
    print("\n" + "=" * 60)
    print(f"{test_name.upper()} SUMMARY:")
    print(f"Passed: {passed}/{total} ({passed/total*100:.1f}%)")
    print(f"Time consumed: {elapsed_time:.2f} seconds ({elapsed_time/total:.3f} sec/molecule)")
    if failed_tests:
        print(f"Failed tests: {failed_tests}")
    else:
        print("All tests passed! ✓")
    print("=" * 60)

if __name__ == "__main__":
    # Test linear peptides
    test_peptides(
        filename='HELM_LINEAR.csv',
        test_name='Linear Peptides',
        molfile_column='molfile(HELM)',
        helm_column='HELM'
    )
    
    # # Test cyclic peptides
    # test_peptides(
    #     filename='HELM_cyclic.csv',
    #     test_name='Cyclic Peptides',
    #     molfile_column='molfile(sequence)',
    #     helm_column='sequence'
    # )
    
    # # Test BILN peptides
    # test_peptides(
    #     filename='BILN_W_HELM.csv',
    #     test_name='BILN Peptides',
    #     molfile_column='molfile(BILN)',
    #     helm_column='helm(BILN)',
    #     extra_column='BILN'
    # )