import pandas as pd
import os
import time
from pipeline import convert_molecules_batch


def normalize_helm_for_comparison(helm_str):
    """
    Normalize HELM string for comparison by handling known library duplicates.
    
    Known issues in HELMCoreLibrary.json:
    - Bmt and Bmt_E have identical SMILES (both lack proper E/Z stereochemistry)
    - Lys_Ac is chemically equivalent to ac.K (acetyl cap + lysine)
      When Lys_Ac splits into ac.K, it affects cyclic notation (residue count and positions)
    
    Args:
        helm_str: HELM notation string
    
    Returns:
        Normalized HELM string
    """
    if not helm_str:
        return helm_str
    
    # Treat Bmt_E as Bmt since they're identical in the library
    normalized = helm_str.replace('[Bmt_E]', '[Bmt]')
    normalized = normalized.replace('.Bmt_E.', '.Bmt.')
    
    # Lys_Ac (acetylated lysine) = ac (acetyl cap) + K (lysine)
    # When at position 1 in cyclic peptides, this affects both:
    # 1. The sequence: [Lys_Ac] → [ac].K
    # 2. The cyclic connection: N:R2-1:R1 → (N+1):R2-2:R1
    import re
    
    # Pattern 1: [Lys_Ac] at start with cyclic notation
    # e.g., {[Lys_Ac]...}$PEPTIDE1,PEPTIDE1,16:R2-1:R1$ → {[ac].K...}$PEPTIDE1,PEPTIDE1,17:R2-2:R1$
    pattern1 = r'(\{)(\[Lys_Ac\])(\.[^\}]+\})(\$PEPTIDE1,PEPTIDE1,)(\d+)(:R2-1:R1)'
    def replace1(match):
        n = int(match.group(5))
        return f"{match.group(1)}[ac].K{match.group(3)}{match.group(4)}{n+1}:R2-2:R1"
    normalized = re.sub(pattern1, replace1, normalized)
    
    # Pattern 2: General [Lys_Ac] → ac.K (for non-position-1 cases or linear)
    normalized = normalized.replace('[Lys_Ac]', 'ac.K')
    normalized = normalized.replace('.Lys_Ac.', '.ac.K.')
    
    return normalized


def test_peptides(filename, test_name, molfile_column, helm_column, extra_column=None, library_path=None):
    """
    Generic function to test peptide conversion from molecule to HELM notation.
    
    Args:
        filename: Name of the CSV file in test-sets folder
        test_name: Display name for the test (e.g., "LINEAR PEPTIDES")
        molfile_column: Name of the column containing molfile data
        helm_column: Name of the column containing reference HELM notation
        extra_column: Optional column name to display additional info (e.g., BILN notation)
        library_path: Optional path to custom monomer library JSON file
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
    
    # Load custom library if path provided
    library_json = None
    if library_path:
        if not os.path.exists(library_path):
            print(f"ERROR: Library file not found: {library_path}")
            return
        
        print(f"Reading custom library from: {library_path}")
        with open(library_path, 'r', encoding='utf-8') as f:
            library_json = f.read()
    
    # Convert all molecules in batch
    print(f"Converting {len(molfiles)} molecules...")
    start_time = time.time()
    results = convert_molecules_batch(molfiles, library_json=library_json)
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
            
            # Normalize both HELM strings for comparison (handles library duplicates like Bmt/Bmt_E)
            normalized_true = normalize_helm_for_comparison(true_helm)
            normalized_predicted = normalize_helm_for_comparison(predicted_helm)
            
            if normalized_true == normalized_predicted:
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
        print("All tests passed!")
    print("=" * 60)

if __name__ == "__main__":

    # # Test problems
    # test_peptides(
    #     filename='PROBLEMS.csv',
    #     test_name='problems',
    #     molfile_column='molfile(sequence)',
    #     helm_column='sequence'
    # )

    # # Test linear peptides
    # test_peptides(
    #     filename='HELM_LINEAR.csv',
    #     test_name='Linear Peptides',
    #     molfile_column='molfile(HELM)',
    #     helm_column='HELM'
    # )
    
    # Test cyclic peptides
    test_peptides(
        filename='HELM_cyclic.csv',
        test_name='Cyclic Peptides',
        molfile_column='molfile(sequence)',
        helm_column='sequence'
    )
    
    # # Test BILN peptides
    # test_peptides(
    #     filename='BILN_W_HELM_2.csv',
    #     test_name='BILN Peptides',
    #     molfile_column='molfile(helm(BILN))',
    #     helm_column='helm(BILN)',
    #     extra_column='BILN'
    # )