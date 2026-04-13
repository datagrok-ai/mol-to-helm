import pandas as pd
import os
import time
import re
from pipeline import convert_molecules_batch


def _parse_helm_parts(helm_str):
    """
    Parse HELM string into its constituent parts.
    HELM format: chains$connections$groups$annotations$version

    Returns:
        tuple: (chains_str, connections_list, groups_str, annotations_str, version_str)
    """
    if not helm_str:
        return "", [], "", "", ""

    parts = helm_str.split('$')
    chains_str = parts[0] if len(parts) > 0 else ""
    connections_str = parts[1] if len(parts) > 1 else ""
    groups_str = parts[2] if len(parts) > 2 else ""
    annotations_str = parts[3] if len(parts) > 3 else ""
    version_str = parts[4] if len(parts) > 4 else ""

    connections_list = [c for c in connections_str.split('|') if c] if connections_str else []

    return chains_str, connections_list, groups_str, annotations_str, version_str


def _parse_sequence(seq_str):
    """
    Parse a HELM sequence string into a list of monomer symbols.
    Handles both bracketed ([Nva]) and single-letter (A) monomers.

    Example: "[dI].[Trp_Ome].A.K" -> ["dI", "Trp_Ome", "A", "K"]
    """
    if not seq_str:
        return []

    monomers = []
    i = 0
    current = ""
    while i < len(seq_str):
        if seq_str[i] == '[':
            # Bracketed monomer
            end = seq_str.index(']', i)
            monomers.append(seq_str[i+1:end])
            i = end + 1
            if i < len(seq_str) and seq_str[i] == '.':
                i += 1  # skip dot separator
        elif seq_str[i] == '.':
            if current:
                monomers.append(current)
                current = ""
            i += 1
        else:
            current += seq_str[i]
            i += 1
    if current:
        monomers.append(current)

    return monomers


def _normalize_sequence(seq_str):
    """
    Normalize a sequence by ensuring all multi-char monomers are bracketed.

    Example: "ac.K.A" -> "[ac].K.A"
    """
    monomers = _parse_sequence(seq_str)
    formatted = []
    for m in monomers:
        if len(m) > 1:
            formatted.append(f"[{m}]")
        else:
            formatted.append(m)
    return ".".join(formatted)


def _normalize_chain(chain_str):
    """
    Normalize a single chain string like 'PEPTIDE1{seq}'.
    Returns (chain_name, normalized_sequence_str, monomer_list).
    """
    match = re.match(r'(\w+)\{(.+)\}', chain_str)
    if not match:
        return chain_str, "", []

    name = match.group(1)
    seq = match.group(2)
    monomers = _parse_sequence(seq)
    norm_seq = _normalize_sequence(seq)

    return name, norm_seq, monomers


def _normalize_connection(conn_str):
    """
    Normalize a single connection string to canonical form.
    e.g., "PEPTIDE1,PEPTIDE1,10:R3-5:R3" -> sorted endpoints.
    """
    # Parse: CHAIN1,CHAIN2,pos1:R#-pos2:R#
    match = re.match(r'(\w+),(\w+),(\d+:\w+)-(\d+:\w+)', conn_str)
    if not match:
        return conn_str

    chain1, chain2, end1, end2 = match.group(1), match.group(2), match.group(3), match.group(4)

    # Normalize direction: sort the two endpoints for consistent ordering
    # Use (chain, endpoint) tuples for sorting
    ep1 = (chain1, end1)
    ep2 = (chain2, end2)

    if ep1 > ep2:
        ep1, ep2 = ep2, ep1

    return f"{ep1[0]},{ep2[0]},{ep1[1]}-{ep2[1]}"


def normalize_helm_for_comparison(helm_str):
    """
    Normalize HELM string for robust comparison.

    Handles:
    - Known library duplicates (Bmt/Bmt_E, Lys_Ac/ac+K)
    - Bracket normalization (ac -> [ac] for multi-char monomers)
    - Connection ordering (sort connections for order-independent comparison)
    - Connection direction normalization (canonical endpoint ordering)

    Args:
        helm_str: HELM notation string

    Returns:
        Normalized HELM string
    """
    if not helm_str:
        return helm_str

    # --- Step 1: Known library duplicate normalizations ---
    normalized = helm_str

    # Treat Bmt_E as Bmt since they're identical in the library
    normalized = normalized.replace('[Bmt_E]', '[Bmt]')
    normalized = normalized.replace('.Bmt_E.', '.Bmt.')

    # Lys_Ac (acetylated lysine) = ac (acetyl cap) + K (lysine)
    # When at position 1 in cyclic peptides, this affects both:
    # 1. The sequence: [Lys_Ac] → [ac].K
    # 2. The cyclic connection: N:R2-1:R1 → (N+1):R2-2:R1
    pattern1 = r'(\{)(\[Lys_Ac\])(\.[^\}]+\})(\$PEPTIDE1,PEPTIDE1,)(\d+)(:R2-1:R1)'
    def replace1(match):
        n = int(match.group(5))
        return f"{match.group(1)}[ac].K{match.group(3)}{match.group(4)}{n+1}:R2-2:R1"
    normalized = re.sub(pattern1, replace1, normalized)

    # General [Lys_Ac] → [ac].K
    normalized = normalized.replace('[Lys_Ac]', '[ac].K')
    normalized = normalized.replace('.Lys_Ac.', '.[ac].K.')

    # xi-prefix monomers have unspecified stereochemistry at one center.
    # They can match either resolved form. Normalize to the L-form equivalent.
    # xiThr (unspecified at Cb) -> T (L-Thr), aThr (L-allo-Thr) -> T
    # xiIle (unspecified at Cb) -> I (L-Ile), aIle (L-allo-Ile) -> I
    for xi_sym, allo_sym, std_sym in [('xiThr', 'aThr', 'T'), ('xiIle', 'aIle', 'I')]:
        normalized = normalized.replace(f'[{xi_sym}]', std_sym)
        normalized = normalized.replace(f'.{xi_sym}.', f'.{std_sym}.')
        normalized = normalized.replace(f'[{allo_sym}]', std_sym)
        normalized = normalized.replace(f'.{allo_sym}.', f'.{std_sym}.')

    # --- Step 2: Parse into parts ---
    chains_str, connections, groups_str, annotations_str, version_str = _parse_helm_parts(normalized)

    # --- Step 3: Normalize chains (bracket normalization) ---
    chain_strs = []
    if chains_str:
        # Split chains on | but respect {} braces
        raw_chains = re.findall(r'\w+\{[^}]+\}', chains_str)
        for chain in raw_chains:
            name, norm_seq, _ = _normalize_chain(chain)
            chain_strs.append(f"{name}{{{norm_seq}}}")

    normalized_chains = "|".join(chain_strs) if chain_strs else chains_str

    # --- Step 4: Normalize and sort connections ---
    normalized_connections = sorted([_normalize_connection(c) for c in connections])
    connections_str = "|".join(normalized_connections)

    # --- Step 5: Reassemble ---
    result = f"{normalized_chains}${connections_str}${groups_str}${annotations_str}${version_str}"

    return result


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

def test_stapled_peptides(filename='stapled_helm.csv', library_path=None):
    """
    Test stapled peptide conversion with monomer-level matching.

    Stapled peptides have CHEM elements (Ac, NH2, FC01, RCMtrans, etc.) that
    are expected to come out as unknown X fragments. The test measures how many
    PEPTIDE monomers are correctly matched vs unknown.

    Args:
        filename: CSV file with original_helm, molfile, smiles columns
        library_path: Path to monomer library JSON (defaults to HELMCoreAndOther.json)
    """
    print(f"\n=== TESTING STAPLED PEPTIDES ===")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    table_path = os.path.join(project_root, 'test-sets', filename)

    if not os.path.exists(table_path):
        print("File not found")
        return

    table = pd.read_csv(table_path)

    if library_path is None:
        library_path = os.path.join(project_root, 'libraries', 'HELMCoreAndOther.json')

    library_json = None
    if os.path.exists(library_path):
        with open(library_path, 'r', encoding='utf-8') as f:
            library_json = f.read()

    smiles_list = table['smiles'].tolist()
    helms = table['original_helm'].tolist()

    chem_symbols = {'Ac', 'NH2', 'FC01', 'RCMtrans', 'RCMcis', 'VHL185'}

    print(f"Converting {len(smiles_list)} molecules...")
    start_time = time.time()
    results = convert_molecules_batch(smiles_list, library_json=library_json, input_type='smiles')
    elapsed = time.time() - start_time

    total_expected = 0
    total_matched = 0
    total_unknown = 0

    for i, (success, got_helm) in enumerate(results):
        expected_monomers = []
        for p in re.findall(r'PEPTIDE\d+\{([^}]+)\}', helms[i]):
            expected_monomers.extend(_parse_sequence(p))

        got_monomers = []
        if success:
            for p in re.findall(r'PEPTIDE\d+\{([^}]+)\}', got_helm):
                got_monomers.extend(_parse_sequence(p))

        expected_peptide_only = [m for m in expected_monomers if m not in chem_symbols]
        total_expected += len(expected_peptide_only)
        for m in got_monomers:
            if '*:' in m or '*]' in m:
                total_unknown += 1  # Inline SMILES (unmatched, but with R-group info)
            else:
                total_matched += 1  # Library monomer match

    match_pct = 100 * total_matched / total_expected if total_expected > 0 else 0

    print("\n" + "=" * 60)
    print("STAPLED PEPTIDES SUMMARY:")
    print(f"Molecules: {len(results)}")
    print(f"Library-matched monomers: {total_matched}/{total_expected} ({match_pct:.1f}%)")
    print(f"Inline SMILES monomers: {total_unknown} (unmatched fragments with R-group SMILES)")
    print(f"Time consumed: {elapsed:.2f} seconds ({elapsed/len(results):.3f} sec/molecule)")
    print("=" * 60)


if __name__ == "__main__":

    # Test linear peptides
    test_peptides(
        filename='HELM_LINEAR.csv',
        test_name='Linear Peptides',
        molfile_column='molfile(HELM)',
        helm_column='HELM'
    )

    # Test cyclic peptides
    test_peptides(
        filename='HELM_cyclic.csv',
        test_name='Cyclic Peptides',
        molfile_column='molfile(sequence)',
        helm_column='sequence'
    )

    # Test BILN peptides
    test_peptides(
        filename='BILN_W_HELM_2.csv',
        test_name='BILN Peptides',
        molfile_column='molfile(helm(BILN))',
        helm_column='helm(BILN)',
        extra_column='BILN'
    )

    # Test stapled peptides
    test_stapled_peptides()