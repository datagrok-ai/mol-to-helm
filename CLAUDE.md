# mol-to-helm

Converts peptide molecules (molfile/SMILES) to HELM notation using the HELM Core Library.

## Mandatory Workflow (CRITICAL)

After **every code iteration**, you MUST complete these steps in order:

1. **Run tests** - Ensure no regressions:
   ```bash
   cd logics && py -c "
   from main import test_peptides
   test_peptides('HELM_LINEAR.csv', 'Linear', 'molfile(HELM)', 'HELM')
   test_peptides('HELM_cyclic.csv', 'Cyclic', 'molfile(sequence)', 'sequence')
   test_peptides('BILN_W_HELM_2.csv', 'BILN', 'molfile(helm(BILN))', 'helm(BILN)', extra_column='BILN')
   "
   ```
   Expected: Linear 20/20, Cyclic 41/41, BILN 1/5.

2. **Update documentation** if behavior changed:
   - `README.md` - capabilities, API, test results
   - `DESCRIPTION.md` - technical details, data structures, algorithms
   - `PROJECT_CONTEXT.md` - architecture, fixes history, design decisions
   - This `CLAUDE.md` file - if workflow or rules change

3. **Run aggregation** if any file in `logics/` changed:
   ```bash
   cd aggregated && py aggregate_logics.py
   ```
   Then verify syntax: `py -c "compile(open('mol_to_helm_aggregated.py').read(), 'f', 'exec'); print('OK')"`

## Commands

- Use `py` not `python` (Windows)
- Tests: `cd logics && py main.py`
- Aggregate: `cd aggregated && py aggregate_logics.py`

## Project Structure

```
logics/              # Core source code
  pipeline.py        # Entry point: convert_molecules_batch()
  monomer_library.py # Library loading + R-group matching
  fragment_processor.py # Bond detection, fragmentation, recovery
  fragment_graph.py  # Graph data structures
  monomer_matcher.py # Fragment-to-monomer matching
  helm_generator.py  # HELM notation generation
  main.py            # Test runner with HELM comparison
libraries/           # HELMCoreLibrary.json (322 peptide monomers)
test-sets/           # CSV test data (linear, cyclic, BILN)
aggregated/          # Single-file output for Datagrok
```

## Architecture

Pipeline: parse mol -> detect bonds (SMARTS) -> fragment -> match fragments to library (R-group combinatorics) -> recover unmatched (3-phase) -> generate HELM.

Recovery phases: exact merge -> stereo-agnostic individual -> stereo-agnostic merge (both-unmatched pairs only).

## Key Conventions

- Monomer matching uses lazy-cached R-group removal combinations
- `remove_stereochemistry_from_smiles()` preserves brackets for non-organic atoms (Se, Te)
- HELM comparison normalizes: connection order, brackets, Bmt/Bmt_E, Lys_Ac/ac+K, xi-stereo
- Single-cycle + branch peptides route through `_generate_simple_helm` (not multi-chain)
