# MolToHelm
Converts peptide molecule structures (RDKit molecules) to HELM (Hierarchical Editing Language for Macromolecules) notation.

## Capabilities
- Linear peptides (100% accuracy on test set)
- Cyclic peptides with head-to-tail connections (100% accuracy on test set)
- Disulfide bridges (S-S bonds between Cys residues)
- Stapled peptides (RCMtrans/RCMcis olefin metathesis, FC01 thioether linkers)
- Branch/side-chain modifications (acetyl caps, etc.)
- Multi-chain structures (BILN peptides - partial support)
- Modified and non-standard amino acids (1,919 PEPTIDE monomers from HELMCore+Other Library)
- Accepts both molfile and SMILES input (auto-detection)
- Custom monomer library support (pass JSON at runtime)

## Setup

### Prerequisites
- Python 3.8+
- pip

### Installation

**Windows:**
```bash
py -m pip install -r requirements.txt
```

**Linux/Mac:**
```bash
pip install -r requirements.txt
```

Or use the setup scripts: `setup.bat` (Windows) / `./setup.sh` (Linux/Mac).

Verify with:
```bash
py verify_setup.py   # Windows
python verify_setup.py  # Linux/Mac
```

## Project Structure
```
mol-to-helm/
├── libraries/
│   └── HELMCoreLibrary.json       # HELM monomer library (322 peptide monomers)
├── logics/
│   ├── main.py                     # Test runner with HELM comparison
│   ├── pipeline.py                 # Main conversion pipeline (entry point)
│   ├── monomer_library.py          # Monomer library + R-group matching
│   ├── monomer_matcher.py          # Fragment-to-monomer matching
│   ├── fragment_processor.py       # Bond detection, fragmentation, recovery
│   ├── fragment_graph.py           # Graph data structures (nodes, links, cycles)
│   └── helm_generator.py           # HELM notation generation
├── aggregated/
│   ├── aggregate_logics.py         # Aggregation script for Datagrok
│   └── mol_to_helm_aggregated.py   # Single-file output for Datagrok
├── test-sets/
│   ├── HELM_LINEAR.csv             # Linear peptide tests (20 molecules)
│   ├── HELM_cyclic.csv             # Cyclic peptide tests (41 molecules)
│   ├── BILN_W_HELM.csv             # BILN peptide tests
│   └── BILN_W_HELM_2.csv           # BILN peptide tests (5 molecules)
├── requirements.txt
├── DESCRIPTION.md                  # Detailed module documentation
├── PROJECT_CONTEXT.md              # Architecture and design decisions
└── README.md
```

## Usage

### Running Tests
```bash
cd logics
py main.py        # Windows
python main.py    # Linux/Mac
```

### API
```python
from pipeline import convert_molecules_batch

# From molfiles
results = convert_molecules_batch([molfile_string], input_type="molfile")

# From SMILES
results = convert_molecules_batch(["CC(N)C(=O)NCC(=O)O"], input_type="smiles")

# Auto-detect format
results = convert_molecules_batch([mol_string])

# With custom library
with open("my_library.json") as f:
    results = convert_molecules_batch([mol_string], library_json=f.read())

# Each result is (success: bool, helm_notation: str)
for success, helm in results:
    if success:
        print(helm)
```

### Pipeline Flow
```
Input (molfile/SMILES)
  -> Parse molecule (RDKit)
  -> Detect cleavable bonds (peptide + disulfide + staple sidechain)
  -> Fragment molecule at bonds
  -> Build FragmentGraph (nodes + links)
  -> Match each fragment to library monomer (R-group combinatorics)
  -> Recovery: merge unmatched fragments with neighbors
  -> Recovery: stereo-agnostic matching for poor-quality input
  -> Recovery: merge pairs of unmatched fragments (split monomers)
  -> Generate HELM notation from graph
Output: HELM string
```

## Test Results
| Test Set | Result |
|---|---|
| Linear Peptides | 20/20 (100%) |
| Cyclic Peptides | 41/41 (100%) |
| BILN Peptides | 1/5 (20%) |
| Stapled Peptides | 3678/3756 monomers matched (97.9%) |

## Performance

Matching uses a pre-built SMILES hash index for O(1) lookup per fragment. Per-molecule time is constant regardless of library size.

| Library | Parse | Index Build | Per Molecule |
|---|---|---|---|
| 322 (HELMCore) | 20ms | 100ms | 8ms |
| 1,919 (HELMCore+Other) | 184ms | 2.0s | 51ms |

Stereo-agnostic index uses RDKit re-canonicalization (not string manipulation) to ensure consistent SMILES regardless of how molecules were constructed. Index building deduplicates monomers with identical SMILES+R-groups.

## Requirements
- rdkit >= 2023.3.1
- pandas >= 2.0.0
- numpy >= 1.24.0

## License
See LICENSE file for details.
