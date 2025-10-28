# MolToHelm
Functionality for converting peptide molecules to HELM representation

## Description
This project converts peptide molecule structures (RDKit molecules) to HELM (Hierarchical Editing Language for Macromolecules) notation. It supports both linear and cyclic peptides.

## Setup Environment

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Installation

#### Quick Setup Verification
Before installation, you can verify your environment:

**Windows:**
```bash
py verify_setup.py
```

**Linux/Mac:**
```bash
python verify_setup.py
```

#### Option 1: Using setup scripts (Easiest)

**Windows:**
```bash
setup.bat
```

**Linux/Mac:**
```bash
chmod +x setup.sh
./setup.sh
```

#### Option 2: Manual installation with pip
**Windows:**
```bash
py -m pip install -r requirements.txt
```

**Linux/Mac:**
```bash
pip install -r requirements.txt
```

#### Option 3: Using conda
```bash
conda create -n mol-to-helm python=3.9
conda activate mol-to-helm
conda install -c conda-forge rdkit pandas numpy
```

## Project Structure
```
mol-to-helm/
├── libraries/
│   └── HELMCoreLibrary.json    # HELM monomer library
├── logics/
│   ├── main.py                 # Main test runner
│   ├── pipeline.py             # Main conversion pipeline
│   ├── monomer_library.py      # Monomer library management
│   └── target_processor.py     # Molecular fragment processing
├── test-sets/
│   ├── HELM_LINEAR.csv         # Linear peptide test data
│   ├── HELM_cyclic.csv         # Cyclic peptide test data
│   └── BILN_W_HELM.csv         # Additional test data
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

## Usage

### Running Tests
To test the conversion of linear and cyclic peptides:

**Windows:**
```bash
cd logics
py main.py
```

**Linux/Mac:**
```bash
cd logics
python main.py
```

### Using in Your Code
```python
from rdkit import Chem
from pipeline import mol_to_helm

# Load your molecule
mol = Chem.MolFromSmiles("your_smiles_here")

# Convert to HELM notation
helm_notation = mol_to_helm(mol, is_cyclic=False)
print(helm_notation)
```

## Features
- Convert RDKit molecules to HELM notation
- Support for linear peptides
- Support for cyclic peptides
- Monomer library based on HELM Core Library
- Fragment detection and matching

## Requirements
- rdkit>=2023.3.1
- pandas>=2.0.0
- numpy>=1.24.0

## License
See LICENSE file for details.
