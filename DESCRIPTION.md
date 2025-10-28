# Monomer Library Documentation

## Overview

The Monomer Library is the core reference database used for converting molecular structures to HELM notation. It contains detailed information about amino acids and other monomers, including their structures, symbols, and chemical properties.

## Library Source

The library is loaded from `libraries/HELMCoreLibrary.json`, which is a JSON file containing the HELM Core Library standard monomers.

---

## Data Structures

### 1. MonomerData Class

The `MonomerData` class represents a single monomer entry with the following attributes:

| Field | Type | Description |
|-------|------|-------------|
| `symbol` | `str` | The HELM symbol used in notation (e.g., "A" for Alanine, "C" for Cysteine) |
| `name` | `str` | The full name of the monomer (e.g., "Alanine", "Cysteine") |
| `mol` | `Chem.Mol` | RDKit molecule object representing the monomer structure |
| `smiles` | `str` | Original SMILES string from the library |
| `smiles_canonical` | `str` | Canonical SMILES string with R-groups removed (used for matching) |
| `r_groups` | `dict` | Dictionary storing R-group information (currently not populated) |

**Example:**
```python
monomer = MonomerData()
monomer.symbol = "A"
monomer.name = "Alanine"
monomer.smiles = "C[C@H](N[*:1])C(=O)[*:2]"
monomer.smiles_canonical = "C[C@H](N)C=O"  # R-groups removed
```

---

### 2. MonomerLibrary Class

The `MonomerLibrary` class manages the collection of monomers and provides lookup functionality through multiple indices:

#### Attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| `monomers` | `dict[str, MonomerData]` | Main dictionary mapping symbols to MonomerData objects |
| `symbol_to_monomer` | `dict[str, MonomerData]` | Fast lookup by HELM symbol |
| `smiles_to_monomer` | `dict[str, MonomerData]` | Fast lookup by canonical SMILES |
| `name_to_monomer` | `dict[str, MonomerData]` | Fast lookup by normalized name (lowercase, no spaces/dashes) |

#### Methods:

| Method | Description |
|--------|-------------|
| `load_from_helm_json(json_path)` | Load monomers from a HELM JSON library file |
| `find_monomer_by_smiles(smiles)` | Find a monomer by its SMILES string (exact match) |
| `find_monomer_by_symbol(symbol)` | Find a monomer by its HELM symbol |

---

## JSON Library Structure

Each monomer in `HELMCoreLibrary.json` has the following fields:

### Core Fields:

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `symbol` | `string` | HELM notation symbol | `"A"` |
| `name` | `string` | Full monomer name | `"Alanine"` |
| `smiles` | `string` | SMILES representation with R-group attachments | `"C[C@H](N[*:1])C(=O)[*:2]"` |
| `molfile` | `string` | MOL file format structure | (multi-line string) |
| `polymerType` | `string` | Type of polymer | `"PEPTIDE"`, `"RNA"`, `"DNA"`, `"CHEM"` |
| `monomerType` | `string` | Monomer classification | `"Backbone"`, `"Branch"`, `"Undefined"` |
| `naturalAnalog` | `string` | Natural amino acid analog | `"A"` |

### Additional Fields:

| Field | Type | Description |
|-------|------|-------------|
| `id` | `number` | Internal identifier |
| `rgroups` | `array` | R-group attachment point definitions |
| `meta` | `object` | Additional metadata |
| `author` | `string` | Author/source of the entry |
| `createDate` | `string/null` | Creation date |

### R-Groups Structure:

Each entry in the `rgroups` array defines an attachment point:

```json
{
  "alternateId": "R1-H",
  "capGroupName": "H",
  "capGroupSMILES": "[*:1][H]",
  "label": "R1"
}
```

---

## Loading Process

When the library is loaded:

1. **Read JSON**: The `HELMCoreLibrary.json` file is parsed
2. **Parse Monomers**: Each monomer entry is processed:
   - Extract `symbol` and `name`
   - Parse `smiles` or `molfile` into RDKit molecule object
   - Generate canonical SMILES with R-groups removed
3. **Build Indices**: Multiple lookup dictionaries are created:
   - By symbol (e.g., `"A"` → Alanine)
   - By canonical SMILES (for structure matching)
   - By normalized name (lowercase, no spaces/dashes)
4. **Validation**: Only monomers with valid structures are included

---

## Canonical SMILES Generation

The library creates a special canonical SMILES for matching purposes:

1. **Input**: Original SMILES with R-groups (e.g., `"C[C@H](N[*:1])C(=O)[*:2]"`)
2. **Process**: 
   - Remove all dummy atoms (R-groups, represented as atoms with atomic number 0)
   - Generate canonical SMILES
3. **Output**: Clean canonical SMILES (e.g., `"C[C@H](N)C=O"`)

This allows matching peptide fragments (which have free ends) against library entries (which have R-group placeholders).

---

## Lookup Examples

### By Symbol:
```python
library = MonomerLibrary()
library.load_from_helm_json("libraries/HELMCoreLibrary.json")

# Direct symbol lookup
monomer = library.find_monomer_by_symbol("A")
# Returns MonomerData for Alanine
```

### By SMILES:
```python
# Fragment SMILES from peptide
fragment_smiles = "C[C@H](N)C=O"

# Find matching monomer
monomer = library.find_monomer_by_smiles(fragment_smiles)
# Returns MonomerData for Alanine if match found
```

### By Normalized Name:
```python
# Normalized names (lowercase, no spaces/dashes)
monomer = library.name_to_monomer.get("alanine")
monomer = library.name_to_monomer.get("phe4me")  # "Phe_4Me" normalized
```

---

## Statistics

The HELMCoreLibrary.json contains:
- **200+** peptide backbone monomers
- Standard 20 amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
- D-amino acids (prefixed with 'd' or 'D-')
- Modified amino acids (methylated, halogenated, etc.)
- Non-natural amino acids
- Chemical modifications and linkers

---

## Usage in Pipeline

The monomer library is used in the conversion pipeline as follows:

1. **Initialization**: Library is loaded once at startup (singleton pattern)
2. **Fragmentation**: Input molecule is fragmented at peptide bonds
3. **Matching**: Each fragment is compared against library:
   - First: Exact SMILES match
   - Second: Hardcoded mapping fallback
   - Third: Normalization attempts (terminal variations)
4. **HELM Generation**: Matched monomers' symbols are assembled into HELM notation

---

## Performance Considerations

- **Caching**: Library is loaded once and reused across multiple conversions
- **Multiple Indices**: Three separate dictionaries enable O(1) lookups by different keys
- **Preprocessing**: Canonical SMILES are pre-computed during loading, not at query time
- **Validation**: Invalid structures are filtered out during loading to prevent runtime errors

---

## Extending the Library

To add custom monomers:

1. Edit `libraries/HELMCoreLibrary.json`
2. Add a new entry following the standard format:
```json
{
  "symbol": "CustomAA",
  "name": "Custom Amino Acid",
  "smiles": "your_smiles_here",
  "polymerType": "PEPTIDE",
  "monomerType": "Backbone",
  "naturalAnalog": "A",
  "rgroups": [...]
}
```
3. Reload the library in your application

The system will automatically:
- Parse the structure
- Generate canonical SMILES
- Add to all lookup indices

