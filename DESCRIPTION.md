# Monomer Library Documentation

## Overview

The Monomer Library is the core reference database used for converting molecular structures to HELM notation. It contains detailed information about amino acids and other monomers, including their structures, symbols, and chemical properties.

## Library Source

The library is loaded from `libraries/HELMCoreLibrary.json`, which is a JSON file containing the HELM Core Library standard monomers. Only PEPTIDE-type monomers are loaded (322 entries). RNA, DNA, and CHEM monomers are skipped to avoid symbol collisions (e.g., A, C, G, T overlap between amino acids and nucleotides).

---

## Data Structures

### 1. MonomerData Class (`monomer_library.py`)

| Field | Type | Description |
|-------|------|-------------|
| `symbol` | `str` | HELM symbol (e.g., `"A"`, `"dF"`, `"Phe_4Sdihydroorotamido"`) |
| `name` | `str` | Full name (e.g., `"Alanine"`) |
| `mol` | `Chem.Mol` | RDKit molecule object with R-group dummy atoms |
| `smiles` | `str` | Original SMILES from library (with `[*:1]`, `[*:2]` R-groups) |
| `r_groups` | `dict` | R-group label -> cap SMILES (e.g., `{"R1": "[*:1][H]", "R2": "O[*:2]"}`) |
| `r_group_count` | `int` | Number of R-groups |
| `capped_smiles_cache` | `dict` | Cache: `frozenset(removed_rgroups)` -> canonical SMILES |

**Example:**
```python
monomer.symbol = "A"
monomer.name = "Alanine"
monomer.smiles = "C[C@H](N[*:1])C(=O)[*:2]"
monomer.r_groups = {"R1": "[*:1][H]", "R2": "O[*:2]"}
monomer.r_group_count = 2
```

#### Key Method: `get_capped_smiles_for_removed_rgroups(removed_rgroups)`

Generates canonical SMILES with specific R-groups removed (lazy, cached).

- R-groups in `removed_rgroups`: dummy atom removed (bond was broken during fragmentation)
- R-groups NOT in `removed_rgroups`: capped per library spec (e.g., OH for R2, H for R1)

```python
# Monomer with R1, R2:
monomer.get_capped_smiles_for_removed_rgroups(frozenset({"R1", "R2"}))
# -> "CC(N)C=O"  (both R-groups removed, free ends)

monomer.get_capped_smiles_for_removed_rgroups(frozenset({"R1"}))
# -> "CC(N)C(=O)O"  (R1 removed, R2 capped with OH)
```

---

### 2. MonomerLibrary Class (`monomer_library.py`)

#### Lookup Indices

| Attribute | Key Type | Description |
|-----------|----------|-------------|
| `monomers` | symbol | Main dict: symbol -> MonomerData |
| `symbol_to_monomer` | symbol | Same as monomers |
| `name_to_monomer` | normalized name | Lowercase, no spaces/dashes/underscores |

#### Key Methods

| Method | Description |
|--------|-------------|
| `load_from_helm_json(path)` | Load monomers from HELM JSON library file |
| `find_monomer_by_symbol(symbol)` | O(1) lookup by HELM symbol |
| `find_monomer_by_fragment_smiles(smiles, num_connections)` | Match fragment using R-group combinatorics |
| `find_monomer_by_fragment_smiles_no_stereo(smiles, num_connections)` | Stereo-agnostic matching (fallback for poor input) |

---

## Matching Algorithm

### Graph-Aware R-Group Matching

The matching system uses R-group information from the library to eliminate ambiguity.

**Index Building (one-time, at library load):**
1. For each monomer with M R-groups, generate all 2^M - 1 capped SMILES variants
2. Store in hash index: `(canonical_smiles, num_removed)` -> `MonomerData`
3. Deduplication: monomers with identical SMILES+R-groups share one set of computations
4. Also builds stereo-agnostic index with `remove_stereochemistry_from_smiles`

**Matching (per fragment, O(1)):**
1. Get fragment canonical SMILES and connection count
2. Dict lookup: `_smiles_index[(smiles, num_connections)]` -> monomer

**Why this eliminates N/D confusion:**
- Asparagine (N) has 2 R-groups: R1 (N-term), R2 (C-term)
- Aspartic acid (D) has 3 R-groups: R1, R2, R3 (side chain carboxyl)
- With 2 connections removed: N gives `NC(=O)C[C@H](N)C=O`, D still has `[*:3]` -> different SMILES

### Stereo-Agnostic Matching

For input molecules with missing stereochemistry:
- Strips `@`/`@@` markers from SMILES
- Only strips brackets from organic subset atoms (B,C,N,O,P,S,F,Cl,Br,I)
- Preserves brackets for Se, Te, etc. (required for valid SMILES)
- Uses molecular graph isomorphism as fallback when canonical SMILES differ

---

## Recovery System

Three-phase recovery for unmatched fragments:

1. **Exact merge recovery**: Merge unmatched fragment with neighbor, try exact match
2. **Stereo-agnostic individual recovery**: Match individual unmatched fragments without stereo
3. **Stereo-agnostic merge recovery**: Merge pairs of BOTH-unmatched neighbors, try stereo-agnostic match. Handles split monomers like `Phe_4Sdihydroorotamido` whose internal amide bonds get incorrectly cleaved.

---

## Fragment Graph System (`fragment_graph.py`)

### LinkageType Enum
- `PEPTIDE` - C(=O)-N peptide bond
- `DISULFIDE` - S-S disulfide bridge
- `ESTER`, `ETHER`, `THIOETHER`, `UNKNOWN`

### FragmentNode
- `id`, `mol`, `smiles`, `monomer`, `is_c_terminal`, `is_n_terminal`

### FragmentLink
- `from_node_id`, `to_node_id`, `linkage_type`, `from_atom_idx`, `to_atom_idx`

### FragmentGraph
- `nodes`, `links`
- `get_ordered_nodes()` - sequential traversal from N-terminus
- `is_cyclic()` - detects head-to-tail peptide bonds
- `find_all_cycles()` - DFS cycle detection for multi-chain structures
- `get_connected_components()` - for multi-chain separation

---

## HELM Generation (`helm_generator.py`)

### Simple HELM (linear + single-cycle)
- Backbone traversal from N-terminus
- Branch detection (nodes without R1, like `ac`)
- Branch R-group detection from library (R1 or R2)
- Cyclic connection: `PEPTIDE1,PEPTIDE1,N:R2-1:R1`
- Disulfide: `PEPTIDE1,PEPTIDE1,X:R3-Y:R3`

### Multi-Chain HELM (multiple cycles)
- Each cycle becomes a separate PEPTIDE chain
- Cross-chain connections via R3

---

## HELM Comparison (`main.py`)

The `normalize_helm_for_comparison()` function handles notation equivalences:

- **Connection ordering**: connections can appear in any order
- **Connection direction**: canonical endpoint ordering
- **Bracket normalization**: `ac` -> `[ac]` for multi-char monomers
- **Library duplicates**: `Bmt_E` -> `Bmt`, `Lys_Ac` -> `[ac].K`
- **Unspecified stereo**: `xiThr`/`aThr` -> `T`, `xiIle`/`aIle` -> `I`

---

## JSON Library Entry Format

```json
{
  "symbol": "A",
  "name": "Alanine",
  "smiles": "C[C@H](N[*:1])C(=O)[*:2]",
  "polymerType": "PEPTIDE",
  "monomerType": "Backbone",
  "naturalAnalog": "A",
  "rgroups": [
    {"label": "R1", "capGroupName": "H", "capGroupSMILES": "[*:1][H]"},
    {"label": "R2", "capGroupName": "OH", "capGroupSMILES": "O[*:2]"}
  ]
}
```

---

## Bond Detection SMARTS (`fragment_processor.py`)

```
Peptide bond: [#6]-[C;X3;!r5;!r6](=[O;X1])-[N;X2,X3]~[#6;X3,X4]
  - First carbon: any (#6), including aromatic (for NMe2Abz)
  - Carbonyl: sp2 (X3), NOT in 5/6-membered ring (preserves lactams like Pyr)
  - Nitrogen: X2 (proline) or X3 (standard), any bond type (~) to next carbon
  - Alpha-C: sp3 (X4) or sp2 (X3, dehydro amino acids)

Disulfide bond: [C;X4]-[S;X2]-[S;X2]-[C;X4]

Primary amine (N-terminus): [N;H2,H3;X3,X4]-[C;X3,X4]
```
