# Fixes Applied to R-Group Matching System

## Issues Found and Fixed

### Issue 1: Library Loading Wrong Monomers ❌ → ✅

**Problem:**
- Alanine (A) was loading as "Adenine" (RNA base) with only R1
- Cysteine (C) was loading as "Cytosine" (RNA base) with only R1
- The library contains multiple polymer types (PEPTIDE, RNA, CHEM, etc.) with overlapping symbols

**Root Cause:**
- The `_parse_monomer` method wasn't filtering by `polymerType`
- RNA bases (A, C, G, T, U) were overwriting amino acids with the same symbols

**Fix:**
```python
def _parse_monomer(self, monomer_dict: dict):
    # IMPORTANT: Only load PEPTIDE monomers (amino acids)
    # The library contains RNA, CHEM, etc. with overlapping symbols (A, C, G, T, U)
    polymer_type = monomer_dict.get('polymerType', '')
    if polymer_type != 'PEPTIDE':
        return None
```

**Result:**
- ✅ Alanine (A) now loads correctly with R1, R2
- ✅ Cysteine (C) now loads correctly with R1, R2, R3
- ✅ All 322 PEPTIDE monomers loaded (down from 700 mixed types)

---

### Issue 2: R-Groups Not Being Removed ❌ → ✅

**Problem:**
- When removing R-groups, all results were identical
- R1 removed, R2 removed, both removed - all showed the same SMILES
- The dummy atoms weren't actually being removed

**Root Cause:**
- Code was checking `atom.GetIsotope()` but SMILES notation `[*:1]` uses **atom map numbers**, not isotopes
- All dummy atoms had `Isotope: 0`
- The `:1`, `:2`, `:3` are atom map numbers accessed via `GetAtomMapNum()`

**Before (wrong):**
```python
isotope = atom.GetIsotope()  # Always 0!
if isotope > 0:
    r_label = f"R{isotope}"
```

**After (correct):**
```python
map_num = atom.GetAtomMapNum()  # Gets 1, 2, 3, etc.
if map_num > 0:
    r_label = f"R{map_num}"
```

**Result:**
- ✅ R1 removed: `C[C@H](N)C(=O)[*:2]` - Different!
- ✅ R2 removed: `C[C@@H](C=O)N[*:1]` - Different!
- ✅ Both removed: `C[C@H](N)C=O` - Different!

---

## Verification Test Results

### Alanine (A)
```
Symbol: A
Name: Alanine
SMILES: C[C@H](N[*:1])C(=O)[*:2]
R-groups: {'R1': '[*:1][H]', 'R2': 'O[*:2]'}
R-group count: 2

Lazy generation:
- R1 removed: C[C@H](N)C(=O)[*:2]        ← N-terminus open
- R2 removed: C[C@@H](C=O)N[*:1]         ← C-terminus open
- R1,R2 removed: C[C@H](N)C=O            ← Both open (middle)
```

### Cysteine (C)
```
Symbol: C
Name: Cysteine
SMILES: O=C([C@H](CS[*:3])N[*:1])[*:2]
R-groups: {'R1': '[*:1][H]', 'R2': 'O[*:2]', 'R3': '[*:3][H]'}
R-group count: 3
```

### Aspartic acid (D)
```
Symbol: D
Name: Aspartic acid
SMILES: O=C([C@H](CC(=O)[*:3])N[*:1])[*:2]
R-groups: {'R1': '[*:1][H]', 'R2': 'O[*:2]', 'R3': 'O[*:3]'}
R-group count: 3
```
**Note:** R3 is the side chain carboxyl! This is why D is never confused with N.

### Asparagine (N)
```
Symbol: N
Name: Asparagine
SMILES: NC(=O)C[C@H](N[*:1])C(=O)[*:2]
R-groups: {'R1': '[*:1][H]', 'R2': 'O[*:2]'}
R-group count: 2
```
**Note:** Only 2 R-groups (no R3), so side chain `-CH2-C(=O)NH2` stays intact.

---

## Impact on Matching

### Before Fixes
```python
# Fragment from middle of Ala-Xxx-Yyy
# Expected: 2 connections (left and right)

# But library loaded RNA bases:
library['A'] = Adenine (1 R-group)  ❌

# And R-group removal didn't work:
# All attempts returned same SMILES ❌
```

### After Fixes
```python
# Fragment from middle of Ala-Xxx-Yyy
fragment_smiles = "C[C@H](N)C=O"
num_connections = 2

# Library has correct amino acid:
library['A'] = Alanine (2 R-groups)  ✅

# Try combinations where 2 R-groups removed:
combinations(['R1', 'R2'], 2) = [('R1', 'R2')]

# Generate SMILES with R1, R2 removed:
smiles = alanine.get_capped_smiles_for_removed_rgroups({'R1', 'R2'})
# Returns: "C[C@H](N)C=O"  ✅

# Match!
if smiles == fragment_smiles:
    return alanine  ✅
```

---

## Summary of Changes

### Files Modified
- `logics/monomer_library.py`:
  - Added polymer type filter (line 125-129)
  - Fixed R-group removal to use `GetAtomMapNum()` instead of `GetIsotope()` (line 70)

### Lines Changed
- **2 critical fixes** in monomer_library.py
- Total changes: ~10 lines
- Impact: Massive!

### Before
- ❌ Loading RNA/CHEM monomers instead of amino acids
- ❌ R-groups not being removed
- ❌ All lazy-generated SMILES identical
- ❌ Matching would fail completely

### After  
- ✅ Loading only PEPTIDE monomers
- ✅ R-groups correctly removed using atom map numbers
- ✅ Each lazy-generated SMILES is unique
- ✅ Matching will work correctly!

---

## What's Working Now

1. ✅ **Library loads 322 PEPTIDE monomers** (not 700 mixed types)
2. ✅ **Alanine has R1, R2** (not just R1)
3. ✅ **R-group removal works** (different SMILES for different removals)
4. ✅ **Lazy generation works** (generates on-demand, caches results)
5. ✅ **Asparagine vs Aspartic distinction** (2 vs 3 R-groups)
6. ✅ **No linting errors**

---

## Ready for Testing

The system is now ready to test with actual peptide molecules! The revolutionary R-group based matching should work correctly:

1. Fragment molecules
2. Count connections in graph
3. Try all C(M,N) combinations of R-group removals
4. Match fragment SMILES against generated library SMILES
5. Return correct monomer!

🎉 **Both critical bugs fixed!** 🎉

