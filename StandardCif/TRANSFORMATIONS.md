# mmCIF Standardisation — Transformation Reference

This document describes every data transformation applied by the
`standard-cif` pipeline when converting a depositor-supplied mmCIF file
into the OneDep upload-compatible format (`D_XYZ_model-upload-convert_P1.cif`).

---

## Overview

The pipeline reads a single mmCIF data block, runs six transformers in a
fixed order, and writes the result to a new file.  Each transformer is
responsible for exactly one mmCIF category (Single Responsibility
Principle).  Transformers are stateless except where they must pass a
derived mapping to the next step.

```
input.cif
    │
    ├─ 1. EntityBuilder          → _entity, _struct_asym.entity_id
    ├─ 2. ChainRemapper          → _struct_asym.id, _atom_site.label_asym_id
    ├─ 3. EntityPolySeqBuilder   → _entity_poly_seq, _pdbx_prerelease_seq
    ├─ 4. AtomSiteTransformer    → _atom_site (multiple columns)
    ├─ 5. ChemCompEnricher       → _chem_comp
    └─ 6. AtomTypeBuilder        → _atom_type
    │
output.cif
```

The order is significant: EntityBuilder must run before ChainRemapper
(it needs original chain IDs), and both must finish before
AtomSiteTransformer (which needs the final chain→entity mapping).

---

## 1. EntityBuilder — `transformers/entity.py`

**Categories written:** `_entity`, `_struct_asym.entity_id`

### What it does

Scans `_atom_site` to discover the chemical composition of every chain
and assigns entity IDs following OneDep conventions:

| Priority | Chain type | Detection rule | `_entity.type` | `src_method` |
|---|---|---|---|---|
| 1 (lowest ID) | Polymer | Contains at least one standard amino acid | `polymer` | `man` |
| 2 | Non-polymer | All residues are non-standard small molecules | `non-polymer` | `syn` |
| 3 (highest ID) | Water | All residues are HOH / WAT | `water` | `nat` |

**Grouping rule:** multiple chains that contain exactly the same set of
residue codes share one entity ID (e.g. nine DMS molecules all map to a
single DMS entity).  Each polymer chain always gets its own entity.

**`pdbx_description`** is carried over from the original `_entity` table
by matching on `type`, so the polymer chain description (e.g.
"RNA-directed RNA polymerase NS5") is preserved even when entity IDs
are renumbered.

### Orphan chains

Some input files produced by Gemmi have atoms whose `label_asym_id` is
`.` (a CIF null) but whose `auth_asym_id` is a valid letter (e.g. `W`
for waters, `B` for a stray ligand).  These atoms are not registered in
`_struct_asym`.

The EntityBuilder handles this by:
1. Collecting all atoms with `label_asym_id='.'` and grouping them into
   virtual chains keyed as `orphan:<auth_asym_id>` (e.g. `orphan:W`).
2. Assigning entity IDs to those virtual chains using the same
   composition rules above.
3. Appending new rows to `_struct_asym` for each virtual chain, sorted
   by entity ID so that non-polymer orphans appear before water orphans.
4. Updating the `label_asym_id` of the affected `_atom_site` rows to
   the virtual key, so that the ChainRemapper can later assign them
   proper single-letter IDs.

---

## 2. ChainRemapper — `transformers/chain_remapper.py`

**Categories written:** `_struct_asym.id`, `_atom_site.label_asym_id`,
`_struct_conn.ptnr1/2_label_asym_id`

### What it does

Gemmi-generated mmCIF files use multi-character, non-standard chain IDs
such as `Axp`, `Ax1`, `Ax02`.  OneDep requires single-letter uppercase
identifiers starting from `A`.

The remapper:
1. Reads all `_struct_asym.id` values in document order.
2. Assigns sequential labels `A, B, C, … Z, AA, AB, …` (the same
   sequence used by the PDB).
3. Replaces the IDs in `_struct_asym.id` in place.
4. Iterates over `_atom_site.label_asym_id` and substitutes every
   occurrence using the mapping built in step 2.
5. Applies the same substitution to `_struct_conn` partner chain columns
   if that category is present.

The mapping is exposed as `ChainRemapper.mapping` (a `dict[str, str]`)
so that the pipeline can re-key the `chain → entity_id` dictionary
produced by the EntityBuilder.

**Example**

| Before | After |
|---|---|
| `Axp` | `A` |
| `Ax1` | `B` |
| `Ax2` | `C` |
| `orphan:B` | `P` |
| `orphan:W` | `Q` |

---

## 3. EntityPolySeqBuilder — `transformers/entity_poly_seq.py`

**Categories written:** `_entity_poly_seq`, `_pdbx_prerelease_seq`

### What it does

Reads the canonical one-letter-code polymer sequence from
`_entity_poly.pdbx_seq_one_letter_code` and derives two new categories.

#### `_entity_poly_seq`

A row-per-residue enumeration of the canonical sequence:

| Column | Value |
|---|---|
| `entity_id` | from `_entity_poly` |
| `num` | 1-based position in the canonical sequence |
| `mon_id` | three-letter residue code (decoded from one-letter code) |
| `hetero` | always `n` |

The one-letter → three-letter decoding covers all 20 standard amino
acids plus selenocysteine (`U → SEC`), pyrrolysine (`O → PYL`), and
unknown (`X → UNK`).

#### `_pdbx_prerelease_seq`

A convenience pair of key-value items storing the raw one-letter-code
sequence for the first (or only) polymer entity:

```
_pdbx_prerelease_seq.entity_id            1
_pdbx_prerelease_seq.seq_one_letter_code  YHGSYEAPT…
```

---

## 4. AtomSiteTransformer — `transformers/atom_site.py`

**Category written:** `_atom_site` (multiple columns)

This is the most complex transformer.  It enriches the atom coordinate
table with several items required by OneDep.

### Column-by-column changes

| Column | Action | Detail |
|---|---|---|
| `id` | **Renumbered** | Sequential integers starting from 1, in document order |
| `label_entity_id` | **Assigned** | Looked up from the `chain → entity_id` mapping produced by EntityBuilder and re-keyed by ChainRemapper |
| `label_seq_id` | **Assigned** | For polymer atoms: canonical sequence position (see below). For HETATM / non-polymer / water: `.` |
| `auth_comp_id` | **Filled** | Copied from `label_comp_id` (source files often leave this blank) |
| `auth_atom_id` | **Filled** | Copied from `label_atom_id` when the original value is `?` or `.` |
| `pdbx_auth_seq_id` | **Added** | Copy of `auth_seq_id` |
| `pdbx_auth_comp_id` | **Added** | Copy of `auth_comp_id` (after the fill above) |
| `pdbx_auth_asym_id` | **Added** | Copy of `auth_asym_id` |
| `pdbx_auth_atom_name` | **Added** | Copy of `auth_atom_id` |

Columns `auth_comp_id` and `auth_atom_id` are created in the loop if
they are absent from the input (using `loop.add_columns`).

### `label_seq_id` derivation

The transformer tries two strategies in order:

1. **`_pdbx_poly_seq_scheme`** (preferred): if this category is present
   *and non-empty*, each polymer residue's `auth_seq_id` is mapped to
   its canonical sequence position (`seq_id`) via the scheme table.
   This yields the authoritative OneDep position (which may differ from
   the author numbering).

2. **Fallback — sequential per-chain numbering**: if the scheme is
   absent or empty, the transformer scans `_atom_site` and assigns
   `label_seq_id = 1, 2, 3, …` based on the order of first appearance
   of unique `auth_seq_id` values within each polymer chain.  This is
   a best-effort approximation; it will not reproduce the canonical
   sequence offset produced by OneDep when residues at the N-terminus
   are missing from the model.

---

## 5. ChemCompEnricher — `transformers/chem_comp.py`

**Category written:** `_chem_comp`

### What it does

Collects the unique component IDs already present in `_chem_comp.id`
and rewrites the loop with authoritative CCD (Chemical Component
Dictionary) data for the following columns:

| Column | Content |
|---|---|
| `id` | Three-letter code (unchanged) |
| `type` | Linking type, e.g. `L-peptide linking`, `non-polymer` |
| `mon_nstd_flag` | `y` for standard amino acids, `.` otherwise |
| `name` | Full chemical name |
| `pdbx_synonyms` | Alternative names (or `?`) |
| `formula` | Hill-notation molecular formula |
| `formula_weight` | Monoisotopic molecular weight |

### Data sources

Data are looked up in priority order:

1. **User-supplied CCD file** (`--ccd-path` CLI option): a standard
   wwPDB `components.cif` file.  Every block in the file is parsed; the
   block name is the component ID.  These entries override the built-in
   table.

2. **Built-in table**: a hardcoded dictionary covering the 20 standard
   amino acids and ~15 common solvents/ions (HOH, DMS, GOL, ZN, CL,
   MG, NA, CA, SO4, PO4, EDO, PEG, MPD, IMD, ACT).

For any component absent from both sources all data columns are set
to `?`.

---

## 6. AtomTypeBuilder — `transformers/atom_type.py`

**Category written:** `_atom_type`

### What it does

Collects the unique set of element symbols from
`_atom_site.type_symbol`, sorts them alphabetically, and rebuilds the
`_atom_type` loop with a single column:

```
loop_
_atom_type.symbol
C
CL
N
O
P
S
ZN
```

This matches the OneDep format which requires `_atom_type` to enumerate
every element present in the model exactly once, in alphabetical order.

---

## Known Limitations

### `label_seq_id` canonical offset

When `_pdbx_poly_seq_scheme` is absent or empty (common in Gemmi-generated
files), `label_seq_id` is assigned as a sequential counter starting
from 1 for the first *observed* residue.  The OneDep system derives this
number from a full sequence alignment against the depositor-provided
canonical sequence, so the value may differ for structures where
N-terminal or other residues are unmodelled.  Providing a pre-populated
`_pdbx_poly_seq_scheme` in the input resolves this.

### Non-standard ligand `_chem_comp` data

Ligands not in the built-in table and not covered by a supplied CCD
file will have all `_chem_comp` columns set to `?`.  Passing the full
wwPDB `components.cif` via `--ccd-path` resolves this.

### `_atom_site` row ordering

OneDep expects ATOM records before HETATM records, and waters last.
The current pipeline preserves the row order from the input file.
Reordering is explicitly deferred ("IGNORE FOR NOW" in the
specification).

---

## Installation and usage

```bash
# Install
cd StandardCif
pip install -e .

# Basic usage
standard-cif --input model.cif --output standard_model.cif

# With a full CCD file for complete _chem_comp data
standard-cif --input model.cif --output standard_model.cif \
             --ccd-path /data/components.cif

# Verbose logging to stderr
standard-cif --input model.cif --output standard_model.cif --verbose
```

Dependencies: **gemmi ≥ 0.6**, **click ≥ 8.0**, Python ≥ 3.10.
