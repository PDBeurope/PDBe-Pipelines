"""
Enriches the _chem_comp category with data from the Chemical Component
Dictionary (CCD).

The CCD provides authoritative values for:
  type, mon_nstd_flag, name, pdbx_synonyms, formula, formula_weight

Data sources (tried in priority order):
  1. A user-supplied CCD mmCIF file (--ccd-path CLI option).
  2. A hardcoded table of the 20 standard amino acids + common solvent molecules.

For non-standard residues that are absent from both sources the values
remain '?' (unknown).
"""
from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import gemmi

# ── static CCD data for common residues ──────────────────────────────────────

class _CompData(NamedTuple):
    type: str
    mon_nstd_flag: str  # 'y' for standard, '.' for non-polymer
    name: str
    pdbx_synonyms: str
    formula: str
    formula_weight: str


# Standard 20 amino acids + common small molecules found in PDB structures
_BUILTIN_CCD: dict[str, _CompData] = {
    "ALA": _CompData("L-peptide linking", "y", "ALANINE", "?", "C3 H7 N O2", "89.093"),
    "ARG": _CompData("L-peptide linking", "y", "ARGININE", "?", "C6 H15 N4 O2 1", "175.209"),
    "ASN": _CompData("L-peptide linking", "y", "ASPARAGINE", "?", "C4 H8 N2 O3", "132.118"),
    "ASP": _CompData("L-peptide linking", "y", "ASPARTIC ACID", "?", "C4 H7 N O4", "133.103"),
    "CYS": _CompData("L-peptide linking", "y", "CYSTEINE", "?", "C3 H7 N O2 S", "121.158"),
    "GLN": _CompData("L-peptide linking", "y", "GLUTAMINE", "?", "C5 H10 N2 O3", "146.144"),
    "GLU": _CompData("L-peptide linking", "y", "GLUTAMIC ACID", "?", "C5 H9 N O4", "147.129"),
    "GLY": _CompData("peptide linking", "y", "GLYCINE", "?", "C2 H5 N O2", "75.067"),
    "HIS": _CompData("L-peptide linking", "y", "HISTIDINE", "?", "C6 H10 N3 O2 1", "156.162"),
    "ILE": _CompData("L-peptide linking", "y", "ISOLEUCINE", "?", "C6 H13 N O2", "131.173"),
    "LEU": _CompData("L-peptide linking", "y", "LEUCINE", "?", "C6 H13 N O2", "131.173"),
    "LYS": _CompData("L-peptide linking", "y", "LYSINE", "?", "C6 H15 N2 O2 1", "147.195"),
    "MET": _CompData("L-peptide linking", "y", "METHIONINE", "?", "C5 H11 N O2 S", "149.211"),
    "PHE": _CompData("L-peptide linking", "y", "PHENYLALANINE", "?", "C9 H11 N O2", "165.189"),
    "PRO": _CompData("L-peptide linking", "y", "PROLINE", "?", "C5 H9 N O2", "115.130"),
    "SER": _CompData("L-peptide linking", "y", "SERINE", "?", "C3 H7 N O3", "105.093"),
    "THR": _CompData("L-peptide linking", "y", "THREONINE", "?", "C4 H9 N O3", "119.119"),
    "TRP": _CompData("L-peptide linking", "y", "TRYPTOPHAN", "?", "C11 H12 N2 O2", "204.225"),
    "TYR": _CompData("L-peptide linking", "y", "TYROSINE", "?", "C9 H11 N O3", "181.189"),
    "VAL": _CompData("L-peptide linking", "y", "VALINE", "?", "C5 H11 N O2", "117.146"),
    # Common solvent / ions
    "HOH": _CompData("non-polymer", ".", "WATER", "?", "H2 O", "18.015"),
    "WAT": _CompData("non-polymer", ".", "WATER", "?", "H2 O", "18.015"),
    "GOL": _CompData("non-polymer", ".", "GLYCEROL", "GLYCERIN; PROPANE-1,2,3-TRIOL", "C3 H8 O3", "92.094"),
    "DMS": _CompData("non-polymer", ".", "DIMETHYL SULFOXIDE", "?", "C2 H6 O S", "78.133"),
    "ZN":  _CompData("non-polymer", ".", "ZINC ION", "?", "Zn 2", "65.409"),
    "CL":  _CompData("non-polymer", ".", "CHLORIDE ION", "?", "Cl -1", "35.453"),
    "MG":  _CompData("non-polymer", ".", "MAGNESIUM ION", "?", "Mg 2", "24.305"),
    "NA":  _CompData("non-polymer", ".", "SODIUM ION", "?", "Na 1", "22.989"),
    "CA":  _CompData("non-polymer", ".", "CALCIUM ION", "?", "Ca 2", "40.078"),
    "SO4": _CompData("non-polymer", ".", "SULFATE ION", "?", "O4 S", "96.062"),
    "PO4": _CompData("non-polymer", ".", "PHOSPHATE ION", "?", "H2 O4 P", "95.979"),
    "EDO": _CompData("non-polymer", ".", "1,2-ETHANEDIOL", "ETHYLENE GLYCOL", "C2 H6 O2", "62.068"),
    "PEG": _CompData("non-polymer", ".", "DI(HYDROXYETHYL)ETHER", "?", "C4 H10 O3", "106.121"),
    "MPD": _CompData("non-polymer", ".", "(4R)-2-METHYLPENTANE-2,4-DIOL", "?", "C6 H14 O2", "118.174"),
    "IMD": _CompData("non-polymer", ".", "IMIDAZOLE", "?", "C3 H4 N2", "68.077"),
    "ACT": _CompData("non-polymer", ".", "ACETATE ION", "?", "C2 H3 O2", "59.044"),
}


_TARGET_TAGS = [
    "id",
    "type",
    "mon_nstd_flag",
    "name",
    "pdbx_synonyms",
    "formula",
    "formula_weight",
]


class ChemCompEnricher:
    """Transformer: enrich _chem_comp with CCD data.

    Args:
        ccd_path: Optional path to a CCD mmCIF file.  When supplied every
                  component in the file overrides the built-in table.
    """

    def __init__(self, ccd_path: Path | None = None) -> None:
        self._ccd: dict[str, _CompData] = dict(_BUILTIN_CCD)
        if ccd_path is not None:
            self._load_ccd(ccd_path)

    def transform(self, block: gemmi.cif.Block) -> None:
        """Rewrite _chem_comp loop with enriched CCD data."""
        ids = self._collect_ids(block)
        if not ids:
            return

        # Build loop items
        rows: list[list[str]] = []
        for comp_id in ids:
            data = self._ccd.get(comp_id)
            if data:
                rows.append([
                    comp_id,
                    data.type,
                    data.mon_nstd_flag,
                    data.name,
                    data.pdbx_synonyms,
                    data.formula,
                    data.formula_weight,
                ])
            else:
                rows.append([comp_id, "?", "?", "?", "?", "?", "?"])

        # Remove existing _chem_comp items and rebuild
        block.find_mmcif_category("_chem_comp.")  # ensure category is present
        self._write_loop(block, rows)

    # ------------------------------------------------------------------
    def _collect_ids(self, block: gemmi.cif.Block) -> list[str]:
        """Return unique comp_ids in document order from _chem_comp."""
        table = block.find("_chem_comp.", ["id"])
        if table is None:
            return []
        seen: set[str] = set()
        ids: list[str] = []
        for row in table:
            cid = row[0]
            if cid not in seen:
                seen.add(cid)
                ids.append(cid)
        return ids

    def _write_loop(
        self, block: gemmi.cif.Block, rows: list[list[str]]
    ) -> None:
        """Recreate the _chem_comp loop with enriched columns."""
        col_data: dict[str, list[str]] = {t: [] for t in _TARGET_TAGS}
        for row in rows:
            for i, tag in enumerate(_TARGET_TAGS):
                col_data[tag].append(row[i])
        block.set_mmcif_category("_chem_comp.", col_data)

    def _load_ccd(self, path: Path) -> None:
        """Load CCD data from an mmCIF file and merge into the lookup table."""
        doc = gemmi.cif.read(str(path))
        for ccd_block in doc:
            comp_id = ccd_block.name
            table = ccd_block.find(
                "_chem_comp.",
                ["type", "mon_nstd_flag", "name", "pdbx_synonyms",
                 "formula", "formula_weight"],
            )
            if table is None:
                continue
            for row in table:
                self._ccd[comp_id] = _CompData(
                    type=gemmi.cif.as_string(row[0]),
                    mon_nstd_flag=gemmi.cif.as_string(row[1]),
                    name=gemmi.cif.as_string(row[2]),
                    pdbx_synonyms=gemmi.cif.as_string(row[3]),
                    formula=gemmi.cif.as_string(row[4]),
                    formula_weight=gemmi.cif.as_string(row[5]),
                )
                break  # one row per block


