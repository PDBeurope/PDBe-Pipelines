"""
Transforms _atom_site to OneDep standard format.

Changes applied:
  - label_entity_id: assigned from chain→entity mapping.
  - label_seq_id:    assigned for polymer atoms from _pdbx_poly_seq_scheme
                     (if present) or computed sequentially; HETATM gets '.'.
  - auth_comp_id:    copied from label_comp_id.
  - auth_atom_id:    copied from label_atom_id when not present.
  - pdbx_auth_seq_id / pdbx_auth_comp_id / pdbx_auth_asym_id /
    pdbx_auth_atom_name: copied from corresponding auth_* columns.
  - id:              renumbered sequentially starting from 1.

Depends on ChainRemapper having already remapped label_asym_id.
"""
from __future__ import annotations

from collections import defaultdict

import gemmi

_POLYMER_RESIDUES: frozenset[str] = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS "
    "MET PHE PRO SER THR TRP TYR VAL".split()
)


class AtomSiteTransformer:
    """Enrich _atom_site with standard OneDep fields.

    Args:
        chain_to_entity: mapping of (new) label_asym_id → entity_id string.
    """

    def __init__(self, chain_to_entity: dict[str, str]) -> None:
        self._chain_to_entity = chain_to_entity

    def transform(self, block: gemmi.cif.Block) -> None:
        # First ensure optional columns exist in the loop
        loop_col = block.find_loop("_atom_site.id")
        if loop_col is None:
            return
        loop = loop_col.get_loop()
        self._ensure_columns(loop, ["auth_comp_id", "auth_atom_id"])

        seq_map = self._build_seq_mapping(block)

        # Use the block.find() table which gives mutable row proxies
        table = block.find("_atom_site.", [
            "id",             # 0
            "label_asym_id",  # 1
            "label_entity_id",# 2
            "label_seq_id",   # 3
            "label_comp_id",  # 4
            "label_atom_id",  # 5
            "auth_seq_id",    # 6
            "auth_comp_id",   # 7
            "auth_atom_id",   # 8
            "auth_asym_id",   # 9 (read-only, for pdbx_auth_asym_id)
        ])
        if table is None:
            return

        for i, row in enumerate(table):
            asym = row[1]
            comp = row[4]
            atom = row[5]
            auth_seq = row[6]

            row[0] = str(i + 1)  # renumber id
            row[2] = self._chain_to_entity.get(asym, "?")  # label_entity_id

            # label_seq_id
            if comp in _POLYMER_RESIDUES:
                row[3] = seq_map.get((asym, auth_seq), ".")
            else:
                row[3] = "."

            # auth_comp_id ← label_comp_id
            row[7] = comp

            # auth_atom_id ← label_atom_id if missing
            if row[8] in ("?", "."):
                row[8] = atom

        # Add pdbx_auth_* columns after the core modifications
        self._add_pdbx_auth_columns(block, loop)

    # ------------------------------------------------------------------
    def _ensure_columns(
        self, loop: gemmi.cif.Loop, names: list[str]
    ) -> None:
        """Add missing columns to the loop with default value '?'."""
        existing = set(loop.tags)
        for name in names:
            full = f"_atom_site.{name}"
            if full not in existing:
                loop.add_columns([full], "?", -1)

    def _build_seq_mapping(
        self, block: gemmi.cif.Block
    ) -> dict[tuple[str, str], str]:
        """Return (label_asym_id, auth_seq_id) → label_seq_id."""
        mapping: dict[tuple[str, str], str] = {}

        scheme = block.find(
            "_pdbx_poly_seq_scheme.",
            ["asym_id", "seq_id", "pdb_seq_num"],
        )
        if scheme is not None:
            for row in scheme:
                asym, seq_id, auth_seq = row[0], row[1], row[2]
                if auth_seq not in ("?", "."):
                    mapping[(asym, auth_seq)] = seq_id
            if mapping:
                return mapping
            # Empty scheme — fall through to fallback

        # Fallback: sequential per-chain numbering for polymer atoms
        atom_table = block.find(
            "_atom_site.", ["label_asym_id", "auth_seq_id", "label_comp_id"]
        )
        if atom_table is None:
            return mapping

        seen: dict[str, list[str]] = defaultdict(list)
        for row in atom_table:
            asym, auth_seq, comp = row[0], row[1], row[2]
            if comp in _POLYMER_RESIDUES and auth_seq not in seen[asym]:
                seen[asym].append(auth_seq)

        for asym, seqs in seen.items():
            for idx, auth_seq in enumerate(seqs, start=1):
                mapping[(asym, auth_seq)] = str(idx)

        return mapping

    def _add_pdbx_auth_columns(
        self,
        block: gemmi.cif.Block,
        loop: gemmi.cif.Loop,
    ) -> None:
        """Append pdbx_auth_* columns mirroring the auth_* columns."""
        existing_tags = set(loop.tags)
        pairs = [
            # (new pdbx_auth tag,        source tag)
            ("_atom_site.pdbx_auth_seq_id",    "_atom_site.auth_seq_id"),
            ("_atom_site.pdbx_auth_comp_id",   "_atom_site.auth_comp_id"),
            ("_atom_site.pdbx_auth_asym_id",   "_atom_site.auth_asym_id"),
            ("_atom_site.pdbx_auth_atom_name", "_atom_site.auth_atom_id"),
        ]
        to_add = [(dst, src) for dst, src in pairs if dst not in existing_tags]
        if not to_add:
            return

        # Add placeholder columns to the loop
        for dst, _src in to_add:
            loop.add_columns([dst], "?", -1)

        # Now fill each new column by reading source via block.find()
        for dst, src in to_add:
            src_name = src.replace("_atom_site.", "")
            dst_name = dst.replace("_atom_site.", "")
            src_col = block.find("_atom_site.", [src_name])
            dst_col = block.find("_atom_site.", [dst_name])
            if src_col is None or dst_col is None:
                continue
            for src_row, dst_row in zip(src_col, dst_col):
                dst_row[0] = src_row[0]
