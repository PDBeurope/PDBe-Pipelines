"""
Builds _atom_type from the unique element symbols in _atom_site.type_symbol.

The resulting loop contains a single column (_atom_type.symbol) with
symbols sorted alphabetically, matching the OneDep standard.
"""
from __future__ import annotations

import gemmi


class AtomTypeBuilder:
    """Transformer: rebuild _atom_type from _atom_site element symbols."""

    def transform(self, block: gemmi.cif.Block) -> None:
        symbols = self._collect_symbols(block)
        if not symbols:
            return
        block.set_mmcif_category("_atom_type.", {"symbol": sorted(symbols)})

    def _collect_symbols(self, block: gemmi.cif.Block) -> set[str]:
        table = block.find("_atom_site.", ["type_symbol"])
        if table is None:
            return set()
        return {row[0] for row in table if row[0] not in (".", "?")}
