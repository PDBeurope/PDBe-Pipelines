"""
Remaps non-standard Gemmi chain IDs (e.g. Axp, Ax1) in _struct_asym
to sequential single-letter uppercase identifiers (A, B, C, …, Z, AA, AB, …).

Also updates the chain ID column inside _atom_site.label_asym_id so
the two categories remain consistent.  A mapping dict is stored on the
instance so that other transformers can look up old→new IDs.
"""
import string
from itertools import product

import gemmi


def _label_generator():
    """Yields A, B, … Z, AA, AB, … in order."""
    for length in range(1, 10):
        for combo in product(string.ascii_uppercase, repeat=length):
            yield "".join(combo)


class ChainRemapper:
    """Transformer: standardise _struct_asym chain IDs.

    After calling :meth:`transform`, :attr:`mapping` holds the
    ``old_id → new_id`` mapping for all chains.
    """

    def __init__(self) -> None:
        self.mapping: dict[str, str] = {}

    def transform(self, block: gemmi.cif.Block) -> None:
        struct_asym = block.find("_struct_asym.", ["id"])
        if struct_asym is None:
            return

        gen = _label_generator()
        old_ids = [row[0] for row in struct_asym]
        self.mapping = {old: next(gen) for old in old_ids}

        # Update _struct_asym.id
        for row in struct_asym:
            row[0] = self.mapping[row[0]]

        # Update _atom_site.label_asym_id
        atom_site = block.find("_atom_site.", ["label_asym_id"])
        if atom_site is not None:
            for row in atom_site:
                old = row[0]
                if old in self.mapping:
                    row[0] = self.mapping[old]

        # Update _struct_conn partner chain references
        for col in ["ptnr1_label_asym_id", "ptnr2_label_asym_id"]:
            table = block.find("_struct_conn.", [col])
            if table is not None:
                for row in table:
                    if row[0] in self.mapping:
                        row[0] = self.mapping[row[0]]
