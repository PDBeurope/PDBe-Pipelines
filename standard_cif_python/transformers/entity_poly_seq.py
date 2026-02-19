"""
Builds _entity_poly_seq from _entity_poly.pdbx_seq_one_letter_code.

The _entity_poly_seq category lists each residue of the canonical polymer
sequence:
  entity_id | num | mon_id | hetero

The one-letter-code sequence is decoded to three-letter codes using a
standard mapping.  Non-standard residues not in the mapping are written
as 'UNK'.

Also populates _pdbx_prerelease_seq with the same sequence if the category
is not already present (or updates it).
"""
from __future__ import annotations

import gemmi

# One-letter → three-letter amino acid mapping
_ONE_TO_THREE: dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    "U": "SEC", "O": "PYL", "X": "UNK",
}


class EntityPolySeqBuilder:
    """Transformer: build _entity_poly_seq and refresh _pdbx_prerelease_seq."""

    def transform(self, block: gemmi.cif.Block) -> None:
        poly_seqs = self._read_sequences(block)
        if not poly_seqs:
            return
        self._write_entity_poly_seq(block, poly_seqs)
        self._write_prerelease_seq(block, poly_seqs)

    # ------------------------------------------------------------------
    def _read_sequences(
        self, block: gemmi.cif.Block
    ) -> dict[str, str]:
        """Return {entity_id: one_letter_sequence} from _entity_poly."""
        result: dict[str, str] = {}
        table = block.find(
            "_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code"]
        )
        if table is None:
            return result
        for row in table:
            entity_id = row[0]
            raw_seq = gemmi.cif.as_string(row[1])
            # Strip newlines and spaces (sequences can be multi-line)
            seq = "".join(raw_seq.split())
            result[entity_id] = seq
        return result

    def _write_entity_poly_seq(
        self,
        block: gemmi.cif.Block,
        poly_seqs: dict[str, str],
    ) -> None:
        """Write _entity_poly_seq loop."""
        ids: list[str] = []
        nums: list[str] = []
        mons: list[str] = []
        heteros: list[str] = []

        for entity_id, seq in sorted(poly_seqs.items()):
            for pos, letter in enumerate(seq, start=1):
                mon = _ONE_TO_THREE.get(letter, "UNK")
                ids.append(entity_id)
                nums.append(str(pos))
                mons.append(mon)
                heteros.append("n")

        block.set_mmcif_category("_entity_poly_seq.", {
            "entity_id": ids,
            "num": nums,
            "mon_id": mons,
            "hetero": heteros,
        })

    def _write_prerelease_seq(
        self,
        block: gemmi.cif.Block,
        poly_seqs: dict[str, str],
    ) -> None:
        """Write or update _pdbx_prerelease_seq."""
        # Only write for the first (or single) polymer entity
        if not poly_seqs:
            return

        # If multiple entities, write one row per entity
        for entity_id, seq in sorted(poly_seqs.items()):
            # Format sequence in 80-char wrapped lines (mmCIF semicolon style
            # is handled by gemmi when writing; we just store the plain seq)
            block.set_pair(
                "_pdbx_prerelease_seq.entity_id", entity_id
            )
            block.set_pair(
                "_pdbx_prerelease_seq.seq_one_letter_code",
                gemmi.cif.quote(seq),
            )
            # Only write first entity here (multi-entity would need a loop)
            break
