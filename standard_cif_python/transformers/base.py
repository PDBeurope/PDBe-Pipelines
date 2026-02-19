"""Base protocol for CIF category transformers."""
from typing import Protocol
import gemmi


class CifTransformer(Protocol):
    """Each transformer operates on a gemmi.cif.Block in place."""

    def transform(self, block: gemmi.cif.Block) -> None:
        """Apply transformation to the mmCIF data block."""
        ...
