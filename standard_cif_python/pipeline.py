"""
Orchestrates the full mmCIF standardisation pipeline.

Transformers run in the following order so each step can rely on the
results of the previous ones:

  1. EntityBuilder       – assigns entity IDs to chains; also updates
                           _struct_asym.entity_id.
  2. ChainRemapper       – renames Gemmi non-standard chain IDs to A/B/C/…
                           in both _struct_asym and _atom_site.
  3. EntityPolySeqBuilder– derives _entity_poly_seq and _pdbx_prerelease_seq
                           from the canonical polymer sequence.
  4. AtomSiteTransformer – enriches _atom_site (entity_id, seq_id, auth
                           columns, sequential ids).
  5. ChemCompEnricher    – enriches _chem_comp with CCD data.
  6. AtomTypeBuilder     – rebuilds _atom_type from element symbols.
"""
from __future__ import annotations

import logging
from pathlib import Path

import gemmi

from .transformers.chain_remapper import ChainRemapper
from .transformers.atom_site import AtomSiteTransformer
from .transformers.chem_comp import ChemCompEnricher
from .transformers.entity import EntityBuilder
from .transformers.entity_poly_seq import EntityPolySeqBuilder
from .transformers.atom_type import AtomTypeBuilder

log = logging.getLogger(__name__)


def standardise(
    input_path: Path,
    output_path: Path,
    ccd_path: Path | None = None,
) -> None:
    """Run the full standardisation pipeline on one mmCIF file.

    Args:
        input_path:  Path to the input mmCIF file.
        output_path: Destination path for the standardised mmCIF file.
        ccd_path:    Optional path to a CCD components mmCIF file used
                     to enrich _chem_comp entries.
    """
    log.info("Reading %s", input_path)
    doc = gemmi.cif.read(str(input_path))
    block = doc.sole_block()

    # 1. Build entities first (needed by ChainRemapper for _struct_asym)
    log.info("Building entity table …")
    entity_builder = EntityBuilder()
    entity_builder.transform(block)

    # 2. Remap chain IDs; exposes chain_to_entity after this call
    log.info("Remapping chain IDs …")
    chain_remapper = ChainRemapper()
    chain_remapper.transform(block)

    # Re-key entity mapping with the *new* chain IDs from the remapper
    new_chain_to_entity = {
        chain_remapper.mapping.get(old, old): eid
        for old, eid in entity_builder.chain_to_entity.items()
    }

    # 3. Derive polymer sequences
    log.info("Building entity poly-seq …")
    EntityPolySeqBuilder().transform(block)

    # 4. Enrich _atom_site
    log.info("Transforming atom_site …")
    AtomSiteTransformer(new_chain_to_entity).transform(block)

    # 5. Enrich _chem_comp from CCD
    log.info("Enriching chem_comp …")
    ChemCompEnricher(ccd_path=ccd_path).transform(block)

    # 6. Rebuild _atom_type
    log.info("Rebuilding atom_type …")
    AtomTypeBuilder().transform(block)

    log.info("Writing %s", output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.write_file(str(output_path))
    log.info("Done.")
