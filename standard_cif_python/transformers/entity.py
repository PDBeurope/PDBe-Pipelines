"""
Builds and updates the _entity and _struct_asym.entity_id categories.

Entity assignment rules:
  - Polymer chains (contain standard amino acids): type = "polymer",
    src_method = copied from existing _entity if available.
  - Water chains (contain HOH/WAT): type = "water", src_method = "nat".
  - Non-polymer chains: type = "non-polymer", src_method = "syn".

Chains sharing the same unique set of component IDs are grouped into
the same entity (e.g. multiple DMS molecules are all entity "DMS").

Polymer entities are always assigned first (lowest entity IDs), followed
by non-polymer, then water.

The resulting chain→entity_id mapping is exposed as
`EntityBuilder.chain_to_entity` for downstream transformers.
"""
from __future__ import annotations

from collections import defaultdict

import gemmi

_POLYMER_RESIDUES: frozenset[str] = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS "
    "MET PHE PRO SER THR TRP TYR VAL".split()
)
_WATER_RESIDUES: frozenset[str] = frozenset({"HOH", "WAT", "H2O"})


def _as_str(value: str) -> str:
    """Strip CIF quoting from a raw gemmi value."""
    return gemmi.cif.as_string(value)


class EntityBuilder:
    """Transformer: generate _entity rows and assign _struct_asym.entity_id.

    After calling :meth:`transform` the attribute :attr:`chain_to_entity`
    holds the mapping ``label_asym_id → entity_id`` for use by
    :class:`AtomSiteTransformer`.
    """

    def __init__(self) -> None:
        self.chain_to_entity: dict[str, str] = {}

    def transform(self, block: gemmi.cif.Block) -> None:
        # Only operate on chains registered in _struct_asym
        struct_asym_ids = self._get_struct_asym_ids(block)
        chain_comps = self._gather_chain_compositions(block, struct_asym_ids)

        # Detect orphan atoms (label_asym_id='.') and create virtual chains
        orphan_comps = self._gather_orphan_compositions(block)
        chain_comps.update(orphan_comps)

        entity_map = self._assign_entities(chain_comps)
        self.chain_to_entity = {
            chain: str(eid) for chain, eid in entity_map.items()
        }

        # Add orphan chains to _struct_asym if not present
        if orphan_comps:
            self._add_orphan_chains_to_struct_asym(block, orphan_comps, entity_map)

        # Map orphan atom rows using their comp_id → entity
        self._map_orphan_atoms(block)

        self._update_struct_asym(block)
        self._write_entity_table(block, chain_comps, entity_map)

    # ------------------------------------------------------------------
    def _gather_orphan_compositions(
        self, block: gemmi.cif.Block
    ) -> dict[str, list[str]]:
        """
        Find atoms with label_asym_id='.' and group them into virtual
        chains keyed by their auth_asym_id.

        Returns {virtual_chain_key: sorted_comp_ids}.
        The virtual chain key uses the format 'orphan:<auth_asym_id>'.
        """
        table = block.find(
            "_atom_site.", ["label_asym_id", "label_comp_id", "auth_asym_id"]
        )
        if table is None:
            return {}

        raw: dict[str, set[str]] = defaultdict(set)
        for row in table:
            if row[0] == ".":
                key = f"orphan:{row[2]}"
                raw[key].add(row[1])

        return {key: sorted(comps) for key, comps in raw.items()}

    def _add_orphan_chains_to_struct_asym(
        self,
        block: gemmi.cif.Block,
        orphan_comps: dict[str, list[str]],
        entity_map: dict[str, int],
    ) -> None:
        """
        Append orphan chain entries to the _struct_asym loop so that
        the ChainRemapper will assign proper single-letter IDs to them.
        The entity_id is set from entity_map.
        Orphan chains are added in entity_id order (lower entity IDs first).
        """
        loop_col = block.find_loop("_struct_asym.id")
        if loop_col is None:
            return
        loop = loop_col.get_loop()
        tags = list(loop.tags)
        n_cols = loop.width()

        # Sort by entity_id so non-polymer orphans come before water
        sorted_keys = sorted(
            orphan_comps.keys(),
            key=lambda k: entity_map.get(k, 9999),
        )

        for key in sorted_keys:
            eid = str(entity_map.get(key, "?"))
            new_row = ["?"] * n_cols
            if "_struct_asym.id" in tags:
                new_row[tags.index("_struct_asym.id")] = key
            if "_struct_asym.entity_id" in tags:
                new_row[tags.index("_struct_asym.entity_id")] = eid
            loop.add_row(new_row, -1)

    def _map_orphan_atoms(self, block: gemmi.cif.Block) -> None:
        """
        Assign label_asym_id and label_entity_id for atoms whose
        label_asym_id is '.', using the orphan chain key
        'orphan:<auth_asym_id>'.
        """
        table = block.find(
            "_atom_site.",
            ["label_asym_id", "label_comp_id", "auth_asym_id",
             "label_entity_id"],
        )
        if table is None:
            return

        for row in table:
            if row[0] != ".":
                continue
            key = f"orphan:{row[2]}"
            eid = self.chain_to_entity.get(key, "?")
            row[0] = key   # set label_asym_id to the orphan key temporarily
            row[3] = eid   # set label_entity_id

    def _get_struct_asym_ids(self, block: gemmi.cif.Block) -> list[str]:
        """Return all chain IDs from _struct_asym in document order."""
        table = block.find("_struct_asym.", ["id"])
        if table is None:
            return []
        return [row[0] for row in table]

    def _gather_chain_compositions(
        self,
        block: gemmi.cif.Block,
        allowed_chains: list[str],
    ) -> dict[str, list[str]]:
        """
        Returns {label_asym_id: sorted_unique_comp_ids} for chains in
        allowed_chains only.  Atoms with '.' asym_id are skipped here
        (they share a chain with other atoms via auth_asym_id already).
        """
        allowed = set(allowed_chains)
        table = block.find(
            "_atom_site.", ["label_asym_id", "label_comp_id"]
        )
        if table is None:
            return {}

        raw: dict[str, set[str]] = defaultdict(set)
        for row in table:
            asym, comp = row[0], row[1]
            if asym in allowed:
                raw[asym].add(comp)

        # Preserve _struct_asym order
        return {chain: sorted(raw[chain]) for chain in allowed_chains if chain in raw}

    def _assign_entities(
        self,
        chain_comps: dict[str, list[str]],
    ) -> dict[str, int]:
        """
        Assign entity IDs.
        - Polymer chains get the lowest IDs (one per polymer chain).
        - Non-polymer chains sharing composition share an entity ID.
        - Water chains share one entity ID.
        """
        # Categorise chains
        polymer_chains: list[str] = []
        non_polymer_chains: list[str] = []
        water_chains: list[str] = []

        for chain, comps in chain_comps.items():
            comp_set = set(comps)
            if comp_set & _POLYMER_RESIDUES:
                polymer_chains.append(chain)
            elif comp_set <= _WATER_RESIDUES:
                water_chains.append(chain)
            else:
                non_polymer_chains.append(chain)

        chain_to_entity: dict[str, int] = {}
        next_id = 1

        # Polymers first
        for chain in polymer_chains:
            chain_to_entity[chain] = next_id
            next_id += 1

        # Non-polymer: group by composition
        comp_key_to_id: dict[str, int] = {}
        for chain in non_polymer_chains:
            key = "|".join(chain_comps[chain])
            if key not in comp_key_to_id:
                comp_key_to_id[key] = next_id
                next_id += 1
            chain_to_entity[chain] = comp_key_to_id[key]

        # Water: all share one entity
        if water_chains:
            water_id = next_id
            next_id += 1
            for chain in water_chains:
                chain_to_entity[chain] = water_id

        return chain_to_entity

    def _update_struct_asym(self, block: gemmi.cif.Block) -> None:
        """Write entity_id values into _struct_asym."""
        table = block.find("_struct_asym.", ["id", "entity_id"])
        if table is None:
            return
        for row in table:
            chain = row[0]
            row[1] = self.chain_to_entity.get(chain, "?")

    def _write_entity_table(
        self,
        block: gemmi.cif.Block,
        chain_comps: dict[str, list[str]],
        entity_map: dict[str, int],
    ) -> None:
        """Rebuild _entity with type, src_method, pdbx_description columns."""
        # Invert: entity_id → list of chains (preserve order)
        entity_chains: dict[int, list[str]] = defaultdict(list)
        for chain, eid in entity_map.items():
            entity_chains[eid].append(chain)

        # Read existing entity descriptions keyed by OLD entity type
        # (so we can preserve polymer entity descriptions regardless of new ID)
        existing_desc_by_type = self._read_existing_entity_descriptions_by_type(block)

        rows: dict[str, list[str]] = {
            "id": [], "type": [], "src_method": [], "pdbx_description": [],
            "formula_weight": [], "pdbx_number_of_molecules": [],
            "pdbx_ec": [], "pdbx_mutation": [], "pdbx_fragment": [],
            "details": [],
        }

        for eid in sorted(entity_chains.keys()):
            chains = entity_chains[eid]
            comps = set(chain_comps.get(chains[0], []))
            is_polymer = bool(comps & _POLYMER_RESIDUES)
            is_water = comps <= _WATER_RESIDUES

            if is_polymer:
                etype, src = "polymer", "man"
            elif is_water:
                etype, src = "water", "nat"
            else:
                etype, src = "non-polymer", "syn"

            desc = existing_desc_by_type.get(etype, "?")

            rows["id"].append(str(eid))
            rows["type"].append(etype)
            rows["src_method"].append(src)
            rows["pdbx_description"].append(desc)
            rows["formula_weight"].append("?")
            rows["pdbx_number_of_molecules"].append(str(len(chains)))
            rows["pdbx_ec"].append("?")
            rows["pdbx_mutation"].append("?")
            rows["pdbx_fragment"].append("?")
            rows["details"].append("?")

        block.set_mmcif_category("_entity.", rows)

    def _read_existing_entity_descriptions_by_type(
        self, block: gemmi.cif.Block
    ) -> dict[str, str]:
        """Return {entity_type: pdbx_description} from the existing _entity table.

        Useful for transferring descriptions (e.g. the polymer chain name)
        even when entity IDs change during standardisation.
        """
        result: dict[str, str] = {}
        table = block.find("_entity.", ["type", "pdbx_description"])
        if table is None:
            return result
        for row in table:
            etype = _as_str(row[0])
            desc = _as_str(row[1])
            if desc and desc != "?":
                result.setdefault(etype, desc)
        return result
