#!/usr/bin/env python

from pdbeccdutils.core.boundmolecule import BoundMolecule


from arpeggio.core import InteractionComplex
from arpeggio.core import exceptions as ar_exceptions

import os
import json
import argparse
import pickle as pck


def get_orig_auth_and_operator(
    contacts: dict[str, dict[str, str]]
) -> dict[str, dict[str, str]]:
    """Adds orig_auth_asym_id and operator fields to contacts

    Args:
        contacts: A dictionary containing intreactions information from arpeggio

    Returns:
        A dictionary of contacts from arpeggio after adding orig_auth_asym_id and
        operator fields
    """

    if "_" in contacts["bgn"]["auth_asym_id"]:
        contacts["bgn"]["orig_auth_asym_id"] = contacts["bgn"]["auth_asym_id"].split(
            "_"
        )[0]
        contacts["bgn"]["operator"] = (
            "_" + contacts["bgn"]["auth_asym_id"].split("_")[1]
        )
    else:
        contacts["bgn"]["orig_auth_asym_id"] = contacts["bgn"]["auth_asym_id"]
        contacts["bgn"]["operator"] = ""

    if "_" in contacts["end"]["auth_asym_id"]:
        contacts["end"]["orig_auth_asym_id"] = contacts["end"]["auth_asym_id"].split(
            "_"
        )[0]
        contacts["end"]["operator"] = (
            "_" + contacts["end"]["auth_asym_id"].split("_")[1]
        )
    else:
        contacts["end"]["orig_auth_asym_id"] = contacts["end"]["auth_asym_id"]
        contacts["end"]["operator"] = ""
    return contacts

class ProtLigInteractions:
    """Wrapper for the Arpeggio computation. Parses bound molecules and
    calculates interactions.
    """

    def __init__(self, structure, to_discard=None):
        """Create protein - ligand interaction object.

        Args:
            structure (str): Path to the structure to be processed
            to_discard (list, optional): Defaults to []. list of residue
                names to be discarded prior to protein-ligand interaction
                lookup.
        """
        self.compl = None
        self.path = structure

        to_discard = [] if to_discard is None else to_discard

    def initialize(self):
        """Set up Arpeggio basics. Create all the internal representation"""
        self.compl = InteractionComplex(self.path)
        self.compl.structure_checks()
        self.compl.address_ambiguities()
        self.compl.initialize()

    def get_all_interactions(
        self,
        bound_molecules,
        interaction_cutoff=5.0,
        compensation_factor=0.1,
        include_neighbors=False,
    ):
        """Retrieve interactions for all bound molecules found in the
        entry.

        Args:
            interaction_cutoff (float, optional): Defaults to 5.0. Distance
                cutoff for grid points to be `interacting` with the entity.
            compensation_factor (float, optional): Defaults to 0.1.
                Compensation factor for VdW radii dependent interaction types.
            include_neighbors (bool, optional): Defaults to False. Include
                non-bonding interactions between residues that are next to
                each other in sequence.

        Returns:
            list of dict of str: list of interactions in a dictionary like schema.
        """
        results = {}

        for i, bm in enumerate(bound_molecules, start=1):
            results[f"bm{i}"] = self.get_interaction(
                bm, interaction_cutoff, compensation_factor, include_neighbors
            )

        return results

    def get_interaction(
        self,
        bm,
        interaction_cutoff=5.0,
        compensation_factor=0.1,
        include_neighbors=False,
    ):
        """Retrieve interactions in the protein within a given selection.

        Args:
            bm (BoundMolecule): Bound molecule
            interaction_cutoff (float, optional): Defaults to 5.0. Distance
                cutoff for grid points to be `interacting` with the entity.
            compensation_factor (float, optional): Defaults to 0.1.
                Compensation factor for VdW radii dependent interaction types.
            include_neighbors (bool, optional): Defaults to False. Include
                non-bonding interactions between residues that are next to
                each other in sequence.

        Returns:
            dict of str: Interactions in a dictionary like schema.
        """
        selection = bm.to_arpeggio()
        self.compl.run_arpeggio(
            selection, interaction_cutoff, compensation_factor, include_neighbors
        )

        return self.compl.get_contacts()

def generate_interactions(
        pdb_id:str, 
        structure: str,
        bound_molecules: list[BoundMolecule],
        discarded_ligands: list[str],
        out: str
    ) -> None:
        """Generates interactions among residues within 5A⁰ of each bound-molecule
        in the input structure

        Args:
            structure: Path to protonated assembly file
            bound_molecules: List of BoundMolecule instances

        Raises:
            SelectionError: If the selection of BoundMolecule to arpeggio is invalid

        """

        # interactions.json
        result_interactions = {"entry": pdb_id, "boundMoleculeInteractions": []}

        interactions = ProtLigInteractions(structure, discarded_ligands)

        interactions.initialize()

        for i, bm in enumerate(bound_molecules, start=1):

            try:
                contacts = interactions.get_interaction(bm)
                contacts_filtered = [
                    c
                    for c in contacts
                    if c["interacting_entities"]
                    in ("INTER", "INTRA_SELECTION", "SELECTION_WATER")
                ]
                for c in contacts_filtered:
                    c = get_orig_auth_and_operator(c)

                result_interactions["boundMoleculeInteractions"].append(
                    {
                        "id": f"bm{i}",
                        "contacts": contacts_filtered,
                    }
                )

            except ar_exceptions.SelectionError as e:
                if str(e) != "entity not found":
                    raise e

        with open(os.path.join(out, "interactions.json"), "w") as f:
            json.dump(result_interactions, f, sort_keys=True, indent=4)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute interactions for a structure")
    parser.add_argument(
        "structure_id",
        type=str,
        help="PDB ID of the structure to be processed"
    )

    parser.add_argument(
        "cif_structure",
        type=str,
        help="Path to the cif structure to be processed"
    )

    parser.add_argument(
        "--is_assembly",
        action="store_true",
        default=False,
        help="Whether the input structure is a biological assembly"
    )
    args = parser.parse_args()
    # discard water molecules from interactions, 
    # TODO: this should be a parameter in the future
    discarded_ligands = ["HOH"]
    
    from generate_bm import generate_boundmolecules

    bound_molecules, clc_reader_results = generate_boundmolecules(args.cif_structure,
                                                                    discarded_ligands,
                                                                    args.is_assembly)

    generate_interactions(
        args.structure_id,
        args.cif_structure,
        bound_molecules,
        discarded_ligands,
        "./"
    )