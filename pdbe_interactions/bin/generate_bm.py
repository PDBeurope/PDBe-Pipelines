#!/usr/bin/env python3

import os
import json
import argparse
from pdbeccdutils.core import clc_reader
from pdbeccutils.core.clc_reader import CLCReaderResult
from pdbeccdutils.core.boundmolecule import infer_bound_molecules, BoundMolecule

def generate_boundmolecules(
        structure: str,
        discarded_ligands: list[str],
        is_assembly: bool
        ) -> list[BoundMolecule]:
    """
    Identifies bound-molecules and writes them to bound_molecules.json
    Args:
        structure: Path to preferred assembly
        discarded_ligands: Lit of ligands to ignore
        is_assembly: Indicates whether the input is assembly or asymmetric unit
    Returns:
        List of BoundMolecule instances
    """
    bound_molecules = infer_bound_molecules(
        structure, discarded_ligands, assembly=is_assembly
    )
    
    clc_reader_results = []
    for i, bm in enumerate(bound_molecules, start=1):
        bm_id = f"bm{i}"
        reader_result = clc_reader.infer_multiple_chem_comp(structure, bm, bm_id, sanitize=True)
        if reader_result:
            clc_reader_results.append(reader_result)

    return (bound_molecules, clc_reader_results)


def write_out_bm(
        input_id: str,
        bound_molecules: list[BoundMolecule],
        clc_reader_results: list[CLCReaderResult],
        is_assembly: bool,
        out_dir: str
) -> None:
    """Writes details of bound-molecules to json file
    Args:
        input_id: 
        bound_molecules: List of BoundMolecule instances
        clc_reader_results: List of CLCReaderResult 
        is_assembly: Indicates whether the assembly file was asymmetric unit or not
    """
    result_bag = {
        "entry": input_id,
        "boundMolecules": [],
        "is_assembly": is_assembly,
    }
    for i, bm in enumerate(bound_molecules, start=1):
        bm_id = f"bm{i}"
        inchi = ''
        inchikey = ''
        for clc_result in clc_reader_results:
            component = clc_result.component
            if bm.is_equivalent(clc_result.bound_molecule):
                inchi = component.inchi
                inchikey = component.inchikey
        result_bag["boundMolecules"].append(
            {"id": bm_id, "composition": bm.to_dict(),
             "inchi": inchi, "inchikey":inchikey}
        )
    bm_file = os.path.join(out_dir, f"{input_id}_bound_molecules.json")
    with open(bm_file, "w") as f:
        json.dump(result_bag, f, sort_keys=True, indent=4)


def create_parser():
    
    parser = argparse.ArgumentParser(
        description="Fix updated mmCIF file and remove alternate conformations."
    )
    parser.add_argument(
        "--input_id",
        help="Input ID (e.g., PDB ID)"
    )
    parser.add_argument(
        "--input_file",
        help="Path to input mmCIF file"
    )

    parser.add_argument(
        "--is_assembly",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--discard",
        help="Comma separated string of ligands to ignore",
        default="HOH,UNK",
        required=False
    )

    parser.add_argument(
        "--output_dir",
        help="Directory where processed file will be written"
    )

def main():
    
    parser = create_parser()
    args = parser.parse_args()
    discarded_ligands = args.discard.split(",")

    (bound_molecules, clc_reader_results) = generate_boundmolecules(args.input_file,
                                                                    discarded_ligands,
                                                                    args.is_assembly)

    if len(bound_molecules) > 0:
        write_out_bm(args.input_id, 
                    bound_molecules,
                    clc_reader_results,
                    args.is_assembly,
                    args.output_dir)