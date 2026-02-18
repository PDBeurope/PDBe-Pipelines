#!/usr/bin/env python3

import os
import argparse
from pdbeccdutils.helpers import cif_tools

def remove_alt_conf(input_id, output_dir):
    fixed_mmcif_file = os.path.join(output_dir, f"{input_id}_processed.cif")
    cif_tools.fix_updated_mmcif(input_id, fixed_mmcif_file)
    if not os.path.isfile(fixed_mmcif_file):
        raise Exception(
            f"Preprocessing of {input_id} failed"
            )
    return fixed_mmcif_file


def main():
    parser = argparse.ArgumentParser(
        description="Fix updated mmCIF file and remove alternate conformations."
    )
    parser.add_argument(
        "input_id",
        help="Input ID (e.g., PDB ID)"
    )
    parser.add_argument(
        "output_dir",
        help="Directory where processed file will be written"
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    output_file = remove_alt_conf(args.input_id, args.output_dir)
    print(f"Processed file written to: {output_file}")
    