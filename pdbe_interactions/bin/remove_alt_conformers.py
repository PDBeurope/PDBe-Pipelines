import os
from pdbeccdutils.helpers import cif_tools

def remove_alt_conf(input_id, output_dir):
    fixed_mmcif_file = os.path.join(output_dir, f"{input_id}_processed.cif")
    cif_tools.fix_updated_mmcif(input_id, fixed_mmcif_file)
    if not os.path.isfile(fixed_mmcif_file):
        raise Exception(
            f"Preprocessing of {input_id} failed"
            )
    return fixed_mmcif_file

