from pdbeccdutils.helpers import cif_tools
import os

def remove_alt_conformers(input_id, input_file, out_dir):
    fixed_mmcif_file = os.path.join(out_dir, f"{input_id}_processed.cif")
    cif_tools.fix_updated_mmcif(input_file, fixed_mmcif_file)
    if not os.path.isfile(fixed_mmcif_file):
        raise Exception(
                f"Preprocessing of {input_id} failed"
            )
    return fixed_mmcif_file
