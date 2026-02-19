#!/usr/bin/env bash -C -e -u -o pipefail
mkdir out_dir
copy_file.py --input structure_zika.cif --output out_dir/test.cif
