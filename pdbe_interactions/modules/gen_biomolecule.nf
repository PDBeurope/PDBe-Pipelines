process gen_biomolecule {

    //container 'community.wave.seqera.io/library/pip_pdbeccdutils:2b3a61e65d7d9d21'
    container 'community.wave.seqera.io/library/xorg-libxrender_pip_pdbeccdutils_rdkit:5021caa5f0b867f2'
    // https://wave.seqera.io/view/builds/bd-4fe84714d5258193_1
    input:
    tuple val(meta), file(bioassembly_file)

    output:
    tuple val(meta), file("${meta.id}_bound_molecules.json"), emit: biomolecule_jsons

    script:
    """
    generate_bm.py \
        --input_id ${meta.id} \
        --input_file ${bioassembly_file} \
        --output_dir ./
    """
}