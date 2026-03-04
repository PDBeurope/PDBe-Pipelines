process fix_assembly {

    container "community.wave.seqera.io/library/pip_gemmi:aa372a360d8aedc4"

    input:
    tuple val(meta), path(assembly_cif)
    output:
    tuple val(meta), path("fixed_${meta.id}.cif"), emit: fixed_cifs

    script:
    """
    fix_assembly.py \
        --input_file ${assembly_cif} \
        --output_file fixed_${meta.id}.cif
    """
}