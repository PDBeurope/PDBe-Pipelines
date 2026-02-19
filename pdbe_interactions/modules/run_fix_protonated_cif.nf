process run_fix_protonated_cif {

    container 'community.wave.seqera.io/library/pip_gemmi:f669d80622e17dca'

    input:
    tuple val(meta), file(protonated_file), file(assembly_fixed_cif)

    output:
    tuple val(meta), file("${meta.id}_protonated_fixed.cif"), emit: protonated_fixed_cifs

    script:
    """
    fix_protonated.py \
        --assembly_file ${assembly_fixed_cif} \
        --input_protonated ${protonated_file} \
        --output_protonated ./${meta.id}_protonated_fixed.cif
    """
}