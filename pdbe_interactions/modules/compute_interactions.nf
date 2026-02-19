process compute_interactions {

    container 'community.wave.seqera.io/library/openbabel_pip_pdbe-arpeggio_pdbeccdutils:a733d295728b5b70'
    // https://wave.seqera.io/view/builds/bd-4fe84714d5258193_1
    
    input:
    tuple val(meta), file(protonated_file)

    output:
    tuple val(meta), file("${meta.id}_interactions.json"), emit: interactions_jsons

    script:
    """
    compute_interactions.py ${meta.id} ${protonated_file}
    mv interactions.json ${meta.id}_interactions.json
    """
}