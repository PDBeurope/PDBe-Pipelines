process remove_alt_conformations {

    //container 'community.wave.seqera.io/library/pip_pdbeccdutils:de9d3afe23c1656a'
    container 'community.wave.seqera.io/library/pip_pdbeccdutils:2b3a61e65d7d9d21'

    input:
    tuple val(meta), path(cif_file)

    output:
    tuple val(meta), path("${meta.id}_processed.cif"), emit: no_alt_conf_cifs

    script:
    """
    remove_alt_conformers.py ${meta.id} ${cif_file} ./
    """
}
