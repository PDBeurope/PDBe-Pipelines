process run_protonation {

    //container "dockerhub.ebi.ac.uk/pdbe/containers/chimerax-protonation:1.11.2"

    container "chimerax:1.0"

    input:
    tuple val(meta), path(bound_molecule_file)

    output:
    tuple val(meta), path("protonated_${meta.id}.cif"), emit: protonated_files

    script:
    def out_name = "protonated_${meta.id}.cif"
    """
    /ChimeraX/ChimeraX.app/bin/ChimeraX \
        --nogui \
        --cmd 'open "${bound_molecule_file}"; addh; save "${out_name}" format mmcif' \
        --silent
    """
}