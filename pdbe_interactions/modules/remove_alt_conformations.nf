process remove_alt_conformations {

    container 'docker://community.wave.seqera.io/library/pip_pdbeccdutils:de9d3afe23c1656a'

    input:
    tuple val(meta), file(mmcif_file)

    output:
    tuple val(meta), file("cleaned_${mmcif_file}"), emit: cleaned_pdb_files

    script:
    """
    # Command to remove alternative conformations from the PDB file
    # For example, using a hypothetical tool 'remove_alt_conf' 
    
    #ccdutils ${mmcif_file} > cleaned_${mmcif_file}
    mv ${mmcif_file} cleaned_${mmcif_file}
    """
}
