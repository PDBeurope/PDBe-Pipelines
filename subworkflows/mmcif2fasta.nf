#!/usr/bin/env nextflow


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.mmcif_file   = params.mmcif_file        // required (file, dir, or glob)
params.merged_fasta = params.merged_fasta ?: "merged.fasta"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MMCIF_TO_FASTA {

    tag { mmcif.simpleName }

    input:
    path mmcif

    output:
    path "${mmcif.simpleName}.fasta"

    script:
    """
    python -m pdbe_sifts.sifts_mmcif2fasta \
        --cif ${mmcif} \
        --out ${mmcif.simpleName}.fasta
    """
}



process MERGE_FASTA {

    tag "merge_fasta"

    input:
    path fasta_files

    output:
    path params.merged_fasta

    script:
    """
    cat ${fasta_files.join(' ')} > ${params.merged_fasta}
    """
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_MMCIF {

    main:

    if( !params.mmcif_file ) {
        error "--mmcif_file is required"
    }

    mmcif_ch = Channel.fromPath(params.mmcif_file, checkIfExists: true)

    fasta_ch = MMCIF_TO_FASTA(mmcif_ch)

    MERGE_FASTA(fasta_ch.collect())
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ENTRYPOINT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    PROCESS_MMCIF()
}