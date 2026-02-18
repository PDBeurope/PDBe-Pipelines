#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PDBEurope/pdbe_sifts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/PDBEurope/pdbe_sifts
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MMSEQS_CREATEDB } from './modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CREATETAXDB } from './modules/nf-core/mmseqs/createtaxdb/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// (Removed unused PDBEUROPE_PDBE_SIFTS workflow)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_UNIPROT {
    tag "download_uniprot"
    label 'process_low'

    input:
    val url

    output:
    path "uniprot_sprot.fasta.gz"

    script:
    """
    set -euo pipefail
    curl -L --fail -o uniprot_sprot.fasta.gz "${url}"
    """
}

process EXTRACT_FASTA {
    tag "extract_uniprot"
    label 'process_low'

    input:
    path gz

    output:
    path "uniprot.fasta"

    script:
    """
    set -euo pipefail
    gzip -c -d ${gz} > uniprot.fasta
    """
}

workflow UNIPROT_DB_BUILD {

    main:
    // Parameters with sensible defaults
    def url     = params.uniprot_url
    def threads = params.create_db_threads
    def prefix  = params.output_db_name  

    // Create channels
    Channel.of(url).set { ch_url }

    // Download and extract
    gz_ch     = DOWNLOAD_UNIPROT(ch_url)
    fasta_ch  = EXTRACT_FASTA(gz_ch)

    // Prepare meta and pair with FASTA for the module
    meta_ch = Channel.value([ id: prefix ])
    createdb_in = meta_ch.combine(fasta_ch)

    MMSEQS_CREATEDB(createdb_in)
    
    /*
    MMSEQS_CREATETAXDB(
        MMSEQS_CREATEDB.out.db,
        Channel.value([ id: 'tmp' ]).combine(Channel.fromPath(params.taxdump_dir)),
        //Channel.value([ id: 'taxdump' ]).combine(Channel.fromPath(params.taxdump_dir)),
        Channel.value([ id: 'mapping' ]).combine(Channel.fromPath(params.tax_mapping_file))
    )
    */
}


workflow {
    main:
    UNIPROT_DB_BUILD()
}


