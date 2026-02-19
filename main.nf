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
include { MMSEQS_CREATEINDEX } from './modules/nf-core/mmseqs/createindex/main'

// Generic Download
process DOWNLOAD {

    input:
    val url

    output:
    path "*"

    script:
    """
    fname=\$(basename $url)
    echo "Downloading $url"
    curl -L $url -o \$fname
    """
}

//Generic Extract
process EXTRACT {

    input:
    path archive

    output:
    path "*"

    script:
    """
    echo "Extracting $archive"
        case "$archive" in
      *.tar.gz|*.tgz)
        tar -xzf "$archive"
        ;;
      *.gz)
        gunzip -c "$archive" > \$(basename "$archive" .gz)
        ;;
      *)
        echo "Unsupported format: $archive"
        exit 1
        ;;
    esac
    """
}

workflow SEQENCE_DB_BUILD {

    //main:
    // Parameters with sensible defaults
    def url     = params.uniprot_url    
    def threads = params.create_db_threads
    def prefix  = params.output_db_name

    def fasta_gz = params.input_fasta ? Channel.fromPath(params.input_fasta) : DOWNLOAD(Channel.value(url))
    fasta_ch = EXTRACT(fasta_gz)

    createdb_ch = Channel.value([ id: prefix ]).combine(fasta_ch)
    //createdb_ch.view()

    MMSEQS_CREATEDB(createdb_ch)
    emit:
        db = MMSEQS_CREATEDB.out.db
}

workflow TAX_DB_BUILD {

    take:
    seq_db_ch

    main:
    //def tax_url = params.tax_url
    def mapping_file = params.mapping_file
    def taxdump_dir = params.taxdump_dir

    taxdump_ch = Channel.fromPath(params.taxdump_dir, checkIfExists: true)
    mapping_ch = Channel.fromPath(params.mapping_file, checkIfExists: true)
    MMSEQS_CREATETAXDB(
        seq_db_ch,
        Channel.value([id: 'taxdump']).combine(taxdump_ch),
        Channel.value([id: 'mapping']).combine(mapping_ch)
    )
    emit:
        bam = "bam"
    
}

workflow INDEX_BUILD {
    take:
    seq_db_ch
    tax_db_ch

    main:
    //tmp_path_ch = Channel.fromPath(params.idx_tmp_path)
    //seq_db_ch.view()
    MMSEQS_CREATEINDEX(seq_db_ch)

}

workflow {
    main:
    SEQENCE_DB_BUILD()
    TAX_DB_BUILD(SEQENCE_DB_BUILD.out.db)
    INDEX_BUILD(SEQENCE_DB_BUILD.out.db, TAX_DB_BUILD.out.bam)
    
}


