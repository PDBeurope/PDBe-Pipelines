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

workflow TAXDB_BUILD {

    take:
    seq_db_ch

    main:
    def tax_url = params.tax_url
    def mapping_file = [params.mapping_file] ?: [null]
    def taxdump_dir = [params.taxdump_dir] ?: DOWNLOAD(Channel.value(tax_url))
    println "Mapping file" + mapping_file
    println "taxdump dir " + taxdump_dir
    //taxdump_ch = Channel.fromPath(taxdump_dir, checkIfExists: true)
    //mapping_ch = Channel.fromPath(mapping_file, checkIfExists: true)
    MMSEQS_CREATETAXDB(
        seq_db_ch,
        Channel.value([id: 'taxdump']).combine(taxdump_dir),
        Channel.value([id: 'mapping']).combine(mapping_file)
    )
}

workflow SEQENCEDB_BUILD {

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


workflow {
    SEQENCEDB_BUILD()
    TAXDB_BUILD(SEQENCEDB_BUILD.out.db)
}


