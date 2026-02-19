/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include {STANDARDCIF} from "../modules/local/standardcif"
include {DENSITY_FITNESS} from "../modules/local/densityfitness"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STANDARDIZECIF {

    take:
    ch_samplesheet // channel: samplesheet read in from --input


    main:
    ch_input = channel.fromPath(ch_samplesheet)
                      .splitCsv(header: true)
                      .view { row -> "CSV row: $row" }   // <-- prints each parsed row

    cif_in = ch_input.map { row -> row.cif_path }.view { p -> "cif_in: $p" }
    mtz_in = ch_input.map { row -> row.mtz_path }.view { p -> "mtz_in: $p" }
    

    STANDARDCIF(cif_in)
    // DENSITY_FITNESS(STANDARDCIF.out.cif_out, mtz_in)

    // emit:
    // DENSITY_FITNESS.out

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
