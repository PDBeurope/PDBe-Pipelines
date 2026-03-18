include { remove_alt_conformations } from './modules/remove_alt_conformations.nf'
include { gen_biological_assembly } from './modules/gen_biological_assembly.nf'
include { gen_bound_molecules } from './modules/gen_bound_molecules.nf'
include { run_protonation } from './modules/run_protonation.nf'
include { compute_interactions } from './modules/compute_interactions.nf'
include { run_fix_protonated_cif } from './modules/run_fix_protonated_cif.nf'
include { fix_assembly } from './modules/fix_assembly.nf'

workflow {

    main:
        input_ch = parse_manifest(params.manifest)
        input_ch.view()
        // RUN PREPROCESING WORKFLOW
        // run remove alternative conformations
        remove_alt_conformations(input_ch)
        remove_alt_conformations.out.no_alt_conf_cifs.view()

        // generate biological assembly
        // run_model_server(remove_alt_conformations.out.no_alt_conf_cifs)

        // generate bound molecules
        //gen_bound_molecules(remove_alt_conformations.out.no_alt_conf_cifs)

        // fix assemnbly
        fix_assembly(remove_alt_conformations.out.no_alt_conf_cifs)

        // protonate structures
        run_protonation(remove_alt_conformations.out.no_alt_conf_cifs)

        // fix protonated cif
        //run_fix_protonated_cif(run_protonation.out.protonated_files)

        // run compute interactions
        compute_interactions(fix_assembly.out.fixed_cifs)

    publish:
        intx_jsons = compute_interactions.out.interactions_jsons

}

output {
    intx_jsons {
       path {"interactions_jsons/"}
       mode 'copy'
        index {
            path "interactions_jsons/index.json"
       }
    }
}

def parse_manifest(mnf) {
    // Function to parse the manifest file and return a list of PDB files

    def mnf_rows = channel.fromPath(mnf).splitCsv(header: true, sep: ',')
            .map { row ->
            // set meta
                def meta = [
                      id: row.id,
                      // sample_id is explictily used on the
                      // publishing of files paths
                    ]
                def cif_path = row.cif_path.replace('projectDir', projectDir.toString())
                tuple(meta, cif_path)
            }

    return mnf_rows
}