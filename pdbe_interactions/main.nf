include { remove_alt_conformations } from './modules/remove_alt_conformations.nf'
include { gen_biological_assembly } from './modules/gen_biological_assembly.nf'
include { gen_bound_molecules } from './modules/gen_bound_molecules.nf'
include { run_protonation } from './modules/run_protonation.nf'
workflow {

    main:
        input_ch = parse_manifest(params.manifest)
        input_ch.view()
        // RUN PREPROCESING WORKFLOW
        // run remove alternative conformations
        remove_alt_conformations(input_ch)
        remove_alt_conformations.out.no_alt_conf_cifs.view()
        // generate biological assembly
        //gen_biological_assembly(remove_alt_conformations.out.cleaned_pdb_files)

        // gen biomolecule
        //gen_biological_assembly.out.biological_assembly_files
        gen_bound_molecules(remove_alt_conformations.out.no_alt_conf_cifs)
    
        // protonate structures
        //run_protonation(gen_biomolecule.out.biomolecule_jsons)

        // run chimerax to generate interactions
        // run_chimerax(run_protonation.out.protonated_files)

        //
    
    //publish:
    //    protonated_files = run_protonation.out.protonated_files

}

//output {
//    protonated_files {
//        path "protonated_structures"
//       mode 'copy'
//        index {
//            path "protonated_structures/index.json"
//       }
//    }
//}

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
                tuple(meta, row.cif_path)
            }

    return mnf_rows
}