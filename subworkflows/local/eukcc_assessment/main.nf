/*
 * Taxonomy + QC for every genome in the experiment, from one EukCC run.
 */

include { EUKCC_SINGLE } from '../../../modules/local/eukcc_single'

workflow EUKCC_ASSESSMENT {

    take:
    genomes  // channel: [ val(meta), path(fasta) ]  uncompressed, meta.id == fasta.baseName

    main:
    EUKCC_SINGLE( genomes, file(params.eukcc_db, checkIfExists: true) )

    taxonomy = EUKCC_SINGLE.out.eukcc_result
        .splitCsv(header: true, elem: 1)
        .map { meta, row ->
            [ meta, row + [taxid: taxid_of(row.ncbi_lng)] ]
        }

    emit:
    taxonomy = taxonomy                           // [ val(meta), [genome, completeness, contamination, ncbi_lng, taxid] ]
    eukcc    = EUKCC_SINGLE.out.eukcc_result      // [ val(meta), path(csv) ]
}

// "1-131567-2759-33154-4751-451864-5204" -> "5204"
def taxid_of( ncbi_lng ) {
    def lineage = ncbi_lng?.trim()
    if ( !lineage || lineage == "NA" ) {
        return "NA"
    }
    return lineage.tokenize('-').last()
}
