/*
 * Taxonomy + QC for every genome in the experiment, from one EukCC run.
 *
 * The reference genome and the MAGs of both cycles all go through this subworkflow, so
 * both sides of the comparison are classified by identical code. Two genomes are "the
 * same organism" when they carry the same NCBI taxid, and that taxid is the last element
 * of EukCC's `ncbi_lng` chain. Genomes EukCC cannot place get taxid `NA`, which matches
 * nothing (proposal.md section 5).
 *
 * The prokaryotic counterpart (GTDB-Tk) would plug in here, see proposal.md section 6.
 */

include { EUKCC_SINGLE } from '../../../modules/local/eukcc_single'

workflow ASSIGN_TAXONOMY {

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
