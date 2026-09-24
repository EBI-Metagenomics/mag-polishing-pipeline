/*
 * Cycle 1 target MAG -> Branchwater hits -> ENA fastqs -> one co-assembly dataset per
 * value of --n_concat_samples.
 *
 * Branchwater runs on one MAG per sample (the one carrying the target taxid), so the
 * datasets are "this sample's reads plus the reads of the N public runs that contain the
 * most of this organism".
 */

include { BRANCHWATER              } from '../branchwater'
include { SELECT_BRANCHWATER_HITS  } from '../../../modules/local/select_branchwater_hits'
include { ENA_FETCH_FASTQ          } from '../../../modules/local/ena_fetch_fastq'
include { CONCAT_FASTQ             } from '../../../modules/local/concat_fastq'

workflow BUILD_CONCAT_DATASETS {

    take:
    target_mags   // channel: [ val(meta), path(mag) ]        meta.sample links back to the sample
    sample_reads  // channel: [ val(sample), [ fastq_1, fastq_2 ] ]

    main:
    def n_values = params.n_concat_samples.toString().tokenize(',').collect { it.trim() as int }

    BRANCHWATER(
        target_mags,
        file(params.branchwater_index, checkIfExists: true),
        file(params.branchwater_metadata_db, checkIfExists: true)
    )

    SELECT_BRANCHWATER_HITS(
        BRANCHWATER.out.hits.join( BRANCHWATER.out.metadata ),
        n_values.max()
    )

    // ranked runs, kept as an ordered list per sample so "the first N" stays meaningful
    selected = SELECT_BRANCHWATER_HITS.out.runs.map { meta, csv ->
        [ meta.sample, csv.splitCsv(header: true) ]
    }

    // every accession is downloaded once (storeDir), whatever the N values ask for
    ENA_FETCH_FASTQ(
        selected
            .flatMap { _sample, rows ->
                rows.collect { [ it.run_accession, it.fastq_1, it.fastq_2, it.md5_1, it.md5_2 ] }
            }
            .unique()
    )

    fetched = ENA_FETCH_FASTQ.out.reads.toList()

    CONCAT_FASTQ(
        selected
            .join( sample_reads )
            .combine( fetched.map { [it] } )
            .flatMap { sample, rows, reads, downloads ->
                def by_accession = downloads.collectEntries { accession, r1, r2 -> [accession, [r1, r2]] }

                def selected_branchwater_hits = rows.collect { it.run_accession }
                    .findAll { it != sample && by_accession[it] }

                if ( !selected_branchwater_hits ) {
                    log.warn "${sample}: no usable Branchwater hit - the only match is the sample's own run, or nothing was downloadable from ENA. No co-assembly, so no cycle 2 for this sample."
                    return []
                }

                def depths = concat_depths( n_values, selected_branchwater_hits.size() )
                if ( depths != n_values.unique().sort() ) {
                    log.warn "${sample}: --n_concat_samples asks for ${n_values.join(',')} but only ${selected_branchwater_hits.size()} hit(s) are usable; building ${depths.join(',')} instead."
                }

                depths.collect { n ->
                    def accessions = selected_branchwater_hits.take(n)
                    def meta = [
                        id            : "${sample}_n${n}",
                        sample        : sample,
                        n_concat      : n,
                        sources       : [sample] + accessions,
                        single_end    : false,
                        assembler     : params.assembler,
                    ]
                    [
                        meta,
                        [reads[0]] + accessions.collect { by_accession[it][0] },
                        [reads[1]] + accessions.collect { by_accession[it][1] },
                    ]
                }
            }
    )

    emit:
    reads      = CONCAT_FASTQ.out.reads        // [ val(meta), [ fastq_1, fastq_2 ] ]
    provenance = CONCAT_FASTQ.out.provenance
    hits       = BRANCHWATER.out.hits
    versions   = BRANCHWATER.out.versions.mix( SELECT_BRANCHWATER_HITS.out.versions )
}

/*
 * The co-assembly depths to build given what Branchwater actually delivered: every
 * requested depth that fits, plus one capped dataset when the deepest request does not -
 * capped and deduplicated, so a run with 3 usable hits and `--n_concat_samples 1,5,10`
 * builds n1 and n3 rather than three datasets of which two are identical.
 */
def concat_depths( requested, usable ) {
    if ( usable < 1 ) {
        return []
    }
    def fits = requested.findAll { it <= usable }
    return ( fits + [ Math.min( requested.max(), usable ) ] ).unique().sort()
}
