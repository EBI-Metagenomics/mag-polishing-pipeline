/*
 * Local reimplementation of EBI-Metagenomics/branchwater-nf.
 *
 * The upstream pipeline is 3 processes but uses the workflow output DSL
 * (`publish:` + `output {}`), which needs Nextflow >= 25.04, while miassembler and GGP
 * pin 24.04.
 */

include { SOURMASH_SKETCH      } from '../../../modules/local/sourmash_sketch'
include { SOURMASH_MANYSEARCH  } from '../../../modules/local/sourmash_manysearch'
include { BRANCHWATER_METADATA } from '../../../modules/local/branchwater_metadata'

workflow BRANCHWATER {

    take:
    genomes      // channel: [ val(meta), path(fasta) ]
    index        // path: sourmash rocksdb index
    metadata_db  // path: branchwater metadata duckdb

    main:
    ch_versions = Channel.empty()

    SOURMASH_SKETCH( genomes )
    ch_versions = ch_versions.mix( SOURMASH_SKETCH.out.versions.first() )

    SOURMASH_MANYSEARCH( SOURMASH_SKETCH.out.signature, index )
    ch_versions = ch_versions.mix( SOURMASH_MANYSEARCH.out.versions.first() )

    BRANCHWATER_METADATA( SOURMASH_MANYSEARCH.out.hits, metadata_db )
    ch_versions = ch_versions.mix( BRANCHWATER_METADATA.out.versions.first() )

    emit:
    signatures = SOURMASH_SKETCH.out.signature       // [ val(meta), path(sig) ]
    hits       = SOURMASH_MANYSEARCH.out.hits        // [ val(meta), path(csv) ]  containment + cANI
    metadata   = BRANCHWATER_METADATA.out.metadata   // [ val(meta), path(csv) ]  acc + librarylayout
    versions   = ch_versions
}
