/*
 * Adapted from EBI-Metagenomics/branchwater-nf (process GET_METADATA), rewritten for
 * Nextflow 24.04.
 *
 * Joins the manysearch hits with the Branchwater metadata DuckDB. The hits CSV carries
 * containment / query_containment_ani but no sample metadata; the metadata carries
 * `acc` and `librarylayout` but no containment. Downstream filtering needs both, so the
 * hits CSV is emitted alongside the join.
 */
process BRANCHWATER_METADATA {

    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/duckdb-cli:1.5.5--5ae76bf8b4435b22' :
        'community.wave.seqera.io/library/duckdb-cli:1.5.5--76cefe0dde3b99d5' }"

    input:
    tuple val(meta), path(hits_csv)
    path metadata_db

    output:
    tuple val(meta), path("${prefix}_hits_metadata.csv"), emit: metadata
    path "versions.yml"                                 , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    duckdb ${metadata_db} -readonly -c "
    COPY (
      SELECT metadata.*
      FROM metadata
      INNER JOIN read_csv_auto('${hits_csv}') hits
        ON metadata.acc = hits.match_name
    ) TO '${prefix}_hits_metadata.csv' WITH (FORMAT csv, HEADER);
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duckdb: \$(duckdb --version | sed 's/ .*//; s/^v//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hits_metadata.csv
    touch versions.yml
    """
}
