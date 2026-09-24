/*
 * The sample's own reads plus the first N Branchwater hits, one co-assembly dataset.
 *
 * Plain `cat` of the gzip members: valid gzip, and metaSPAdes tolerates the duplicated
 * read names across runs. If GGP's back-mapping ever misbehaves, prefix the read names
 * with the run accession here (proposal.md section 9).
 */
process CONCAT_FASTQ {

    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.12'

    input:
    tuple val(meta), path(reads_1, stageAs: "r1/*"), path(reads_2, stageAs: "r2/*")

    output:
    tuple val(meta), path("${meta.id}_{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_provenance.tsv"), emit: provenance

    script:
    """
    cat ${reads_1} > ${meta.id}_1.fastq.gz
    cat ${reads_2} > ${meta.id}_2.fastq.gz

    printf 'dataset\\tn_concat_samples\\tsource\\n' > ${meta.id}_provenance.tsv
    for source in ${meta.sources.join(' ')}; do
        printf '${meta.id}\\t${meta.n_concat}\\t%s\\n' "\$source" >> ${meta.id}_provenance.tsv
    done
    """

    stub:
    """
    echo | gzip > ${meta.id}_1.fastq.gz
    echo | gzip > ${meta.id}_2.fastq.gz
    printf 'dataset\\tn_concat_samples\\tsource\\n${meta.id}\\t${meta.n_concat}\\tstub\\n' > ${meta.id}_provenance.tsv
    """
}
