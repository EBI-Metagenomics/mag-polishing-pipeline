/*
 * One run accession -> its two fastqs, md5 verified.
 *
 * storeDir is the whole caching story: a run picked for several samples, or for several
 * --n_concat_samples values, is downloaded once (proposal.md section 9, download volume).
 */
process ENA_FETCH_FASTQ {

    tag "${run_accession}"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.12'

    storeDir "${params.ena_cache_dir}"

    input:
    tuple val(run_accession), val(fastq_1), val(fastq_2), val(md5_1), val(md5_2)

    output:
    tuple val(run_accession), path("${run_accession}_1.fastq.gz"), path("${run_accession}_2.fastq.gz"), emit: reads

    script:
    """
    ena_fetch_fastq.py --url ${fastq_1} --md5 ${md5_1} --output ${run_accession}_1.fastq.gz
    ena_fetch_fastq.py --url ${fastq_2} --md5 ${md5_2} --output ${run_accession}_2.fastq.gz
    """

    stub:
    """
    echo | gzip > ${run_accession}_1.fastq.gz
    echo | gzip > ${run_accession}_2.fastq.gz
    """
}
