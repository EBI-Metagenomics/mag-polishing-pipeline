process SELECT_BRANCHWATER_HITS {

    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.12'

    input:
    tuple val(meta), path(hits_csv), path(metadata_csv)
    val n_runs

    output:
    tuple val(meta), path("${prefix}_selected_runs.csv"), emit: runs
    path "versions.yml"                                 , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_branchwater_hits.py \\
        --hits ${hits_csv} \\
        --metadata ${metadata_csv} \\
        --n-runs ${n_runs} \\
        --exclude ${meta.sample} \\
        --output ${prefix}_selected_runs.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "run_accession,fastq_1,fastq_2,md5_1,md5_2,containment" > ${prefix}_selected_runs.csv
    for i in \$(seq 1 ${n_runs}); do
        echo "SRRSTUB\${i},ftp://x/SRRSTUB\${i}_1.fastq.gz,ftp://x/SRRSTUB\${i}_2.fastq.gz,md5a,md5b,0.9" >> ${prefix}_selected_runs.csv
    done
    touch versions.yml
    """
}
