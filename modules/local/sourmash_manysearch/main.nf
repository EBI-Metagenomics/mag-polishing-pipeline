/*
 * Adapted from EBI-Metagenomics/branchwater-nf (process MANY_SEARCH), rewritten for
 * Nextflow 24.04. The upstream process also stages the genome fasta, which its script
 * never uses; dropped here.
 *
 * manysearch is the sourmash branchwater plugin:
 * https://github.com/sourmash-bio/sourmash_plugin_branchwater
 *
 * --cores is ours, upstream does not pass it: the plugin defaults to a rayon pool sized
 * to the whole machine (128 threads were spawned inside a 4 cpu SLURM allocation).
 */
process SOURMASH_MANYSEARCH {

    tag "${meta.id}"

    container 'quay.io/microbiome-informatics/sourmash:4.8.14_branchwater-plugin_0.9.13'

    input:
    tuple val(meta), path(signature)
    path index

    output:
    tuple val(meta), path("${prefix}.csv"), emit: hits
    path "versions.yml"                   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash scripts manysearch \\
        ${signature} \\
        ${index} \\
        -o ${prefix}.csv \\
        -k ${params.branchwater_k} \\
        -s ${params.branchwater_scaled} \\
        -t ${params.branchwater_threshold} \\
        --cores ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/^sourmash //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "query_name,query_md5,match_name,containment,intersect_hashes,ksize,scaled,moltype,match_md5,jaccard,max_containment,query_containment_ani" > ${prefix}.csv
    touch versions.yml
    """
}
