/*
 * Adapted from EBI-Metagenomics/branchwater-nf (process SKETCH), rewritten for
 * Nextflow 24.04: publishing moved to conf/modules.config and versions.yml added.
 *
 * Takes the genome as it comes out of GGP, i.e. `.fa.gz` / `.fna.gz`; `sourmash sketch`
 * reads gzip natively, so there is no reason to sketch the decompressed copy that EukCC
 * and BUSCO need.
 */
process SOURMASH_SKETCH {

    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/microbiome-informatics/sourmash:4.8.14_branchwater-plugin_0.9.13'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${prefix}.sig"), emit: signature
    path "versions.yml"                   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash sketch dna \\
        -p k=${params.branchwater_k},scaled=${params.branchwater_scaled} \\
        -o ${prefix}.sig \\
        ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/^sourmash //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sig
    touch versions.yml
    """
}
