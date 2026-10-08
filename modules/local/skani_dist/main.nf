/*
 * `skani dist` of every MAG of one assembly (cycle 1, or one cycle-2 N) against the
 * sample's reference genome.
 *
 * This is how each target MAG is chosen: the bin that aligns to the reference, not
 * the one that shares its EukCC taxid. EukCC often stops at genus, and one sample can
 * carry several populations of the same genus - SRR26991367 has three under taxid 70447,
 * and picking by taxid + completeness swapped one for another in cycle 2 n1.
 *
 * skani only reports pairs above ~80% ANI, so a candidate with no row did not align at
 * all. `Align_fraction_ref` is the fraction of the reference the candidate covers.
 */
process SKANI_DIST {

    tag "${meta.id}"

    container 'quay.io/biocontainers/skani:0.3.2--h79ce301_0'

    input:
    tuple val(meta), path(reference), path(candidates, stageAs: 'candidates/*')

    output:
    tuple val(meta), path("${prefix}_skani.tsv"), emit: hits
    path "versions.yml"                         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    skani dist \\
        -t ${task.cpus} \\
        -r ${reference} \\
        -q ${candidates} \\
        -o ${prefix}_skani.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^skani //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n' > ${prefix}_skani.tsv
    for query in ${candidates}; do
        printf '%s\\t%s\\t99.00\\t80.00\\t60.00\\tref\\tquery\\n' ${reference} "\$query" >> ${prefix}_skani.tsv
    done
    touch versions.yml
    """
}
