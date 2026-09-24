/*
 * MAGs leave GGP as .fa.gz and reference genomes arrive as .fa.gz / .fna.gz, while EukCC,
 * BUSCO and the assembly stats all want a plain fasta. Decompressing once up front also
 * makes `fasta.baseName` the single key every downstream join uses.
 */
process GUNZIP_GENOME {

    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.12'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: genome

    script:
    """
    case "${genome}" in
        *.gz) gunzip -c ${genome} > ${meta.id}.fa ;;
        *)    cp -L ${genome} genome.tmp && mv genome.tmp ${meta.id}.fa ;;
    esac
    """

    stub:
    """
    touch ${meta.id}.fa
    """
}
