process CALCULATE_ASSEMBLY_STATS {

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    label 'process_light'

    input:
    path genomes_fnas, stageAs: "new_genomes_dir"

    output:
    path "new_genome_stats.tsv", emit: stats_file

    script:
    """
    # Compute length, N50, GC content
    precompute_assembly_stats.py -i new_genomes_dir -o new_genome_stats.tsv
    """
    stub:
    """
    printf 'Genome\\tLength\\tN50\\tGC_content\\tN_contigs\\n' > new_genome_stats.tsv
    for genome in new_genomes_dir/*; do
        printf '%s\\t1000\\t500\\t50.0\\t10\\n' "\$(basename \${genome%.*})" >> new_genome_stats.tsv
    done
    """
}
