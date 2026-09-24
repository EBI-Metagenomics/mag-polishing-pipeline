/*
 * `eukcc single` on one genome.
 *
 * One run serves two purposes (proposal.md section 5): completeness/contamination for the
 * comparison table, and `ncbi_lng` - a dash separated chain of NCBI taxids whose last
 * element is the assigned taxid - for the "same organism" match. Named EUKCC_SINGLE and
 * not EUKCC because genomes-generation ships its own EUKCC process whose config selector
 * would otherwise apply to this one.
 */
process EUKCC_SINGLE {

    container 'community.wave.seqera.io/library/python_metaeuk_pplacer_epa-ng_pruned:0b7ea587ebcad440'
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)
    path eukcc_db

    output:
    tuple val(meta), path("${meta.id}_eukcc.csv"), emit: eukcc_result

    script:
    """
    # When EukCC does not find any marker genes it exits with status code 201; allow that
    # exit code so the genome is reported with an NA taxid instead of killing the run
    eukcc single \\
        --out ${meta.id}_eukcc_results \\
        --threads ${task.cpus} \\
        --db ${eukcc_db} \\
        ${fasta} || [ \$? -eq 201 ]

    # eukcc.tsv is tab separated and looks like:
    #     fasta                      completeness  contamination  ncbi_lng
    #     /path/to/MGYG000000001.fa  95.24         1.19           1-131567-2759-33154-4751
    # comma separate it, replace the path with the genome id and set our own header
    if [ -s ${meta.id}_eukcc_results/eukcc.tsv ]; then
        awk '{gsub(".*/", "", \$1); \$1=\$1; OFS=","; print}' ${meta.id}_eukcc_results/eukcc.tsv |\\
            cut -d',' -f1,2,3,4 |\\
            awk -F',' 'NR==1 {print "genome,completeness,contamination,ncbi_lng"; next} {print "${meta.id}," \$2 "," \$3 "," (\$4 == "" ? "NA" : \$4)}' \\
            > ${meta.id}_eukcc.csv
    else
        printf 'genome,completeness,contamination,ncbi_lng\\n${meta.id},NA,NA,NA\\n' > ${meta.id}_eukcc.csv
    fi
    """

    stub:
    """
    printf 'genome,completeness,contamination,ncbi_lng\\n${meta.id},90.0,1.0,1-131567-2759-33154-4751\\n' > ${meta.id}_eukcc.csv
    """
}
