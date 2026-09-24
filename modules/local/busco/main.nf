process BUSCO {

    tag "${fasta.baseName}"

    container 'quay.io/biocontainers/busco:5.8.0--pyhdfd78af_0'
    
    beforeScript "rm -rf *.fa* || true"

    input:
    path fasta
    path busco_db
    val busco_mode

    output:
    path "short_summary.specific_${fasta.baseName}.txt", emit: busco_summary
    path "${fasta.baseName}", emit: busco_folder

    script:
    """
    if [ -d "${fasta.baseName}" ]; then
        rm -rf "${fasta.baseName}"
    fi

    busco  --offline \
            -i ${fasta} \
            -m '${busco_mode}' \
            -o ${fasta.baseName} \
            --auto-lineage-euk \
            --download_path ${busco_db} \
            -c ${task.cpus}

    #   parse and output genomes id and busco scores as csv
    result_file=\$(ls ${fasta.baseName}/short_summary.specific*.${fasta.baseName}.txt | head -n 1)

    if [ -f "\${result_file}" ]; then
        result=\$(grep 'C:' "\${result_file}" | head -n 1 | sed 's?\t??g')
        echo "${fasta.name}\t\${result}" > "short_summary.specific_${fasta.baseName}.txt"
    else
        echo "No result file found starting short_summary.specific..."
        exit 1
    fi
    """

    stub:
    """
    mkdir ${fasta.baseName}
    echo "${fasta.name}\tC:90.0%[S:89.0%,D:1.0%],F:2.0%,M:8.0%,n:255" > short_summary.specific_${fasta.baseName}.txt
    """
}
