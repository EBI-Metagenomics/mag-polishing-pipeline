process COLLECT_ASSEMBLIES {

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    label 'process_light'

    input:
    path fastas, stageAs: "input/*"

    output:
    path "assemblies", emit: assemblies_dir

    script:
    """
    mkdir assemblies
    cp -L input/* assemblies/
    """
    stub:
    """
    mkdir assemblies
    cp -L input/* assemblies/
    """
}
