//
// Subworkflow with functionality specific to the EBI-Metagenomics/mag-polishing-pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, genome ->
            [ meta + [ single_end: false ], [ fastq_1, fastq_2 ], genome ]
        }
        .set { ch_samplesheet }

    validateInputSamplesheet( ch_samplesheet )

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    // what miassembler's SHORT_READS_ASSEMBLER accepts for paired-end reads.
    // Declared inside the function on purpose: a script-level `def` in a Nextflow file is
    // a local of the implicit run() method, so functions cannot see it.
    def ASSEMBLERS = ["metaspades", "megahit"]

    /*
     * miassembler picks the assembler from `meta.assembler` and errors with a message
     * about the read layout when it does not recognise the value, so check it here where
     * the message can name the typo. metaSPAdes is the default; megahit is the escape
     * hatch for cycle-2 datasets too deep to assemble in the available memory.
     */
    if ( !ASSEMBLERS.contains(params.assembler) ) {
        error("--assembler must be one of ${ASSEMBLERS.join(', ')}, got '${params.assembler}'")
    }

    /*
     * GGP's uploader submits genomes to ENA for real. It is off by default and this
     * pipeline has no reason to turn it on, but if someone does, fail here rather than
     * after two assembly cycles: `genome_upload -u` takes
     * --ena_assembly_study_accession verbatim and ours defaults to a filename prefix,
     * and the uploader also refuses without --metagenome.
     */
    if ( params.upload_mags || params.upload_bins ) {
        if ( !(params.ena_assembly_study_accession ==~ /^(PRJ|[EDS]RP)\w+$/) ) {
            error("--upload_mags/--upload_bins need --ena_assembly_study_accession to be a registered ENA study accession, got '${params.ena_assembly_study_accession}'")
        }
        if ( !params.metagenome?.trim() ) {
            error("--upload_mags/--upload_bins need --metagenome (e.g. 'soil metagenome')")
        }
        log.warn "GGP will SUBMIT genomes to ENA${params.test_upload ? ' (test service)' : ' (live)'} as study ${params.ena_assembly_study_accession}."
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(ch_samplesheet) {
    // every metric in the comparison is keyed back to the genome basename, so a repeated
    // reference genome file name would silently merge two samples' rows
    ch_samplesheet
        .map { _meta, _reads, genome -> genome.name }
        .toList()
        .subscribe { names ->
            def dups = names.countBy { it }.findAll { _name, count -> count > 1 }.keySet()
            if ( dups ) {
                error("Please check input samplesheet -> reference genome file names must be unique across the run, duplicated: ${dups.join(', ')}")
            }
        }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "miassembler (https://github.com/EBI-Metagenomics/miassembler),",
            "metaSPAdes (Nurk et al. 2017)" + (params.assembler == "megahit" ? " or MEGAHIT (Li et al. 2015)," : ","),
            "genomes-generation (https://github.com/EBI-Metagenomics/genomes-generation),",
            "sourmash branchwater (Irber et al. 2022),",
            "EukCC (Saary et al. 2020),",
            "BUSCO (Manni et al. 2021)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome Research, 27(5), 824-834. doi: 10.1101/gr.213959.116</li>",
            "<li>Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674-1676. doi: 10.1093/bioinformatics/btv033</li>",
            "<li>Irber, L., et al. (2022). Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. bioRxiv. doi: 10.1101/2022.01.11.475838</li>",
            "<li>Saary, P., Mitchell, A. L., & Finn, R. D. (2020). Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC. Genome Biology, 21(1), 244. doi: 10.1186/s13059-020-02155-4</li>",
            "<li>Manni, M., Berkeley, M. R., Seppey, M., Simao, F. A., & Zdobnov, E. M. (2021). BUSCO update. Molecular Biology and Evolution, 38(10), 4647-4654. doi: 10.1093/molbev/msab199</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

