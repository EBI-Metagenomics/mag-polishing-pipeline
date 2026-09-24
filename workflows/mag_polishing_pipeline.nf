/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP_GENOME as GUNZIP_REFERENCE     } from '../modules/local/gunzip_genome'
include { GUNZIP_GENOME as GUNZIP_CYCLE1        } from '../modules/local/gunzip_genome'
include { GUNZIP_GENOME as GUNZIP_CYCLE2        } from '../modules/local/gunzip_genome'

include { ASSIGN_TAXONOMY as TAXONOMY_REFERENCE } from '../subworkflows/local/assign_taxonomy'
include { ASSIGN_TAXONOMY as TAXONOMY_CYCLE1    } from '../subworkflows/local/assign_taxonomy'
include { ASSIGN_TAXONOMY as TAXONOMY_CYCLE2    } from '../subworkflows/local/assign_taxonomy'

include { ASSEMBLE_AND_BIN as CYCLE1            } from '../subworkflows/local/assemble_and_bin'
include { ASSEMBLE_AND_BIN as CYCLE2            } from '../subworkflows/local/assemble_and_bin'

include { BUILD_CONCAT_DATASETS                 } from '../subworkflows/local/build_concat_datasets'
include { COMPARE                               } from '../subworkflows/local/compare'

include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// constant batch tags, >= 7 characters because miassembler substrings them for its paths
def CYCLE1_TAG = "MAGCYC1"
def CYCLE2_TAG = "MAGCYC2"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MPP {

    take:
    ch_samplesheet // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ], path(genome) ]

    main:

    ch_versions = Channel.empty()

    sample_ids = ch_samplesheet.map { meta, _reads, _genome -> meta.id }.collect()

    /*
     * The reference genome: the target taxid, and the run's fail-fast input validation.
     */
    GUNZIP_REFERENCE(
        ch_samplesheet.map { meta, _reads, genome ->
            def sample = meta.id
            [ [id: genome_name(genome), sample: sample, origin: "reference", slot: "reference"], genome ]
        }
    )
    TAXONOMY_REFERENCE( GUNZIP_REFERENCE.out.genome )

    target_taxid = TAXONOMY_REFERENCE.out.taxonomy.map { meta, row ->
        def lineage = row.ncbi_lng?.trim()
        if ( !lineage || lineage == "NA" || !lineage.tokenize('-').contains("2759") ) {
            error "Reference genome ${meta.id} (sample ${meta.sample}) is not eukaryotic: EukCC lineage '${lineage}'. This pipeline is eukaryote only, see proposal.md section 1."
        }
        [ meta.sample, row.taxid ]
    }

    /*
     * Cycle 1 - the sample on its own, assembled and binned, and the MAG carrying the
     * target taxid is what Branchwater then searches with.
     *
     * --skip_first_assembly drops that whole cycle and searches with the reference genome
     * from the samplesheet instead. The run is then reference + cycle 2 only: no cycle-1
     * assembly, no cycle-1 GGP, no `cycle1` slot in any table. Use it when the sample has
     * already been assembled elsewhere, or when the question is only what the public runs
     * add to a genome you already have.
     */
    if ( params.skip_first_assembly ) {
        cycle1_taxonomy = Channel.empty()
        cycle1_target   = Channel.empty()
        cycle1_versions = Channel.empty()

        // the samplesheet genome as it is, not the GUNZIP copy: sourmash reads gzip natively
        branchwater_query = ch_samplesheet.map { meta, _reads, genome ->
            def sample = meta.id
            [ [id: genome_name(genome), sample: sample], genome ]
        }
    }
    else {
        CYCLE1(
            ch_samplesheet.map { meta, reads, _genome ->
            def sample = meta.id [ assembler_meta(sample, CYCLE1_TAG), reads ] }
        )

        cycle1_genomes = CYCLE1.out.bins.combine( sample_ids.map { [it] } ).map { bin, ids ->
            [ [id: genome_name(bin), sample: owner_of(bin.name, ids), origin: "cycle1", slot: "cycle1"], bin ]
        }

        GUNZIP_CYCLE1( cycle1_genomes )
        TAXONOMY_CYCLE1( GUNZIP_CYCLE1.out.genome )

        cycle1_taxonomy = TAXONOMY_CYCLE1.out.taxonomy
        cycle1_target   = pick_target( cycle1_taxonomy, target_taxid, GUNZIP_CYCLE1.out.genome )
        cycle1_versions = CYCLE1.out.versions

        // the compressed bin, not the copy GUNZIP made: sourmash reads gzip natively
        branchwater_query = cycle1_target.map { meta, _row, _fasta -> [meta.id, meta] }
            .join( cycle1_genomes.map { meta, bin -> [meta.id, bin] } )
            .map { _id, meta, bin -> [meta, bin] }
    }

    /*
     * Branchwater -> ENA -> one co-assembly dataset per N.
     */
    BUILD_CONCAT_DATASETS(
        branchwater_query,
        ch_samplesheet.map { meta, reads, _genome ->
            def sample = meta.id [sample, reads] }
    )

    /*
     * Cycle 2 - the same sample co-assembled with the hits.
     */
    CYCLE2(
        BUILD_CONCAT_DATASETS.out.reads.map { meta, reads ->
            [ assembler_meta(meta.id, CYCLE2_TAG) + [sample: meta.sample, n_concat: meta.n_concat], reads ]
        }
    )

    dataset_ids = BUILD_CONCAT_DATASETS.out.reads.map { meta, _reads -> meta.id }.collect()

    GUNZIP_CYCLE2(
        CYCLE2.out.bins.combine( dataset_ids.map { [it] } ).map { bin, ids ->
            def dataset = owner_of(bin.name, ids)
            [
                [
                    id    : genome_name(bin),
                    sample: dataset.replaceFirst(/_n\d+$/, ''),
                    origin: "cycle2",
                    slot  : "n" + (dataset =~ /_n(\d+)$/)[0][1]
                ],
                bin
            ]
        }
    )
    TAXONOMY_CYCLE2( GUNZIP_CYCLE2.out.genome )

    cycle2_target = pick_target( TAXONOMY_CYCLE2.out.taxonomy, target_taxid, GUNZIP_CYCLE2.out.genome, true )

    /*
     * Tables.
     */
    all_taxonomy = TAXONOMY_REFERENCE.out.taxonomy
        .mix( cycle1_taxonomy, TAXONOMY_CYCLE2.out.taxonomy )

    all_taxonomy
        .map { meta, row ->
            [meta.id, meta.sample, meta.origin, row.taxid, row.completeness, row.contamination, row.ncbi_lng].join('\t') + '\n'
        }
        .collectFile(
            name: "eukcc_taxonomy.tsv",
            storeDir: "${params.outdir}/taxonomy",
            sort: true,
            seed: "genome\tsample\torigin\ttaxid\tcompleteness\tcontamination\tncbi_lng\n"
        )

    cycle1_taxonomy.mix( TAXONOMY_CYCLE2.out.taxonomy )
        .map { meta, row ->
            [meta.id, meta.sample, meta.origin, meta.slot, row.taxid, row.completeness, row.contamination].join('\t') + '\n'
        }
        .collectFile(
            name: "all_mags.tsv",
            storeDir: "${params.outdir}/mags",
            sort: true,
            seed: "genome\tsample\torigin\tslot\ttaxid\tcompleteness\tcontamination\n"
        )

    // one row per expected slot, NA when the cycle recovered no MAG with the target taxid
    found_targets = cycle1_target.mix( cycle2_target )
        .map { meta, row, _fasta ->
            [ "${meta.sample}\t${meta.slot}".toString(), [meta.id, row.taxid, row.completeness, row.contamination] ]
        }
        .toList()

    /*
     * One row per slot that could be built: cycle 1 for every sample, plus the
     * co-assembly datasets BUILD_CONCAT_DATASETS actually assembled - a sample whose
     * only Branchwater hit is its own run has no cycle-2 slot at all, and one with fewer
     * usable hits than --n_concat_samples asks for has its depths capped (it warns).
     */
    cycle1_slots = params.skip_first_assembly
        ? Channel.empty()
        : ch_samplesheet.map { meta, _reads, _genome -> "${meta.id}\tcycle1".toString() }

    expected_slots = cycle1_slots
        .mix(
            BUILD_CONCAT_DATASETS.out.reads.map { meta, _reads ->
                "${meta.sample}\tn${meta.n_concat}".toString()
            }
        )

    expected_slots
        .combine( found_targets.map { [it] } )
        .map { key, found ->
            def values = found.collectEntries()[key] ?: ["NA", "NA", "NA", "NA"]
            ([key] + values).join('\t') + '\n'
        }
        .collectFile(
            name: "target_mags.tsv",
            storeDir: "${params.outdir}/mags",
            sort: true,
            seed: "sample\tslot\tgenome\ttaxid\tcompleteness\tcontamination\n"
        )

    /*
     * The comparison set: reference, cycle 1, one cycle 2 MAG per N, in that order.
     */
    comparison_set = GUNZIP_REFERENCE.out.genome.map { meta, fasta -> [meta.sample, -1, fasta] }
        .mix(
            cycle1_target.map { meta, _row, fasta -> [meta.sample, 0, fasta] },
            cycle2_target.map { meta, _row, fasta -> [meta.sample, meta.slot.substring(1) as int, fasta] }
        )
        .map { sample, rank, fasta -> [sample, [rank, fasta]] }
        .groupTuple()
        .map { sample, entries -> [ sample, entries.sort { it[0] }.collect { it[1] } ] }

    COMPARE(
        comparison_set,
        all_taxonomy.map { meta, row -> [meta.id, row] }
    )

    //
    // Collate and save software versions
    //
    ch_versions = ch_versions.mix(
        cycle1_versions,
        CYCLE2.out.versions,
        BUILD_CONCAT_DATASETS.out.versions,
    )

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'mag_polishing_pipeline_software_versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    comparison     = COMPARE.out.metrics   // channel: path(assembly_qc_metrics.tsv)
    versions       = ch_collated_versions  // channel: path(versions.yml)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * The meta SHORT_READS_ASSEMBLER and GGP need. `study_accession` is the cycle tag, which
 * is what separates the two cycles in miassembler's publish paths.
 */
def assembler_meta( id, study_accession ) {
    [
        id                    : id,
        study_accession       : study_accession,
        single_end            : false,
        assembler             : params.assembler,
        human_reference       : params.human_reference,
        phix_reference        : params.phix_reference,
        contaminant_reference : params.contaminant_reference,
        lambdaphage_reference : null,
    ]
}

// "SAMPLE_001.fna.gz" -> "SAMPLE_001", "SRR123456_concoct_3.fa.gz" -> "SRR123456_concoct_3"
def genome_name( genome ) {
    genome.name.replaceFirst(/\.(fa|fna|fasta)(\.gz)?$/, '')
}

// bins are named <dataset id>_<binner>_<n>.fa.gz; longest match wins so that SAMPLE_1
// never steals a bin belonging to SAMPLE_11
def owner_of( bin_name, ids ) {
    def owner = ids.findAll { bin_name.startsWith("${it}_") }.max { it.size() }
    if ( !owner ) {
        error "Cannot tell which dataset produced ${bin_name}; known datasets: ${ids.join(', ')}"
    }
    return owner
}

def as_number( value ) {
    try {
        return value as Double
    }
    catch ( Exception ignored ) {
        return -1d
    }
}

/*
 * The MAGs carrying the sample's target taxid, best completeness first, one per slot
 * (cycle 1 has one slot per sample, cycle 2 one per sample and N).
 */
def pick_target( taxonomy, target_taxid, genomes, per_slot = false ) {
    return taxonomy
        .map { meta, row -> [meta.sample, meta, row] }
        .combine( target_taxid, by: 0 )
        .filter { _sample, _meta, row, taxid -> row.taxid != "NA" && row.taxid == taxid }
        .map { sample, meta, row, _taxid -> [ per_slot ? "${sample}\t${meta.slot}".toString() : sample, [meta, row] ] }
        .groupTuple()
        .map { _key, candidates -> candidates.max { as_number(it[1].completeness) } }
        .map { meta, row -> [meta.id, meta, row] }
        .join( genomes.map { meta, fasta -> [meta.id, fasta] } )
        .map { _id, meta, row, fasta -> [meta, row, fasta] }
}
