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

include { SKANI_DIST as SKANI_CYCLE1           } from '../modules/local/skani_dist'
include { SKANI_DIST as SKANI_CYCLE2           } from '../modules/local/skani_dist'

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
     * The reference genome: what every target MAG is aligned against, and the run's
     * fail-fast input validation.
     */
    GUNZIP_REFERENCE(
        ch_samplesheet.map { meta, _reads, genome ->
            [ [id: genome_name(genome), sample: meta.id, origin: "reference", source: "reference"], genome ]
        }
    )
    TAXONOMY_REFERENCE( GUNZIP_REFERENCE.out.genome )

    eukaryotic_reference = TAXONOMY_REFERENCE.out.taxonomy.map { meta, row ->
        def lineage = row.ncbi_lng?.trim()
        if ( !lineage || lineage == "NA" || !lineage.tokenize('-').contains("2759") ) {
            error "Reference genome ${meta.id} (sample ${meta.sample}) is not eukaryotic: EukCC lineage '${lineage}'. This pipeline is eukaryote only, see proposal.md section 1."
        }
        [ meta.sample ]
    }

    // joined on the check above, so no MAG is ever aligned against a reference that failed it
    anchors = GUNZIP_REFERENCE.out.genome
        .map { meta, fasta -> [meta.sample, fasta] }
        .join( eukaryotic_reference )

    /*
     * Cycle 1 - the sample on its own, assembled and binned, and the MAG closest to the
     * reference is what Branchwater then searches with.
     *
     * --skip_first_assembly drops that whole cycle and searches with the reference genome
     * from the samplesheet instead. The run is then reference + cycle 2 only: no cycle-1
     * assembly, no cycle-1 GGP, no `cycle1` source in any table. Use it when the sample has
     * already been assembled elsewhere, or when the question is only what the public runs
     * add to a genome you already have.
     */
    if ( params.skip_first_assembly ) {
        cycle1_taxonomy = Channel.empty()
        cycle1_skani    = Channel.empty()
        cycle1_target   = Channel.empty()
        cycle1_versions = Channel.empty()

        // the samplesheet genome as it is, not the GUNZIP copy: sourmash reads gzip natively
        branchwater_query = ch_samplesheet.map { meta, _reads, genome ->
            [ [id: genome_name(genome), sample: meta.id], genome ]
        }
    }
    else {
        CYCLE1(
            ch_samplesheet.map { meta, reads, _genome -> [ assembler_meta(meta.id, CYCLE1_TAG), reads ] }
        )

        cycle1_genomes = CYCLE1.out.bins.combine( sample_ids.map { [it] } ).map { bin, ids ->
            [ [id: genome_name(bin), sample: owner_of(bin.name, ids), origin: "cycle1", source: "cycle1"], bin ]
        }

        GUNZIP_CYCLE1( cycle1_genomes )
        TAXONOMY_CYCLE1( GUNZIP_CYCLE1.out.genome )

        SKANI_CYCLE1( source_alignments( GUNZIP_CYCLE1.out.genome, anchors ) )

        cycle1_taxonomy = TAXONOMY_CYCLE1.out.taxonomy
        cycle1_skani    = SKANI_CYCLE1.out.hits
        cycle1_target   = pick_target( cycle1_skani, cycle1_taxonomy, GUNZIP_CYCLE1.out.genome )
        cycle1_versions = CYCLE1.out.versions.mix( SKANI_CYCLE1.out.versions )

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
        ch_samplesheet.map { meta, reads, _genome -> [meta.id, reads] }
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
                    source: "n" + (dataset =~ /_n(\d+)$/)[0][1]
                ],
                bin
            ]
        }
    )
    TAXONOMY_CYCLE2( GUNZIP_CYCLE2.out.genome )

    SKANI_CYCLE2( source_alignments( GUNZIP_CYCLE2.out.genome, anchors ) )

    cycle2_target = pick_target( SKANI_CYCLE2.out.hits, TAXONOMY_CYCLE2.out.taxonomy, GUNZIP_CYCLE2.out.genome )

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
            [meta.id, meta.sample, meta.origin, meta.source, row.taxid, row.completeness, row.contamination].join('\t') + '\n'
        }
        .collectFile(
            name: "all_mags.tsv",
            storeDir: "${params.outdir}/mags",
            sort: true,
            seed: "genome\tsample\torigin\tsource\ttaxid\tcompleteness\tcontamination\n"
        )

    // every skani row of every MAG, so the choice in target_mags.tsv can be checked
    cycle1_skani.mix( SKANI_CYCLE2.out.hits )
        .flatMap { meta, tsv ->
            tsv.splitCsv(header: true, sep: '\t').collect { hit ->
                [meta.sample, meta.source, genome_name(file(hit.Query_file)), hit.ANI, hit.Align_fraction_ref, hit.Align_fraction_query].join('\t') + '\n'
            }
        }
        .collectFile(
            name: "skani_vs_reference.tsv",
            storeDir: "${params.outdir}/mags",
            sort: true,
            seed: "sample\tsource\tgenome\tani\taf_reference\taf_genome\n"
        )

    // one target MAG per sample and cycle/N, NA when none aligned to the reference
    found_targets = cycle1_target.mix( cycle2_target )
        .map { meta, row, _fasta ->
            [ "${meta.sample}\t${meta.source}".toString(), [meta.id, row.ani, row.af_reference, row.taxid, row.completeness, row.contamination] ]
        }
        .toList()

    /*
     * One row per source that could be built: cycle 1 for every sample, plus the
     * co-assembly datasets BUILD_CONCAT_DATASETS actually assembled - a sample whose
     * only Branchwater hit is its own run has no cycle-2 source at all, and one with fewer
     * usable hits than --n_concat_samples asks for has its depths capped (it warns).
     */
    cycle1_sources = params.skip_first_assembly
        ? Channel.empty()
        : ch_samplesheet.map { meta, _reads, _genome -> "${meta.id}\tcycle1".toString() }

    expected_sources = cycle1_sources
        .mix(
            BUILD_CONCAT_DATASETS.out.reads.map { meta, _reads ->
                "${meta.sample}\tn${meta.n_concat}".toString()
            }
        )

    expected_sources
        .combine( found_targets.map { [it] } )
        .map { key, found ->
            def values = found.collectEntries()[key] ?: ["NA"] * 6
            ([key] + values).join('\t') + '\n'
        }
        .collectFile(
            name: "target_mags.tsv",
            storeDir: "${params.outdir}/mags",
            sort: true,
            seed: "sample\tsource\tgenome\tani\taf_reference\ttaxid\tcompleteness\tcontamination\n"
        )

    /*
     * The comparison set: reference, cycle 1, one cycle 2 MAG per N - COMPARE orders them.
     * Each target carries its skani numbers against the reference and the runs whose reads
     * were assembled into it: the sample alone in cycle 1, the sample plus the Branchwater
     * hits in cycle 2. The reference has neither.
     */
    concatenated_runs = BUILD_CONCAT_DATASETS.out.reads.map { meta, _reads ->
        [ [meta.sample, "n${meta.n_concat}".toString()], meta.sources.join(',') ]
    }

    comparison_set = GUNZIP_REFERENCE.out.genome
        .map { meta, fasta -> [ meta + [ani: "NA", af_reference: "NA", runs: "NA"], fasta ] }
        .mix(
            cycle1_target.map { meta, row, fasta ->
                [ meta + [ani: row.ani, af_reference: row.af_reference, runs: meta.sample], fasta ]
            },
            cycle2_target
                .map { meta, row, fasta -> [ [meta.sample, meta.source], meta, row, fasta ] }
                .join( concatenated_runs )
                .map { _key, meta, row, fasta, runs ->
                    [ meta + [ani: row.ani, af_reference: row.af_reference, runs: runs], fasta ]
                }
        )

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
        SKANI_CYCLE2.out.versions,
        BUILD_CONCAT_DATASETS.out.versions,
    )

    ch_versions = ch_versions.filter { versions_file ->
        if ( parses_as_versions_map( versions_file ) ) {
            return true
        }
        log.warn "Skipping unparseable versions.yml: ${versions_file}"
        return false
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'mag_polishing_pipeline_software_versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    comparison     = COMPARE.out.metrics   // channel: path(output.tsv)
    versions       = ch_collated_versions  // channel: path(versions.yml)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parses_as_versions_map( yaml_file ) {
    try {
        return new org.yaml.snakeyaml.Yaml().load( yaml_file ) instanceof Map
    }
    catch ( Exception ignored ) {  // noqa - an unparseable version string is not a failure
        return false
    }
}

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
 * One SKANI_DIST input per assembly: all of its MAGs against the sample's reference.
 * Cycle 1 has one assembly per sample, cycle 2 one per sample and N.
 */
def source_alignments( genomes, anchors ) {
    return genomes
        .map { meta, fasta -> [ [meta.sample, meta.source], fasta ] }
        .groupTuple()
        .map { key, fastas -> [ key[0], key[1], fastas ] }
        .combine( anchors, by: 0 )
        .map { sample, source, fastas, reference ->
            [ [id: "${sample}_${source}".toString(), sample: sample, source: source], reference, fastas ]
        }
}

/*
 * The hit that is the reference organism: at least `min_ani` ANI, then the one covering
 * most of the reference, then the most complete. 95% ANI is the usual species boundary.
 * null when nothing qualifies - the target MAG is then an NA row, never a failure.
 */
def best_match( hits, completeness, min_ani = 95 ) {
    return hits
        .findAll { as_number(it.ani) >= min_ani }
        .max { a, b ->
            as_number(a.af_reference) <=> as_number(b.af_reference) ?:
                as_number(completeness[a.genome]) <=> as_number(completeness[b.genome])
        }
}

/*
 * The MAG closest to the reference, one per sample and cycle/N, as [meta, taxonomy row + ani/af_reference, fasta].
 * The EukCC taxid stays in the row as information only; it no longer decides anything.
 */
def pick_target( skani_hits, taxonomy, genomes ) {
    return skani_hits
        .map { meta, tsv -> [ [meta.sample, meta.source], tsv.splitCsv(header: true, sep: '\t') ] }
        .join( taxonomy.map { meta, row -> [ [meta.sample, meta.source], [meta, row] ] }.groupTuple() )
        .flatMap { _key, hits, candidates ->
            // a genome EukCC gave up on has no taxonomy row, and so cannot be a target
            def by_id = candidates.collectEntries { meta, row -> [ (meta.id): [meta, row] ] }
            def best  = best_match(
                hits.collect { hit ->
                    [ genome: genome_name(file(hit.Query_file)), ani: hit.ANI, af_reference: hit.Align_fraction_ref ]
                }.findAll { by_id.containsKey(it.genome) },
                by_id.collectEntries { id, entry -> [ (id): entry[1].completeness ] }
            )
            if ( !best ) {
                return []
            }
            def (meta, row) = by_id[best.genome]
            [ [ meta.id, meta, row + [ani: best.ani, af_reference: best.af_reference] ] ]
        }
        .join( genomes.map { meta, fasta -> [meta.id, fasta] } )
        .map { _id, meta, row, fasta -> [meta, row, fasta] }
}
