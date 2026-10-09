/*
 * The comparison table: the reference genome against the MAG each cycle recovered, one
 * row per genome.
 *
 * EukCC is taken as an input rather than run here, because assign_taxonomy.nf already
 * runs it on every genome in the experiment to assign the taxid - one EukCC run per
 * genome, whichever side of the comparison it is on (proposal.md section 5).
 */

include { BUSCO                    } from '../../../modules/local/busco'
include { CALCULATE_ASSEMBLY_STATS } from '../../../modules/local/calculate_assembly_stats'
include { COLLECT_ASSEMBLIES       } from '../../../modules/local/collect_assemblies'

workflow COMPARE {

    take:
    genomes     // channel: [ val(meta), path(fasta) ]  meta: sample, source (reference | cycle1 | nN), runs, ani, af_reference
    eukcc_rows  // channel: [ val(genome), [completeness: .., contamination: ..] ]  genome == fasta.baseName

    main:
    assemblies = genomes.map { _meta, fasta -> fasta }

    // every metric is keyed back to the assembly basename, so they have to be unique
    assemblies.map { it.baseName }.toList().subscribe { names ->
        def dups = names.countBy { it }.findAll { _name, count -> count > 1 }.keySet()
        if ( dups ) {
            error "Assembly file names must be unique across the comparison, duplicated: ${dups.join(', ')}"
        }
    }

    BUSCO(
        assemblies,
        file(params.busco_db, checkIfExists: true),
        params.busco_mode
    )

    COLLECT_ASSEMBLIES( assemblies.collect() )
    CALCULATE_ASSEMBLY_STATS( COLLECT_ASSEMBLIES.out.assemblies_dir )

    // "<assembly file name>\tC:..%[S:..%,D:..%],F:..%,M:..%,n:.."
    busco_scores = BUSCO.out.busco_summary
        .splitText()
        .map { line -> def fields = line.trim().split('\t', 2); [file(fields[0]).baseName, fields[1]] }

    assembly_stats = CALCULATE_ASSEMBLY_STATS.out.stats_file
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.Genome, row] }

    header = [
        "sample", "assembly", "source", "concatenated_runs", "length", "n50", "gc_content", "n_contigs",
        "eukcc_completeness", "eukcc_contamination", "busco", "ani_reference", "af_reference"
    ]

    // each side is wrapped in a list because combine() spreads list items into the tuple
    metrics = genomes.toList().map { [it] }
        .combine( busco_scores.toList().map { [it] } )
        .combine( eukcc_rows.toList().map { [it] } )
        .combine( assembly_stats.toList().map { [it] } )
        .map { rows, busco_list, eukcc_list, stats_list ->
            def busco = busco_list.collectEntries()
            def eukcc = eukcc_list.collectEntries()
            def stats = stats_list.collectEntries()
            def lines = rows
                .sort { a, b -> a[0].sample <=> b[0].sample ?: source_rank(a[0].source) <=> source_rank(b[0].source) }
                .collect { meta, fasta ->
                    def name = fasta.baseName
                    [
                        meta.sample,
                        fasta.name,
                        source_label(meta.source),
                        meta.runs,
                        stats[name]?.Length,
                        stats[name]?.N50,
                        stats[name]?.GC_content,
                        stats[name]?.N_contigs,
                        eukcc[name]?.completeness,
                        eukcc[name]?.contamination,
                        busco[name],
                        meta.ani,
                        meta.af_reference
                    ].collect { it ?: "NA" }.join('\t')
                }
            ([header.join('\t')] + lines).join('\n') + '\n'
        }
        .collectFile(name: "output.tsv", storeDir: "${params.outdir}")

    emit:
    metrics = metrics
}

// reference first, then cycle 1, then the cycle-2 depths in ascending N
def source_rank( source ) {
    source == "reference" ? -1 : source == "cycle1" ? 0 : source.substring(1) as int
}

// "reference" -> "reference", "cycle1" -> "1st cycle", "n5" -> "2nd cycle n5"
def source_label( source ) {
    source == "reference" ? "reference" : source == "cycle1" ? "1st cycle" : "2nd cycle ${source}".toString()
}
