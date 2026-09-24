/*
 * The comparison table: the reference genome against the MAG each cycle recovered.
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
    samples     // channel: [ val(sample), [ path(fasta), ... ] ]
    eukcc_rows  // channel: [ val(genome), [completeness: .., contamination: ..] ]  genome == fasta.baseName

    main:
    assemblies = samples.flatMap { _sample, fastas -> fastas }

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
        "sample", "assembly", "length", "n50", "gc_content", "n_contigs",
        "busco", "eukcc_completeness", "eukcc_contamination"
    ]

    // each side is wrapped in a list because combine() spreads list items into the tuple
    metrics = samples.toList().map { [it] }
        .combine( busco_scores.toList().map { [it] } )
        .combine( eukcc_rows.toList().map { [it] } )
        .combine( assembly_stats.toList().map { [it] } )
        .map { rows, busco_list, eukcc_list, stats_list ->
            def busco = busco_list.collectEntries()
            def eukcc = eukcc_list.collectEntries()
            def stats = stats_list.collectEntries()
            def lines = rows.collect { sample, fastas ->
                [
                    sample,
                    cell( fastas.collect { it.name } ),
                    cell( fastas.collect { stats[it.baseName]?.Length } ),
                    cell( fastas.collect { stats[it.baseName]?.N50 } ),
                    cell( fastas.collect { stats[it.baseName]?.GC_content } ),
                    cell( fastas.collect { stats[it.baseName]?.N_contigs } ),
                    cell( fastas.collect { busco[it.baseName] } ),
                    cell( fastas.collect { eukcc[it.baseName]?.completeness } ),
                    cell( fastas.collect { eukcc[it.baseName]?.contamination } )
                ].join('\t')
            }
            ([header.join('\t')] + lines).join('\n') + '\n'
        }
        .collectFile(name: "assembly_qc_metrics.tsv", storeDir: "${params.outdir}/compare")

    emit:
    metrics = metrics
}

// one value per assembly, "/" separated, in the order the comparison set was built
def cell(values) {
    values.collect { it ?: "NA" }.join("/")
}
