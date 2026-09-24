# EBI-Metagenomics/mag-polishing-pipeline: Usage

## Samplesheet input

`--input`, CSV with a header:

```csv
sample,fastq_1,fastq_2,genome
SAMPLE_001,/data/SAMPLE_001_1.fastq.gz,/data/SAMPLE_001_2.fastq.gz,/data/SAMPLE_001.fna.gz
```

| column                | notes                                                                           |
| --------------------- | ------------------------------------------------------------------------------- |
| `sample`              | **at least 7 characters** — miassembler substrings it to build its output paths |
| `fastq_1` / `fastq_2` | paired-end only, `.fastq.gz`                                                    |
| `genome`              | the reference genome/MAG, **must be eukaryotic**, `.fa(.gz)` / `.fna(.gz)`      |

`assets/samplesheet_example.csv` is a copy of the above.

### Running

#### Locally (Docker)

```bash
export NXF_VER=24.04.3

nextflow run . -profile docker \
  --input assets/samplesheet_example.csv \
  --n_concat_samples 5 \
  --eukcc_db /path/to/eukcc2_db \
  --busco_db /path/to/busco \
  --cat_db_folder /path/to/cat_db \
  --cat_taxonomy_db /path/to/cat_taxonomy \
  --ref_genome /path/to/human.fna \
  --branchwater_index /path/to/branchwater_db/ \
  --branchwater_metadata_db /path/to/metadata.duckdb \
  --ena_cache_dir /path/to/ena_cache \
  --outdir results
```

#### On Codon

```bash
export NXF_VER=24.04.3

nextflow run . -profile codon \
  --input /path/to/samplesheet.csv \
  --n_concat_samples 5 \
  --outdir results
```

### Parameters

Every param the pipeline declares, grouped the way `nextflow.config` groups them. Anything
marked _from the profile_ is set by `-profile codon` and only needs a value when running
elsewhere. Params the composed submodules read but this pipeline never exercises are in
[Inert params](#inert-params) at the end.

#### Pipeline

| param                   | default     | meaning                                                                                                                                                               |
| ----------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`               | —           | required, see [Input](#input)                                                                                                                                         |
| `--n_concat_samples`    | `1,5,10`    | comma separated; one co-assembly dataset per value, capped at the number of usable Branchwater hits. See [Co-assembly depths](#co-assembly-depths---n_concat_samples) |
| `--skip_first_assembly` | `false`     | skip cycle 1 entirely and search Branchwater with the samplesheet genome; no cycle-1 assembly, no cycle-1 GGP, no `cycle1` slot in any table                          |
| `--outdir`              | `results`   |                                                                                                                                                                       |
| `--ena_cache_dir`       | `ena_cache` | shared fastq cache (`storeDir`); every run accession is downloaded once, across samples and N values. Expect hundreds of GB                                           |
| `--publish_dir_mode`    | `copy`      | any Nextflow `publishDir` mode; `symlink` saves space but breaks once `work/` is deleted                                                                              |

#### Branchwater

| param                       | default          | meaning                                                       |
| --------------------------- | ---------------- | ------------------------------------------------------------- |
| `--branchwater_index`       | from the profile | sourmash rocksdb index searched by `manysearch`               |
| `--branchwater_metadata_db` | from the profile | Branchwater metadata DuckDB, used to resolve hits to ENA runs |
| `--branchwater_k`           | `21`             | k-mer size; must match the index (`bw_k21`)                   |
| `--branchwater_scaled`      | `1000`           | sourmash sketch scaling factor                                |
| `--branchwater_threshold`   | `0.1`            | minimum containment for a hit to be considered                |

#### Assembly (miassembler)

| param                     | default      | meaning                                                                                                                                               |
| ------------------------- | ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--assembler`             | `metaspades` | `metaspades` or `megahit`. Anything else is rejected before the run starts. Applies to **both** cycles                                                |
| `--assembly_memory`       | `100`        | GB for the assembler, +50% per retry. Per sample it can be overridden with `meta.assembly_memory`                                                     |
| `--spades_only_assembler` | `true`       | passes `--only-assembler` to metaSPAdes, skipping read error correction. Big time and memory saving; turn off only if you want SPAdes' own correction |
| `--max_spades_retries`    | `3`          | retries for metaSPAdes, each with 50% more memory                                                                                                     |
| `--max_megahit_retries`   | `3`          | same for MEGAHIT                                                                                                                                      |
| `--spades_version`        | `3.15.5`     | recorded in `meta.assembler_version`, and part of the assembly publish path                                                                           |
| `--megahit_version`       | `1.2.9`      | same for MEGAHIT                                                                                                                                      |

#### Assembly QC (miassembler)

These gate what leaves the assembler and reaches GGP.

| param                                     | default       | meaning                                                                                                   |
| ----------------------------------------- | ------------- | --------------------------------------------------------------------------------------------------------- |
| `--short_reads_min_contig_length`         | `500`         | contigs shorter than this are dropped (`seqkit seq`), and it is QUAST's `--min-contig`                    |
| `--short_reads_contig_threshold`          | `2`           | an assembly with fewer contigs than this fails QC and produces no MAG                                     |
| `--short_reads_low_reads_count_threshold` | `1000`        | fewer reads than this after fastp and the run is marked QC-failed                                         |
| `--short_reads_filter_ratio_threshold`    | `0.1`         | fastp keeping this fraction or less of the reads is also a QC failure                                     |
| `--min_qcov` / `--min_pid`                | `0.3` / `0.4` | query coverage and identity floors for `filterpaf`, which decides what counts as a contaminant contig hit |

#### Decontamination references (miassembler)

Resolved as `<reference_genomes_folder>/<name>/<name>.fna`, so the names are directory
names, not paths.

| param                        | default                                      | meaning                                     |
| ---------------------------- | -------------------------------------------- | ------------------------------------------- |
| `--reference_genomes_folder` | from the profile                             | root holding the decontamination references |
| `--human_reference`          | from the profile (`human_GCF_000001405.40`)  | human decontamination of the reads          |
| `--phix_reference`           | from the profile (`phiX174_GCF_000819615.1`) | PhiX removal from the contigs               |
| `--contaminant_reference`    | `null`                                       | optional extra host/contaminant genome      |

#### Binning (genomes-generation)

| param                                                                                                                                           | default                                                                             | meaning                                                                                                |
| ----------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `--skip_prok`                                                                                                                                   | `true`                                                                              | **keep it true.** This is a eukaryotic workflow at the moment; the prokaryotic branch is not exercised |
| `--skip_euk`                                                                                                                                    | `false`                                                                             | the eukaryotic branch is the one that produces the MAGs being compared                                 |
| `--skip_decontamination`                                                                                                                        | `false`                                                                             | GGP's own read decontamination, against `--ref_genome`                                                 |
| `--skip_preprocessing_input`                                                                                                                    | `true`                                                                              | **keep it true.** GGP's preprocessing expects ENA-derived inputs; ours are not                         |
| `--merge_pairs`                                                                                                                                 | `false`                                                                             | merge overlapping read pairs in GGP's fastp step                                                       |
| `--min_contig_size`                                                                                                                             | `1500`                                                                              | minimum contig length fed to the binners                                                               |
| `--metabat2_rng_seed`                                                                                                                           | `1`                                                                                 | MetaBAT2's `--seed`; fixed so binning is reproducible                                                  |
| `--publish_all`                                                                                                                                 | `false`                                                                             | GGP's own switch for publishing its intermediates. The genomes are published regardless                |
| `--ena_assembly_study_accession`                                                                                                                | `mag-polishing`                                                                     | prefix on GGP's report file names — **not** an accession, unless uploads are on (see below)            |
| `--subdir_euks`, `--subdir_proks`, `--subdir_bins`, `--subdir_mags`, `--subdir_stats`, `--subdir_taxonomy`, `--subdir_coverage`, `--subdir_rna` | `eukaryotes`, `prokaryotes`, `bins`, `mags`, `stats`, `taxonomy`, `coverage`, `rna` | names of GGP's output subdirectories                                                                   |

#### ENA upload (genomes-generation)

All off. See [Uploading to ENA](#uploading-to-ena) — these submit for real.

| param                             | default                       | meaning                                              |
| --------------------------------- | ----------------------------- | ---------------------------------------------------- |
| `--upload_mags` / `--upload_bins` | `false`                       | submit the MAGs / the bins to ENA                    |
| `--test_upload`                   | `false`                       | `true` submits to ENA's test service instead of live |
| `--upload_tpa`                    | `false`                       | `--tpa`, third party annotation                      |
| `--upload_force`                  | `false`                       | `--force`, overwrite an existing submission          |
| `--metagenome`                    | `""`                          | required by `genome_upload`, e.g. `soil metagenome`  |
| `--biomes`                        | `""`                          |                                                      |
| `--centre_name`                   | `""` (`EMG` from the profile) | submitting centre                                    |

#### Reference databases

| param                                                        | default          | meaning                                                                                                         |
| ------------------------------------------------------------ | ---------------- | --------------------------------------------------------------------------------------------------------------- |
| `--eukcc_db`                                                 | from the profile | EukCC reference data — taxonomy _and_ completeness/contamination                                                |
| `--busco_db`                                                 | from the profile | BUSCO lineage downloads                                                                                         |
| `--busco_mode`                                               | `genome`         | BUSCO `-m`; `genome`, `transcriptome` or `proteins`                                                             |
| `--cat_db_folder` / `--cat_diamond_db` / `--cat_taxonomy_db` | from the profile | CAT/BAT taxonomy assignment inside GGP                                                                          |
| `--ref_genome`                                               | from the profile | GGP's decontamination reference. Read even when `--skip_decontamination` is on, so it must point at a real file |

#### Environment

| param                    | default          | meaning                                                                            |
| ------------------------ | ---------------- | ---------------------------------------------------------------------------------- |
| `--singularity_cachedir` | from the profile | shared Singularity image cache, used by `-profile codon` as `singularity.cacheDir` |

#### Inert params

Declared because the composed submodules resolve them at config time, but never reached by
this pipeline's path.

| param                                                                                                                               | default                  | why it is inert                                                                                                                |
| ----------------------------------------------------------------------------------------------------------------------------------- | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| `--checkm2_db`, `--gunc_db`, `--gtdbtk_db`, `--rfam_rrna_models`                                                                    | from the profile         | prokaryotic branch only, and `--skip_prok true`                                                                                |
| `--long_reads_min_read_length`, `--long_reads_ont_quality_threshold`, `--long_reads_pacbio_quality_threshold`, `--max_flye_retries` | `200`, `0.8`, `0.9`, `3` | miassembler's long-read path; we only call `SHORT_READS_ASSEMBLER`                                                             |
| `--short_reads_min_contig_length_metat`                                                                                             | `200`                    | metatranscriptome path                                                                                                         |
| `--diamond_db`                                                                                                                      | from the profile         | miassembler's frameshift correction, long reads only                                                                           |
| `--private_study`, `--study_accession`, `--reads_accession`                                                                         | `false`, `null`, `null`  | miassembler's ENA fetch path. We build `meta` ourselves in `main.nf`, so `study_accession` is only a fallback that never fires |
| `--download_data`                                                                                                                   | `false`                  | GGP's ENA download path, replaced by our own samplesheet                                                                       |
| `--multiqc_title`                                                                                                                   | `null`                   | GGP v1.3.3 has its `MULTIQC` import commented out, so no report is produced                                                    |

### Co-assembly depths: `--n_concat_samples`

This is the knob the experiment turns. Everything else in the pipeline is fixed for a
given sample; `--n_concat_samples` is what you sweep to ask _how many public runs are
worth adding before the MAG stops improving_.

It takes a comma-separated list of depths, and **each value becomes its own co-assembly
dataset, its own de novo assembly and its own MAG**. With `--n_concat_samples 1,5,10` a
sample produces:

| value | dataset                                  | assembled as   |
| ----- | ---------------------------------------- | -------------- |
| `1`   | sample reads + the top-1 Branchwater hit | `<sample>_n1`  |
| `5`   | sample reads + the top-5 hits            | `<sample>_n5`  |
| `10`  | sample reads + the top-10 hits           | `<sample>_n10` |

The hits are the same ranked list every time — `n5` is `n1` plus the next four — so the
depths are nested and the comparison reads as a curve. Each one lands in its own column
slot of `compare/assembly_qc_metrics.tsv`, in ascending order, next to the reference
genome and the cycle-1 MAG. That is the whole point: a single value tells you what one
co-assembly produced, a list tells you whether more data kept helping.

Reasonable sweeps:

- `--n_concat_samples 1` — not a benchmark, just "co-assemble with the best hit". Cheapest
  useful run.
- `--n_concat_samples 1,5,10` (the default) — the standard curve, coarse enough to see the
  shape without paying for every step.
- `--n_concat_samples 1,2,3,5,10,20` — a fine sweep for one or two samples when you are
  actually looking for the point of diminishing returns. Expensive; see below.

#### This is storage hungry

Every value in the list writes a **full pair of concatenated fastqs** into the work
directory and then feeds them to metaSPAdes. The reads are not shared between depths —
`n5` does not reuse `n1`'s file, it is a new concatenation — so the fastqs written for one
sample are roughly:

```
sum over N of (1 sample + N runs)
```

For `1,5,10` that is 2 + 6 + 11 = **19 run-equivalents** of fastq for one sample, on top of
the ENA cache holding the 10 distinct runs themselves. For `1,2,3,5,10,20` it is
2 + 3 + 4 + 6 + 11 + 21 = **47**. Ten samples on the default list is ~190. At a few GB per
paired run — typical for a metagenome in ENA — that is hundreds of GB of intermediates
before a single assembly starts.

Then each dataset is assembled: the `n=10` assembly is metaSPAdes on 11x the depth of the
original sample, and its work directory is sized accordingly. Budget for the assemblies,
not just the reads.

What softens it:

- `--ena_cache_dir` is a `storeDir` shared across samples and across N values, so every run
  accession is downloaded exactly once no matter how many datasets use it. Point it
  somewhere with room and keep it between runs.
- The concatenated fastqs live in `work/` and are intermediates — once the run is done and
  you no longer need `-resume`, `nextflow clean` or deleting `work/` reclaims all of it.
  Nothing in `results/` depends on them; `concat_datasets/<sample>_nN/` keeps only the
  provenance tsv saying which accessions went in.
- Depths are capped and deduplicated against what Branchwater actually delivered, so
  asking for more than exists costs nothing extra (see Gotchas).

If storage is the binding constraint, `--assembler megahit` cuts the assembly side
sharply, and a short list on many samples usually answers more than a long list on one.

### Uploading to ENA

GGP ships an ENA uploader (`genomes_uploader` + `ena-webin-cli`). It is **off**, and it
should stay off for now: this first version of pipeline builds MAGs to benchmark them, and the
cycle-2 genomes come out of co-assemblies of other people's public runs. Submitting them
is a deliberate decision, not a side effect of running the workflow.

`--upload_mags true` and/or `--upload_bins true` turn it on. Both submit **for real** —
there is no dry run, only ENA's test service.

| param                            | default                             | meaning                                                                                                        |
| -------------------------------- | ----------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `--upload_mags`                  | `false`                             | submit the dereplicated MAGs                                                                                   |
| `--upload_bins`                  | `false`                             | submit the bins                                                                                                |
| `--test_upload`                  | `false`                             | `true` submits to ENA's **test** service (`ena-webin-cli -test`). Always do this first                         |
| `--ena_assembly_study_accession` | `mag-polishing`                     | passed verbatim to `genome_upload -u`, so with uploads on it has to be a registered ENA study (`PRJ…`, `ERP…`) |
| `--metagenome`                   | `""`                                | required, e.g. `soil metagenome`                                                                               |
| `--biomes`                       | `""`                                |                                                                                                                |
| `--centre_name`                  | `""` (`EMG` under `-profile codon`) |                                                                                                                |
| `--upload_tpa`                   | `false`                             | `--tpa`, third party annotation                                                                                |
| `--upload_force`                 | `false`                             | `--force`, overwrite an existing submission                                                                    |

Credentials are Nextflow secrets, not params:

```bash
nextflow secrets set WEBIN_ACCOUNT '...'
nextflow secrets set WEBIN_PASSWORD '...'
```

`main.nf` fails fast when uploads are on but `--ena_assembly_study_accession` is not an
accession or `--metagenome` is empty — better than finding out after two assembly cycles.
It also logs a warning naming the study and whether the submission is live or test.

Results land in `${outdir}/upload/{mags,bins}/ena_submission_summary.txt`.

Note that `--ena_assembly_study_accession` does double duty: with uploads off it is only a
filename prefix on GGP's report files, which is why the default is a name and not an
accession.

### Gotchas

- **Sample ids under 7 characters** abort the run with an explicit message.
- **Genome file names must be unique** across a comparison — every metric is keyed back to
  the basename.
- **EukCC exits 201 when it finds no marker genes.** Tolerated: the genome gets taxid `NA`,
  matches nothing, and the run continues. The slot is reported as `NA` in
  `mags/target_mags.tsv`.
- **A cycle that recovers no MAG with the target taxid** is never a failure; it is an `NA`
  row and a shorter comparison set for that sample.
- **Fewer usable hits than `--n_concat_samples` asks for** is not a failure either. The
  requested depths are capped at what Branchwater delivered and deduplicated, with a
  warning: 3 usable hits and `--n_concat_samples 1,5,10` builds `n1` and `n3`, not three
  datasets of which two would be identical. A sample whose only hit is its own run gets
  no cycle-2 dataset at all - there is nothing to co-assemble with - and
  `mags/target_mags.tsv` simply has no row for it.
- **`-resume` works per process** across the whole graph, including both submodules — the
  point of composing everything in one runtime.

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start requesting a large amount of
memory. We recommend adding the following line to your environment to limit this (typically
in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
