# mag-polishing-pipeline

Improves MAG quality by recruiting additional related samples through sequence similarity search. Assembles an initial sample via [miassembler](https://github.com/EBI-Metagenomics/miassembler) (optional), bin MAGs via [genomes-generation pipeline](https://github.com/EBI-Metagenomics/genomes-generation) (GGP), then queries branchwater for matching samples.

For each sample (paired reads + a reference genome) the pipeline assembles and bins the sample on its own, finds the public runs containing the target organism, and re-assembles the sample together with the top N of them (`--n_concat_samples`). The output is one table comparing the provided reference genome, the cycle-1 MAG (optional) and one cycle-2 MAG per co-assembly depth (per `--n_concat_samples` value).

With `--skip_first_assembly` the first cycle is dropped: Branchwater searches with the parsed reference genome and the run is reference vs cycle 2 only.

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/EBI-Metagenomics/mag-polishing-pipeline)
[![GitHub Actions CI Status](https://github.com/EBI-Metagenomics/mag-polishing-pipeline/actions/workflows/nf_tests.yml/badge.svg)](https://github.com/EBI-Metagenomics/mag-polishing-pipeline/actions/workflows/nf_tests.yml)
[![GitHub Actions Linting Status](https://github.com/EBI-Metagenomics/mag-polishing-pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/EBI-Metagenomics/mag-polishing-pipeline/actions/workflows/linting.yml)
[![Python Tests Status](https://github.com/EBI-Metagenomics/mag-polishing-pipeline/actions/workflows/python_tests.yml/badge.svg)](https://github.com/EBI-Metagenomics/mag-polishing-pipeline/actions/workflows/python_tests.yml)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.2.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.2.0)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Requirements

- **Nextflow 24.04.3**
- Singularity or Docker, plus the reference databases listed below.
- Java 17+, git.

## Configure the environment

```bash
git clone https://github.com/EBI-Metagenomics/mag-polishing-pipeline.git
cd mag-polishing-pipeline
git submodule update --init --recursive
```

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation)
> on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with your input data:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,genome
SAMPLE_001,/data/SAMPLE_001_1.fastq.gz,/data/SAMPLE_001_2.fastq.gz,/data/SAMPLE_001.fna.gz
```

Each row is one sample: paired-end reads plus the eukaryotic reference genome it is
compared against. Then:

```bash
NXF_VER=24.04.3 nextflow run EBI-Metagenomics/mag-polishing-pipeline \
   -profile <docker/singularity/codon> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config
> files including those provided by the `-c` Nextflow option can be used to provide any
> configuration _**except for parameters**_; see
> [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](docs/usage.md)
and the [parameter documentation](docs/usage.md#parameters).

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Checks

```bash
# python unit tests: bin/*.py
pip install -r requirements-dev.txt
pytest

# nextflow unit tests: functions, local modules, local subworkflows
NXF_VER=24.04.3 nf-test test

export NXF_VER=24.04.3

# builds the whole DAG, both cycles, without running anything
nextflow run . -preview -profile test --outdir /tmp/mag-polishing

# the merged config, including both upstream modules.config files
nextflow config . -profile codon

# the nf-core template contract. The version MUST match `nf_core_version` in .nf-core.yml
NXF_SYNTAX_PARSER=v1 uvx --from 'nf-core==3.2.0' nf-core pipelines lint
```

`tests/` follows the layout of
[genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline/tree/master/tests):

```
tests/
├── scripts/       pytest for bin/*.py, with fixtures/
├── functions/     nf-test for the Groovy helpers (genome_name, owner_of, concat_depths)
├── modules/       nf-test for modules/local/*, with fixtures/
├── subworkflows/  nf-test for subworkflows/local/*, with fixtures/
└── fixtures/      the samplesheet and empty reads that -profile test points at
```

Every nf-test drives a `stub:` block and runs with no container engine
(`tests/nextflow.config`), so the whole suite is seconds and needs neither Docker nor the
reference databases. The two exceptions are `GUNZIP_GENOME`'s tests, which run for real
because `gunzip` is on any host.

A full run of `main.nf` is not testable here: several miassembler and GGP processes
(`SPADES` among them) have no `stub` block, so they would actually run. That is why the
pipeline-level check is `-preview` rather than an nf-test.

> [!IMPORTANT]
> **Pin the nf-core tools version.** This repo is on template 3.2.0 (`nf_core_version` in
> `.nf-core.yml`), and nf-core's linter is not compatible across template generations:
>
> | tools | on this repo | why |
> |---|---|---|
> | **3.2.0** | 191 passed, 0 failed | matches the template |
> | 4.0.2 | `KeyError: 'manifest.name'` | broken upstream — it crashes on a pipeline it generated itself |
> | 4.1.0 | `IndexError: string index out of range` | 4.x linter cannot read a 3.2.0 template |
>
> The CI workflow gets this right automatically: `.github/workflows/linting.yml` reads
> `nf_core_version` out of `.nf-core.yml` and installs exactly that.

> [!NOTE]
> `NXF_SYNTAX_PARSER=v1` is only needed when the Nextflow on your `PATH` is 25.04 or newer,
> which `nf-core` shells out to. miassembler's `conf/modules.config` calls
> `study_reads_folder(meta)` from its `publishDir` closures, and the strict config parser
> introduced in Nextflow 25 rejects function definitions in `nextflow.config`. The pipeline
> itself runs on 24.04.3, where the default parser handles it.

## Credits

EBI-Metagenomics/mag-polishing-pipeline was originally written by Filipe Dezordi.

This pipeline is an orchestrator. It composes, in a single Nextflow runtime:

| Component | Where | Pinned at |
|---|---|---|
| [miassembler](https://github.com/EBI-Metagenomics/miassembler) | git submodule, `pipelines/miassembler` | `v3.1.16` |
| [genomes-generation](https://github.com/EBI-Metagenomics/genomes-generation) | git submodule, `pipelines/genomes-generation` | `v1.3.3` |
| [branchwater-nf](https://github.com/EBI-Metagenomics/branchwater-nf) | reimplemented locally, `subworkflows/local/branchwater` | — |
| `compare` | local, `subworkflows/local/compare` | — |

## Contributions and Support

If you would like to contribute to this pipeline, please see the
[contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the
[`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the
[nf-core](https://nf-co.re) community, reused here under the
[MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg,
> Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
