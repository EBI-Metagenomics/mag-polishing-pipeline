# `bin/`

Scripts the pipeline's processes call by name.

| file | origin | how |
| --- | --- | --- |
| `assembly_stats.py` | genomes-catalogue-pipeline | copy |
| `precompute_assembly_stats.py` | genomes-catalogue-pipeline | copy |
| `bat_taxo_process.py` | genomes-generation | symlink |
| `cov_recycler.py` | genomes-generation | symlink |
| `create_qc_table.py` | genomes-generation | symlink |
| `detect_contamination.py` | genomes-generation | symlink |
| `ggp_binlinks.py` | genomes-generation | symlink |
| `logging_stats.py` | genomes-generation | symlink |
| `propagate_taxonomy_to_bins.py` | genomes-generation | symlink |
| `select_seqs_notin_ids.py` | genomes-generation | symlink |
| `calculate_assembly_coverage.py` | miassembler | symlink |

## Symlinks

Scripts that belong to the two pipelines under `pipelines/`, which call them while running
as part of this one. They have to be reachable from here, so each is a link to the real
file in its submodule.

Only the top-level pipeline's `bin/` is visible to running processes, so a script left in a
submodule's own `bin/` is never found.
