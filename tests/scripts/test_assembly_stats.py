"""Unit tests for bin/assembly_stats.py and bin/precompute_assembly_stats.py.

These two produce the length / N50 / GC / contig-count columns of
`compare/assembly_qc_metrics.tsv`, which is the pipeline's deliverable, so the arithmetic
is worth pinning down.
"""

import csv
from pathlib import Path

import pytest

from assembly_stats import calculate_stats, fasta_iter, read_genome, run_assembly_stats

FIXTURES = Path(__file__).parent / "fixtures" / "assembly_stats"
TINY = FIXTURES / "tiny.fa"

# contig_1 ACGTACGTAA   len 10, 4 G/C
# contig_2 GGGGCCCCGGGGCCCC  len 16, 16 G/C
# contig_3 ATATATAT     len  8, 0 G/C
# total 34 bp, 20 G/C -> 58.82%


class TestFastaIter:
    def test_yields_header_and_sequence(self):
        records = list(fasta_iter(TINY))
        assert [h for h, _ in records] == ["contig_1", "contig_2", "contig_3"]
        assert records[0][1] == "ACGTACGTAA"

    def test_sequence_is_uppercased(self, tmp_path):
        path = tmp_path / "lower.fa"
        path.write_text(">c\nacgt\n")
        assert list(fasta_iter(path))[0][1] == "ACGT"

    def test_multiline_sequence_is_joined(self, tmp_path):
        path = tmp_path / "wrapped.fa"
        path.write_text(">c\nACGT\nACGT\n")
        assert list(fasta_iter(path))[0][1] == "ACGTACGT"


class TestReadGenome:
    def test_contig_lengths_and_gc(self):
        contig_lens, gc = read_genome(TINY)
        assert contig_lens == [10, 16, 8]
        assert gc == pytest.approx(58.82)

    def test_whole_number_gc_is_returned_as_int(self, tmp_path):
        path = tmp_path / "half.fa"
        path.write_text(">c\nACGT\n")  # 2 of 4 -> exactly 50%
        _, gc = read_genome(path)
        assert gc == 50
        assert isinstance(gc, int)


class TestCalculateStats:
    def test_tiny_assembly(self):
        stats = calculate_stats([10, 16, 8], 58.82)
        assert stats["N_contigs"] == 3
        assert stats["Length"] == 34
        assert stats["GC_content"] == 58.82
        # sorted desc 16,10,8; cumsum 16,26,34; half of 34 is 17 -> first cumsum >= 17 is 26
        assert stats["N50"] == 10

    def test_single_contig_n50_is_its_length(self):
        assert calculate_stats([100], 50)["N50"] == 100

    def test_equal_contigs(self):
        stats = calculate_stats([50, 50], 40)
        assert stats["Length"] == 100
        assert stats["N50"] == 50
        assert stats["N_contigs"] == 2


class TestRunAssemblyStats:
    def test_end_to_end_on_the_fixture(self):
        stats = run_assembly_stats(TINY)
        assert stats == {
            "N_contigs": 3,
            "GC_content": pytest.approx(58.82),
            "Length": 34,
            "N50": 10,
        }


class TestPrecomputeAssemblyStats:
    def test_writes_one_row_per_genome_with_the_accession_stripped(self, tmp_path):
        from precompute_assembly_stats import main

        genomes = tmp_path / "genomes"
        genomes.mkdir()
        (genomes / "SAMPLE_001.fa").write_text(TINY.read_text())
        (genomes / "SAMPLE_001_n5_concoct_3.fna").write_text(">c\nACGT\n")

        out = tmp_path / "stats.tsv"
        main(str(genomes), str(out))

        rows = {r["Genome"]: r for r in csv.DictReader(out.open(), delimiter="\t")}
        assert set(rows) == {"SAMPLE_001", "SAMPLE_001_n5_concoct_3"}
        assert rows["SAMPLE_001"]["Length"] == "34"
        assert rows["SAMPLE_001"]["N_contigs"] == "3"
        # the basename is what compare/ keys every metric back to
        assert rows["SAMPLE_001_n5_concoct_3"]["Length"] == "4"
