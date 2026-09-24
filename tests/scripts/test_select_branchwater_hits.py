"""Unit tests for bin/select_branchwater_hits.py."""

import pytest

from select_branchwater_hits import rank_candidates, read_csv, select, write_csv

# SRR000001 scores 1.0/1.0 - the shape of a real self hit.
# SRR000003 is SINGLE, so it must never be a candidate.
HITS = [
    {"match_name": "SRR000001", "containment": "1.0", "query_containment_ani": "1.0"},
    {"match_name": "SRR000002", "containment": "0.5", "query_containment_ani": "0.98"},
    {"match_name": "SRR000003", "containment": "0.9", "query_containment_ani": "0.99"},
    {"match_name": "SRR000004", "containment": "0.7", "query_containment_ani": "0.99"},
]
METADATA = [
    {"acc": "SRR000001", "librarylayout": "PAIRED"},
    {"acc": "SRR000002", "librarylayout": "PAIRED"},
    {"acc": "SRR000003", "librarylayout": "SINGLE"},
    {"acc": "SRR000004", "librarylayout": "paired"},  # case must not matter
]


def accs(candidates):
    return [c["acc"] for c in candidates]


class TestRankCandidates:
    def test_drops_self_hit_by_score_when_no_exclude(self):
        """With no --exclude the 1.0/1.0 hit is scored out, and SINGLE is dropped."""
        assert accs(rank_candidates(HITS, METADATA)) == ["SRR000004", "SRR000002"]

    def test_exclude_not_among_hits_falls_back_to_scoring(self):
        assert accs(rank_candidates(HITS, METADATA, exclude="SRRABSENT")) == [
            "SRR000004",
            "SRR000002",
        ]

    def test_exclude_among_hits_drops_by_accession_only(self):
        """Once --exclude matches a hit, the score heuristic is not applied at all: the
        only row dropped is that accession. SRR000001 therefore survives despite scoring
        1.0/1.0, which is what this pins down.

        The accession takes over because the score test cannot find the real self hit: a
        MAG is an assembly of its run, not the run, so it contains ~0.999 of it rather
        than 1.0. The old threshold never fired and the self hit ate one of the N slots.
        """
        assert accs(rank_candidates(HITS, METADATA, exclude="SRR000002")) == [
            "SRR000001",
            "SRR000004",
        ]

    def test_score_fallback_drops_only_the_first_hit_over_the_threshold(self):
        """There is exactly one self hit. A second run scoring >= 0.99 on both is a
        different run that happens to be near identical, and is a real candidate - it
        must survive. manysearch emits hits highest containment first, so the first row
        over the threshold is the self hit."""
        hits = [
            {"match_name": "SELF", "containment": "1.0", "query_containment_ani": "1.0"},
            {"match_name": "TWIN", "containment": "0.995", "query_containment_ani": "0.999"},
            {"match_name": "NEAR", "containment": "0.99", "query_containment_ani": "0.99"},
            {"match_name": "FAR", "containment": "0.7", "query_containment_ani": "0.99"},
        ]
        metadata = [{"acc": a, "librarylayout": "PAIRED"} for a in ("SELF", "TWIN", "NEAR", "FAR")]
        assert accs(rank_candidates(hits, metadata)) == ["TWIN", "NEAR", "FAR"]

    def test_score_fallback_drops_nothing_when_no_hit_reaches_the_threshold(self):
        hits = [
            {"match_name": "SRR000002", "containment": "0.5", "query_containment_ani": "0.98"},
            {"match_name": "SRR000004", "containment": "0.7", "query_containment_ani": "0.99"},
        ]
        assert accs(rank_candidates(hits, METADATA)) == ["SRR000004", "SRR000002"]

    def test_score_fallback_needs_both_containment_and_ani_over_the_threshold(self):
        """0.99 containment with a low ANI is not a self hit, so nothing is dropped."""
        hits = [{"match_name": "SRR000002", "containment": "1.0", "query_containment_ani": "0.5"}]
        assert accs(rank_candidates(hits, METADATA)) == ["SRR000002"]

    def test_ranked_by_containment_descending(self):
        containments = [c["containment"] for c in rank_candidates(HITS, METADATA)]
        assert containments == sorted(containments, reverse=True)

    def test_library_layout_is_case_insensitive(self):
        assert "SRR000004" in accs(rank_candidates(HITS, METADATA))

    def test_missing_metadata_row_drops_the_hit(self):
        assert accs(rank_candidates(HITS, [])) == []

    def test_absent_containment_is_treated_as_zero(self):
        hits = [{"match_name": "SRR000002", "containment": "", "query_containment_ani": ""}]
        assert rank_candidates(hits, METADATA) == [{"acc": "SRR000002", "containment": 0.0}]

    def test_match_name_is_split_on_whitespace(self):
        hits = [{"match_name": "SRR000002 some description", "containment": "0.5"}]
        assert accs(rank_candidates(hits, METADATA)) == ["SRR000002"]


class TestSelect:
    def test_skips_runs_missing_from_ena_and_takes_the_next(self):
        def lookup(acc):
            if acc == "SRR000004":
                return None  # ranked first, but not downloadable
            return {"fastq_1": "ftp://x_1.gz", "fastq_2": "ftp://x_2.gz", "md5_1": "a", "md5_2": "b"}

        picked = select(HITS, METADATA, 1, lookup=lookup)
        assert [row["run_accession"] for row in picked] == ["SRR000002"]
        assert picked[0]["containment"] == 0.5

    def test_zero_runs_selects_nothing(self):
        assert select(HITS, METADATA, 0, lookup=lambda _: {}) == []

    def test_stops_at_n_runs(self):
        lookup = lambda _: {"fastq_1": "a", "fastq_2": "b", "md5_1": "c", "md5_2": "d"}
        assert len(select(HITS, METADATA, 1, lookup=lookup)) == 1

    def test_returns_fewer_than_requested_when_nothing_is_available(self):
        assert select(HITS, METADATA, 5, lookup=lambda _: None) == []


class TestCsvRoundTrip:
    def test_write_then_read(self, tmp_path):
        rows = [
            {
                "run_accession": "SRR000002",
                "fastq_1": "ftp://a_1.gz",
                "fastq_2": "ftp://a_2.gz",
                "md5_1": "aaa",
                "md5_2": "bbb",
                "containment": 0.5,
            }
        ]
        path = tmp_path / "runs.csv"
        write_csv(rows, path)
        assert read_csv(path) == [{k: str(v) for k, v in rows[0].items()}]
