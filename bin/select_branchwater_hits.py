#!/usr/bin/env python3
"""Pick the public runs to co-assemble with a sample.

Joins the Branchwater `manysearch` hits (containment, query_containment_ani) with the
Branchwater metadata dump (acc, librarylayout) on hits.match_name == metadata.acc, drops
the self hit, keeps paired-end runs, ranks by containment and then verifies - through the
ENA portal API - that the fastqs are actually downloadable before emitting the top N.

The self hit is dropped by accession (--exclude, the sample id). When that accession is
not among the hits - the samplesheet id is not the SRA run, or the run is not in the
index - it falls back to dropping whatever scores >= 0.99 on both containment and ANI. A
MAG never contains 100% of its own run's hashes (it is an assembly of it, not it), so the
old 0.9999-on-both test never fired and the self hit ate one of the N slots.

The Branchwater index is SRA derived, so a hit is not necessarily in ENA; that is why the
availability check is here and not in the download step.
"""

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
SELF_HIT_THRESHOLD = 0.99
OUTPUT_FIELDS = ["run_accession", "fastq_1", "fastq_2", "md5_1", "md5_2", "containment"]


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def rank_candidates(hits, metadata, exclude=None, self_hit_threshold=SELF_HIT_THRESHOLD):
    """[{acc, containment}] ranked by containment desc, self hit and non-paired dropped."""
    layout_by_acc = {row["acc"]: (row.get("librarylayout") or "").upper() for row in metadata}
    accessions = [hit["match_name"].split()[0] for hit in hits]
    by_accession = exclude in accessions

    candidates = []
    for acc, hit in zip(accessions, hits):
        containment = float(hit.get("containment") or 0)
        ani = float(hit.get("query_containment_ani") or 0)
        if by_accession:
            if acc == exclude:
                continue  # the sample's own run
        elif containment >= self_hit_threshold and ani >= self_hit_threshold:
            continue  # no usable accession, so score it instead
        if layout_by_acc.get(acc) != "PAIRED":
            continue
        candidates.append({"acc": acc, "containment": containment})

    candidates.sort(key=lambda c: c["containment"], reverse=True)
    return candidates


def ena_read_run(accession):
    """ENA portal record for a run, or None when it has no paired fastqs on the FTP."""
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_md5",
            "format": "json",
        }
    )
    try:
        with urllib.request.urlopen(f"{ENA_PORTAL}?{query}", timeout=60) as response:
            records = json.load(response)
    except Exception as error:  # noqa: BLE001 - an unavailable run is not a failure
        print(f"{accession}: ENA lookup failed ({error})", file=sys.stderr)
        return None

    if not records:
        return None

    ftp = [url for url in (records[0].get("fastq_ftp") or "").split(";") if url]
    md5 = [value for value in (records[0].get("fastq_md5") or "").split(";") if value]
    if len(ftp) != len(md5) or len(ftp) < 2:
        return None

    # a run with 3 files is _1/_2 plus the unpaired leftovers, which we do not want
    paired = [(url, value) for url, value in zip(ftp, md5) if url.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    if len(paired) != 2:
        return None
    paired.sort()
    return {
        "fastq_1": f"ftp://{paired[0][0]}" if not paired[0][0].startswith("ftp") else paired[0][0],
        "fastq_2": f"ftp://{paired[1][0]}" if not paired[1][0].startswith("ftp") else paired[1][0],
        "md5_1": paired[0][1],
        "md5_2": paired[1][1],
    }


def select(hits, metadata, n_runs, exclude=None, lookup=ena_read_run):
    selected = []
    for candidate in rank_candidates(hits, metadata, exclude=exclude):
        if len(selected) == n_runs:
            break
        record = lookup(candidate["acc"])
        if record is None:
            print(f"{candidate['acc']}: no paired fastqs in ENA, skipped", file=sys.stderr)
            continue
        selected.append(
            dict(record, run_accession=candidate["acc"], containment=candidate["containment"])
        )
    return selected


def write_csv(rows, path):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", help="manysearch csv")
    parser.add_argument("--metadata", help="branchwater metadata csv")
    parser.add_argument("--n-runs", type=int, help="how many runs to emit")
    parser.add_argument("--exclude", help="the sample's own run accession, dropped from the hits")
    parser.add_argument("--output", help="output csv")
    args = parser.parse_args()

    for required in ("hits", "metadata", "n_runs", "output"):
        if getattr(args, required) is None:
            parser.error(f"--{required.replace('_', '-')} is required")

    rows = select(read_csv(args.hits), read_csv(args.metadata), args.n_runs, exclude=args.exclude)
    if len(rows) < args.n_runs:
        print(
            f"warning: only {len(rows)} usable runs for {args.n_runs} requested",
            file=sys.stderr,
        )
    write_csv(rows, args.output)


if __name__ == "__main__":
    main()
