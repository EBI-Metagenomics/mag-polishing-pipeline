#!/usr/bin/env python3
"""Download one fastq from ENA and verify its md5.

Used with a `storeDir` cache, so every run accession is fetched once and reused across
every --n_concat_samples value and every sample.
"""

import argparse
import hashlib
import os
import shutil
import sys
import tempfile
import urllib.request

CHUNK = 1 << 20


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, destination, expected_md5=None, attempts=3):
    """Download to a temp file next to the destination, verify, then move into place."""
    if "://" not in url:
        url = f"ftp://{url}"

    last_error = None
    for attempt in range(1, attempts + 1):
        handle, temporary = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(destination)))
        os.close(handle)
        try:
            with urllib.request.urlopen(url, timeout=300) as response, open(temporary, "wb") as out:
                shutil.copyfileobj(response, out, CHUNK)
            if expected_md5:
                observed = md5sum(temporary)
                if observed != expected_md5:
                    raise ValueError(f"md5 mismatch: expected {expected_md5}, got {observed}")
            os.replace(temporary, destination)
            return destination
        except Exception as error:
            last_error = error
            print(f"attempt {attempt}/{attempts} failed for {url}: {error}", file=sys.stderr)
            if os.path.exists(temporary):
                os.remove(temporary)

    raise SystemExit(f"could not download {url}: {last_error}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url")
    parser.add_argument("--output")
    parser.add_argument("--md5", default=None)
    args = parser.parse_args()

    if not args.url or not args.output:
        parser.error("--url and --output are required")

    download(args.url, args.output, args.md5)


if __name__ == "__main__":
    main()
