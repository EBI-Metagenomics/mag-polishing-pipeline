"""Unit tests for bin/ena_fetch_fastq.py."""

import pytest

from ena_fetch_fastq import download, md5sum

CONTENT = b"mag\n"
CONTENT_MD5 = "8b74d652835268d0413b696e4f080729"


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "source.txt"
    path.write_bytes(CONTENT)
    return path


class TestMd5sum:
    def test_known_digest(self, source):
        assert md5sum(source) == CONTENT_MD5

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty"
        path.write_bytes(b"")
        assert md5sum(path) == "d41d8cd98f00b204e9800998ecf8427e"


class TestDownload:
    def test_writes_the_file_and_verifies_the_md5(self, tmp_path, source):
        destination = tmp_path / "copy.txt"
        download(f"file://{source}", destination, CONTENT_MD5)
        assert destination.read_bytes() == CONTENT

    def test_md5_mismatch_aborts(self, tmp_path, source):
        with pytest.raises(SystemExit, match="md5 mismatch"):
            download(f"file://{source}", tmp_path / "bad.txt", "0" * 32, attempts=1)

    def test_md5_mismatch_leaves_no_partial_file(self, tmp_path, source):
        """The download goes to a temp file and is only moved into place once verified,
        so a failure must not leave a half-written cache entry behind - a storeDir would
        treat it as a completed download forever after."""
        destination = tmp_path / "bad.txt"
        with pytest.raises(SystemExit):
            download(f"file://{source}", destination, "0" * 32, attempts=1)
        assert not destination.exists()

    def test_md5_mismatch_cleans_up_every_temp_file(self, tmp_path, source):
        destination = tmp_path / "bad.txt"
        with pytest.raises(SystemExit):
            download(f"file://{source}", destination, "0" * 32, attempts=3)
        # only the fixture remains, no leftover mkstemp files
        assert [p.name for p in tmp_path.iterdir()] == ["source.txt"]

    def test_no_md5_skips_verification(self, tmp_path, source):
        destination = tmp_path / "copy.txt"
        download(f"file://{source}", destination)
        assert destination.read_bytes() == CONTENT

    def test_retries_then_gives_up(self, tmp_path, monkeypatch):
        attempts = []

        def failing_urlopen(*args, **kwargs):
            attempts.append(1)
            raise OSError("boom")

        monkeypatch.setattr("ena_fetch_fastq.urllib.request.urlopen", failing_urlopen)
        with pytest.raises(SystemExit, match="could not download"):
            download("ftp://nowhere/x.gz", tmp_path / "out.gz", attempts=3)
        assert len(attempts) == 3

    def test_bare_host_path_is_given_an_ftp_scheme(self, tmp_path, monkeypatch):
        """ENA's fastq_ftp field has no scheme, so the script has to add one."""
        seen = {}

        def fake_urlopen(url, *args, **kwargs):
            seen["url"] = url
            raise OSError("stop here")

        monkeypatch.setattr("ena_fetch_fastq.urllib.request.urlopen", fake_urlopen)
        with pytest.raises(SystemExit):
            download("ftp.sra.ebi.ac.uk/vol1/x_1.fastq.gz", tmp_path / "out.gz", attempts=1)
        assert seen["url"] == "ftp://ftp.sra.ebi.ac.uk/vol1/x_1.fastq.gz"
