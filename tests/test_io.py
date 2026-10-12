"""Input layer: content-based decompression, stdin and failure modes."""
import bz2
import gzip
import lzma

import pytest

from fastxkit import io as fxio

FASTA = ">a\nACGT\n>b\nGG\n"


def read_all(path):
    with fxio.open_text(path) as handle:
        return handle.read()


@pytest.mark.parametrize("compress", [gzip.compress, bz2.compress, lzma.compress])
def test_compression_detected_from_content_not_name(make_file, compress):
    path = make_file("reads.txt", compress(FASTA.encode()))
    assert read_all(path) == FASTA


def test_plain_text_with_compressed_looking_name(make_file):
    path = make_file("reads.fa.gz", FASTA)
    assert read_all(path) == FASTA


def test_concatenated_gzip_members(make_file):
    data = gzip.compress(b">a\nAC\n") + gzip.compress(b">b\nGT\n")
    assert read_all(make_file("x.bin", data)) == ">a\nAC\n>b\nGT\n"


def test_empty_file(make_file):
    assert read_all(make_file("empty.fa", "")) == ""


def test_crlf_is_translated(make_file):
    assert read_all(make_file("w.fa", ">a\r\nAC\r\n")) == ">a\nAC\n"


def test_stdin_dash(run_cli):
    code, out, _ = run_cli("seq", "-s", "-", stdin=b">a\nACGT\n")
    assert (code, out) == (0, "ACGT\n")


def test_stdin_when_no_file_given(run_cli):
    code, out, _ = run_cli("seq", "-s", stdin=b">a\nACGT\n")
    assert (code, out) == (0, "ACGT\n")


def test_compressed_stdin(run_cli):
    code, out, _ = run_cli("seq", "-s", stdin=gzip.compress(b">a\nACGT\n"))
    assert (code, out) == (0, "ACGT\n")


def test_missing_file_is_reported_with_exit_code_1(run_cli):
    code, out, err = run_cli("seq", "/nonexistent/in.fa")
    assert (code, out) == (1, "")
    assert "/nonexistent/in.fa" in err


def test_truncated_gzip_fails_cleanly(run_cli, make_file):
    path = make_file("t.gz", gzip.compress(b">a\n" + b"ACGT" * 5000)[:-20])
    code, _, err = run_cli("seq", "-s", path)
    assert code == 1 and err.startswith("fastxkit:")


def test_unrecognised_content_exit_code_1(run_cli, make_file):
    code, out, err = run_cli("stats", make_file("x.txt", "hello\n"))
    assert code == 1 and "fastxkit:" in err


def test_malformed_fastq_exit_code_1_and_message(run_cli, make_file):
    path = make_file("bad.fq", "@r\nACGT\n+\nII\n")
    code, _, err = run_cli("seq", path)
    assert code == 1
    assert "bad.fq" in err and "line 1" in err


def test_records_before_error_are_still_written(run_cli, make_file):
    path = make_file("bad.fq", "@a\nA\n+\nI\n@b\nAC\n+\nI\n")
    code, out, _ = run_cli("seq", path)
    assert code == 1 and out == "@a\nA\n+\nI\n"


@pytest.mark.parametrize("compress", [lambda data: data, gzip.compress])
def test_utf8_bom_is_skipped(make_file, compress):
    path = make_file("bom.fa", compress(b"\xef\xbb\xbf>a\nAC\n"))
    assert read_all(path) == ">a\nAC\n"


def test_utf8_bom_on_stdin_for_fasta_and_fastq(run_cli):
    code, out, _ = run_cli("seq", "-s", stdin=b"\xef\xbb\xbf>a\nACGT\n")
    assert (code, out) == (0, "ACGT\n")
    code, out, _ = run_cli("seq", "-s", stdin=b"\xef\xbb\xbf@a\nACGT\n+\nIIII\n")
    assert (code, out) == (0, "ACGT\n")


@pytest.mark.parametrize("name, data", [
    ("trunc.fa.gz", gzip.compress(b">a\n" + b"ACGT" * 5000 + b"\n")[:-30]),
    ("bad.fa.bz2", b"BZh9" + b"\x00" * 40),
    ("bad.fa.xz", b"\xfd7zXZ\x00" + b"\x01" * 40),
])
def test_corrupt_compressed_input_names_the_file(run_cli, make_file, name, data):
    path = make_file(name, data)
    code, out, err = run_cli("seq", "-s", path)
    assert code == 1 and path in err


def test_corrupt_stdin_is_reported_as_stdin(run_cli):
    code, _, err = run_cli("seq", "-s", stdin=gzip.compress(b">a\nAC\n")[:-8])
    assert code == 1 and "stdin" in err
