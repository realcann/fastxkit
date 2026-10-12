"""Parser behaviour: FASTA/FASTQ detection, edge cases and errors."""
import io

import pytest

from fastxkit import records


def parse(text):
    reader = records.FastxReader(io.StringIO(text))
    return reader.format, list(reader)


def test_multiline_fasta_is_joined():
    fmt, recs = parse(">a desc\nACGT\nTT\n>b\nGG\n")
    assert fmt == "FASTA"
    assert [(r.name, r.seq, r.qual) for r in recs] == [
        ("a desc", "ACGTTT", None), ("b", "GG", None)]


def test_fastq_four_line_records():
    fmt, recs = parse("@r1 x\nACGT\n+\nIIII\n@r2\nGG\n+r2\nI#\n")
    assert fmt == "FASTQ"
    assert [(r.name, r.seq, r.qual) for r in recs] == [
        ("r1 x", "ACGT", "IIII"), ("r2", "GG", "I#")]


def test_id_is_first_word_of_name():
    _, recs = parse(">seq1 sample record\nA\n")
    assert recs[0].id == "seq1"


def test_empty_input_has_no_format_and_no_records():
    assert parse("") == (None, [])
    assert parse("\n\n  \n") == (None, [])


def test_crlf_line_endings():
    _, fasta = parse(">a\r\nAC\r\nGT\r\n")
    assert (fasta[0].name, fasta[0].seq) == ("a", "ACGT")
    _, fastq = parse("@r\r\nAC\r\n+\r\nII\r\n")
    assert (fastq[0].name, fastq[0].seq, fastq[0].qual) == ("r", "AC", "II")


def test_blank_lines_are_skipped():
    _, recs = parse("\n>a\n\nAC\n\n>b\nG\n\n")
    assert [r.seq for r in recs] == ["AC", "G"]
    _, recs = parse("\n@r\nAC\n+\nII\n\n@s\nG\n+\nI\n\n")
    assert [r.seq for r in recs] == ["AC", "G"]


def test_no_trailing_newline():
    _, recs = parse(">a\nACGT")
    assert recs[0].seq == "ACGT"
    _, recs = parse("@r\nAC\n+\nII")
    assert recs[0].qual == "II"


def test_fasta_header_without_sequence():
    _, recs = parse(">a\n>b\nAC\n")
    assert [(r.name, r.seq) for r in recs] == [("a", ""), ("b", "AC")]


def test_fastq_empty_sequence_record():
    _, recs = parse("@r\n\n+\n\n@s\nA\n+\nI\n")
    assert [(r.seq, r.qual) for r in recs] == [("", ""), ("A", "I")]


def test_fastq_quality_may_start_with_at_sign():
    _, recs = parse("@r\nAC\n+\n@@\n")
    assert recs[0].qual == "@@"


def test_fastq_length_mismatch_is_an_error():
    reader = records.FastxReader(io.StringIO("@r\nACGT\n+\nIII\n"), "x.fq")
    with pytest.raises(records.FormatError, match=r"x\.fq.*line 1.*4.*3"):
        list(reader)


def test_error_reports_record_after_valid_ones():
    reader = records.FastxReader(
        io.StringIO("@a\nA\n+\nI\n@b\nAC\n+\nI\n"))
    it = iter(reader)
    assert next(it).name == "a"
    with pytest.raises(records.FormatError, match="line 5"):
        next(it)


def test_data_before_any_header_is_an_error():
    with pytest.raises(records.FormatError, match="line 1"):
        parse("ACGT\n>a\nAC\n")


def test_fastq_missing_header_marker_is_an_error():
    with pytest.raises(records.FormatError, match="'@'"):
        parse("@a\nA\n+\nI\nB\nA\n+\nI\n")


def test_fastq_missing_plus_line_is_an_error():
    with pytest.raises(records.FormatError, match="'\\+'"):
        parse("@a\nA\nX\nI\n")


def test_truncated_fastq_is_an_error():
    with pytest.raises(records.FormatError, match="truncated"):
        parse("@a\nACGT\n+\n")


def test_parser_is_lazy():
    def lines():
        yield ">a\n"
        yield "AC\n"
        yield ">b\n"
        raise AssertionError("read past the first record")

    it = iter(records.FastxReader(lines()))
    assert next(it).seq == "AC"
