"""faidx command."""
import os

import pytest

from fastxkit.commands import faidx

FA = ">s1 first\nACGTAC\nGTAC\n>s2\nTTTT\n"


def fai(path):
    with open(path + ".fai", encoding="utf-8") as stream:
        return stream.read()


def test_index_columns_are_samtools_style(run_cli, make_file):
    path = make_file("a.fa", FA)
    assert run_cli("faidx", path)[:2] == (0, "")
    assert fai(path) == "s1\t10\t10\t6\t7\ns2\t4\t26\t4\t5\n"


def test_index_with_crlf_counts_the_longer_line_width(run_cli, make_file):
    path = make_file("a.fa", FA.replace("\n", "\r\n"))
    run_cli("faidx", path)
    assert fai(path) == "s1\t10\t11\t6\t8\ns2\t4\t30\t4\t6\n"


def test_index_without_final_newline(run_cli, make_file):
    path = make_file("a.fa", ">a\nACGT\nAC")
    run_cli("faidx", path)
    assert fai(path) == "a\t6\t3\t4\t5\n"


def test_index_of_single_unterminated_line(run_cli, make_file):
    path = make_file("a.fa", ">a\nACGT")
    run_cli("faidx", path)
    assert fai(path) == "a\t4\t3\t4\t4\n"


def test_empty_sequence_record(run_cli, make_file):
    path = make_file("a.fa", ">e\n>f\nAC\n")
    run_cli("faidx", path)
    assert fai(path) == "e\t0\t3\t0\t0\nf\t2\t6\t2\t3\n"


def test_empty_file_gives_empty_index(run_cli, make_file):
    path = make_file("a.fa", "")
    assert run_cli("faidx", path)[0] == 0
    assert fai(path) == ""


def test_trailing_blank_lines_are_tolerated(run_cli, make_file):
    path = make_file("a.fa", ">a\nACGT\n\n\n")
    run_cli("faidx", path)
    assert fai(path) == "a\t4\t3\t4\t5\n"


@pytest.mark.parametrize("content", [
    ">a\nACG\nACGT\n",      # a later line is longer than the first
    ">a\nACGT\nAC\nACGT\n",  # a short line before the end
    ">a\nAC\n\nAC\n",       # blank line inside the sequence
])
def test_uneven_line_lengths_are_an_error(run_cli, make_file, content):
    code, _, err = run_cli("faidx", make_file("a.fa", content))
    assert code == 1
    assert "line" in err


def test_duplicate_names_are_an_error(run_cli, make_file):
    code, _, err = run_cli("faidx", make_file("a.fa", ">a\nA\n>a x\nC\n"))
    assert code == 1
    assert "duplicate" in err


def test_non_fasta_input_is_an_error(run_cli, make_file):
    assert run_cli("faidx", make_file("a.fq", "@r\nA\n+\nI\n"))[0] == 1
    assert run_cli("faidx", make_file("a.gz", b"\x1f\x8b\x08\x00junk"))[0] == 1


def test_missing_file_is_an_error(run_cli):
    assert run_cli("faidx", "/nonexistent/a.fa")[0] == 1


def test_extract_whole_sequence_and_region(run_cli, make_file):
    path = make_file("a.fa", FA)
    out = run_cli("faidx", path, "s1", "s1:5-8", "s2:1-1")[1]
    assert out == ">s1\nACGTACGTAC\n>s1:5-8\nACGT\n>s2:1-1\nT\n"


def test_extract_across_crlf_lines(run_cli, make_file):
    path = make_file("a.fa", FA.replace("\n", "\r\n"))
    assert run_cli("faidx", path, "s1:6-8")[1] == ">s1:6-8\nCGT\n"


def test_extract_edge_coordinates(run_cli, make_file):
    path = make_file("a.fa", FA)
    assert run_cli("faidx", path, "s1:1-10")[1] == ">s1:1-10\nACGTACGTAC\n"
    assert run_cli("faidx", path, "s1:10-10")[1] == ">s1:10-10\nC\n"
    assert run_cli("faidx", path, "s1:7-6")[0] == 1


def test_extract_end_beyond_length_is_clamped(run_cli, make_file):
    path = make_file("a.fa", FA)
    assert run_cli("faidx", path, "s1:9-99")[1] == ">s1:9-99\nAC\n"


def test_extract_builds_a_missing_index(run_cli, make_file):
    path = make_file("a.fa", FA)
    run_cli("faidx", path, "s2")
    assert fai(path) == "s1\t10\t10\t6\t7\ns2\t4\t26\t4\t5\n"


def set_mtime(path, seconds):
    os.utime(path, (seconds, seconds))


def test_extract_uses_a_fresh_existing_index(run_cli, make_file):
    path = make_file("a.fa", FA)
    index = make_file("a.fa.fai", "x\t4\t26\t4\t5\n")
    set_mtime(path, 1_000_000)
    set_mtime(index, 2_000_000)
    assert run_cli("faidx", path, "x")[1] == ">x\nTTTT\n"
    assert fai(path) == "x\t4\t26\t4\t5\n"


def test_older_index_is_rebuilt_silently(run_cli, make_file):
    path = make_file("a.fa", FA)
    index = make_file("a.fa.fai", "x\t4\t26\t4\t5\n")
    set_mtime(index, 1_000_000)
    set_mtime(path, 2_000_000)
    code, out, err = run_cli("faidx", path, "s2")
    assert (code, out, err) == (0, ">s2\nTTTT\n", "")
    assert fai(path) == "s1\t10\t10\t6\t7\ns2\t4\t26\t4\t5\n"


def test_index_reaching_beyond_the_file_is_rebuilt(run_cli, make_file):
    path = make_file("a.fa", FA)
    index = make_file("a.fa.fai", "s2\t4000\t26\t4\t5\n")
    set_mtime(path, 1_000_000)
    set_mtime(index, 2_000_000)
    assert run_cli("faidx", path, "s2")[1] == ">s2\nTTTT\n"
    assert fai(path) == "s1\t10\t10\t6\t7\ns2\t4\t26\t4\t5\n"


def test_symlinked_index_is_an_error_and_target_is_untouched(
        run_cli, make_file, tmp_path):
    path = make_file("a.fa", FA)
    victim = make_file("victim.txt", "keep me")
    (tmp_path / "a.fa.fai").symlink_to(victim)
    code, _, err = run_cli("faidx", path)
    assert code == 1 and "symbolic link" in err
    assert (tmp_path / "victim.txt").read_text() == "keep me"


def test_regular_index_is_rewritten_in_place(run_cli, make_file):
    path = make_file("a.fa", FA)
    make_file("a.fa.fai", "junk that is much longer than the real index " * 5)
    assert run_cli("faidx", path)[0] == 0
    assert fai(path) == "s1\t10\t10\t6\t7\ns2\t4\t26\t4\t5\n"


def test_name_containing_a_colon(run_cli, make_file):
    path = make_file("a.fa", ">chr1:1-5\nACGTA\n")
    assert run_cli("faidx", path, "chr1:1-5")[1] == ">chr1:1-5\nACGTA\n"
    assert run_cli("faidx", path, "chr1:1-5:2-3")[1] == ">chr1:1-5:2-3\nCG\n"


def test_extract_wraps_at_line_width(run_cli, make_file):
    path = make_file("a.fa", ">a\n" + "A" * 70 + "\n")
    out = run_cli("faidx", "-w", "30", path, "a")[1]
    assert [len(line) for line in out.splitlines()] == [2, 30, 30, 10]


@pytest.mark.parametrize("region", [
    "nope", "nope:1-2", "s1:0-5", "s1:11-12", "s1:5", "s1:a-b", "s1:",
])
def test_unknown_name_or_bad_region_is_an_error(run_cli, make_file, region):
    code, out, err = run_cli("faidx", make_file("a.fa", FA), region)
    assert (code, out) == (1, "")
    assert "not found" in err or "region" in err


def test_malformed_index_is_an_error(run_cli, make_file):
    path = make_file("a.fa", FA)
    make_file("a.fa.fai", "s1\t10\n")
    code, _, err = run_cli("faidx", path, "s1")
    assert code == 1
    assert "index" in err


def test_parse_region():
    index = {"a": faidx.IndexEntry("a", 100, 3, 60, 61)}
    assert faidx.parse_region("a", index) == (index["a"], 1, 100)
    assert faidx.parse_region("a:5-200", index) == (index["a"], 5, 100)
