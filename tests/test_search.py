"""grep and locate commands."""
import pytest

FA = ">s1 alpha one\nACGTNN\n>s10 beta\nacgt\n>s2 alpha two\nTTTT\n"


def ids(out):
    return [line[1:].split()[0] for line in out.splitlines() if line[0] == ">"]


def locate_rows(out):
    return [line.split("\t") for line in out.splitlines()[1:]]


def test_grep_id_is_an_exact_match(run_cli):
    assert ids(run_cli("grep", "-p", "s1", stdin=FA.encode())[1]) == ["s1"]


def test_grep_multiple_patterns_match_any(run_cli):
    out = run_cli("grep", "-p", "s1", "-p", "s2", stdin=FA.encode())[1]
    assert ids(out) == ["s1", "s2"]


def test_grep_by_name_uses_the_full_name(run_cli):
    out = run_cli("grep", "-n", "-p", "s2 alpha two", stdin=FA.encode())[1]
    assert ids(out) == ["s2"]
    assert run_cli("grep", "-n", "-p", "s2", stdin=FA.encode())[1] == ""


def test_grep_by_sequence_is_a_substring_match(run_cli):
    out = run_cli("grep", "-s", "-p", "GTNN", stdin=FA.encode())[1]
    assert ids(out) == ["s1"]


def test_grep_ignore_case(run_cli):
    out = run_cli("grep", "-s", "-i", "-p", "ACGT", stdin=FA.encode())[1]
    assert ids(out) == ["s1", "s10"]
    assert ids(run_cli("grep", "-i", "-p", "S2", stdin=FA.encode())[1]) == ["s2"]


def test_grep_regex_is_a_partial_match(run_cli):
    out = run_cli("grep", "-r", "-p", "^s1", stdin=FA.encode())[1]
    assert ids(out) == ["s1", "s10"]


def test_grep_special_characters_are_literal_without_regex(run_cli):
    data = b">a.b\nA\n>aXb\nA\n"
    assert ids(run_cli("grep", "-p", "a.b", stdin=data)[1]) == ["a.b"]
    assert ids(run_cli("grep", "-r", "-p", "a.b", stdin=data)[1]) == ["a.b", "aXb"]


def test_grep_invert(run_cli):
    out = run_cli("grep", "-v", "-p", "s1", stdin=FA.encode())[1]
    assert ids(out) == ["s10", "s2"]


def test_grep_no_match_and_empty_input(run_cli):
    assert run_cli("grep", "-p", "zzz", stdin=FA.encode())[:2] == (0, "")
    assert run_cli("grep", "-p", "zzz")[:2] == (0, "")


def test_grep_crlf_and_multiline_sequence(run_cli):
    data = b">a\r\nAC\r\nGT\r\n>b\r\nTT\r\n"
    assert run_cli("grep", "-s", "-p", "CG", stdin=data)[1] == ">a\nACGT\n"


def test_grep_fastq_keeps_quality(run_cli):
    data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n"
    assert run_cli("grep", "-p", "r2", stdin=data)[1] == "@r2\nTTTT\n+\nIIII\n"


def test_grep_requires_a_pattern(run_cli):
    assert run_cli("grep", stdin=FA.encode())[0] == 2


def test_grep_rejects_name_with_seq(run_cli):
    assert run_cli("grep", "-n", "-s", "-p", "a", stdin=FA.encode())[0] == 2


def test_grep_invalid_regex_is_an_error(run_cli):
    code, _, err = run_cli("grep", "-r", "-p", "(", stdin=FA.encode())
    assert code == 1
    assert "invalid regular expression" in err


def test_locate_both_strands_in_forward_coordinates(run_cli):
    code, out, _ = run_cli("locate", "-p", "AAC", stdin=b">a\nAACGTT\n")
    assert code == 0
    assert out.splitlines()[0] == (
        "seqID\tpatternName\tpattern\tstrand\tstart\tend\tmatched")
    assert locate_rows(out) == [["a", "AAC", "AAC", "+", "1", "3", "AAC"],
                                ["a", "AAC", "AAC", "-", "4", "6", "AAC"]]


def test_locate_fixed_matches_overlap(run_cli):
    out = run_cli("locate", "-P", "-p", "AA", stdin=b">a\nAAAA\n")[1]
    assert [row[4:6] for row in locate_rows(out)] == [
        ["1", "2"], ["2", "3"], ["3", "4"]]


def test_locate_only_positive_strand(run_cli):
    out = run_cli("locate", "-P", "-p", "AAC", stdin=b">a\nAACGTT\n")[1]
    assert [row[3] for row in locate_rows(out)] == ["+"]


def test_locate_special_characters_are_literal_without_regex(run_cli):
    out = run_cli("locate", "-P", "-p", "A+", stdin=b">a\nAA+A\n")[1]
    assert [row[4:6] for row in locate_rows(out)] == [["2", "3"]]


def test_locate_is_case_sensitive_unless_asked(run_cli):
    data = b">a\nacgtACGT\n"
    assert locate_rows(run_cli("locate", "-P", "-p", "ACGT", stdin=data)[1]) == [
        ["a", "ACGT", "ACGT", "+", "5", "8", "ACGT"]]
    out = run_cli("locate", "-P", "-i", "-p", "ACGT", stdin=data)[1]
    assert [row[4:6] for row in locate_rows(out)] == [["1", "4"], ["5", "8"]]


def test_locate_no_match_prints_only_the_header(run_cli):
    out = run_cli("locate", "-p", "GGGG", stdin=b">a\nAAAA\n")[1]
    assert len(out.splitlines()) == 1


def test_locate_empty_input_prints_only_the_header(run_cli):
    assert len(run_cli("locate", "-p", "A")[1].splitlines()) == 1


def test_locate_regex_skips_empty_matches(run_cli):
    out = run_cli("locate", "-P", "-r", "-p", "G*", stdin=b">a\nAGGA\n")[1]
    assert [row[4:6] for row in locate_rows(out)] == [["2", "3"], ["3", "3"]]


def test_locate_uses_the_id_and_handles_multiline_crlf(run_cli):
    out = run_cli("locate", "-P", "-p", "CG", stdin=b">x desc\r\nAC\r\nGT\r\n")[1]
    assert locate_rows(out) == [["x", "CG", "CG", "+", "2", "3", "CG"]]


@pytest.mark.parametrize("pattern", ["A", "AAC"])
def test_locate_accepts_several_patterns(run_cli, pattern):
    out = run_cli("locate", "-P", "-p", pattern, "-p", "GT",
                  stdin=b">a\nAACGTT\n")[1]
    assert {row[1] for row in locate_rows(out)} == {pattern, "GT"}


STRAND_FA = b">fwd\nAAACCCGGG\n>rev\nCCCGGGTTT\n>none\nAAAAAAAAA\n"


def test_grep_seq_searches_both_strands_by_default(run_cli):
    out = run_cli("grep", "-s", "-p", "AAACC", stdin=STRAND_FA)[1]
    assert ids(out) == ["fwd", "rev"]


def test_grep_seq_only_positive_strand(run_cli):
    out = run_cli("grep", "-s", "-P", "-p", "AAACC", stdin=STRAND_FA)[1]
    assert ids(out) == ["fwd"]


def test_grep_seq_both_strands_ignore_case_and_invert(run_cli):
    out = run_cli("grep", "-s", "-i", "-p", "aaacc", stdin=STRAND_FA)[1]
    assert ids(out) == ["fwd", "rev"]
    out = run_cli("grep", "-s", "-v", "-p", "AAACC", stdin=STRAND_FA)[1]
    assert ids(out) == ["none"]


def test_grep_seq_regex_also_searches_reverse_complement_of_record(run_cli):
    out = run_cli("grep", "-s", "-r", "-p", "A{3}C+", stdin=STRAND_FA)[1]
    assert ids(out) == ["fwd", "rev"]
    out = run_cli("grep", "-s", "-r", "-P", "-p", "A{3}C+", stdin=STRAND_FA)[1]
    assert ids(out) == ["fwd"]


def test_grep_seq_protein_searches_forward_strand_only(run_cli):
    data = b">p\nMKLVFFA\n"
    assert ids(run_cli("grep", "-s", "-p", "MKL", stdin=data)[1]) == ["p"]
    # Reverse complement of MKL would be LKM... only DNA/RNA gets two strands.
    assert run_cli("grep", "-s", "-p", "LKM", stdin=data)[1] == ""
    assert run_cli("grep", "-s", "-p", "AAAC", stdin=b">p\nMKLVFFA\n")[1] == ""


def test_grep_seq_empty_input(run_cli):
    assert run_cli("grep", "-s", "-p", "ACGT", stdin=b"") == (0, "", "")


def test_grep_P_does_not_change_id_matching(run_cli):
    assert ids(run_cli("grep", "-P", "-p", "fwd", stdin=STRAND_FA)[1]) == ["fwd"]


def test_locate_ignore_case_reports_pattern_letters(run_cli):
    out = run_cli("locate", "-i", "-P", "-p", "acg", stdin=b">a\nACGT\n")[1]
    assert locate_rows(out) == [["a", "acg", "acg", "+", "1", "3", "acg"]]


def test_locate_without_ignore_case_reports_sequence_letters(run_cli):
    rows = locate_rows(run_cli("locate", "-p", "ACG", stdin=b">a\nACGT\n")[1])
    assert rows[0][-1] == "ACG"
