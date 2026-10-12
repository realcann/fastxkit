"""replace command."""
import pytest

FA = ">seq1 first\nACGTN\n>seq2 second\nNNAC\n"


def test_replaces_in_the_name_by_default(run_cli):
    out = run_cli("replace", "-p", "seq", "-r", "S", stdin=FA.encode())[1]
    assert out == ">S1 first\nACGTN\n>S2 second\nNNAC\n"


def test_replaces_every_occurrence(run_cli):
    out = run_cli("replace", "-s", "-p", "N", "-r", "-", stdin=FA.encode())[1]
    assert out == ">seq1 first\nACGT-\n>seq2 second\n--AC\n"


def test_pattern_is_a_regular_expression(run_cli):
    out = run_cli("replace", "-p", r"^seq(\d)", "-r", "id_$1",
                  stdin=FA.encode())[1]
    assert out.splitlines()[0] == ">id_1 first"


def test_braced_group_reference_and_whole_match(run_cli):
    out = run_cli("replace", "-p", r"seq(\d)", "-r", "${1}x$0", stdin=FA.encode())
    assert out[1].splitlines()[0] == ">1xseq1 first"


def test_dollar_escape_and_backslash_are_literal(run_cli):
    out = run_cli("replace", "-p", "seq", "-r", r"$$\1", stdin=FA.encode())
    assert out[1].splitlines()[0] == r">$\11 first"


def test_unmatched_optional_group_becomes_empty(run_cli):
    out = run_cli("replace", "-p", "(x)?seq", "-r", "<$1>", stdin=FA.encode())
    assert out[1].splitlines()[0] == "><>1 first"


def test_ignore_case(run_cli):
    out = run_cli("replace", "-i", "-p", "SEQ", "-r", "S", stdin=FA.encode())[1]
    assert out.splitlines()[0] == ">S1 first"


def test_no_match_leaves_input_unchanged(run_cli):
    assert run_cli("replace", "-p", "zzz", "-r", "y", stdin=FA.encode())[1] == FA


def test_empty_replacement_deletes(run_cli):
    out = run_cli("replace", "-p", " first", "-r", "", stdin=FA.encode())[1]
    assert out.splitlines()[0] == ">seq1"


def test_replacing_sequence_of_fastq_with_same_length(run_cli):
    data = b"@r\nACNT\n+\nIIII\n"
    assert run_cli("replace", "-s", "-p", "N", "-r", "G", stdin=data)[1] == (
        "@r\nACGT\n+\nIIII\n")


def test_fastq_sequence_length_change_is_an_error(run_cli):
    data = b"@r\nACNT\n+\nIIII\n"
    code, _, err = run_cli("replace", "-s", "-p", "N", "-r", "GG", stdin=data)
    assert code == 1
    assert "quality" in err


def test_multiline_and_crlf(run_cli):
    data = b">a\r\nAC\r\nNN\r\n"
    assert run_cli("replace", "-s", "-p", "CN", "-r", "xx", stdin=data)[1] == (
        ">a\nAxxN\n")


def test_empty_input(run_cli):
    assert run_cli("replace", "-p", "a", "-r", "b")[:2] == (0, "")


def test_group_reference_beyond_pattern_is_an_error(run_cli):
    code, _, err = run_cli("replace", "-p", "seq", "-r", "$1", stdin=FA.encode())
    assert code == 1
    assert "group" in err


def test_invalid_regex_is_an_error(run_cli):
    code, _, err = run_cli("replace", "-p", "(", "-r", "x", stdin=FA.encode())
    assert code == 1
    assert "invalid regular expression" in err


@pytest.mark.parametrize("missing", [["-p", "a"], ["-r", "b"]])
def test_pattern_and_replacement_are_required(run_cli, missing):
    assert run_cli("replace", *missing, stdin=FA.encode())[0] == 2
