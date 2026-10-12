"""fx2tab and tab2fx."""
import pytest

FA = ">s1 d\nACGTNN\n>s2\nGGCC\n"
FQ = "@r1 x\nACGT\n+\nIIII\n"


def test_fx2tab_fasta_has_empty_quality_column(run_cli):
    assert run_cli("fx2tab", stdin=FA.encode())[1] == "s1 d\tACGTNN\t\ns2\tGGCC\t\n"


def test_fx2tab_fastq(run_cli):
    assert run_cli("fx2tab", stdin=FQ.encode())[1] == "r1 x\tACGT\tIIII\n"


def test_fx2tab_name_and_id(run_cli):
    assert run_cli("fx2tab", "-n", stdin=FA.encode())[1] == "s1 d\ns2\n"
    assert run_cli("fx2tab", "-n", "-i", stdin=FA.encode())[1] == "s1\ns2\n"
    assert run_cli("fx2tab", "-i", stdin=FA.encode())[1] == "s1\tACGTNN\t\ns2\tGGCC\t\n"


def test_fx2tab_length_gc_and_base_content(run_cli):
    out = run_cli("fx2tab", "-l", "-g", "-B", "AT", stdin=FA.encode())[1]
    assert out == "s1 d\tACGTNN\t\t6\t50.00\t33.33\ns2\tGGCC\t\t4\t100.00\t0.00\n"


def test_fx2tab_multiple_base_contents(run_cli):
    out = run_cli("fx2tab", "-B", "A", "-B", "C", stdin=b">a\nAACG\n")[1]
    assert out == "a\tAACG\t\t50.00\t25.00\n"


def test_fx2tab_avg_quality(run_cli):
    out = run_cli("fx2tab", "-q", stdin=b"@r\nAC\n+\n!I\n")[1]
    assert out == "r\tAC\t!I\t3.01\n"


def test_fx2tab_empty_sequence(run_cli):
    assert run_cli("fx2tab", "-l", "-g", stdin=b">a\n")[1] == "a\t\t\t0\t0.00\n"


def test_fx2tab_empty_input(run_cli):
    assert run_cli("fx2tab")[:2] == (0, "")


def test_fx2tab_multiline_and_crlf(run_cli):
    assert run_cli("fx2tab", stdin=b">a\r\nAC\r\nGT\r\n")[1] == "a\tACGT\t\n"


def test_tab2fx_fasta_wraps_at_60(run_cli):
    out = run_cli("tab2fx", stdin=f"a\t{'A' * 70}\t\n".encode())[1]
    assert out == f">a\n{'A' * 60}\n{'A' * 10}\n"


def test_tab2fx_two_columns_is_fasta(run_cli):
    assert run_cli("tab2fx", stdin=b"a\tAC\n")[1] == ">a\nAC\n"


def test_tab2fx_fastq_when_quality_present(run_cli):
    assert run_cli("tab2fx", stdin=b"r 1\tAC\tII\n")[1] == "@r 1\nAC\n+\nII\n"


def test_tab2fx_ignores_blank_lines_and_crlf(run_cli):
    assert run_cli("tab2fx", stdin=b"\na\tAC\r\n\n")[1] == ">a\nAC\n"


def test_tab2fx_empty_input(run_cli):
    assert run_cli("tab2fx")[:2] == (0, "")


@pytest.mark.parametrize("line, fragment", [
    (b"onlyname\n", "line 1"),
    (b"a\tACGT\tII\n", "length"),
])
def test_tab2fx_bad_rows(run_cli, line, fragment):
    code, _, err = run_cli("tab2fx", stdin=line)
    assert code == 1 and fragment in err


@pytest.mark.parametrize("content", [FA, FQ])
def test_roundtrip_through_pipe(run_cli, content):
    table = run_cli("fx2tab", stdin=content.encode())[1]
    assert run_cli("tab2fx", stdin=table.encode())[1] == content


def test_fx2tab_exact_halves_use_standard_formatting(run_cli):
    # 51 G/C of 96 bases is exactly 53.125, which Python formats as 53.12.
    seq = "G" * 51 + "A" * 45
    assert run_cli("fx2tab", "-g", stdin=f">a\n{seq}\n".encode())[1] == (
        f"a\t{seq}\t\t53.12\n")


def test_stats_still_rounds_exact_halves_up(run_cli):
    seq = "G" * 51 + "A" * 45
    out = run_cli("stats", "-a", stdin=f">a\n{seq}\n".encode())[1]
    assert "53.13" in out
