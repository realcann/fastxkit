"""seq command."""
import pytest

FA = ">s1 first\nACGTNN\n>s2\nacgt\n"
FQ = "@r1 x\nACGT\n+\nABCD\n"


def seq(run_cli, content, *flags, make_file=None, name="in.fa"):
    return run_cli("seq", *flags, stdin=content.encode())


def test_reverse_only(run_cli):
    assert seq(run_cli, ">a\nAACG\n", "-r")[1] == ">a\nGCAA\n"


def test_complement_only(run_cli):
    assert seq(run_cli, ">a\nAACG\n", "-p")[1] == ">a\nTTGC\n"


def test_reverse_complement_keeps_case_and_n(run_cli):
    assert seq(run_cli, ">a\nAacgN\n", "-r", "-p")[1] == ">a\nNcgtT\n"


def test_iupac_complement(run_cli):
    assert seq(run_cli, ">a\nRYKMBDHV\n", "-p")[1] == ">a\nYRMKVHDB\n"


def test_fastq_reverse_reverses_quality_but_complement_does_not(run_cli):
    assert seq(run_cli, FQ, "-r")[1] == "@r1 x\nTGCA\n+\nDCBA\n"
    assert seq(run_cli, FQ, "-p")[1] == "@r1 x\nTGCA\n+\nABCD\n"


def test_upper_and_lower(run_cli):
    assert seq(run_cli, ">a\nAcGt\n", "-u")[1] == ">a\nACGT\n"
    assert seq(run_cli, ">a\nAcGt\n", "-l")[1] == ">a\nacgt\n"


def test_upper_and_lower_conflict(run_cli):
    code, out, err = seq(run_cli, ">a\nA\n", "-u", "-l")
    assert code == 2 and out == "" and err


def test_name_and_id_only(run_cli):
    assert seq(run_cli, FA, "-n")[1] == "s1 first\ns2\n"
    assert seq(run_cli, FA, "-n", "-i")[1] == "s1\ns2\n"


def test_id_alone_shortens_header(run_cli):
    assert seq(run_cli, FA, "-i")[1] == ">s1\nACGTNN\n>s2\nacgt\n"


def test_seq_only(run_cli):
    assert seq(run_cli, FA, "-s")[1] == "ACGTNN\nacgt\n"


def test_wrapping_width(run_cli):
    assert seq(run_cli, ">a\nACGTACGTAC\n", "-w", "4")[1] == ">a\nACGT\nACGT\nAC\n"
    assert seq(run_cli, ">a\nACGT\n", "-w", "4")[1] == ">a\nACGT\n"


def test_width_does_not_affect_fastq(run_cli):
    assert seq(run_cli, FQ, "-w", "2")[1] == FQ


def test_multiline_fasta_is_rewrapped_to_default_60(run_cli):
    long = "A" * 130
    out = seq(run_cli, f">a\n{long[:50]}\n{long[50:]}\n")[1]
    assert out.splitlines()[1:] == ["A" * 60, "A" * 60, "A" * 10]


def test_length_bounds_are_inclusive(run_cli):
    data = ">a\nAAAA\n>b\nAAAAA\n>c\nAAAAAA\n"
    assert seq(run_cli, data, "-m", "5", "-M", "5")[1] == ">b\nAAAAA\n"
    assert seq(run_cli, data, "-m", "5")[1] == ">b\nAAAAA\n>c\nAAAAAA\n"
    assert seq(run_cli, data, "-M", "4")[1] == ">a\nAAAA\n"


def test_remove_gaps(run_cli):
    assert seq(run_cli, ">a\nA-C.G T\n", "-g")[1] == ">a\nACGT\n"


def test_remove_gaps_keeps_fastq_quality_aligned(run_cli):
    assert seq(run_cli, "@r\nA-CG\n+\nabcd\n", "-g")[1] == "@r\nACG\n+\nacd\n"


def test_empty_input_gives_empty_output(run_cli):
    assert run_cli("seq")[:2] == (0, "")


def test_crlf_input_gives_lf_output(run_cli):
    assert seq(run_cli, ">a\r\nAC\r\n", "-u")[1] == ">a\nAC\n"


def test_empty_sequence_record_survives(run_cli):
    assert seq(run_cli, ">a\n>b\nAC\n")[1] == ">a\n>b\nAC\n"


def test_reads_multiple_files_in_order(run_cli, make_file):
    a = make_file("a.fa", ">a\nA\n")
    b = make_file("b.fq", "@b\nC\n+\nI\n")
    assert run_cli("seq", a, b)[1] == ">a\nA\n@b\nC\n+\nI\n"


@pytest.mark.parametrize("flag", ["-w", "-m", "-M"])
def test_negative_numbers_rejected(run_cli, flag):
    assert seq(run_cli, ">a\nA\n", flag, "-1")[0] == 2
