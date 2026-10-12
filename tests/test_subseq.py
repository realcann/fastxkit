"""subseq command."""
FA = ">a desc\nACGTACGTAC\n>b\nTT\n"
FQ = "@r1 x\nACGTAC\n+\nIJKLMN\n"


def test_positive_range_is_one_based_inclusive(run_cli):
    assert run_cli("subseq", "-r", "2:4", stdin=FA.encode())[1] == (
        ">a desc\nCGT\n>b\nT\n")


def test_negative_range_counts_from_the_end(run_cli):
    assert run_cli("subseq", "-r", "-3:-1", stdin=FA.encode())[1] == (
        ">a desc\nTAC\n>b\nTT\n")


def test_whole_sequence_edges(run_cli):
    out = run_cli("subseq", "-r", "1:-1", stdin=FA.encode())[1]
    assert out == FA


def test_end_beyond_length_is_clamped(run_cli):
    assert run_cli("subseq", "-r", "9:99", stdin=FA.encode())[1] == (
        ">a desc\nAC\n>b\n\n")


def test_start_beyond_length_gives_an_empty_record(run_cli):
    assert run_cli("subseq", "-r", "11:20", stdin=FA.encode())[1] == (
        ">a desc\n\n>b\n\n")


def test_start_before_the_beginning_is_clamped(run_cli):
    assert run_cli("subseq", "-r", "-50:2", stdin=FA.encode())[1] == (
        ">a desc\nAC\n>b\nTT\n")


def test_start_after_end_gives_an_empty_record(run_cli):
    assert run_cli("subseq", "-r", "5:2", stdin=b">a\nACGTACGT\n")[1] == ">a\n\n"


def test_fastq_slices_quality_too(run_cli):
    assert run_cli("subseq", "-r", "2:3", stdin=FQ.encode())[1] == (
        "@r1 x\nCG\n+\nJK\n")


def test_fastq_empty_slice(run_cli):
    assert run_cli("subseq", "-r", "9:9", stdin=FQ.encode())[1] == (
        "@r1 x\n\n+\n\n")


def test_multiline_and_crlf_input(run_cli):
    data = b">a\r\nACG\r\nTAC\r\n"
    assert run_cli("subseq", "-r", "3:5", stdin=data)[1] == ">a\nGTA\n"


def test_line_width_wraps_long_output(run_cli):
    data = b">a\n" + b"A" * 70 + b"\n"
    out = run_cli("subseq", "-r", "1:-1", stdin=data)[1]
    assert [len(line) for line in out.splitlines()] == [2, 60, 10]


def test_range_is_required_and_zero_rejected(run_cli):
    assert run_cli("subseq", stdin=FA.encode())[0] == 2
    assert run_cli("subseq", "-r", "0:3", stdin=FA.encode())[0] == 2


def test_empty_input(run_cli):
    assert run_cli("subseq", "-r", "1:2")[:2] == (0, "")
