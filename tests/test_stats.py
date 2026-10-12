"""stats command."""
FQ = "@r1\nACGT\n+\n!!II\n"


def rows(out):
    return [line.split() for line in out.splitlines()]


def test_default_columns(run_cli, make_file):
    path = make_file("a.fa", ">a\nACGT\n>b\nAC\n")
    code, out, _ = run_cli("stats", path)
    assert code == 0
    header, row = rows(out)
    assert header == ["file", "format", "type", "num_seqs", "sum_len",
                      "min_len", "avg_len", "max_len"]
    assert row[1:] == ["FASTA", "DNA", "2", "6", "2", "3", "4"]


def test_all_columns(run_cli, make_file):
    path = make_file("a.fa", ">a\nACGT\n")
    _, out, _ = run_cli("stats", "-a", path)
    assert rows(out)[0][8:] == ["Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num",
                                "Q20(%)", "Q30(%)", "AvgQual", "GC(%)", "sum_n"]


def test_fastq_quality_columns(run_cli, make_file):
    _, out, _ = run_cli("stats", "-a", make_file("a.fq", FQ))
    row = dict(zip(*rows(out)))
    assert row["format"] == "FASTQ"
    assert (row["Q20(%)"], row["Q30(%)"], row["AvgQual"]) == ("50", "50", "3.01")


def test_n_gap_and_gc_counts(run_cli, make_file):
    path = make_file("a.fa", ">a\nGCNN-.AT\n")
    _, out, _ = run_cli("stats", "-a", path)
    row = dict(zip(*rows(out)))
    assert (row["sum_n"], row["sum_gap"], row["GC(%)"]) == ("2", "2", "25")


def test_n50_and_quartiles(run_cli, make_file):
    seqs = "".join(f">s{i}\n{'A' * n}\n" for i, n in enumerate([10, 60, 75, 130, 200]))
    _, out, _ = run_cli("stats", "-a", make_file("a.fa", seqs))
    row = dict(zip(*rows(out)))
    assert [row[k] for k in ("Q1", "Q2", "Q3", "N50", "N50_num")] == [
        "60", "75", "130", "130", "2"]


def test_large_numbers_are_grouped(run_cli, make_file):
    path = make_file("a.fa", f">a\n{'A' * 1500}\n")
    _, out, _ = run_cli("stats", path)
    assert rows(out)[1][4] == "1,500"


def test_type_detection(run_cli, make_file):
    for content, kind in [(">a\nACGU\n", "RNA"), (">a\nMKVLAAG\n", "Protein"),
                          (">a\nAC1!\n", "Unlimit")]:
        _, out, _ = run_cli("stats", make_file("t.fa", content))
        assert rows(out)[1][2] == kind


def test_empty_file_has_zero_counts(run_cli, make_file):
    code, out, _ = run_cli("stats", "-a", make_file("e.fa", ""))
    row = dict(zip(*rows(out)))
    assert code == 0
    assert (row["num_seqs"], row["sum_len"], row["avg_len"], row["N50"]) == (
        "0", "0", "0", "0")


def test_single_record(run_cli, make_file):
    _, out, _ = run_cli("stats", "-a", make_file("s.fa", ">a\nACGT\n"))
    row = dict(zip(*rows(out)))
    assert (row["Q1"], row["Q2"], row["Q3"], row["N50_num"]) == ("4", "4", "4", "1")


def test_columns_are_aligned_for_multiple_files(run_cli, make_file):
    a = make_file("a.fa", ">a\nACGT\n")
    b = make_file("long_name.fq", FQ)
    _, out, _ = run_cli("stats", a, b)
    lines = out.splitlines()
    assert len({len(line) for line in lines}) == 1


def test_stdin_row_is_named_dash(run_cli):
    _, out, _ = run_cli("stats", stdin=b">a\nAC\n")
    assert rows(out)[1][0] == "-"


def test_crlf_input(run_cli, make_file):
    _, out, _ = run_cli("stats", make_file("w.fa", ">a\r\nACGT\r\n"))
    assert rows(out)[1][3:5] == ["1", "4"]


def test_protein_has_zero_gc_and_n(run_cli, make_file):
    path = make_file("p.fa", ">p1\nMKVLAAGIVGLLLAQPSEA\n>p2\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n")
    _, out, _ = run_cli("stats", "-a", path)
    row = dict(zip(*rows(out)))
    assert row["type"] == "Protein"
    assert (row["GC(%)"], row["sum_n"]) == ("0", "0")
