"""qc command."""
import gzip
import hashlib
import json

from fastxkit import __version__, numfmt

SMALL_FQ = "tests/data/small.fq"
FQ = "@r1\nACGN\n+\n!!II\n"


def test_matches_stats_all_golden(run_cli):
    _, stats_out, _ = run_cli("stats", "-a", SMALL_FQ)
    with open("tests/golden/stats_all_fq.out") as golden:
        expected = [line.split() for line in golden.read().splitlines()]
    actual = [line.split() for line in stats_out.splitlines()]
    assert actual[0] == expected[0]
    assert actual[1][1:] == expected[1][1:]  # the file name differs
    stats = dict(zip(*(line.split() for line in stats_out.splitlines())))
    _, out, _ = run_cli("qc", "--json", SMALL_FQ)
    [entry] = json.loads(out)
    assert entry["reads"] == int(stats["num_seqs"])
    assert entry["total_bases"] == int(stats["sum_len"])
    assert entry["min_len"] == int(stats["min_len"])
    assert entry["max_len"] == int(stats["max_len"])
    assert numfmt.trimmed(entry["mean_len"], 1) == stats["avg_len"]
    assert numfmt.trimmed(entry["gc_percent"], 2) == stats["GC(%)"]
    assert numfmt.trimmed(entry["q20_percent"], 2) == stats["Q20(%)"]
    assert numfmt.trimmed(entry["q30_percent"], 2) == stats["Q30(%)"]
    assert numfmt.trimmed(entry["mean_quality"], 2) == stats["AvgQual"]
    assert entry["n_bases"] == int(stats["sum_n"])


def test_golden_row_values(run_cli):
    _, out, _ = run_cli("qc", "--json", SMALL_FQ)
    [entry] = json.loads(out)
    assert (entry["reads"], entry["total_bases"]) == (12, 975)
    assert entry["tool_version"] == __version__
    with open(SMALL_FQ, "rb") as handle:
        assert entry["sha256"] == hashlib.sha256(handle.read()).hexdigest()


def test_gzip_gives_same_counts_but_own_checksum(run_cli, make_file):
    plain = make_file("a.fq", FQ)
    packed = make_file("b.fq.gz", gzip.compress(FQ.encode()))
    [one] = json.loads(run_cli("qc", "--json", plain)[1])
    [two] = json.loads(run_cli("qc", "--json", packed)[1])
    assert one["reads"] == two["reads"] == 1
    assert one["sha256"] != two["sha256"]


def test_single_read_values(run_cli, make_file):
    [entry] = json.loads(run_cli("qc", "--json", make_file("a.fq", FQ))[1])
    assert entry["n_bases"] == 1
    assert entry["q20_percent"] == 50
    assert entry["min_len"] == entry["max_len"] == 4
    assert entry["gc_percent"] == 50


def test_empty_fastq(run_cli, make_file):
    code, out, _ = run_cli("qc", "--json", make_file("e.fq", ""))
    assert code == 0
    [entry] = json.loads(out)
    assert (entry["reads"], entry["total_bases"], entry["min_len"]) == (0, 0, 0)
    assert entry["gc_percent"] == 0


def test_fasta_has_zero_quality(run_cli, make_file):
    [entry] = json.loads(run_cli("qc", "--json", make_file("a.fa", ">a\nGC\n"))[1])
    assert entry["q30_percent"] == 0
    assert entry["gc_percent"] == 100


def test_json_rounds_half_up_to_six_decimals(run_cli, make_file):
    # 1/3 repeating needs rounding at 6 decimals.
    path = make_file("a.fa", ">a\nGAA\n")
    [entry] = json.loads(run_cli("qc", "--json", path)[1])
    assert entry["gc_percent"] == 33.333333


def test_table_output(run_cli, make_file):
    code, out, _ = run_cli("qc", make_file("a.fq", FQ))
    assert code == 0
    header, row = (line.split() for line in out.splitlines())
    assert header[:3] == ["file", "reads", "total_bases"]
    assert "sha256" in header and row[1] == "1"


def test_olgu_output(run_cli, make_file):
    path = make_file("S1_1.fastq.gz", gzip.compress(FQ.encode()))
    code, out, _ = run_cli("qc", "--olgu", path)
    assert code == 0
    olgu = json.loads(out)
    assert sorted(olgu) == sorted(
        f"{key}:S1_1" for key in ("fastq_okuma_sayisi", "fastq_toplam_baz",
                                  "fastq_gc_yuzde", "fastq_q30_yuzde",
                                  "fastq_ort_uzunluk"))
    fact = olgu["fastq_okuma_sayisi:S1_1"]
    assert fact["deger"] == 1
    assert fact["kaynak"] == path
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    assert fact["konum"] == f"fastxkit qc {__version__}; sha256 {digest}"


def test_olgu_duplicate_stem_is_error(run_cli, make_file, tmp_path):
    first = make_file("S1.fq", FQ)
    (tmp_path / "sub").mkdir()
    second = make_file("sub/S1.fq", FQ)
    code, _, err = run_cli("qc", "--olgu", first, second)
    assert code == 1
    assert "S1" in err


def test_malformed_fastq_exits_one(run_cli, make_file):
    code, out, err = run_cli("qc", make_file("bad.fq", "@r\nACGT\n+\n!!\n"))
    assert code == 1
    assert "differs" in err and out == ""


def test_missing_file_exits_one(run_cli):
    code, _, err = run_cli("qc", "/nonexistent/x.fq")
    assert code == 1 and "x.fq" in err


def test_qc_rejects_stdin_as_a_usage_error(run_cli):
    code, out, err = run_cli("qc", "-", stdin=b"@r\nA\n+\nI\n")
    assert (code, out) == (2, "")
    assert "stdin" in err


def test_qc_help_says_stdin_is_unsupported(run_cli):
    assert "stdin is not supported" in " ".join(run_cli("qc", "-h")[1].split())


def test_invalid_quality_character_names_the_file(run_cli, make_file):
    path = make_file("bad.fq", "@r1\nACGT\n+\nII I\n")
    for command in ("qc", "stats"):
        code, _, err = run_cli(command, path)
        assert code == 1 and path in err and "invalid quality" in err
    code, _, err = run_cli("fx2tab", "-q", path)
    assert code == 1 and path in err and "invalid quality" in err
