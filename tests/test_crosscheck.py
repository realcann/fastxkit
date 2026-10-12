"""crosscheck command."""
import json
import pathlib
import subprocess
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DEMO_OLGU = str(pathlib.Path(__file__).resolve().parent / "data4" / "olgu_demo.json")


def fastq(reads):
    return "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(reads))


def olgu_file(make_file, table, wrap=True):
    fact = {"deger": table, "kaynak": "multiqc", "konum": "x"}
    return make_file("olgu.json", json.dumps({"ornek_tablosu": fact} if wrap else table))


def sample_rows(*pairs):
    return [{"ornek": name, "grup": "g", "toplam_okuma": count,
             "essiz_hizalanan_yuzde": 85.1} for name, count in pairs]


def status_of(out):
    return {cells[0]: cells[-1] for cells in
            (line.split() for line in out.splitlines()[1:])}


def test_all_ok(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 3), ("B", 2)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(3))
                    + "\nB," + make_file("b.fq", fastq(2)) + "\n")
    code, out, _ = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 0
    assert status_of(out) == {"A": "OK", "B": "OK"}


def test_float_value_from_multiqc_is_accepted(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 3.0)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(3)) + "\n")
    assert run_cli("crosscheck", olgu, "--map", csv)[0] == 0


def test_mismatch_and_tolerance(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 5)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(4)) + "\n")
    code, out, _ = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 1 and status_of(out) == {"A": "MISMATCH"}
    assert "-1" in out.splitlines()[1].split()
    code, out, _ = run_cli("crosscheck", olgu, "--map", csv, "--tolerance", "1")
    assert code == 0 and status_of(out) == {"A": "OK"}


def test_paired_all_vs_r1(run_cli, make_file):
    r1, r2 = make_file("a1.fq", fastq(3)), make_file("a2.fq", fastq(3))
    csv = make_file("map.csv", f"sample,fastq,read\nA,{r1},r1\nA,{r2},r2\n")
    pairs = olgu_file(make_file, sample_rows(("A", 3)))
    assert run_cli("crosscheck", pairs, "--map", csv)[0] == 1
    assert run_cli("crosscheck", pairs, "--map", csv, "--count", "r1")[0] == 0
    both = olgu_file(make_file, sample_rows(("A", 6)))
    assert run_cli("crosscheck", both, "--map", csv)[0] == 0


def test_crlf_blank_lines_and_missing_statuses(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 1), ("B", 1)))
    csv = make_file("map.csv", "sample,fastq\r\n\r\nA," + make_file("a.fq", fastq(1))
                    + "\r\nC," + make_file("c.fq", fastq(1)) + "\r\n\r\n")
    code, out, _ = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 1
    assert status_of(out) == {"A": "OK", "B": "MISSING_IN_MAP",
                              "C": "MISSING_IN_OLGU"}


def test_bad_fastq_and_missing_file_are_file_error(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 1), ("B", 1)))
    bad = make_file("bad.fq", "@r\nACGT\n+\n!!\n")
    csv = make_file("map.csv", f"sample,fastq\nA,{bad}\nB,/nonexistent/b.fq\n")
    code, out, err = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 3
    assert status_of(out) == {"A": "FILE_ERROR", "B": "FILE_ERROR"}
    assert "ACGT" not in out + err


def test_json_output(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 2)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(1)) + "\n")
    code, out, _ = run_cli("crosscheck", olgu, "--map", csv, "--json")
    assert code == 1
    [row] = json.loads(out)
    assert row == {"sample": "A", "olgu": 2, "fastq": 1, "diff": -1,
                   "status": "MISMATCH"}


@pytest.mark.parametrize("content, message", [
    ("{not json", "olgu"),
    ("{}", "ornek_tablosu"),
    ('{"ornek_tablosu": {"deger": "x"}}', "ornek_tablosu"),
    ('{"ornek_tablosu": {"deger": [{"ornek": "A", "toplam_okuma": 1.5}]}}',
     "toplam_okuma"),
    ('{"ornek_tablosu": {"deger": [{"ornek": "A", "toplam_okuma": -1}]}}',
     "toplam_okuma"),
    ('{"ornek_tablosu": {"deger": [{"ornek": "A", "toplam_okuma": 1},'
     ' {"ornek": "A ", "toplam_okuma": 1}]}}', "duplicate"),
])
def test_bad_olgu(run_cli, make_file, content, message):
    olgu = make_file("olgu.json", content)
    csv = make_file("map.csv", "sample,fastq\n")
    code, out, err = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 3 and message in err and out == ""


def test_missing_olgu(run_cli, make_file):
    csv = make_file("map.csv", "sample,fastq\n")
    code, _, err = run_cli("crosscheck", "/nonexistent/olgu.json", "--map", csv)
    assert code == 3 and "olgu.json" in err


@pytest.mark.parametrize("content", [
    "A,x.fq\n",                     # header missing
    "name,path\nA,x.fq\n",          # wrong header
    "sample,fastq\nA\n",            # short row
    "sample,fastq,read\nA,x.fq,r3\n",
    "sample,fastq,reads\nA,x.fq,r1\n",   # third column is not 'read'
    "sample,fastq\nA,x.fq\nA,x.fq\n",     # same (sample, fastq) twice
    "",
])
def test_bad_map(run_cli, make_file, content):
    olgu = olgu_file(make_file, sample_rows(("A", 1)))
    code, out, err = run_cli("crosscheck", olgu, "--map",
                             make_file("map.csv", content))
    assert code == 3 and "map" in err and out == ""


def test_missing_map_file(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 1)))
    code, _, err = run_cli("crosscheck", olgu, "--map", "/nonexistent/map.csv")
    assert code == 3 and "map.csv" in err


def test_count_r1_requires_read_column(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 1)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(1)) + "\n")
    code, _, err = run_cli("crosscheck", olgu, "--map", csv, "--count", "r1")
    assert code == 3 and "read" in err


def test_relative_fastq_path_is_relative_to_map(run_cli, make_file):
    make_file("a.fq", fastq(1))
    olgu = olgu_file(make_file, sample_rows(("A", 1)))
    csv = make_file("map.csv", "sample,fastq\nA,a.fq\n")
    assert run_cli("crosscheck", olgu, "--map", csv)[0] == 0


def test_demo_end_to_end(tmp_path):
    """Three samples; sample C has one read fewer than the olgu claims."""
    for name, reads in (("A", 30), ("B", 20), ("C", 9)):
        (tmp_path / f"{name}.fq").write_text(fastq(reads))
    rows = [{"ornek": "A", "grup": "control", "toplam_okuma": 30.0},
            {"ornek": "B", "grup": "control", "toplam_okuma": 20},
            {"ornek": "C", "grup": "treated", "toplam_okuma": 10}]
    demo = json.load(open(DEMO_OLGU))
    assert demo["ornek_tablosu"]["deger"] == rows
    (tmp_path / "map.csv").write_text(
        "sample,fastq\n" + "".join(f"{n},{n}.fq\n" for n in "ABC"))
    result = subprocess.run(
        [sys.executable, "-m", "fastxkit", "crosscheck", DEMO_OLGU, "--map",
         str(tmp_path / "map.csv")], capture_output=True, text=True,
        cwd=REPO_ROOT)
    assert result.returncode == 1
    assert status_of(result.stdout) == {"A": "OK", "B": "OK", "C": "MISMATCH"}


def test_nothing_to_compare_is_an_input_error(run_cli, make_file):
    olgu = olgu_file(make_file, [])
    csv = make_file("map.csv", "sample,fastq\n")
    code, out, err = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 3 and "nothing to compare" in err and out == ""


def test_olgu_sample_name_is_stripped(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows((" A ", 1)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(1)) + "\n")
    code, out, _ = run_cli("crosscheck", olgu, "--map", csv)
    assert code == 0 and status_of(out) == {"A": "OK"}


def test_file_error_wins_over_mismatch(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 5), ("B", 1)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(1))
                    + "\nB,/nonexistent/b.fq\n")
    assert run_cli("crosscheck", olgu, "--map", csv)[0] == 3


def test_usage_error_stays_2(run_cli):
    assert run_cli("crosscheck")[0] == 2


def test_sample_and_status_columns_are_left_aligned(run_cli, make_file):
    olgu = olgu_file(make_file, sample_rows(("A", 1), ("LONGNAME", 1)))
    csv = make_file("map.csv", "sample,fastq\nA," + make_file("a.fq", fastq(1))
                    + "\nLONGNAME," + make_file("b.fq", fastq(1)) + "\n")
    lines = run_cli("crosscheck", olgu, "--map", csv)[1].splitlines()
    assert lines[1].startswith("A         ")
    assert lines[1].endswith("OK") and lines[2].endswith("OK")
    assert lines[0].split("  ")[-1] == "status"
