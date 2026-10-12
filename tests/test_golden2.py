"""Byte-exact comparison against reference outputs in tests/golden2."""
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

TESTS_DIR = pathlib.Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"
GOLDEN_DIR = TESTS_DIR / "golden2"
PROJECT_ROOT = TESTS_DIR.parent

# Mirrors the run lines of tests/golden2/make_golden.sh.
GOLDEN_ARGV = {
    "grep_id_fa": ["grep", "-p", "seq2", "small.fa"],
    "grep_name_re_fa": ["grep", "-n", "-r", "-p", "seq[24]", "small.fa"],
    "grep_invert_fa": ["grep", "-n", "-v", "-r", "-p", "seq[24]", "small.fa"],
    "grep_seq_fa": ["grep", "-s", "-i", "-p", "ACGT", "dup.fa"],
    "grep_seq_regex_fa": ["grep", "-s", "-r", "-p", "A{3}G{3}", "dup.fa"],
    "grep_desc_fa": ["grep", "-n", "-r", "-p", "record [13]", "small.fa"],
    "grep_seq_fq": ["grep", "-s", "-p", "GGGCCC", "motif.fq"],
    "locate_fa": ["locate", "-p", "ACGT", "dup.fa"],
    "locate_nonover_fa": ["locate", "-P", "-p", "ACGT", "dup.fa"],
    "locate_regex_fq": ["locate", "-r", "-p", "A+G+", "motif.fq"],
    "subseq_region_fa": ["subseq", "-r", "2:10", "small.fa"],
    "subseq_neg_fa": ["subseq", "-r", "-5:-1", "small.fa"],
    "subseq_flank_fa": ["subseq", "-r", "11:20", "small.fa"],
    "sort_name_fa": ["sort", "-n", "uniq.fa"],
    "sort_len_fa": ["sort", "-l", "uniq.fa"],
    "sort_len_rev_fa": ["sort", "-l", "-r", "uniq.fa"],
    "sort_seq_fa": ["sort", "-s", "uniq.fa"],
    "rmdup_seq_fa": ["rmdup", "-s", "dup.fa"],
    "rmdup_seq_i_fa": ["rmdup", "-s", "-i", "dup.fa"],
    "rmdup_name_fa": ["rmdup", "-n", "dup.fa"],
    "rmdup_id_fa": ["rmdup", "dup.fa"],
    "rename_fa": ["rename", "dup.fa"],
    "replace_name_fa": ["replace", "-p", "seq", "-r", "S", "small.fa"],
    "replace_re_name_fa": ["replace", "-p", r"^seq(\d)", "-r", "id_$1", "small.fa"],
    "replace_seq_fa": ["replace", "-s", "-p", "N", "-r", "-", "small.fa"],
    "replace_pat_fq": ["replace", "-p", "^r", "-r", "read_", "motif.fq"],
}


def _run(argv, cwd=DATA_DIR):
    env = dict(os.environ, PYTHONPATH=str(PROJECT_ROOT))
    return subprocess.run(
        [sys.executable, "-m", "fastxkit", *argv],
        cwd=cwd, input=b"", capture_output=True, env=env, check=True,
    ).stdout


@pytest.mark.parametrize("name", sorted(GOLDEN_ARGV))
def test_golden2(name):
    expected = (GOLDEN_DIR / f"{name}.out").read_bytes()
    assert _run(GOLDEN_ARGV[name]) == expected


def test_golden2_faidx(tmp_path):
    shutil.copy(DATA_DIR / "small.fa", tmp_path / "sf.fa")
    assert _run(["faidx", "sf.fa"], cwd=tmp_path) == b""
    assert ((tmp_path / "sf.fa.fai").read_bytes()
            == (GOLDEN_DIR / "faidx_small_fa.fai").read_bytes())
    extracted = _run(["faidx", "sf.fa", "seq1", "seq5:10-20"], cwd=tmp_path)
    assert extracted == (GOLDEN_DIR / "faidx_extract.out").read_bytes()
