"""Byte-exact comparison against reference outputs in tests/golden3."""
import os
import pathlib
import subprocess
import sys

import pytest

TESTS_DIR = pathlib.Path(__file__).parent
DATA_DIR = TESTS_DIR / "data3"
GOLDEN_DIR = TESTS_DIR / "golden3"
PROJECT_ROOT = TESTS_DIR.parent

# Mirrors CASES and SPLITS in tools/make_golden3.py.
GOLDEN_ARGV = {
    "grep_seq_both_fa": ["grep", "-s", "-p", "ACGTA", "big.fa"],
    "grep_seq_pos_fa": ["grep", "-s", "-P", "-p", "ACGTA", "big.fa"],
    "grep_seq_i_fa": ["grep", "-s", "-i", "-p", "acgta", "big.fa"],
    "grep_seq_re_fq": ["grep", "-s", "-r", "-p", "A{5}", "big.fq"],
    "grep_seq_re_pos_fq": ["grep", "-s", "-r", "-P", "-p", "A{5}", "big.fq"],
    "grep_seq_inv_fa": ["grep", "-s", "-v", "-p", "ACGTA", "big.fa"],
    "locate_i_fq": ["locate", "-i", "-p", "acgt", "big.fq"],
    "locate_i_pos_fq": ["locate", "-i", "-P", "-p", "acgt", "big.fq"],
    "sort_len_fq": ["sort", "-l", "big.fq"],
    "sort_len_rev_fq": ["sort", "-l", "-r", "big.fq"],
    "sort_len_fa": ["sort", "-l", "big.fa"],
    "sort_seq_fq": ["sort", "-s", "big.fq"],
    "fx2tab_lgq_fq": ["fx2tab", "-l", "-g", "-q", "big.fq"],
    "fx2tab_lg_fa": ["fx2tab", "-l", "-g", "big.fa"],
    "fq2fa": ["fq2fa", "big.fq"],
    "fq2fa_small": ["fq2fa", "-o", "/dev/stdout", "dup.fa"],
}
SPLIT_ARGV = {
    "split_p3_fa": ["split", "-p", "3", "big.fa"],
    "split_s15_fq": ["split", "-s", "15", "big.fq"],
    "split_p2_fq": ["split", "-p", "2", "big.fq"],
    "split_s7_fa": ["split", "-s", "7", "big.fa"],
}


def _run(argv):
    env = dict(os.environ, PYTHONPATH=str(PROJECT_ROOT))
    return subprocess.run(
        [sys.executable, "-m", "fastxkit", *argv],
        cwd=DATA_DIR, input=b"", capture_output=True, env=env, check=True,
    ).stdout


@pytest.mark.parametrize("name", sorted(GOLDEN_ARGV))
def test_golden3(name):
    expected = (GOLDEN_DIR / f"{name}.out").read_bytes()
    assert _run(GOLDEN_ARGV[name]) == expected


@pytest.mark.parametrize("name", sorted(SPLIT_ARGV))
def test_golden3_split(name, tmp_path):
    out_dir = tmp_path / "parts"
    assert _run([*SPLIT_ARGV[name], "-O", str(out_dir)]) == b""
    combined = "".join(f"== {p.name}\n{p.read_text()}"
                       for p in sorted(out_dir.iterdir()))
    assert combined == (GOLDEN_DIR / f"{name}.out").read_text()
