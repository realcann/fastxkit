"""Byte-exact comparison against reference outputs in tests/golden."""
import os
import pathlib
import subprocess
import sys

import pytest

TESTS_DIR = pathlib.Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"
GOLDEN_DIR = TESTS_DIR / "golden"
PROJECT_ROOT = TESTS_DIR.parent

# Mirrors the run lines of tests/golden/make_golden.sh.
GOLDEN_ARGV = {
    "stats_basic": ["stats", "small.fa", "small.fq", "small.fq.gz"],
    "stats_all_fq": ["stats", "-a", "small.fq"],
    "stats_all_fa": ["stats", "-a", "small.fa"],
    "seq_rc_fa": ["seq", "-r", "-p", "small.fa"],
    "seq_upper_fa": ["seq", "-u", "small.fa"],
    "seq_lower_fa": ["seq", "-l", "small.fa"],
    "seq_name_fa": ["seq", "-n", "small.fa"],
    "seq_id_fa": ["seq", "-n", "-i", "small.fa"],
    "seq_seq_fa": ["seq", "-s", "small.fa"],
    "seq_w0_fa": ["seq", "-w", "0", "small.fa"],
    "seq_w30_fa": ["seq", "-w", "30", "small.fa"],
    "seq_minlen_fa": ["seq", "-m", "60", "small.fa"],
    "seq_maxlen_fq": ["seq", "-M", "100", "small.fq"],
    "seq_rc_fq": ["seq", "-r", "-p", "small.fq"],
    "head_fa": ["head", "-n", "2", "small.fa"],
    "head_fq": ["head", "-n", "3", "small.fq"],
    "range_fa": ["range", "-r", "2:4", "small.fa"],
    "range_neg_fa": ["range", "-r", "-2:-1", "small.fa"],
    "range_fq": ["range", "-r", "1:2", "small.fq"],
    "fx2tab_fa": ["fx2tab", "small.fa"],
    "fx2tab_fq": ["fx2tab", "small.fq"],
    "fx2tab_lg_fq": ["fx2tab", "-l", "-g", "small.fq"],
    "fx2tab_name_fq": ["fx2tab", "-n", "small.fq"],
    "fx2tab_id_fa": ["fx2tab", "-n", "-i", "small.fa"],
    "fx2tab_q_fq": ["fx2tab", "-q", "small.fq"],
    "fx2tab_lg_fa": ["fx2tab", "-l", "-g", "-B", "AT", "small.fa"],
}
# name -> input file piped through fx2tab, then tab2fx.
ROUNDTRIP = {"tab2fx_fa": "small.fa", "tab2fx_fq": "small.fq"}


def _run(argv, stdin=None):
    env = dict(os.environ, PYTHONPATH=str(PROJECT_ROOT))
    return subprocess.run(
        [sys.executable, "-m", "fastxkit", *argv],
        cwd=DATA_DIR, input=stdin, capture_output=True, env=env, check=True,
    ).stdout


@pytest.mark.parametrize("name", sorted(GOLDEN_ARGV))
def test_golden(name):
    expected = (GOLDEN_DIR / f"{name}.out").read_bytes()
    assert _run(GOLDEN_ARGV[name]) == expected


@pytest.mark.parametrize("name", sorted(ROUNDTRIP))
def test_golden_tab2fx_roundtrip(name):
    expected = (GOLDEN_DIR / f"{name}.out").read_bytes()
    table = _run(["fx2tab", ROUNDTRIP[name]])
    assert _run(["tab2fx"], stdin=table) == expected
