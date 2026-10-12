"""Generate data3/ and golden3/ reference files with a reference seqkit binary.

Usage: python tools/make_golden3.py <reference-binary>
Synthetic data only (tools/diff_reference.make_data, seed 1). Outputs are stored so tests need no binary.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import diff_reference as d  # noqa: E402

TESTS = HERE.parent / "tests"
CASES = {
    "grep_seq_both_fa": "grep -s -p ACGTA big.fa",
    "grep_seq_pos_fa": "grep -s -P -p ACGTA big.fa",
    "grep_seq_i_fa": "grep -s -i -p acgta big.fa",
    "grep_seq_re_fq": "grep -s -r -p A{5} big.fq",
    "grep_seq_re_pos_fq": "grep -s -r -P -p A{5} big.fq",
    "grep_seq_inv_fa": "grep -s -v -p ACGTA big.fa",
    "locate_i_fq": "locate -i -p acgt big.fq",
    "locate_i_pos_fq": "locate -i -P -p acgt big.fq",
    "sort_len_fq": "sort -l big.fq",
    "sort_len_rev_fq": "sort -l -r big.fq",
    "sort_len_fa": "sort -l big.fa",
    "sort_seq_fq": "sort -s big.fq",
    "fx2tab_lgq_fq": "fx2tab -l -g -q big.fq",
    "fx2tab_lg_fa": "fx2tab -l -g big.fa",
    "fq2fa": "fq2fa big.fq",
    "fq2fa_small": "fq2fa -o /dev/stdout dup.fa",
}
SPLITS = {
    "split_p3_fa": "split -p 3 big.fa",
    "split_s15_fq": "split -s 15 big.fq",
    "split_p2_fq": "split -p 2 big.fq",
    "split_s7_fa": "split -s 7 big.fa",
}


def main() -> int:
    ref = sys.argv[1]
    (TESTS / "data3").mkdir(exist_ok=True)
    (TESTS / "golden3").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as t:
        w = Path(t)
        d.make_data(w, 1)
        for f in ("big.fa", "big.fq", "dup.fa"):
            (TESTS / "data3" / f).write_bytes((w / f).read_bytes())
        for name, cmd in CASES.items():
            rc, out = d.run([ref, *cmd.split()], w)
            (TESTS / "golden3" / f"{name}.out").write_bytes(out)
            print(name, rc, len(out))
        for name, cmd in SPLITS.items():
            out_dir = w / f"o_{name}"
            subprocess.run([ref, *cmd.split(), "-O", str(out_dir)], cwd=w, stdin=subprocess.DEVNULL,
                           capture_output=True, timeout=60)
            parts = []
            for p in sorted(out_dir.iterdir()):
                parts.append(f"== {p.name}\n" + p.read_text())
            (TESTS / "golden3" / f"{name}.out").write_text("".join(parts))
            print(name, len(parts), "dosya")
    return 0


if __name__ == "__main__":
    sys.exit(main())
