"""Compare fastxkit output with a reference seqkit binary on freshly generated data.

Usage: python tools/diff_reference.py <reference-binary> [--seed N]
Runs only on synthetic data in a temp dir; the reference is never given network-facing commands.
"""
import argparse
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

CASES = [
    "stats big.fa big.fq rna.fa prot.fa",
    "stats -a big.fa big.fq rna.fa prot.fa",
    "fx2tab -l -g -q big.fq",
    "fx2tab -l -g big.fa prot.fa",
    "fx2tab -n -i big.fa",
    "seq -r -p big.fa", "seq -u big.fa", "seq -l -w 40 big.fa", "seq -m 300 -M 600 big.fa",
    "seq -n -i big.fa", "seq -s -w 0 big.fq",
    "head -n 7 big.fa", "range -r 5:20 big.fq", "range -r -5:-1 big.fa",
    "grep -p c1 big.fa", "grep -p c1 -p c2 big.fa", "grep -r -p c1 big.fa", "grep -r -p ^c1 big.fa",
    "grep -n -p c1 big.fa", "grep -n -r -p len1 big.fa",
    "grep -s -p ACGTA big.fa", "grep -s -i -p acgta big.fa", "grep -v -r -p c1 big.fa",
    "grep -s -r -p A{5} big.fq",
    "locate -p ACGT big.fa", "locate -P -p ACGT big.fa", "locate -r -p A+C big.fa",
    "locate -i -p acgt big.fq",
    "subseq -r 5:30 big.fa", "subseq -r -20:-1 big.fq", "subseq -r 1:1 big.fa",
    "sort -n uniq.fa", "sort -l uniq.fa", "sort -l -r uniq.fa", "sort -s uniq.fa", "sort uniq.fa",
    "sort -i -n uniq.fa", "sort -l big.fq",
    "rmdup -s dup.fa", "rmdup -s -i dup.fa", "rmdup -n dup.fa", "rmdup dup.fa",
    "rename dup.fa", "rename big.fa",
    "replace -p c -r X big.fa", "replace -s -p A{3} -r - big.fa", "replace -p ^(\\w)(\\d) -r ${2}_$1 big.fa",
    "replace -s -i -p acg -r NNN big.fa",
]
STDIN_SKIP = ("faidx",)


def make_data(d: Path, seed: int) -> None:
    rnd = random.Random(seed)

    def dna(n, gc=0.5):
        return "".join(rnd.choice("GC" if rnd.random() < gc else "AT") for _ in range(n))
    with (d / "big.fa").open("w") as f:
        for i in range(40):
            n = rnd.randint(50, 900)
            s = dna(n, rnd.random())
            if i % 7 == 0:
                s = s[:10] + "NN" + s[12:]
            if i % 5 == 0:
                s = s.lower()
            f.write(f">c{i} len{n}\n" + "\n".join(s[j:j + 70] for j in range(0, n, 70)) + "\n")
    with (d / "big.fq").open("w") as f:
        for i in range(300):
            n = rnd.randint(30, 250)
            q = "".join(chr(33 + rnd.randint(0, 41)) for _ in range(n))
            f.write(f"@q{i} x\n{dna(n)}\n+\n{q}\n")
    (d / "rna.fa").write_text(">r1\nACGUACGUAGCU\n>r2\nUUUUGGGCCC\n")
    (d / "prot.fa").write_text(">p1\nMKVLAAGIVGLLLAQPSEA\n>p2\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n")
    (d / "dup.fa").write_text(">a one\nACGTACGT\n>b two\nACGTACGT\n>a three\nGGGG\n>d four\nacgtacgt\n>e five\nTTGCA\n")
    (d / "uniq.fa").write_text(">z1 zeta\nCCCCGGGG\n>a1 alpha\nACGTACGTAC\n>m1 mid\nGGGGGGGGGGGGGGGG\n"
                               ">b1 beta\nACGTACGTAC\n>c1 gamma\nTTT\n>B2 upper\nAC\n")


def run(cmd, cwd):
    p = subprocess.run(cmd, cwd=cwd, stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
    return p.returncode, p.stdout


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    bad = 0
    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        make_data(d, a.seed)
        for case in CASES:
            args = case.split()
            rc_ref, out_ref = run([a.reference, *args], d)
            rc_our, out_our = run([sys.executable, "-m", "fastxkit", *args], d)
            ok = out_ref == out_our and (rc_ref == 0) == (rc_our == 0)
            if not ok:
                bad += 1
                print(f"FARKLI: {case}  (rc ref={rc_ref} bizim={rc_our})")
        print(f"{len(CASES) - bad}/{len(CASES)} eşleşti")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
