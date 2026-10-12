"""qc: independently computed, verifiable QC numbers per sequence file.

The counts come straight from the file contents and use the same
definitions as `stats -a`, so the two commands always agree. `--olgu`
emits facts in the `{"deger", "kaynak", "konum"}` shape used by olgu.json.
"""
import argparse
import hashlib
import json
import pathlib
import sys

from fastxkit import __version__, metrics, numfmt
from fastxkit import io as fxio
from fastxkit.commands import stats

_CHUNK_BYTES = 1 << 20
_JSON_DECIMALS = 6
_COMPRESSION_SUFFIXES = frozenset((".gz", ".bz2", ".xz"))
_TABLE_COLUMNS = ("file", "reads", "total_bases", "min_len", "mean_len",
                  "max_len", "GC(%)", "N", "Q20(%)", "Q30(%)", "mean_qual",
                  "sha256", "version")
_OLGU_FACTS = (  # (olgu key, qc field)
    ("fastq_okuma_sayisi", "reads"),
    ("fastq_toplam_baz", "total_bases"),
    ("fastq_gc_yuzde", "gc_percent"),
    ("fastq_q30_yuzde", "q30_percent"),
    ("fastq_ort_uzunluk", "mean_len"),
)


def _file_argument(text: str) -> str:
    """argparse type for a file name; qc needs seekable files for sha256."""
    if text == fxio.STDIN:
        raise argparse.ArgumentTypeError(
            "stdin ('-') is not supported; pass a file name")
    return text


def sha256_of_file(path: str) -> str:
    """SHA-256 of the raw (still compressed) file bytes."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def compute_qc(path: str) -> dict:
    """Returns the unrounded QC numbers of one file."""
    fmt, totals = stats.collect(path)
    lengths = totals.lengths
    total = sum(lengths)
    is_protein = metrics.detect_type(totals.letters) == "Protein"
    # G, C and N are amino acids in proteins, as in `stats`.
    gc = 0 if is_protein else totals.gc
    n_bases = 0 if is_protein else totals.n_bases
    quality_bases = totals.quality_bases
    return {
        "file": path,
        "format": fmt or "-",
        "reads": len(lengths),
        "total_bases": total,
        "min_len": min(lengths, default=0),
        "mean_len": total / len(lengths) if lengths else 0,
        "max_len": max(lengths, default=0),
        "gc_percent": stats.percent(gc, total),
        "n_bases": n_bases,
        "q20_percent": stats.percent(totals.at_least_20, quality_bases),
        "q30_percent": stats.percent(totals.at_least_30, quality_bases),
        "mean_quality": metrics.error_to_quality(
            totals.error_sum, quality_bases),
        "sha256": sha256_of_file(path),
        "tool_version": __version__,
    }


def _rounded(value):
    """Rounds floats half up to the JSON precision; ints pass through."""
    if isinstance(value, float):
        return float(numfmt.fixed(value, _JSON_DECIMALS))
    return value


def sample_name(path: str) -> str:
    """File name without directory, compression and sequence extension."""
    file_path = pathlib.PurePath(path)
    if file_path.suffix in _COMPRESSION_SUFFIXES:
        file_path = file_path.with_suffix("")
    return file_path.stem


def olgu_facts(results: list[dict]) -> dict:
    """Builds olgu.json style facts, keyed `<fact>:<sample name>`."""
    facts: dict = {}
    for result in results:
        name = sample_name(result["file"])
        location = (f"fastxkit qc {result['tool_version']}; "
                    f"sha256 {result['sha256']}")
        for key, field in _OLGU_FACTS:
            full_key = f"{key}:{name}"
            if full_key in facts:
                raise ValueError(
                    f"two input files share the sample name {name!r}; "
                    "rename one so olgu keys stay unique")
            facts[full_key] = {"deger": _rounded(result[field]),
                               "kaynak": result["file"], "konum": location}
    return facts


def _table_row(result: dict) -> dict[str, str]:
    return {
        "file": result["file"],
        "reads": numfmt.grouped(result["reads"]),
        "total_bases": numfmt.grouped(result["total_bases"]),
        "min_len": numfmt.grouped(result["min_len"]),
        "mean_len": numfmt.trimmed(result["mean_len"], 1),
        "max_len": numfmt.grouped(result["max_len"]),
        "GC(%)": numfmt.trimmed(result["gc_percent"], 2),
        "N": numfmt.grouped(result["n_bases"]),
        "Q20(%)": numfmt.trimmed(result["q20_percent"], 2),
        "Q30(%)": numfmt.trimmed(result["q30_percent"], 2),
        "mean_qual": numfmt.trimmed(result["mean_quality"], 2),
        "sha256": result["sha256"],
        "version": result["tool_version"],
    }


def run(args) -> None:
    results = [compute_qc(path) for path in args.files]
    if args.olgu:
        output = json.dumps(olgu_facts(results), indent=2,
                            ensure_ascii=False) + "\n"
    elif args.json:
        output = json.dumps(
            [{key: _rounded(value) for key, value in result.items()}
             for result in results], indent=2) + "\n"
    else:
        output = stats.render_table(
            _TABLE_COLUMNS, [_table_row(result) for result in results])
    sys.stdout.write(output)


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "qc", help="verifiable QC numbers (reads, GC, Q30, sha256) per file")
    parser.add_argument("files", nargs="+", metavar="FILE", type=_file_argument,
                        help="FASTA/FASTQ files, possibly gz/bz2/xz compressed "
                             "(stdin is not supported)")
    output = parser.add_mutually_exclusive_group()
    output.add_argument("--json", action="store_true",
                        help="JSON list, floats rounded half up to 6 decimals")
    output.add_argument("--olgu", action="store_true",
                        help="olgu.json style facts {deger, kaynak, konum}")
    parser.set_defaults(func=run)
