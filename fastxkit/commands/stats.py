"""stats: per-file summary table."""
import collections
import dataclasses
import sys
from collections.abc import Collection

from fastxkit import metrics, numfmt, records
from fastxkit.commands import _common

BASIC_COLUMNS = ("file", "format", "type", "num_seqs", "sum_len",
                 "min_len", "avg_len", "max_len")
EXTRA_COLUMNS = ("Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num",
                 "Q20(%)", "Q30(%)", "AvgQual", "GC(%)", "sum_n")
_TEXT_COLUMNS = frozenset(("file", "format", "type"))
_SEPARATOR = "  "


@dataclasses.dataclass
class Totals:
    """Running sums over the records of one file."""

    lengths: list[int] = dataclasses.field(default_factory=list)
    letters: set[str] = dataclasses.field(default_factory=set)
    gc: int = 0
    n_bases: int = 0
    gaps: int = 0
    error_sum: float = 0.0
    quality_bases: int = 0
    at_least_20: int = 0
    at_least_30: int = 0

    def add(self, record: records.Record) -> None:
        self.lengths.append(len(record.seq))
        counts = collections.Counter(record.seq)
        self.letters.update(counts)
        self.gc += sum(counts[base] for base in "GCgc")
        self.n_bases += counts["N"] + counts["n"]
        self.gaps += sum(counts[char] for char in metrics.GAP_CHARS)
        if record.qual:
            error_sum, at_20, at_30 = metrics.quality_profile(record.qual)
            self.error_sum += error_sum
            self.quality_bases += len(record.qual)
            self.at_least_20 += at_20
            self.at_least_30 += at_30


def _row(path: str, fmt: str | None, totals: Totals) -> dict[str, str]:
    lengths = totals.lengths
    total = sum(lengths)
    count = len(lengths)
    q1, q2, q3 = metrics.quartiles(lengths)
    n50, n50_count = metrics.n50(lengths)
    quality = metrics.error_to_quality(totals.error_sum, totals.quality_bases)
    seq_type = metrics.detect_type(totals.letters)
    # G, C and N are amino acids in proteins, so nucleotide counts are meaningless.
    is_protein = seq_type == "Protein"
    gc = 0 if is_protein else totals.gc
    n_bases = 0 if is_protein else totals.n_bases
    return {
        "file": path,
        "format": fmt or "-",
        "type": seq_type,
        "num_seqs": numfmt.grouped(count),
        "sum_len": numfmt.grouped(total),
        "min_len": numfmt.grouped(min(lengths, default=0)),
        "avg_len": numfmt.trimmed(total / count if count else 0, 1),
        "max_len": numfmt.grouped(max(lengths, default=0)),
        "Q1": numfmt.trimmed(q1, 1),
        "Q2": numfmt.trimmed(q2, 1),
        "Q3": numfmt.trimmed(q3, 1),
        "sum_gap": numfmt.grouped(totals.gaps),
        "N50": numfmt.grouped(n50),
        "N50_num": numfmt.grouped(n50_count),
        "Q20(%)": numfmt.trimmed(percent(totals.at_least_20, totals.quality_bases), 2),
        "Q30(%)": numfmt.trimmed(percent(totals.at_least_30, totals.quality_bases), 2),
        "AvgQual": numfmt.trimmed(quality, 2),
        "GC(%)": numfmt.trimmed(percent(gc, total), 2),
        "sum_n": numfmt.grouped(n_bases),
    }


def collect(path: str) -> tuple[str | None, Totals]:
    """Reads one file and returns (format, running totals)."""
    totals = Totals()
    with records.open_fastx(path) as reader:
        for record in reader:
            try:
                totals.add(record)
            except ValueError as error:
                raise ValueError(
                    f"{path}: record {record.id!r}: {error}") from None
        return reader.format, totals


def percent(part: int, whole: int) -> float:
    """Percentage of `part` in `whole`; 0 when `whole` is 0."""
    return part / whole * 100 if whole else 0


def _summarize(path: str) -> dict[str, str]:
    fmt, totals = collect(path)
    return _row(path, fmt, totals)


def render_table(columns: tuple[str, ...], rows: list[dict[str, str]],
                 text_columns: Collection[str] = _TEXT_COLUMNS) -> str:
    """Aligns columns: text left, numbers right, two spaces apart.

    Args:
        columns: column names, in output order.
        rows: one name -> text mapping per line.
        text_columns: names of the left-aligned columns.
    """
    widths = {column: max([len(column)] + [len(row[column]) for row in rows])
              for column in columns}

    def line(values: dict[str, str]) -> str:
        cells = [values[c].ljust(widths[c]) if c in text_columns
                 else values[c].rjust(widths[c]) for c in columns]
        return _SEPARATOR.join(cells).rstrip() + "\n"

    header = {column: column for column in columns}
    return "".join(line(values) for values in [header, *rows])


def run(args) -> None:
    columns = BASIC_COLUMNS + EXTRA_COLUMNS if args.all else BASIC_COLUMNS
    rows = [_summarize(path) for path in args.files]
    sys.stdout.write(render_table(columns, rows))


def register(subparsers) -> None:
    parser = subparsers.add_parser("stats", help="summarise sequence files")
    _common.add_input_files(parser)
    parser.add_argument("-a", "--all", action="store_true",
                        help="add quartiles, N50, quality and composition columns")
    parser.set_defaults(func=run)
