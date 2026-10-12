"""fx2tab and tab2fx: convert between FASTA/FASTQ and tab-separated rows."""
import sys

from fastxkit import io as fxio
from fastxkit import metrics, records
from fastxkit.commands import _common

_TAB = "\t"


def table_row(record: records.Record, args) -> str:
    """Builds one fx2tab output line for `record`.

    Numbers use Python's own two-decimal formatting (exact halves round to
    even, as in the reference tool), unlike the half-up `stats` columns.
    Known limitation: -g on protein sequences differs from the reference,
    whose protein GC calculation is not understood.
    """
    name = record.id if args.only_id else record.name
    cells = [name]
    if not args.name:
        cells += [record.seq, record.qual or ""]
    if args.length:
        cells.append(str(len(record.seq)))
    if args.gc:
        cells.append(f"{metrics.gc_percent(record.seq):.2f}")
    for bases in args.base_content:
        cells.append(f"{metrics.base_content(record.seq, bases):.2f}")
    if args.avg_qual:
        cells.append(f"{metrics.mean_quality(record.qual or ''):.2f}")
    return _TAB.join(cells) + "\n"


def run_fx2tab(args) -> None:
    for path in args.files:
        for record in records.iter_records([path]):
            try:
                line = table_row(record, args)
            except ValueError as error:
                raise ValueError(
                    f"{path}: record {record.id!r}: {error}") from None
            sys.stdout.write(line)


def _parse_table(path: str):
    with fxio.open_text(path) as text:
        source = "stdin" if path == fxio.STDIN else path
        for number, line in enumerate(text, start=1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            cells = line.split(_TAB)
            if len(cells) < 2:
                raise records.FormatError(
                    f"{source}: line {number}: expected tab-separated name and "
                    "sequence columns")
            qual = cells[2] if len(cells) > 2 and cells[2] else None
            if qual is not None and len(qual) != len(cells[1]):
                raise records.FormatError(
                    f"{source}: line {number}: sequence length "
                    f"({len(cells[1])}) differs from quality length ({len(qual)})")
            yield records.Record(cells[0], cells[1], qual)


def run_tab2fx(args) -> None:
    for path in args.files:
        for record in _parse_table(path):
            sys.stdout.write(records.format_record(record))


def register(subparsers) -> None:
    fx2tab = subparsers.add_parser("fx2tab", help="convert FASTA/FASTQ to a table")
    _common.add_input_files(fx2tab)
    fx2tab.add_argument("-l", "--length", action="store_true",
                        help="add sequence length column")
    fx2tab.add_argument("-g", "--gc", action="store_true",
                        help="add GC percentage column (N excluded)")
    fx2tab.add_argument("-n", "--name", action="store_true",
                        help="print names only (no sequence/quality columns)")
    fx2tab.add_argument("-i", "--only-id", action="store_true",
                        help="use the ID (first word) instead of the full name")
    fx2tab.add_argument("-q", "--avg-qual", action="store_true",
                        help="add average quality column")
    fx2tab.add_argument("-B", "--base-content", action="append", default=[],
                        metavar="BASES",
                        help="add the percentage of these bases, e.g. AT (repeatable)")
    fx2tab.set_defaults(func=run_fx2tab)

    tab2fx = subparsers.add_parser("tab2fx", help="convert a table to FASTA/FASTQ")
    _common.add_input_files(tab2fx)
    tab2fx.set_defaults(func=run_tab2fx)
