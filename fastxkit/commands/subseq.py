"""subseq: extract a region from every record."""
import dataclasses
import sys

from fastxkit import records
from fastxkit.commands import _common, select


def slice_record(record: records.Record, start: int, end: int) -> records.Record:
    """Cuts 1-based inclusive start..end; negatives count from the end.

    Positions outside the sequence are clamped, so the result may be empty.
    """
    length = len(record.seq)
    first = max((start if start > 0 else length + start + 1) - 1, 0)
    last = min(end if end > 0 else length + end + 1, length)
    qual = None if record.qual is None else record.qual[first:last]
    return dataclasses.replace(record, seq=record.seq[first:last], qual=qual)


def _format(record: records.Record, width: int) -> str:
    text = records.format_record(record, width)
    # An empty FASTA sequence keeps its (blank) sequence line.
    return text + "\n" if record.qual is None and not record.seq else text


def run(args) -> None:
    start, end = args.range
    for record in records.iter_records(args.files):
        sys.stdout.write(_format(slice_record(record, start, end), args.line_width))


def register(subparsers) -> None:
    parser = subparsers.add_parser("subseq", help="extract a region from records")
    _common.add_input_files(parser)
    parser.add_argument(
        "-r", "--range", type=select.parse_range, required=True,
        metavar="START:END",
        help="1-based inclusive region; negative values count from the end")
    parser.add_argument("-w", "--line-width", type=_common.non_negative_int,
                        default=records.DEFAULT_FASTA_WIDTH,
                        help="FASTA line width, 0 for none (default: %(default)s)")
    parser.set_defaults(func=run)
