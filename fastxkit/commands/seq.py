"""seq: transform sequences and filter them by length."""
import dataclasses
import sys

from fastxkit import metrics, records
from fastxkit.commands import _common

_COMPLEMENT_PAIRS = ("AT", "CG", "RY", "KM", "BV", "DH")


def _build_complement_table() -> dict[int, str]:
    mapping = {"U": "A"}
    for first, second in _COMPLEMENT_PAIRS:
        mapping[first], mapping[second] = second, first
    mapping.update({k.lower(): v.lower() for k, v in mapping.items()})
    return str.maketrans(mapping)


COMPLEMENT = _build_complement_table()


def _remove_gaps(record: records.Record) -> records.Record:
    keep = [i for i, base in enumerate(record.seq) if base not in metrics.GAP_CHARS]
    if len(keep) == len(record.seq):
        return record
    qual = None if record.qual is None else "".join(record.qual[i] for i in keep)
    return dataclasses.replace(
        record, seq="".join(record.seq[i] for i in keep), qual=qual)


def _reverse(record: records.Record) -> records.Record:
    qual = None if record.qual is None else record.qual[::-1]
    return dataclasses.replace(record, seq=record.seq[::-1], qual=qual)


def transform(record: records.Record, args) -> records.Record:
    """Applies gap removal, reverse, complement and case change, in that order."""
    if args.remove_gaps:
        record = _remove_gaps(record)
    if args.reverse:
        record = _reverse(record)
    if args.complement:
        record = dataclasses.replace(record, seq=record.seq.translate(COMPLEMENT))
    if args.upper_case:
        record = dataclasses.replace(record, seq=record.seq.upper())
    if args.lower_case:
        record = dataclasses.replace(record, seq=record.seq.lower())
    return record


def _in_length_bounds(record: records.Record, args) -> bool:
    length = len(record.seq)
    if args.min_len is not None and length < args.min_len:
        return False
    return args.max_len is None or length <= args.max_len


def _render(record: records.Record, args) -> str:
    if args.name:
        return f"{record.id if args.only_id else record.name}\n"
    if args.seq:
        return f"{record.seq}\n"
    if args.only_id:
        record = dataclasses.replace(record, name=record.id)
    return records.format_record(record, args.line_width)


def run(args) -> None:
    for record in records.iter_records(args.files):
        if _in_length_bounds(record, args):
            sys.stdout.write(_render(transform(record, args), args))


def register(subparsers) -> None:
    parser = subparsers.add_parser("seq", help="transform and filter sequences")
    _common.add_input_files(parser)
    parser.add_argument("-r", "--reverse", action="store_true",
                        help="reverse sequence (and quality)")
    parser.add_argument("-p", "--complement", action="store_true",
                        help="complement sequence")
    case = parser.add_mutually_exclusive_group()
    case.add_argument("-u", "--upper-case", action="store_true")
    case.add_argument("-l", "--lower-case", action="store_true")
    parser.add_argument("-n", "--name", action="store_true",
                        help="print names only")
    parser.add_argument("-i", "--only-id", action="store_true",
                        help="use the ID (first word) instead of the full name")
    parser.add_argument("-s", "--seq", action="store_true",
                        help="print sequences only")
    parser.add_argument("-w", "--line-width", type=_common.non_negative_int,
                        default=records.DEFAULT_FASTA_WIDTH,
                        help="FASTA line width, 0 for none (default: %(default)s)")
    parser.add_argument("-m", "--min-len", type=_common.non_negative_int,
                        help="keep records at least this long (original length, before gap removal)")
    parser.add_argument("-M", "--max-len", type=_common.non_negative_int,
                        help="keep records at most this long (original length, before gap removal)")
    parser.add_argument("-g", "--remove-gaps", action="store_true",
                        help="remove gap characters ('-', '.', ' ')")
    parser.set_defaults(func=run)
