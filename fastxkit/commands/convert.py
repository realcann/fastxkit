"""fq2fa: convert FASTQ to FASTA."""
import contextlib
import sys

from fastxkit import records
from fastxkit.commands import _common


def run_fq2fa(args) -> None:
    if args.out_file:
        _common.refuse_same_file(args.out_file, args.files)
    with contextlib.ExitStack() as stack:
        out = (stack.enter_context(open(args.out_file, "w", encoding="utf-8"))
               if args.out_file else sys.stdout)
        for record in records.iter_records(args.files):
            # Unwrapped FASTA lines, like the reference tool.
            out.write(records.format_record(
                records.Record(record.name, record.seq), width=0))


def register(subparsers) -> None:
    parser = subparsers.add_parser("fq2fa", help="convert FASTQ to FASTA")
    _common.add_input_files(parser)
    parser.add_argument("-o", "--out-file", metavar="FILE",
                        help="write to FILE instead of stdout")
    parser.set_defaults(func=run_fq2fa)
