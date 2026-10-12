"""sum: order-independent content digest of the sequences in a file.

Output line: `fastxkit.v1_<type><structure>_<md5hex>  <file>`, where type is
the detected DNA/RNA/Protein/Unlimit/- and structure is `s` (sequences only)
or `i` (IDs included, -i). Not compatible with seqkit sum.

The digest covers the lower-cased sequences (gaps removed with -g), never
headers or qualities. Each record is hashed with blake2b, the hashes are
sorted, and md5 is taken over the sorted hashes, the sorted lengths and the
record count.
"""
import hashlib
import sys

from fastxkit import metrics, records
from fastxkit.commands import _common

_DIGEST_BYTES = 16
_ID_SEPARATOR = "\0"


def _sequence_text(record: records.Record, args) -> str:
    text = record.seq.lower()
    if args.remove_gaps:
        text = text.translate(str.maketrans("", "", metrics.GAP_CHARS))
    return text


def compute_sum(path: str, args) -> str:
    """Returns the sum line (without newline) for one file."""
    hashes = []
    lengths = []
    letters: set[str] = set()
    for record in records.iter_records([path]):
        text = _sequence_text(record, args)
        letters.update(text)
        lengths.append(len(text))
        hashed = f"{record.id}{_ID_SEPARATOR}{text}" if args.with_id else text
        hashes.append(hashlib.blake2b(
            hashed.encode(), digest_size=_DIGEST_BYTES).hexdigest())
    combined = hashlib.md5(
        "\n".join([*sorted(hashes), *map(str, sorted(lengths)),
                   str(len(hashes))]).encode()).hexdigest()
    structure = "i" if args.with_id else "s"
    seq_type = metrics.detect_type(letters)
    return f"fastxkit.v1_{seq_type}{structure}_{combined}  {path}"


def run(args) -> None:
    for path in args.files:
        sys.stdout.write(compute_sum(path, args) + "\n")


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "sum", help="order-independent checksum of sequence content")
    _common.add_input_files(parser)
    parser.add_argument("-i", "--with-id", action="store_true",
                        help="include record IDs in the digest")
    parser.add_argument("-g", "--remove-gaps", action="store_true",
                        help="remove gap characters ('-', '.', ' ') first")
    parser.set_defaults(func=run)
