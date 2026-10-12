"""faidx: build a samtools-style .fai index and extract regions from FASTA."""
import dataclasses
import os
import re
import sys
from typing import BinaryIO

from fastxkit import records
from fastxkit.commands import _common

INDEX_SUFFIX = ".fai"
_TAB = "\t"
_REGION = re.compile(r"(.+):(\d+)-(\d+)")


@dataclasses.dataclass
class IndexEntry:
    """One .fai line: `offset` is the byte of the first base."""

    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


def build_index(path: str) -> list[IndexEntry]:
    """Scans a plain FASTA file and returns one entry per record.

    Raises:
        ValueError: on non-FASTA input, repeated names, or sequence lines of
            uneven length (which make byte offsets unpredictable).
    """
    entries: list[IndexEntry] = []
    names: set[str] = set()
    entry = None
    closed = False  # the record has seen a short or blank line
    position = 0
    with open(path, "rb") as stream:
        for number, raw in enumerate(stream, start=1):
            position += len(raw)
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                words = line[1:].decode("utf-8").split(None, 1)
                if not words:
                    raise ValueError(f"{path}: line {number}: empty sequence name")
                if words[0] in names:
                    raise ValueError(
                        f"{path}: line {number}: duplicate sequence name "
                        f"{words[0]!r}")
                names.add(words[0])
                entry = IndexEntry(words[0], 0, position, 0, 0)
                entries.append(entry)
                closed = False
            elif entry is None:
                if line.strip():
                    raise ValueError(
                        f"{path}: line {number}: expected a FASTA header "
                        "starting with '>' (compressed or FASTQ input cannot "
                        "be indexed)")
            elif not line:
                closed = True
            else:
                if closed or (entry.line_bases and len(line) > entry.line_bases):
                    raise ValueError(
                        f"{path}: line {number}: sequence lines of "
                        f"{entry.name!r} must all have the same length, "
                        "except the last")
                if not entry.line_bases:
                    entry.line_bases = len(line)
                    entry.line_width = len(raw) if raw != line else len(line)
                closed = len(line) < entry.line_bases
                entry.length += len(line)
    return entries


def write_index(path: str, entries: list[IndexEntry]) -> None:
    """Writes `<path>.fai` in samtools' five-column format.

    An existing regular index is replaced; a symbolic link in its place is an
    error and its target is left alone.
    """
    with _common.open_new_text_file(
            path + INDEX_SUFFIX, replace_regular=True) as stream:
        for entry in entries:
            stream.write(_TAB.join(map(str, dataclasses.astuple(entry))) + "\n")


def read_index(path: str) -> dict[str, IndexEntry]:
    """Reads an index file into a name -> entry mapping."""
    index: dict[str, IndexEntry] = {}
    with open(path, encoding="utf-8") as stream:
        for number, line in enumerate(stream, start=1):
            cells = line.rstrip("\r\n").split(_TAB)
            try:
                if len(cells) < 5:
                    raise ValueError
                entry = IndexEntry(cells[0], *map(int, cells[1:5]))
            except ValueError:
                raise ValueError(
                    f"{path}: line {number}: malformed index line") from None
            index[entry.name] = entry
    return index


def parse_region(text: str, index: dict[str, IndexEntry]
                 ) -> tuple[IndexEntry, int, int]:
    """Resolves `name` or `name:start-end` (1-based, inclusive) to a region.

    An end past the sequence is clamped.

    Raises:
        ValueError: for an unknown name or a region outside the sequence.
    """
    if text in index:
        entry = index[text]
        return entry, 1, entry.length
    match = _REGION.fullmatch(text)
    if not match or match.group(1) not in index:
        name = match.group(1) if match else text
        raise ValueError(f"sequence {name!r} not found in the index")
    entry = index[match.group(1)]
    start, end = int(match.group(2)), int(match.group(3))
    if start < 1 or end < start:
        raise ValueError(f"invalid region {text!r}: need 1 <= start <= end")
    if start > entry.length:
        raise ValueError(
            f"region {text!r} starts beyond the end of {entry.name!r} "
            f"(length {entry.length})")
    return entry, start, min(end, entry.length)


def _byte_position(entry: IndexEntry, base: int) -> int:
    lines, column = divmod(base, entry.line_bases)
    return entry.offset + lines * entry.line_width + column


def fetch(stream: BinaryIO, entry: IndexEntry, start: int, end: int) -> str:
    """Reads bases start..end (1-based, inclusive) of a valid region."""
    begin = _byte_position(entry, start - 1)
    stream.seek(begin)
    raw = stream.read(_byte_position(entry, end - 1) + 1 - begin)
    return raw.translate(None, b"\r\n").decode("utf-8")


def _is_fresh(path: str, index: dict[str, IndexEntry]) -> bool:
    """False if the index predates the FASTA or points beyond its end."""
    fasta = os.stat(path)
    if os.stat(path + INDEX_SUFFIX).st_mtime_ns < fasta.st_mtime_ns:
        return False
    return all(entry.offset + entry.length <= fasta.st_size
               for entry in index.values())


def run(args) -> None:
    index_path = args.file + INDEX_SUFFIX
    index = None
    if args.regions and os.path.exists(index_path):
        index = read_index(index_path)
        if not _is_fresh(args.file, index):
            index = None
    if index is None:
        entries = build_index(args.file)
        write_index(args.file, entries)
        index = {entry.name: entry for entry in entries}
    with open(args.file, "rb") as stream:
        for text in args.regions:
            entry, start, end = parse_region(text, index)
            sequence = fetch(stream, entry, start, end) if entry.length else ""
            sys.stdout.write(records.format_record(
                records.Record(text, sequence), args.line_width))


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "faidx", help="index a FASTA file and extract regions")
    parser.add_argument("file", metavar="FILE",
                        help="plain (uncompressed) FASTA file")
    parser.add_argument(
        "regions", nargs="*", metavar="REGION",
        help="NAME or NAME:START-END (1-based, inclusive); without any, "
             "FILE.fai is (re)written; an existing FILE.fai that is older "
             "than FILE or does not fit it is rebuilt")
    parser.add_argument("-w", "--line-width", type=_common.non_negative_int,
                        default=records.DEFAULT_FASTA_WIDTH,
                        help="FASTA line width, 0 for none (default: %(default)s)")
    parser.set_defaults(func=run)
