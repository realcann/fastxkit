"""FASTA/FASTQ records: streaming parser and serialiser."""
import contextlib
import dataclasses
import itertools
from collections.abc import Iterable, Iterator

from fastxkit import io as fxio

FASTA = "FASTA"
FASTQ = "FASTQ"
DEFAULT_FASTA_WIDTH = 60


class FormatError(ValueError):
    """Raised when input is not valid FASTA/FASTQ."""


@dataclasses.dataclass(frozen=True)
class Record:
    """One sequence record; `qual` is None for FASTA."""

    name: str
    seq: str
    qual: str | None = None

    @property
    def id(self) -> str:
        """First whitespace-delimited word of the name."""
        words = self.name.split(None, 1)
        return words[0] if words else ""


class FastxReader:
    """Lazily parses FASTA or FASTQ text lines; the format comes from line 1.

    Attributes:
        format: "FASTA", "FASTQ", or None when the input has no records.
    """

    def __init__(self, lines: Iterable[str], source: str = "<input>"):
        self._lines = enumerate(lines, start=1)
        self._source = source
        self._first = self._next_nonblank()
        self.format = None
        if self._first is not None:
            self.format = self._detect(self._first)

    def _error(self, number: int, message: str) -> FormatError:
        return FormatError(f"{self._source}: line {number}: {message}")

    def _next_nonblank(self) -> tuple[int, str] | None:
        for number, line in self._lines:
            line = line.strip()
            if line:
                return number, line
        return None

    def _detect(self, first: tuple[int, str]) -> str:
        number, line = first
        if line[0] == ">":
            return FASTA
        if line[0] == "@":
            return FASTQ
        raise self._error(
            number, f"expected a header starting with '>' or '@', got {line[:20]!r}")

    def __iter__(self) -> Iterator[Record]:
        if self._first is None:
            return iter(())
        if self.format == FASTA:
            return self._parse_fasta()
        return self._parse_fastq()

    def _parse_fasta(self) -> Iterator[Record]:
        name = None
        chunks: list[str] = []
        for _, line in itertools.chain([self._first], self._lines):
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield Record(name, "".join(chunks))
                name, chunks = line[1:], []
            else:
                chunks.append(line)
        yield Record(name, "".join(chunks))

    def _parse_fastq(self) -> Iterator[Record]:
        entry = self._first
        while entry is not None:
            number, header = entry
            if header[0] != "@":
                raise self._error(number, "expected '@' at start of FASTQ record")
            body = [next(self._lines, None) for _ in range(3)]
            if None in body:
                raise self._error(number, "truncated FASTQ record")
            seq, plus, qual = (line.strip() for _, line in body)
            if not plus.startswith("+"):
                raise self._error(number + 2, "expected '+' separator line")
            if len(seq) != len(qual):
                raise self._error(
                    number,
                    f"sequence length ({len(seq)}) differs from quality "
                    f"length ({len(qual)}) in record {header[1:]!r}")
            yield Record(header[1:], seq, qual)
            entry = self._next_nonblank()


@contextlib.contextmanager
def open_fastx(path: str) -> Iterator[FastxReader]:
    """Opens a (possibly compressed) file or '-' as a FastxReader."""
    with fxio.open_text(path) as text:
        yield FastxReader(text, "stdin" if path == fxio.STDIN else path)


def iter_records(paths: Iterable[str]) -> Iterator[Record]:
    """Yields records from each path in turn, one file open at a time."""
    for path in paths:
        with open_fastx(path) as reader:
            yield from reader


def format_record(record: Record, width: int = DEFAULT_FASTA_WIDTH) -> str:
    """Serialises a record; `width` wraps FASTA sequence lines (0 = no wrap)."""
    if record.qual is not None:
        return f"@{record.name}\n{record.seq}\n+\n{record.qual}\n"
    if width <= 0 or len(record.seq) <= width:
        body = [record.seq] if record.seq else []
    else:
        body = [record.seq[i:i + width] for i in range(0, len(record.seq), width)]
    return "".join(f"{line}\n" for line in [f">{record.name}", *body])
