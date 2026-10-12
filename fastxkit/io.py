"""Input layer: open possibly compressed files or stdin as UTF-8 text."""
import bz2
import contextlib
import gzip
import io
import lzma
import sys
import zlib
from collections.abc import Iterator
from typing import BinaryIO, TextIO

STDIN = "-"
# ValueError covers FormatError and UnicodeDecodeError; the others are what
# the gzip/bz2/lzma readers raise on truncated or corrupt streams.
DATA_ERRORS = (OSError, EOFError, ValueError, lzma.LZMAError, zlib.error)
# Raised while reading a damaged stream; gzip.BadGzipFile is an OSError.
_CORRUPT_DATA = (EOFError, lzma.LZMAError, zlib.error, UnicodeDecodeError,
                 gzip.BadGzipFile)
_MAGIC_LENGTH = 6
_DECOMPRESSORS = (
    (b"\x1f\x8b", lambda raw: gzip.GzipFile(fileobj=raw)),
    (b"BZh", bz2.BZ2File),
    (b"\xfd7zXZ\x00", lzma.LZMAFile),
)


def _label(path: str) -> str:
    return "stdin" if path == STDIN else path


def _wrap_decompressor(raw: io.BufferedReader) -> BinaryIO:
    # peek() may return fewer bytes than asked on a slow pipe; the magic
    # numbers are short enough that the first chunk practically always has them.
    head = raw.peek(_MAGIC_LENGTH)[:_MAGIC_LENGTH]
    for magic, opener in _DECOMPRESSORS:
        if head.startswith(magic):
            return opener(raw)
    return raw


@contextlib.contextmanager
def open_text(path: str) -> Iterator[TextIO]:
    """Opens `path` ('-' is stdin) as text, detecting compression by content.

    A leading UTF-8 byte order mark is dropped and line endings are
    normalised to '\\n'. Stdin is never closed. Corrupt compressed data and
    bad UTF-8 are re-raised as ValueError naming the file.
    """
    with contextlib.ExitStack() as stack:
        if path == STDIN:
            raw = sys.stdin.buffer
        else:
            raw = stack.enter_context(open(path, "rb"))
        stream = _wrap_decompressor(raw)
        if stream is not raw:
            stack.callback(stream.close)
        text = io.TextIOWrapper(stream, encoding="utf-8-sig", newline=None)
        try:
            yield text
        except _CORRUPT_DATA as error:
            raise ValueError(f"{_label(path)}: {error}") from error
        except OSError as error:
            # The bz2 reader reports bad data as a bare OSError without an
            # errno; real I/O failures (and BrokenPipeError) carry one.
            if error.errno is not None:
                raise
            raise ValueError(f"{_label(path)}: {error}") from error
        finally:
            text.detach()
