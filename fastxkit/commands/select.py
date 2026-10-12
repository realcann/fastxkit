"""head and range: pick records by position."""
import argparse
import collections
import contextlib
import itertools
import re
import sys
from collections.abc import Iterable, Iterator
from typing import TypeVar

from fastxkit import records
from fastxkit.commands import _common

T = TypeVar("T")
_RANGE = re.compile(r"(-?\d+):(-?\d+)")
DEFAULT_HEAD = 10


def head_records(items: Iterable[T], count: int) -> Iterator[T]:
    """Yields the first `count` items without reading further."""
    return itertools.islice(items, count)


def select_range(items: Iterable[T], start: int, end: int) -> Iterator[T]:
    """Yields items start..end (1-based, inclusive); negatives count from the end.

    Memory stays bounded by the number of items that must be held back to
    resolve a negative index, never the whole input.
    """
    if start > 0 and end > 0:
        yield from itertools.islice(items, start - 1, end)
    elif start > 0:
        # Hold back the final -end-1 items; they are excluded.
        pending: collections.deque = collections.deque()
        for index, item in enumerate(items, start=1):
            pending.append((index, item))
            if len(pending) > -end - 1:
                old_index, old_item = pending.popleft()
                if old_index >= start:
                    yield old_item
    elif end < 0:
        tail = list(collections.deque(items, maxlen=-start))
        yield from tail[:max(len(tail) + end + 1, 0)]
    else:
        head = list(itertools.islice(items, end))
        total = len(head) + sum(1 for _ in items)
        yield from head[max(total + start, 0):]


def parse_range(text: str) -> tuple[int, int]:
    """argparse type for 'START:END' with non-zero 1-based indices."""
    match = _RANGE.fullmatch(text)
    if not match:
        raise argparse.ArgumentTypeError(f"expected START:END, got {text!r}")
    start, end = int(match.group(1)), int(match.group(2))
    if start == 0 or end == 0:
        raise argparse.ArgumentTypeError("indices are 1-based; 0 is not allowed")
    return start, end


def _write(selected: Iterable[records.Record]) -> None:
    for record in selected:
        sys.stdout.write(records.format_record(record))


def run_head(args) -> None:
    with contextlib.closing(records.iter_records(args.files)) as source:
        _write(head_records(source, args.number))


def run_range(args) -> None:
    start, end = args.range
    with contextlib.closing(records.iter_records(args.files)) as source:
        _write(select_range(source, start, end))


def register(subparsers) -> None:
    head = subparsers.add_parser("head", help="print the first records")
    _common.add_input_files(head)
    head.add_argument("-n", "--number", type=_common.non_negative_int,
                      default=DEFAULT_HEAD,
                      help="number of records (default: %(default)s)")
    head.set_defaults(func=run_head)

    range_parser = subparsers.add_parser(
        "range", help="print records in a position range")
    _common.add_input_files(range_parser)
    range_parser.add_argument(
        "-r", "--range", type=parse_range, required=True, metavar="START:END",
        help="1-based inclusive range; negative values count from the end")
    range_parser.set_defaults(func=run_range)
