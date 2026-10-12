"""sort, rmdup and rename: reorder records and handle repeated names."""
import hashlib
import operator
import sys
from collections.abc import Callable, Iterable, Iterator

from fastxkit import records
from fastxkit.commands import _common

_DIGEST_BYTES = 16


def _text_field(args) -> Callable[[records.Record], str]:
    """Picks the record text that -s, -n or (by default) the ID selects."""
    if args.by_seq:
        return operator.attrgetter("seq")
    if args.by_name:
        return operator.attrgetter("name")
    return operator.attrgetter("id")


def _sort_key(args) -> Callable[[records.Record], str | int | tuple]:
    if args.by_length:
        # Ties go by ID ascending even when reversed, as in the reference tool.
        sign = -1 if args.reverse else 1
        return lambda record: (sign * len(record.seq), record.id)
    field = _text_field(args)
    if args.ignore_case:
        return lambda record: field(record).casefold()
    return field


def sort_records(items: Iterable[records.Record], key: Callable,
                 reverse: bool) -> list[records.Record]:
    """Stably sorts records; ties keep input order even when reversed.

    Raises:
        ValueError: if two records share an ID.
    """
    ordered = list(items)
    seen: set[str] = set()
    for record in ordered:
        if record.id in seen:
            raise ValueError(f"duplicated sequences found: {record.id!r}")
        seen.add(record.id)
    return sorted(ordered, key=key, reverse=reverse)


def run_sort(args) -> None:
    ordered = sort_records(
        records.iter_records(args.files), _sort_key(args),
        args.reverse and not args.by_length)
    for record in ordered:
        sys.stdout.write(records.format_record(record))


def _dedup_key(args) -> Callable[[records.Record], str]:
    field = _text_field(args)
    if args.ignore_case:
        return lambda record: field(record).casefold()
    return field


def drop_duplicates(items: Iterable[records.Record],
                    key: Callable[[records.Record], str],
                    hash_keys: bool) -> Iterator[records.Record]:
    """Yields the first record of each key.

    Memory grows with the number of distinct keys; `hash_keys` stores a
    fixed-size digest instead of the key, which suits long sequences.
    """
    seen: set = set()
    for record in items:
        value = key(record)
        if hash_keys:
            value = hashlib.blake2b(
                value.encode(), digest_size=_DIGEST_BYTES).digest()
        if value not in seen:
            seen.add(value)
            yield record


def run_rmdup(args) -> None:
    for record in drop_duplicates(records.iter_records(args.files),
                                  _dedup_key(args), hash_keys=args.by_seq):
        sys.stdout.write(records.format_record(record))


def rename_duplicates(items: Iterable[records.Record]) -> Iterator[records.Record]:
    """Gives repeated IDs a `_2`, `_3`, ... suffix; the first one is kept.

    Memory grows with the number of distinct IDs.
    """
    used: set[str] = set()
    occurrences: dict[str, int] = {}
    for record in items:
        identifier = record.id
        count = occurrences.get(identifier, 0) + 1
        new_id = identifier if count == 1 else f"{identifier}_{count}"
        while new_id in used:
            count += 1
            new_id = f"{identifier}_{count}"
        occurrences[identifier] = count
        used.add(new_id)
        if new_id == identifier:
            yield record
        else:
            description = record.name.lstrip()[len(identifier):]
            yield records.Record(new_id + description, record.seq, record.qual)


def run_rename(args) -> None:
    for record in rename_duplicates(records.iter_records(args.files)):
        sys.stdout.write(records.format_record(record))


def register(subparsers) -> None:
    sort = subparsers.add_parser(
        "sort", help="sort records (reads everything into memory)")
    _common.add_input_files(sort)
    key = sort.add_mutually_exclusive_group()
    key.add_argument("-n", "--by-name", action="store_true",
                     help="sort by full name; default is the ID")
    key.add_argument("-l", "--by-length", action="store_true",
                     help="sort by sequence length")
    key.add_argument("-s", "--by-seq", action="store_true",
                     help="sort by sequence")
    sort.add_argument("-r", "--reverse", action="store_true")
    sort.add_argument("-i", "--ignore-case", action="store_true",
                      help="ignore case when sorting by text")
    sort.set_defaults(func=run_sort)

    rmdup = subparsers.add_parser(
        "rmdup", help="remove duplicate records, keeping the first")
    _common.add_input_files(rmdup)
    kind = rmdup.add_mutually_exclusive_group()
    kind.add_argument("-n", "--by-name", action="store_true",
                      help="compare full names; default is the ID")
    kind.add_argument("-s", "--by-seq", action="store_true",
                      help="compare sequences")
    rmdup.add_argument("-i", "--ignore-case", action="store_true")
    rmdup.set_defaults(func=run_rmdup)

    rename = subparsers.add_parser(
        "rename", help="number records that repeat an ID")
    _common.add_input_files(rename)
    rename.set_defaults(func=run_rename)
