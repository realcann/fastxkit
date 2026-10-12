"""replace: regular-expression substitution in names or sequences."""
import dataclasses
import re
import sys
from collections.abc import Callable

from fastxkit import records
from fastxkit.commands import _common

# `$$` is a literal dollar sign, `$1` and `${1}` refer to capture groups.
_REFERENCE = re.compile(r"\$(?:(\$)|(\d+)|\{(\d+)\})")


def parse_replacement(template: str, groups: int) -> Callable[[re.Match], str]:
    """Turns a `$N` template into a function of the match.

    Raises:
        ValueError: if the template refers to a group the pattern lacks.
    """
    parts: list[str | int] = []
    position = 0
    for reference in _REFERENCE.finditer(template):
        parts.append(template[position:reference.start()])
        if reference.group(1):
            parts.append("$")
        else:
            number = int(reference.group(2) or reference.group(3))
            if number > groups:
                raise ValueError(
                    f"replacement refers to group {number}, but the pattern "
                    f"has {groups}")
            parts.append(number)
        position = reference.end()
    parts.append(template[position:])

    def expand(match: re.Match) -> str:
        return "".join(
            part if isinstance(part, str) else match.group(part) or ""
            for part in parts)

    return expand


def replace_in_record(record: records.Record, pattern: re.Pattern,
                      expand: Callable[[re.Match], str],
                      in_sequence: bool) -> records.Record:
    """Returns the record with every match of `pattern` replaced."""
    if not in_sequence:
        return dataclasses.replace(record, name=pattern.sub(expand, record.name))
    new_seq = pattern.sub(expand, record.seq)
    if record.qual is not None and len(new_seq) != len(record.qual):
        raise ValueError(
            f"replacement changes the sequence length of {record.id!r} "
            "but its quality string would not match")
    return dataclasses.replace(record, seq=new_seq)


def run(args) -> None:
    pattern = _common.compile_regex(args.pattern, args.ignore_case)
    expand = parse_replacement(args.replacement, pattern.groups)
    for record in records.iter_records(args.files):
        sys.stdout.write(records.format_record(
            replace_in_record(record, pattern, expand, args.by_seq)))


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "replace", help="replace text in names or sequences")
    _common.add_input_files(parser)
    parser.add_argument("-p", "--pattern", required=True,
                        help="regular expression to look for")
    parser.add_argument("-r", "--replacement", required=True,
                        help="replacement; $1 or ${1} insert a captured group")
    parser.add_argument("-s", "--by-seq", action="store_true",
                        help="edit the sequence instead of the full name")
    parser.add_argument("-i", "--ignore-case", action="store_true")
    parser.set_defaults(func=run)
