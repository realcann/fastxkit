"""split: write records to numbered part files."""
import fnmatch
import glob
import math
import os
from collections.abc import Iterator

from fastxkit import records
from fastxkit.commands import _common

_COMPRESSION_SUFFIXES = (".gz", ".bz2", ".xz")
_PART_DIGITS = 3


def _stem_and_extension(path: str) -> tuple[str, str]:
    name = os.path.basename(path)
    for suffix in _COMPRESSION_SUFFIXES:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return os.path.splitext(name)


def part_file_name(path: str, number: int) -> str:
    """Names part `number` as `<stem>.part_001<ext>`, dropping .gz/.bz2/.xz."""
    stem, extension = _stem_and_extension(path)
    return f"{stem}.part_{number:0{_PART_DIGITS}d}{extension}"


def _existing_parts(path: str, out_dir: str) -> list[str]:
    """Paths in `out_dir` named like parts of `path` (links included)."""
    stem, extension = _stem_and_extension(path)
    pattern = f"{glob.escape(stem)}.part_*{glob.escape(extension)}"
    return [os.path.join(out_dir, name) for name in sorted(os.listdir(out_dir))
            if fnmatch.fnmatchcase(name, pattern)]


def _prepare_out_dir(path: str, out_dir: str, force: bool) -> None:
    """Makes sure no old part can survive or be mixed into the new output.

    A part name is always longer than the input's own name, so --force can
    never match (and delete) the input file itself.

    Raises:
        ValueError: old parts exist and `force` is not set.
    """
    os.makedirs(out_dir, exist_ok=True)
    old_parts = _existing_parts(path, out_dir)
    if not old_parts:
        return
    if not force:
        raise ValueError(
            f"{out_dir} already holds files named like parts of {path} "
            f"(e.g. {os.path.basename(old_parts[0])}); remove them or use "
            "--force to replace them")
    for old in old_parts:
        os.unlink(old)  # unlink never follows a symbolic link


def _part_size(path: str, args) -> int:
    """Records per part; `-p` rounds up so that the last part is the short one."""
    if args.size is not None:
        return args.size
    count = sum(1 for _ in records.iter_records([path]))
    return max(1, math.ceil(count / args.parts))


def _chunks(path: str, size: int) -> Iterator[list[records.Record]]:
    chunk: list[records.Record] = []
    for record in records.iter_records([path]):
        chunk.append(record)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def run(args) -> None:
    path = args.file
    if path == "-":
        raise ValueError("split cannot read from stdin: it needs a file name")
    out_dir = args.out_dir or f"{path}.split"
    _prepare_out_dir(path, out_dir, args.force)
    size = _part_size(path, args)
    for number, chunk in enumerate(_chunks(path, size), start=1):
        target = os.path.join(out_dir, part_file_name(path, number))
        with _common.open_new_text_file(target) as handle:
            handle.writelines(records.format_record(r) for r in chunk)


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "split", help="split a file into parts by part count or record count")
    parser.add_argument("file", metavar="FILE",
                        help="input file (stdin is not supported)")
    amount = parser.add_mutually_exclusive_group(required=True)
    amount.add_argument("-p", "--parts", type=_common.positive_int, metavar="N",
                        help="split into N parts")
    amount.add_argument("-s", "--size", type=_common.positive_int, metavar="N",
                        help="N records per part")
    parser.add_argument("-O", "--out-dir", metavar="DIR",
                        help="output directory, created if missing "
                             "(default: <FILE>.split)")
    parser.add_argument("--force", action="store_true",
                        help="delete existing <name>.part_*<ext> files in the "
                             "output directory first (links are unlinked, "
                             "never followed); without it such files are an "
                             "error")
    parser.set_defaults(func=run)

