"""crosscheck: compares olgu.json read totals with counts from the FASTQ files.

Only paths, numbers and exception class names are printed, never sequence
content or parser messages (which may quote file text).
"""
import csv
import dataclasses
import json
import os
import sys

from fastxkit import io as fxio
from fastxkit import records
from fastxkit.commands import _common, stats

OK = "OK"
MISMATCH = "MISMATCH"
MISSING_IN_MAP = "MISSING_IN_MAP"
MISSING_IN_OLGU = "MISSING_IN_OLGU"
FILE_ERROR = "FILE_ERROR"
_COLUMNS = ("sample", "olgu", "fastq", "diff", "status")
_TEXT_COLUMNS = frozenset(("sample", "status"))
EXIT_NOT_OK = 1
EXIT_INPUT_ERROR = 3
_ABSENT = "-"
_REQUIRED_MAP_HEADER = ("sample", "fastq")
_READ_COLUMN = "read"
_READ_LABELS = ("r1", "r2")
_FACT_KEY = "ornek_tablosu"


@dataclasses.dataclass
class _Row:
    sample: str
    olgu: int | None = None
    fastq: int | None = None
    status: str = OK

    @property
    def diff(self) -> int | None:
        if self.olgu is None or self.fastq is None:
            return None
        return self.fastq - self.olgu


def _load_olgu_counts(path: str) -> dict[str, int]:
    """Returns {sample: toplam_okuma} from the ornek_tablosu fact."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except json.JSONDecodeError as error:
        raise ValueError(f"olgu file {path}: not valid JSON ({error})") from None
    table = document.get(_FACT_KEY) if isinstance(document, dict) else None
    entries = table.get("deger") if isinstance(table, dict) else None
    if not isinstance(entries, list):
        raise ValueError(
            f"olgu file {path}: no {_FACT_KEY}.deger list found")
    counts: dict[str, int] = {}
    for entry in entries:
        sample, total = _parse_entry(path, entry)
        if sample in counts:
            raise ValueError(
                f"olgu file {path}: duplicate sample {sample!r} in {_FACT_KEY}")
        counts[sample] = total
    return counts


def _parse_entry(path: str, entry) -> tuple[str, int]:
    if not isinstance(entry, dict) or not isinstance(entry.get("ornek"), str):
        raise ValueError(f"olgu file {path}: {_FACT_KEY} entry without 'ornek'")
    sample, total = entry["ornek"].strip(), entry.get("toplam_okuma")
    # MultiQC may write whole numbers as floats (20000000.0).
    if isinstance(total, float) and total.is_integer():
        total = int(total)
    if isinstance(total, bool) or not isinstance(total, int):
        raise ValueError(
            f"olgu file {path}: toplam_okuma of {sample!r} is not a whole "
            f"number: {total!r}")
    if total < 0:
        raise ValueError(
            f"olgu file {path}: toplam_okuma of {sample!r} is negative: {total}")
    return sample, total


def _load_map(path: str, count_mode: str) -> dict[str, list[str]]:
    """Returns {sample: [fastq paths]} honouring --count."""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        lines = [[cell.strip() for cell in row]
                 for row in csv.reader(handle) if any(c.strip() for c in row)]
    if not lines:
        raise ValueError(f"map file {path}: empty")
    header = [cell.lower() for cell in lines[0]]
    if tuple(header[:2]) != _REQUIRED_MAP_HEADER:
        raise ValueError(
            f"map file {path}: first line must be the header "
            "sample,fastq[,read]")
    if len(header) > 2 and header[2] != _READ_COLUMN:
        raise ValueError(
            f"map file {path}: third column must be named 'read', "
            f"got {header[2]!r}")
    has_read = len(header) > 2
    if count_mode == "r1" and not has_read:
        raise ValueError(f"map file {path}: --count r1 needs a 'read' column")
    base_dir = os.path.dirname(os.path.abspath(path))
    mapping: dict[str, list[str]] = {}
    seen: set[tuple[str, str]] = set()
    for number, cells in enumerate(lines[1:], start=2):
        sample, fastq, read = _map_cells(path, number, cells, has_read)
        fastq_path = os.path.normpath(os.path.join(base_dir, fastq))
        if (sample, fastq_path) in seen:
            raise ValueError(
                f"map file {path}: row {number}: {fastq!r} is listed twice "
                f"for sample {sample!r}")
        seen.add((sample, fastq_path))
        mapping.setdefault(sample, [])
        if count_mode == "all" or read == "r1":
            mapping[sample].append(fastq_path)
    for sample, files in mapping.items():
        if not files:
            raise ValueError(f"map file {path}: sample {sample!r} has no r1 row")
    return mapping


def _map_cells(path: str, number: int, cells: list[str],
               has_read: bool) -> tuple[str, str, str]:
    if len(cells) < 2 or not cells[0] or not cells[1]:
        raise ValueError(
            f"map file {path}: row {number} needs sample and fastq values")
    read = cells[2].lower() if has_read and len(cells) > 2 else ""
    if has_read and read not in _READ_LABELS:
        raise ValueError(
            f"map file {path}: row {number}: read must be r1 or r2")
    return cells[0], cells[1], read


def _count_reads(files: list[str]) -> int:
    return sum(1 for _ in records.iter_records(files))


def build_rows(olgu: dict[str, int], mapping: dict[str, list[str]],
               tolerance: int) -> list[_Row]:
    """Compares every sample; unreadable files are reported on stderr."""
    rows = []
    for sample in [*olgu, *(s for s in mapping if s not in olgu)]:
        row = _Row(sample, olgu.get(sample))
        rows.append(row)
        if sample not in mapping:
            row.status = MISSING_IN_MAP
            continue
        try:
            row.fastq = _count_reads(mapping[sample])
        except fxio.DATA_ERRORS as error:
            row.status = FILE_ERROR
            print(f"fastxkit: crosscheck: {sample}: {type(error).__name__}",
                  file=sys.stderr)
            continue
        if row.olgu is None:
            row.status = MISSING_IN_OLGU
        elif abs(row.diff) > tolerance:
            row.status = MISMATCH
    return rows


def _cell(value: int | None) -> str:
    return _ABSENT if value is None else str(value)


def _render(rows: list[_Row]) -> str:
    cells = [{"sample": row.sample, "olgu": _cell(row.olgu),
              "fastq": _cell(row.fastq), "diff": _cell(row.diff),
              "status": row.status} for row in rows]
    return stats.render_table(_COLUMNS, cells, _TEXT_COLUMNS)


def _exit_code(rows: list[_Row]) -> int:
    """3 if any file could not be read, else 1 if any row differs, else 0."""
    statuses = {row.status for row in rows}
    if FILE_ERROR in statuses:
        return EXIT_INPUT_ERROR
    return 0 if statuses <= {OK} else EXIT_NOT_OK


def run(args) -> int:
    """Prints the comparison; returns 0 (all OK), 1 (differences) or 3 (bad input)."""
    try:
        olgu = _load_olgu_counts(args.olgu_file)
        mapping = _load_map(args.map, args.count)
        if not olgu and not mapping:
            raise ValueError(
                "nothing to compare: the olgu file lists no samples and the "
                "map file has no rows")
    except (ValueError, OSError) as error:
        print(f"fastxkit: {error}", file=sys.stderr)
        return EXIT_INPUT_ERROR
    rows = build_rows(olgu, mapping, args.tolerance)
    if args.json:
        output = json.dumps(
            [{"sample": row.sample, "olgu": row.olgu, "fastq": row.fastq,
              "diff": row.diff, "status": row.status} for row in rows],
            indent=2) + "\n"
    else:
        output = _render(rows)
    sys.stdout.write(output)
    return _exit_code(rows)


def register(subparsers) -> None:
    parser = subparsers.add_parser(
        "crosscheck",
        help="compare olgu.json read totals with FASTQ read counts")
    parser.add_argument("olgu_file", metavar="OLGU.json")
    parser.add_argument("--map", required=True, metavar="MAP.csv",
                        help="CSV with columns sample,fastq[,read]; relative "
                        "FASTQ paths are relative to the CSV")
    parser.add_argument("--tolerance", type=_common.non_negative_int, default=0,
                        help="allowed absolute difference (default: 0)")
    parser.add_argument("--count", choices=("all", "r1"), default="all",
                        help="count reads of all files of a sample, or only "
                        "its r1 files (default: all)")
    parser.add_argument("--json", action="store_true", help="JSON output")
    parser.set_defaults(func=run)
