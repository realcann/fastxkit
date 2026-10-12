"""grep and locate: find records or match positions by pattern."""
import re
import sys
from collections.abc import Callable, Iterator

from fastxkit import metrics, records
from fastxkit.commands import _common, seq

_TAB = "\t"
LOCATE_COLUMNS = ("seqID", "patternName", "pattern", "strand", "start", "end",
                  "matched")


def build_matcher(patterns: list[str], use_regex: bool, exact: bool,
                  ignore_case: bool) -> Callable[[str], bool]:
    """Returns a test that is true when any pattern matches the text.

    Without `use_regex` the patterns are literal; `exact` then requires the
    whole text to equal the pattern instead of merely containing it.
    """
    compiled = [
        _common.compile_regex(p if use_regex else re.escape(p), ignore_case)
        for p in patterns]
    if exact and not use_regex:
        return lambda text: any(c.fullmatch(text) for c in compiled)
    return lambda text: any(c.search(text) for c in compiled)


def _grep_target(args) -> Callable[[records.Record], str]:
    if args.name:
        return lambda record: record.name
    if args.seq:
        return lambda record: record.seq
    return lambda record: record.id


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(seq.COMPLEMENT)[::-1]


def _is_nucleotide(sequence: str) -> bool:
    return metrics.detect_type(set(sequence)) in ("DNA", "RNA")


def run_grep(args) -> None:
    # ID and name are compared whole; a sequence is searched for the pattern.
    matches = build_matcher(args.pattern, args.regex, not args.seq,
                            args.ignore_case)
    target = _grep_target(args)
    both_strands = args.seq and not args.only_positive_strand
    for record in records.iter_records(args.files):
        found = matches(target(record))
        # Protein has no reverse strand, so only DNA/RNA is searched twice.
        if (not found and both_strands and _is_nucleotide(record.seq)):
            found = matches(_reverse_complement(record.seq))
        if found != args.invert:
            sys.stdout.write(records.format_record(record))


def find_hits(sequence: str, finder: re.Pattern) -> Iterator[tuple[int, int]]:
    """Yields 1-based inclusive (start, end) of non-empty, overlapping hits."""
    for match in finder.finditer(sequence):
        if match.end(1) > match.start(1):
            yield match.start(1) + 1, match.end(1)


def locate_rows(record: records.Record, pattern: str, finder: re.Pattern,
                both_strands: bool,
                report_pattern: bool = False) -> Iterator[str]:
    """Yields tab-separated hits; '-' hits use forward-strand coordinates.

    With `report_pattern` the last column shows the pattern as typed instead
    of the sequence letters (the reference tool does this for literal -i).
    """
    length = len(record.seq)
    strands = [("+", record.seq)]
    if both_strands:
        strands.append(("-", _reverse_complement(record.seq)))
    for strand, sequence in strands:
        for start, end in find_hits(sequence, finder):
            if strand == "-":
                first, last = length - end + 1, length - start + 1
            else:
                first, last = start, end
            yield _TAB.join((record.id, pattern, pattern, strand, str(first),
                             str(last),
                             pattern if report_pattern
                             else sequence[start - 1:end])) + "\n"


def run_locate(args) -> None:
    sys.stdout.write(_TAB.join(LOCATE_COLUMNS) + "\n")
    # A lookahead lets finditer report hits that overlap each other.
    finders = [
        (pattern, _common.compile_regex(
            f"(?=({pattern if args.regex else re.escape(pattern)}))",
            args.ignore_case))
        for pattern in args.pattern]
    for record in records.iter_records(args.files):
        for pattern, finder in finders:
            for row in locate_rows(record, pattern, finder,
                                   not args.only_positive_strand,
                                   args.ignore_case and not args.regex):
                sys.stdout.write(row)


def register(subparsers) -> None:
    grep = subparsers.add_parser("grep", help="select records by pattern")
    _common.add_input_files(grep)
    grep.add_argument("-p", "--pattern", action="append", required=True,
                      metavar="PATTERN",
                      help="pattern to look for (repeatable, any may match)")
    target = grep.add_mutually_exclusive_group()
    target.add_argument("-n", "--by-name", dest="name", action="store_true",
                        help="match the full name; default is the ID (first word)")
    target.add_argument("-s", "--by-seq", dest="seq", action="store_true",
                        help="search the sequence for the pattern")
    grep.add_argument("-r", "--use-regexp", dest="regex", action="store_true",
                      help="treat patterns as regular expressions (partial match)")
    grep.add_argument("-i", "--ignore-case", action="store_true")
    grep.add_argument("-P", "--only-positive-strand", action="store_true",
                      help="with -s, search only the given strand "
                           "(default: DNA/RNA also searched on the reverse "
                           "complement)")
    grep.add_argument("-v", "--invert-match", dest="invert", action="store_true",
                      help="print records that do not match")
    grep.set_defaults(func=run_grep)

    locate = subparsers.add_parser(
        "locate", help="report pattern positions in the sequences")
    _common.add_input_files(locate)
    locate.add_argument("-p", "--pattern", action="append", required=True,
                        metavar="PATTERN",
                        help="pattern to look for (repeatable)")
    locate.add_argument("-r", "--use-regexp", dest="regex", action="store_true",
                        help="treat patterns as regular expressions")
    locate.add_argument("-i", "--ignore-case", action="store_true")
    locate.add_argument("-P", "--only-positive-strand", action="store_true",
                        help="skip the reverse-complement strand")
    locate.set_defaults(func=run_locate)
