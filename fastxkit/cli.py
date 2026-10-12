"""Command-line entry point: argument parsing, dispatch, error reporting."""
import argparse
import os
import re
import sys
from collections.abc import Sequence

from fastxkit import io as fxio
from fastxkit.commands import (checksum, convert, crosscheck, faidx, organize,
                               qc, replace, sampling, search, select, seq,
                               split, stats, subseq, tabular)

_COMMAND_MODULES = (stats, seq, select, sampling, tabular, search, subseq,
                    faidx, organize, replace, convert, split, checksum, qc,
                    crosscheck)
_RANGE_FLAGS = ("-r", "--range")
_NEGATIVE_RANGE = re.compile(r"-\d+:-?\d+")


def build_parser() -> argparse.ArgumentParser:
    """Builds the top-level parser with one sub-parser per command."""
    parser = argparse.ArgumentParser(
        prog="fastxkit", description="FASTA/FASTQ toolkit (standard library only).")
    subparsers = parser.add_subparsers(
        dest="command", required=True, metavar="COMMAND")
    for module in _COMMAND_MODULES:
        module.register(subparsers)
    return parser


def _glue_range_values(argv: Sequence[str]) -> list[str]:
    """Rewrites `-r -2:-1` as `-r-2:-1`, which argparse takes as a value.

    argparse treats a token such as '-2:-1' as an unknown option, so a
    negative start index would otherwise be impossible to pass separately.
    """
    glued: list[str] = []
    for token in argv:
        if (glued and glued[-1] in _RANGE_FLAGS
                and _NEGATIVE_RANGE.fullmatch(token)):
            flag = glued.pop()
            glued.append(f"{flag}{token}" if flag == "-r" else f"{flag}={token}")
        else:
            glued.append(token)
    return glued


def _silence_broken_pipe() -> None:
    # Point stdout at devnull so interpreter shutdown does not raise again.
    os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())


def main(argv: Sequence[str] | None = None) -> int:
    """Runs the CLI and returns the process exit code (0, 1 data error, 2 usage).

    A command's `run` may return an exit code (e.g. crosscheck returns 1 on
    mismatch); None means success.
    """
    arguments = _glue_range_values(sys.argv[1:] if argv is None else argv)
    try:
        args = build_parser().parse_args(arguments)
    except SystemExit as exit_request:
        return exit_request.code or 0
    exit_code = 0
    try:
        exit_code = args.func(args)
        sys.stdout.flush()
    except BrokenPipeError:
        _silence_broken_pipe()
    except fxio.DATA_ERRORS as error:
        sys.stdout.flush()
        print(f"fastxkit: {error}", file=sys.stderr)
        return 1
    return exit_code or 0
