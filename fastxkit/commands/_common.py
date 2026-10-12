"""Argument helpers shared by the commands."""
import argparse
import os
import re
from collections.abc import Iterable
from typing import TextIO

from fastxkit import io as fxio


def add_input_files(parser: argparse.ArgumentParser) -> None:
    """Adds the positional input-file list; '-' or nothing means stdin."""
    parser.add_argument(
        "files", nargs="*", default=[fxio.STDIN], metavar="FILE",
        help="input files, possibly gz/bz2/xz compressed (default: stdin)")


def non_negative_int(text: str) -> int:
    """argparse type for integers >= 0."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"must be >= 0: {value}")
    return value


def positive_int(text: str) -> int:
    """argparse type for integers >= 1."""
    value = non_negative_int(text)
    if value == 0:
        raise argparse.ArgumentTypeError("must be >= 1: 0")
    return value


def compile_regex(pattern: str, ignore_case: bool = False) -> re.Pattern:
    """Compiles a user-supplied regular expression, or raises ValueError."""
    try:
        return re.compile(pattern, re.IGNORECASE if ignore_case else 0)
    except re.error as error:
        raise ValueError(
            f"invalid regular expression {pattern!r}: {error}") from None


def refuse_same_file(output_path: str, input_paths: Iterable[str]) -> None:
    """Raises ValueError if `output_path` is one of the input files.

    Opening the output for writing truncates it, so this must run before the
    output is opened. Symbolic and hard links to an input count as the same
    file; stdin and an output that does not exist yet cannot clash.
    """
    if not os.path.exists(output_path):
        return
    for input_path in input_paths:
        if input_path != fxio.STDIN and os.path.samefile(
                output_path, input_path):
            raise ValueError(
                f"output file {output_path} is the same file as input "
                f"{input_path}; refusing to overwrite the input")


def open_new_text_file(path: str, replace_regular: bool = False) -> TextIO:
    """Opens `path` for writing without ever following a symbolic link.

    Args:
        path: file to write.
        replace_regular: when False the file must not exist yet; when True an
            existing regular file is truncated. A symbolic link at `path` is
            an error either way (its target is never touched).

    Raises:
        FileExistsError: `path` exists and `replace_regular` is False.
        OSError: `path` is a symbolic link or cannot be opened safely.
    """
    if os.path.islink(path):
        raise OSError(f"{path} is a symbolic link; refusing to write through it")
    flags = os.O_WRONLY | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
    flags |= os.O_TRUNC if replace_regular else os.O_EXCL
    descriptor = os.open(path, flags, 0o666)
    return os.fdopen(descriptor, "w", encoding="utf-8")
