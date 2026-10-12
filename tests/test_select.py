"""head and range commands, plus the generic selection helpers."""
import itertools

import pytest

from fastxkit.commands import select

FA = "".join(f">s{i}\nA\n" for i in range(1, 6))


def names(out):
    return [line[1:] for line in out.splitlines() if line.startswith(">")]


def test_head_default_is_ten(run_cli):
    data = "".join(f">s{i}\nA\n" for i in range(12))
    assert len(names(run_cli("head", stdin=data.encode())[1])) == 10


def test_head_more_than_available(run_cli):
    assert names(run_cli("head", "-n", "9", stdin=FA.encode())[1]) == [
        "s1", "s2", "s3", "s4", "s5"]


def test_head_zero_and_empty_input(run_cli):
    assert run_cli("head", "-n", "0", stdin=FA.encode())[:2] == (0, "")
    assert run_cli("head", "-n", "3")[:2] == (0, "")


def test_head_rejects_negative(run_cli):
    assert run_cli("head", "-n", "-1", stdin=FA.encode())[0] == 2


def test_head_does_not_read_whole_input():
    stream = (n for n in itertools.count())
    assert list(select.head_records(stream, 3)) == [0, 1, 2]


def test_head_stops_before_a_later_corrupt_record(run_cli):
    data = b"@a\nA\n+\nI\n@b\nAC\n+\nI\n"
    code, out, _ = run_cli("head", "-n", "1", stdin=data)
    assert (code, out) == (0, "@a\nA\n+\nI\n")


def test_range_basic_and_open_ends_clipped(run_cli):
    assert names(run_cli("range", "-r", "2:3", stdin=FA.encode())[1]) == ["s2", "s3"]
    assert names(run_cli("range", "-r", "4:99", stdin=FA.encode())[1]) == ["s4", "s5"]


def test_range_negative_indices_after_separate_arg(run_cli):
    assert names(run_cli("range", "-r", "-2:-1", stdin=FA.encode())[1]) == ["s4", "s5"]
    assert names(run_cli("range", "-r", "2:-2", stdin=FA.encode())[1]) == ["s2", "s3", "s4"]
    assert names(run_cli("range", "-r", "-3:4", stdin=FA.encode())[1]) == ["s3", "s4"]
    assert names(run_cli("range", "--range", "-5:-4", stdin=FA.encode())[1]) == ["s1", "s2"]


def test_range_empty_when_start_after_end(run_cli):
    assert run_cli("range", "-r", "4:2", stdin=FA.encode())[:2] == (0, "")


@pytest.mark.parametrize("bad", ["0:3", "3", "a:b", "1:0", "1-2", ""])
def test_range_rejects_malformed(run_cli, bad):
    assert run_cli("range", "-r", bad, stdin=FA.encode())[0] == 2


def test_range_requires_the_flag(run_cli):
    assert run_cli("range", stdin=FA.encode())[0] == 2


@pytest.mark.parametrize("total", range(0, 8))
def test_select_range_matches_python_slicing(total):
    items = list(range(1, total + 1))
    for start in [s for s in range(-9, 10) if s != 0]:
        for end in [e for e in range(-9, 10) if e != 0]:
            first = start - 1 if start > 0 else max(total + start, 0)
            last = end if end > 0 else total + end + 1
            expected = items[first:max(last, 0)] if last > 0 else []
            assert list(select.select_range(iter(items), start, end)) == expected, (
                total, start, end)


def test_select_range_positive_stops_early():
    stream = (n for n in itertools.count(1))
    assert list(select.select_range(stream, 2, 4)) == [2, 3, 4]


def test_range_on_fastq(run_cli):
    data = "@a\nA\n+\nI\n@b\nC\n+\nI\n@c\nG\n+\nI\n"
    assert run_cli("range", "-r", "-1:-1", stdin=data.encode())[1] == "@c\nG\n+\nI\n"
