"""sample and shuffle: property tests (no fixed random outputs)."""
import pytest

TOTAL = 40
FA = "".join(f">s{i:02d}\nACGT\n" for i in range(TOTAL))


def names(out):
    return [line[1:] for line in out.splitlines() if line.startswith(">")]


def test_sample_n_returns_exact_count_of_distinct_input_records(run_cli):
    _, out, _ = run_cli("sample", "-n", "7", "-s", "3", stdin=FA.encode())
    got = names(out)
    assert len(got) == 7 and len(set(got)) == 7
    assert set(got) <= {f"s{i:02d}" for i in range(TOTAL)}


def test_sample_n_keeps_input_order(run_cli):
    got = names(run_cli("sample", "-n", "10", "-s", "5", stdin=FA.encode())[1])
    assert got == sorted(got)


def test_sample_is_reproducible_for_same_seed(run_cli):
    args = ("sample", "-n", "10", "-s", "42")
    assert run_cli(*args, stdin=FA.encode()) == run_cli(*args, stdin=FA.encode())


def test_sample_differs_between_seeds(run_cli):
    results = {tuple(names(run_cli("sample", "-n", "10", "-s", str(s),
                                   stdin=FA.encode())[1])) for s in range(6)}
    assert len(results) > 1


def test_sample_n_larger_than_input_returns_everything(run_cli):
    assert len(names(run_cli("sample", "-n", "999", stdin=FA.encode())[1])) == TOTAL


def test_sample_n_zero(run_cli):
    assert run_cli("sample", "-n", "0", stdin=FA.encode())[:2] == (0, "")


def test_sample_proportion_extremes(run_cli):
    assert run_cli("sample", "-p", "0", stdin=FA.encode())[:2] == (0, "")
    assert len(names(run_cli("sample", "-p", "1", stdin=FA.encode())[1])) == TOTAL


def test_sample_proportion_is_roughly_proportional(run_cli):
    big = "".join(f">s{i}\nA\n" for i in range(4000))
    got = len(names(run_cli("sample", "-p", "0.25", "-s", "1", stdin=big.encode())[1]))
    assert 800 < got < 1200


def test_sample_proportion_reproducible(run_cli):
    args = ("sample", "-p", "0.5", "-s", "9")
    assert run_cli(*args, stdin=FA.encode()) == run_cli(*args, stdin=FA.encode())


def test_sample_works_on_fastq_and_keeps_records_intact(run_cli):
    data = "@a\nAC\n+\nII\n@b\nGT\n+\n#I\n"
    out = run_cli("sample", "-n", "1", stdin=data.encode())[1]
    assert out in ("@a\nAC\n+\nII\n", "@b\nGT\n+\n#I\n")


@pytest.mark.parametrize("args", [(), ("-n", "2", "-p", "0.5"), ("-p", "1.5"),
                                  ("-p", "-0.1"), ("-n", "-1")])
def test_sample_argument_errors(run_cli, args):
    assert run_cli("sample", *args, stdin=FA.encode())[0] == 2


def test_sample_empty_input(run_cli):
    assert run_cli("sample", "-n", "3")[:2] == (0, "")


def test_shuffle_is_a_permutation(run_cli):
    _, out, _ = run_cli("shuffle", "-s", "2", stdin=FA.encode())
    assert sorted(names(out)) == [f"s{i:02d}" for i in range(TOTAL)]
    assert names(out) != sorted(names(out))


def test_shuffle_is_reproducible_and_seed_sensitive(run_cli):
    first = run_cli("shuffle", "-s", "7", stdin=FA.encode())
    assert first == run_cli("shuffle", "-s", "7", stdin=FA.encode())
    assert first[1] != run_cli("shuffle", "-s", "8", stdin=FA.encode())[1]


def test_shuffle_keeps_records_intact(run_cli):
    data = "@a\nAC\n+\nII\n@b\nGT\n+\n#I\n"
    out = run_cli("shuffle", stdin=data.encode())[1]
    assert sorted(out.splitlines(keepends=True)) == sorted(data.splitlines(keepends=True))


def test_shuffle_empty_and_single(run_cli):
    assert run_cli("shuffle")[:2] == (0, "")
    assert run_cli("shuffle", stdin=b">a\nA\n")[1] == ">a\nA\n"
