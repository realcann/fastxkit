"""sort, rmdup and rename commands."""
import pytest

from fastxkit import records
from fastxkit.commands import organize


def ids(out):
    return [line[1:].split()[0] for line in out.splitlines() if line[0] == ">"]


def fasta(*items):
    return "".join(f">{name}\n{seq}\n" for name, seq in items).encode()


def test_sort_default_is_by_id(run_cli):
    data = fasta(("b 1", "A"), ("a 2", "A"), ("c", "A"))
    assert ids(run_cli("sort", stdin=data)[1]) == ["a", "b", "c"]


def test_sort_by_name(run_cli):
    data = fasta(("b a", "A"), ("a z", "A"))
    assert ids(run_cli("sort", "-n", stdin=data)[1]) == ["a", "b"]


def test_sort_by_length_is_stable_and_reverse_keeps_ties_in_order(run_cli):
    data = fasta(("x", "AA"), ("y", "A"), ("z", "AA"))
    assert ids(run_cli("sort", "-l", stdin=data)[1]) == ["y", "x", "z"]
    assert ids(run_cli("sort", "-l", "-r", stdin=data)[1]) == ["x", "z", "y"]


def test_sort_by_length_breaks_ties_by_id_in_both_directions(run_cli):
    data = fasta(("z", "AA"), ("y", "A"), ("b", "AA"), ("a", "A"))
    assert ids(run_cli("sort", "-l", stdin=data)[1]) == ["a", "y", "b", "z"]
    assert ids(run_cli("sort", "-l", "-r", stdin=data)[1]) == ["b", "z", "a", "y"]


def test_sort_by_sequence_ignore_case(run_cli):
    data = fasta(("x", "b"), ("y", "A"), ("z", "C"))
    assert ids(run_cli("sort", "-s", stdin=data)[1]) == ["y", "z", "x"]
    assert ids(run_cli("sort", "-s", "-i", stdin=data)[1]) == ["y", "x", "z"]
    assert ids(run_cli("sort", "-n", "-i", stdin=fasta(("b", "A"), ("A", "A")))[1]
               ) == ["A", "b"]


def test_sort_duplicate_id_is_an_error_without_output(run_cli):
    code, out, err = run_cli("sort", stdin=fasta(("a", "A"), ("a", "C")))
    assert (code, out) == (1, "")
    assert "duplicated sequences found" in err


def test_sort_empty_input(run_cli):
    assert run_cli("sort")[:2] == (0, "")


def test_sort_rejects_two_keys(run_cli):
    assert run_cli("sort", "-n", "-l", stdin=fasta(("a", "A")))[0] == 2


def test_sort_fastq_and_crlf(run_cli):
    data = b"@b\r\nAC\r\n+\r\nII\r\n@a\r\nA\r\n+\r\nI\r\n"
    assert run_cli("sort", stdin=data)[1] == "@a\nA\n+\nI\n@b\nAC\n+\nII\n"


DUPS = fasta(("a one", "ACGT"), ("b two", "acgt"), ("a three", "TTTT"),
             ("c one", "ACGT"))


def test_rmdup_by_id_keeps_first(run_cli):
    assert ids(run_cli("rmdup", stdin=DUPS)[1]) == ["a", "b", "c"]


def test_rmdup_by_name(run_cli):
    data = fasta(("a x", "A"), ("a x", "C"), ("a y", "G"))
    assert ids(run_cli("rmdup", "-n", stdin=data)[1]) == ["a", "a"]


def test_rmdup_by_sequence_is_case_sensitive_unless_asked(run_cli):
    assert ids(run_cli("rmdup", "-s", stdin=DUPS)[1]) == ["a", "b", "a"]
    assert ids(run_cli("rmdup", "-s", "-i", stdin=DUPS)[1]) == ["a", "a"]


def test_rmdup_ignore_case_applies_to_ids(run_cli):
    data = fasta(("A", "A"), ("a", "C"))
    assert ids(run_cli("rmdup", stdin=data)[1]) == ["A", "a"]
    assert ids(run_cli("rmdup", "-i", stdin=data)[1]) == ["A"]


def test_rmdup_empty_input_and_fastq(run_cli):
    assert run_cli("rmdup")[:2] == (0, "")
    data = b"@a\nA\n+\nI\n@a\nC\n+\nI\n"
    assert run_cli("rmdup", stdin=data)[1] == "@a\nA\n+\nI\n"


def test_rmdup_multiline_sequences_compare_whole(run_cli):
    data = b">a\nAC\nGT\n>b\nACGT\n"
    assert ids(run_cli("rmdup", "-s", stdin=data)[1]) == ["a"]


def test_rename_numbers_repeated_ids_and_keeps_description(run_cli):
    data = fasta(("a first", "A"), ("a second", "C"), ("a third", "G"))
    assert run_cli("rename", stdin=data)[1] == (
        ">a first\nA\n>a_2 second\nC\n>a_3 third\nG\n")


def test_rename_leaves_unique_ids_alone(run_cli):
    assert run_cli("rename", stdin=fasta(("a", "A"), ("b", "C")))[1] == (
        ">a\nA\n>b\nC\n")


def test_rename_avoids_colliding_with_an_existing_id(run_cli):
    data = fasta(("a", "A"), ("a_2", "C"), ("a", "G"))
    assert ids(run_cli("rename", stdin=data)[1]) == ["a", "a_2", "a_3"]


def test_rename_fastq_and_empty(run_cli):
    data = b"@a\nA\n+\nI\n@a\nC\n+\nI\n"
    assert run_cli("rename", stdin=data)[1] == "@a\nA\n+\nI\n@a_2\nC\n+\nI\n"
    assert run_cli("rename")[:2] == (0, "")


def test_rename_duplicates_is_lazy():
    source = iter([records.Record("a", "A"), records.Record("a", "C")])
    renamed = organize.rename_duplicates(source)
    assert next(renamed).name == "a"
    assert next(renamed).name == "a_2"


@pytest.mark.parametrize("command", ["sort", "rmdup", "rename"])
def test_missing_file_is_a_data_error(run_cli, command):
    assert run_cli(command, "/nonexistent/x.fa")[0] == 1
