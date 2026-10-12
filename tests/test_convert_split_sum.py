"""fq2fa, split and sum commands."""
import os

import pytest

FQ = b"@r1 x\nACGTAC\n+\nIIIIII\n@r2\nGG\n+\nII\n@r3\nTTT\n+\nIII\n"
FA = b">a\nACGT\n>b\nGGGG\n>c\nTT\n>d\nA\n>e\nC\n"


def test_fq2fa_converts_without_wrapping(run_cli):
    long = "A" * 100
    out = run_cli("fq2fa", stdin=f"@r desc\n{long}\n+\n{'I' * 100}\n".encode())[1]
    assert out == f">r desc\n{long}\n"


def test_fq2fa_writes_to_output_file(run_cli, tmp_path):
    target = tmp_path / "o.fa"
    assert run_cli("fq2fa", "-o", str(target), stdin=FQ)[1] == ""
    assert target.read_text() == ">r1 x\nACGTAC\n>r2\nGG\n>r3\nTTT\n"


def test_fq2fa_output_same_file_as_input_is_an_error(run_cli, make_file):
    path = make_file("in.fq", FQ)
    code, out, err = run_cli("fq2fa", "-o", path, path)
    assert (code, out) == (1, "")
    assert "in.fq" in err and "same file" in err
    with open(path, "rb") as stream:
        assert stream.read() == FQ


def test_fq2fa_output_same_file_through_symlink_or_hard_link(
        run_cli, make_file, tmp_path):
    path = make_file("in.fq", FQ)
    (tmp_path / "soft.fq").symlink_to(path)
    os.link(path, tmp_path / "hard.fq")
    for alias in ("soft.fq", "hard.fq"):
        assert run_cli("fq2fa", "-o", str(tmp_path / alias), path)[0] == 1
    with open(path, "rb") as stream:
        assert stream.read() == FQ


def test_fq2fa_fasta_passes_through_and_empty_input(run_cli):
    assert run_cli("fq2fa", stdin=b">a\nAC\n")[1] == ">a\nAC\n"
    assert run_cli("fq2fa", stdin=b"") == (0, "", "")


def split_parts(directory):
    return {p.name: p.read_text() for p in sorted(directory.iterdir())}


def test_split_by_size(run_cli, make_file, tmp_path):
    path = make_file("in.fa", FA)
    out_dir = tmp_path / "out"
    assert run_cli("split", "-s", "2", "-O", str(out_dir), path)[0] == 0
    assert split_parts(out_dir) == {
        "in.part_001.fa": ">a\nACGT\n>b\nGGGG\n",
        "in.part_002.fa": ">c\nTT\n>d\nA\n",
        "in.part_003.fa": ">e\nC\n",
    }


def test_split_by_parts_gives_earlier_parts_the_remainder(run_cli, make_file,
                                                         tmp_path):
    path = make_file("in.fa", FA)
    out_dir = tmp_path / "out"
    run_cli("split", "-p", "3", "-O", str(out_dir), path)
    parts = split_parts(out_dir)
    assert [v.count(">") for v in parts.values()] == [2, 2, 1]


def test_split_fastq_keeps_extension_and_order(run_cli, make_file, tmp_path):
    path = make_file("reads.fq", FQ)
    out_dir = tmp_path / "out"
    run_cli("split", "-s", "2", "-O", str(out_dir), path)
    parts = split_parts(out_dir)
    assert list(parts) == ["reads.part_001.fq", "reads.part_002.fq"]
    assert "".join(parts.values()) == FQ.decode()


def test_split_fewer_records_than_requested(run_cli, make_file, tmp_path):
    path = make_file("in.fa", b">a\nA\n>b\nC\n")
    out_dir = tmp_path / "p"
    run_cli("split", "-p", "5", "-O", str(out_dir), path)
    assert list(split_parts(out_dir)) == ["in.part_001.fa", "in.part_002.fa"]
    big = tmp_path / "s"
    run_cli("split", "-s", "10", "-O", str(big), path)
    assert list(split_parts(big)) == ["in.part_001.fa"]


def test_split_single_record_and_empty_input(run_cli, make_file, tmp_path):
    out_dir = tmp_path / "one"
    run_cli("split", "-p", "3", "-O", str(out_dir),
            make_file("one.fa", b">a\nA\n"))
    assert split_parts(out_dir) == {"one.part_001.fa": ">a\nA\n"}
    empty_dir = tmp_path / "empty"
    assert run_cli("split", "-s", "2", "-O", str(empty_dir),
                   make_file("e.fa", b""))[0] == 0
    assert split_parts(empty_dir) == {}


def test_split_reuses_existing_output_directory(run_cli, make_file, tmp_path):
    out_dir = tmp_path / "exists"
    out_dir.mkdir()
    code = run_cli("split", "-s", "1", "-O", str(out_dir),
                   make_file("in.fa", b">a\nA\n"))[0]
    assert code == 0 and list(split_parts(out_dir)) == ["in.part_001.fa"]


def test_split_default_directory_is_input_name_plus_split(
        run_cli, make_file, tmp_path, monkeypatch):
    path = make_file("in.fa", b">a\nA\n")
    monkeypatch.chdir(tmp_path)
    run_cli("split", "-s", "1", path)
    assert list(split_parts(tmp_path / "in.fa.split")) == ["in.part_001.fa"]


def test_split_strips_compression_suffix_from_name(run_cli, tmp_path):
    import gzip
    source = tmp_path / "in.fa.gz"
    source.write_bytes(gzip.compress(b">a\nA\n"))
    out_dir = tmp_path / "o"
    run_cli("split", "-s", "1", "-O", str(out_dir), str(source))
    assert list(split_parts(out_dir)) == ["in.part_001.fa"]


def test_split_does_not_follow_a_symlink_at_a_part_path(
        run_cli, make_file, tmp_path):
    victim = make_file("victim.txt", "keep me")
    out_dir = tmp_path / "o"
    out_dir.mkdir()
    (out_dir / "in.part_001.fa").symlink_to(victim)
    code, _, err = run_cli("split", "-s", "1", "-O", str(out_dir),
                           make_file("in.fa", b">a\nA\n"))
    assert code == 1 and "in.part_001.fa" in err
    assert (tmp_path / "victim.txt").read_text() == "keep me"


def test_split_refuses_when_old_parts_exist_without_force(
        run_cli, make_file, tmp_path):
    out_dir = tmp_path / "o"
    out_dir.mkdir()
    (out_dir / "in.part_009.fa").write_text(">old\nA\n")
    (out_dir / "unrelated.txt").write_text("x")
    code, _, err = run_cli("split", "-s", "1", "-O", str(out_dir),
                           make_file("in.fa", b">a\nA\n"))
    assert code == 1 and "--force" in err
    assert sorted(p.name for p in out_dir.iterdir()) == [
        "in.part_009.fa", "unrelated.txt"]


def test_split_force_replaces_old_parts_only(run_cli, make_file, tmp_path):
    out_dir = tmp_path / "o"
    out_dir.mkdir()
    (out_dir / "in.part_009.fa").write_text(">old\nA\n")
    (out_dir / "in.part_001.fa").write_text(">old\nA\n")
    (out_dir / "unrelated.txt").write_text("x")
    code = run_cli("split", "--force", "-s", "1", "-O", str(out_dir),
                   make_file("in.fa", b">a\nA\n"))[0]
    assert code == 0
    assert split_parts(out_dir) == {"in.part_001.fa": ">a\nA\n",
                                    "unrelated.txt": "x"}


def test_split_force_unlinks_a_symlink_without_touching_its_target(
        run_cli, make_file, tmp_path):
    victim = make_file("victim.txt", "keep me")
    out_dir = tmp_path / "o"
    out_dir.mkdir()
    (out_dir / "in.part_001.fa").symlink_to(victim)
    code = run_cli("split", "--force", "-s", "1", "-O", str(out_dir),
                   make_file("in.fa", b">a\nA\n"))[0]
    assert code == 0
    assert (tmp_path / "victim.txt").read_text() == "keep me"
    assert not (out_dir / "in.part_001.fa").is_symlink()
    assert (out_dir / "in.part_001.fa").read_text() == ">a\nA\n"


@pytest.mark.parametrize("args", [[], ["-p", "2", "-s", "2"], ["-p", "0"],
                                  ["-s", "0"]])
def test_split_requires_exactly_one_positive_size(run_cli, make_file, args):
    code = run_cli("split", *args, make_file("in.fa", b">a\nA\n"))[0]
    assert code == 2


def test_split_rejects_stdin(run_cli):
    code, _, err = run_cli("split", "-s", "1", "-", stdin=b">a\nA\n")
    assert code == 1 and "stdin" in err


def sum_line(run_cli, *args, stdin):
    return run_cli("sum", *args, stdin=stdin)[1]


def digest(line):
    return line.split("\t")[0].split("  ")[0]


def test_sum_format(run_cli):
    line = sum_line(run_cli, stdin=b">a\nACGT\n")
    assert line.endswith("  -\n")
    head = line.split("  ")[0]
    prefix, hexdigest = head.rsplit("_", 1)
    assert prefix == "fastxkit.v1_DNAs"
    assert len(hexdigest) == 32 and int(hexdigest, 16) >= 0


def test_sum_is_independent_of_record_order(run_cli):
    first = sum_line(run_cli, stdin=b">a\nACGT\n>b\nGG\n")
    second = sum_line(run_cli, stdin=b">b\nGG\n>a\nACGT\n")
    assert first == second


def test_sum_ignores_letter_case_and_headers(run_cli):
    assert (sum_line(run_cli, stdin=b">a\nACGT\n")
            == sum_line(run_cli, stdin=b">zzz other\nacgt\n"))


def test_sum_fasta_and_fastq_with_same_sequences_agree(run_cli):
    fa = digest(sum_line(run_cli, stdin=b">a\nACGT\n"))
    fq = digest(sum_line(run_cli, stdin=b"@a\nACGT\n+\nIIII\n"))
    assert fa.rsplit("_", 1)[1] == fq.rsplit("_", 1)[1]


def test_sum_different_sequences_differ(run_cli):
    assert (digest(sum_line(run_cli, stdin=b">a\nACGT\n"))
            != digest(sum_line(run_cli, stdin=b">a\nACGA\n")))
    assert (digest(sum_line(run_cli, stdin=b">a\nACGT\n"))
            != digest(sum_line(run_cli, stdin=b">a\nACGT\n>b\nACGT\n")))


def test_sum_ids_matter_only_with_flag(run_cli):
    one = b">a\nACGT\n"
    other = b">b\nACGT\n"
    assert sum_line(run_cli, stdin=one) == sum_line(run_cli, stdin=other)
    assert (sum_line(run_cli, "-i", stdin=one)
            != sum_line(run_cli, "-i", stdin=other))


def test_sum_structure_letter_marks_ids_and_type_follows_content(run_cli):
    assert sum_line(run_cli, "-i", stdin=b">a\nACGT\n").startswith(
        "fastxkit.v1_DNAi_")
    assert sum_line(run_cli, stdin=b">a\nMKLV\n").startswith(
        "fastxkit.v1_Proteins_")


def test_sum_gap_removal(run_cli):
    gapped = b">a\nAC-GT.\n"
    plain = b">a\nACGT\n"
    assert sum_line(run_cli, stdin=gapped) != sum_line(run_cli, stdin=plain)
    assert (digest(sum_line(run_cli, "-g", stdin=gapped)).rsplit("_", 1)[1]
            == digest(sum_line(run_cli, stdin=plain)).rsplit("_", 1)[1])


def test_sum_empty_input_and_file_name(run_cli, make_file):
    assert sum_line(run_cli, stdin=b"").startswith("fastxkit.v1_")
    path = make_file("x.fa", b">a\nA\n")
    assert sum_line(run_cli, path, stdin=b"").endswith(f"  {path}\n")
