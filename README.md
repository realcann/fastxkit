# fastxkit

A FASTA/FASTQ toolkit written with the Python standard library only. It covers
the common seqkit-style operations and adds two commands for independently
verifiable QC numbers (`qc`, `crosscheck`).

## Install

Nothing to install: no third-party dependencies. Use Python 3.10+ and run
from the repository root:

    python -m fastxkit COMMAND ...

Input may be plain, gzip, bzip2 or xz compressed (detected by content).

## Commands

- `stats` - summarise sequence files (`-a` adds quartiles, N50, Q20/Q30, GC).
- `seq` - transform and filter sequences.
- `head` - print the first records.
- `range` - print records in a position range.
- `sample` - sample records at random.
- `shuffle` - shuffle record order.
- `fx2tab` - convert FASTA/FASTQ to a table.
- `tab2fx` - convert a table to FASTA/FASTQ.
- `grep` - select records by pattern.
- `locate` - report pattern positions in the sequences.
- `subseq` - extract a region from records.
- `faidx` - index a FASTA file and extract regions. Writes `<file>.fai` next to the FASTA; when regions are requested, an existing `.fai` is reused only if it is not older than the FASTA and fits its size, otherwise it is rebuilt silently. A symbolic link at the `.fai` path is an error.
- `sort` - sort records (reads everything into memory).
- `rmdup` - remove duplicate records, keeping the first.
- `rename` - number records that repeat an ID.
- `replace` - replace text in names or sequences.
- `fq2fa` - convert FASTQ to FASTA.
- `split` - split a file into parts by part count or record count. Never overwrites: an existing `<stem>.part_*<ext>` file in the output directory is an error unless `--force` is given, which unlinks those files (links themselves, not their targets) first.
- `sum` - order-independent checksum of sequence content.
- `qc` - verifiable QC numbers per file, files only, no stdin (reads, bases, lengths, GC, N, Q20, Q30, mean quality, sha256).
- `crosscheck` - compare read totals in an `olgu.json` with counts from the FASTQ files.

## Differences from seqkit and known limitations

- Behaviour was compared with seqkit on synthetic data; outputs match for the
  implemented commands except for the items below.
- No zstd, lz4 or BAM input.
- Output is never compressed.
- `sum` digests are not compatible with `seqkit sum`.
- Protein input with `fx2tab -g` (GC content) differs from seqkit.
- `qc` and `stats` share one implementation, so their numbers always agree.
- `fq2fa -o` refuses to write to a file that is also an input (same file,
  including links); other commands that write files never follow symbolic
  links.

## Known differences / design notes

- `sort` fails on duplicate IDs in every mode (consistent with the reference
  tool).
- `stats` and `qc` count N in the GC denominator; `fx2tab -g` excludes N
  (both consistent with the reference tool), so the three commands can report
  different GC for the same file.
- `seq -m/-M` filter on the original length, before gap removal.
- `locate -r` wraps the user regex in a lookahead group, so back-references
  (`\1`) and a mid-pattern `(?i)` do not work.
- For FASTA input `qc`/`stats` print Q20, Q30 and mean quality as 0: "absent"
  and "zero" cannot be told apart.
- Phred+64 input is accepted without a warning.
- Multi-line FASTQ is not supported (four lines per record).

## Tests

    python -m pytest -q

Optional reference comparison (needs a seqkit binary you provide; uses only
synthetic data in a temporary directory):

- `tools/diff_reference.py <seqkit-binary>` - runs a list of commands in both
  tools and reports differences.
- `tools/make_golden3.py <seqkit-binary>` - regenerates `tests/data3/` and
  its golden outputs, so the test suite itself needs no binary.

## qc and crosscheck

    fastxkit qc sample_R1.fastq.gz sample_R2.fastq.gz          # readable table
    fastxkit qc --json sample_R1.fastq.gz                      # 6-decimal JSON
    fastxkit qc --olgu sample_R1.fastq.gz > fastq_olgu.json    # olgu-style facts

`--olgu` writes `{"<fact>:<file stem>": {"deger", "kaynak", "konum"}}`, where
`konum` records the tool version and the file's sha256. Values are rounded
half up to 6 decimals in JSON modes.

    fastxkit crosscheck olgu.json --map map.csv [--tolerance N] [--count all|r1] [--json]

`map.csv` needs the header `sample,fastq` and optionally `read` (`r1`/`r2`);
paired-end files go on separate rows, relative paths are relative to the CSV.
The `ornek_tablosu` values of `olgu.json` are compared per sample with the
FASTQ read counts. Statuses: `OK`, `MISMATCH`, `MISSING_IN_MAP`,
`MISSING_IN_OLGU`, `FILE_ERROR`. Exit code: 0 if every sample is `OK`, 1 if at least one sample is `MISMATCH`/`MISSING_IN_MAP`/`MISSING_IN_OLGU`, 3 for input errors (bad `olgu.json`, bad map, nothing to compare, or any `FILE_ERROR`); 2 is argparse's usage error.
Use `--count r1` when the MultiQC total counts read pairs. Only paths and
numbers are printed, never sequence content.

## License

MIT (see `LICENSE`). seqkit is a separate MIT-licensed project; fastxkit does
not use any of its code and only compares behaviour against it.
