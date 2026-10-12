"""Pure sequence, quality and length-distribution computations."""
import collections
import math
from collections.abc import Collection, Sequence

PHRED_OFFSET = 33
GAP_CHARS = "- ."
_DNA_LETTERS = frozenset("ACGTRYMKSWHBVDN" + GAP_CHARS)
_RNA_LETTERS = frozenset("ACGURYMKSWHBVDN" + GAP_CHARS)
_PROTEIN_LETTERS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO*" + GAP_CHARS)


def median(ordered: Sequence[float]) -> float:
    """Median of an already sorted, non-empty sequence."""
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def quartiles(lengths: Sequence[int]) -> tuple[float, float, float]:
    """Returns (Q1, Q2, Q3); the median belongs to both halves when n is odd."""
    if not lengths:
        return 0, 0, 0
    ordered = sorted(lengths)
    half = (len(ordered) + 1) // 2
    return median(ordered[:half]), median(ordered), median(ordered[-half:])


def n50(lengths: Sequence[int]) -> tuple[int, int]:
    """Returns (N50, number of sequences needed to reach it)."""
    total = sum(lengths)
    cumulative = 0
    for count, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if 2 * cumulative >= total:
            return length, count
    return 0, 0


def quality_profile(quality: str) -> tuple[float, int, int]:
    """Returns (sum of error probabilities, #bases >= Q20, #bases >= Q30)."""
    error_sum = 0.0
    at_least_20 = at_least_30 = 0
    for char, count in collections.Counter(quality).items():
        score = ord(char) - PHRED_OFFSET
        if score < 0:
            raise ValueError(f"invalid quality character {char!r}")
        error_sum += count * 10 ** (-score / 10)
        if score >= 20:
            at_least_20 += count
        if score >= 30:
            at_least_30 += count
    return error_sum, at_least_20, at_least_30


def error_to_quality(error_sum: float, bases: int) -> float:
    """Converts a summed error probability over `bases` into a Phred score."""
    if bases == 0 or error_sum <= 0:
        return 0
    return -10 * math.log10(error_sum / bases)


def mean_quality(quality: str) -> float:
    """Average quality via mean error probability (not mean of scores)."""
    return error_to_quality(quality_profile(quality)[0], len(quality))


def gc_percent(seq: str) -> float:
    """GC percentage over non-N bases."""
    counts = collections.Counter(seq.upper())
    valid = len(seq) - counts["N"]
    if valid == 0:
        return 0
    return (counts["G"] + counts["C"]) / valid * 100


def base_content(seq: str, bases: str) -> float:
    """Percentage of the sequence length made of any of `bases` (any case)."""
    if not seq:
        return 0
    wanted = set(bases.upper())
    counts = collections.Counter(seq.upper())
    return sum(counts[base] for base in wanted) / len(seq) * 100


def detect_type(letters: Collection[str]) -> str:
    """Classifies a set of sequence letters as DNA, RNA, Protein or Unlimit."""
    upper = {letter.upper() for letter in letters}
    if not upper:
        return "-"
    for name, alphabet in (("DNA", _DNA_LETTERS), ("RNA", _RNA_LETTERS),
                           ("Protein", _PROTEIN_LETTERS)):
        if upper <= alphabet:
            return name
    return "Unlimit"
