"""Pure sequence/quality computations."""
import pytest

from fastxkit import metrics


def test_quartiles_odd_and_even():
    assert metrics.quartiles([10, 60, 75, 130, 200]) == (60, 75, 130)
    assert metrics.quartiles([50, 50, 50, 75, 75, 75, 75, 75, 100, 100, 100, 150]) == (62.5, 75, 100)
    assert metrics.quartiles([7]) == (7, 7, 7)
    assert metrics.quartiles([]) == (0, 0, 0)


def test_n50():
    assert metrics.n50([200, 130, 75, 60, 10]) == (130, 2)
    assert metrics.n50([]) == (0, 0)
    assert metrics.n50([5]) == (5, 1)
    assert metrics.n50([2, 2, 2, 2]) == (2, 2)  # exactly half is enough


def test_quality_profile_and_mean():
    err_sum, at20, at30 = metrics.quality_profile("!!II")
    assert (at20, at30) == (2, 2)
    assert err_sum == pytest.approx(2 + 2e-4)
    assert metrics.mean_quality("IIII") == pytest.approx(40)
    assert metrics.mean_quality("") == 0


def test_gc_percent_ignores_n():
    assert metrics.gc_percent("GCNNAT") == pytest.approx(50)
    assert metrics.gc_percent("NNN") == 0
    assert metrics.gc_percent("gcat") == pytest.approx(50)


def test_base_content_is_case_insensitive_share_of_length():
    assert metrics.base_content("AaTtGC", "AT") == pytest.approx(4 / 6 * 100)
    assert metrics.base_content("", "AT") == 0


@pytest.mark.parametrize("seq, expected", [
    ("ACGTN", "DNA"), ("ACGUN", "RNA"), ("MKVLAAGIVG", "Protein"),
    ("AC1T!", "Unlimit"), ("ACGT-.", "DNA"), ("", "-")])
def test_detect_type(seq, expected):
    assert metrics.detect_type(set(seq)) == expected
