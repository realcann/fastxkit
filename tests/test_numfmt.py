"""Number formatting helpers (half-up rounding, trimming, grouping)."""
import pytest

from fastxkit import numfmt


@pytest.mark.parametrize("value, digits, expected", [
    (81.25, 1, "81.3"),      # banker's rounding would give 81.2
    (0.125, 2, "0.13"),
    (95.0, 1, "95"),
    (9.8987, 2, "9.9"),
    (0, 2, "0"),
    (1234567.891, 2, "1,234,567.89"),
])
def test_trimmed(value, digits, expected):
    assert numfmt.trimmed(value, digits) == expected


def test_grouped_int():
    assert numfmt.grouped(999) == "999"
    assert numfmt.grouped(1000) == "1,000"


@pytest.mark.parametrize("value, expected", [
    (9.095013, "9.10"), (0.125, "0.13"), (70, "70.00"), (0, "0.00")])
def test_fixed2(value, expected):
    assert numfmt.fixed(value, 2) == expected
