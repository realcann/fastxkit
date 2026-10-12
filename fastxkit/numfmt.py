"""Number formatting with half-up rounding (Python's round() is half-even)."""
import decimal


def _quantize(value: float, digits: int) -> decimal.Decimal:
    exponent = decimal.Decimal(1).scaleb(-digits)
    return decimal.Decimal(repr(float(value))).quantize(
        exponent, rounding=decimal.ROUND_HALF_UP)


def fixed(value: float, digits: int) -> str:
    """Formats with exactly `digits` decimals, rounding half up."""
    return f"{_quantize(value, digits):.{digits}f}"


def grouped(number: int) -> str:
    """Formats an integer with thousands separators."""
    return f"{number:,}"


def trimmed(value: float, digits: int) -> str:
    """Rounds half up to `digits` decimals, drops trailing zeros, groups."""
    text = f"{_quantize(value, digits):,.{digits}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text
