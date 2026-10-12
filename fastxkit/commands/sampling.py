"""sample and shuffle: seeded random selection and ordering."""
import argparse
import random
import sys
from collections.abc import Iterable, Iterator

from fastxkit import records
from fastxkit.commands import _common

DEFAULT_SEED = 11


def sample_count(items: Iterable, count: int, rng: random.Random) -> Iterator:
    """Reservoir-samples `count` items, returned in input order."""
    reservoir: list[tuple[int, object]] = []
    for index, item in enumerate(items):
        if index < count:
            reservoir.append((index, item))
            continue
        slot = rng.randint(0, index)
        if slot < count:
            reservoir[slot] = (index, item)
    reservoir.sort(key=lambda pair: pair[0])
    return (item for _, item in reservoir)


def sample_fraction(items: Iterable, fraction: float,
                    rng: random.Random) -> Iterator:
    """Keeps each item independently with probability `fraction`."""
    return (item for item in items if rng.random() < fraction)


def _fraction(text: str) -> float:
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {text!r}") from None
    if not 0 <= value <= 1:
        raise argparse.ArgumentTypeError(f"must be within 0..1: {value}")
    return value


def _write(selected: Iterable[records.Record]) -> None:
    for record in selected:
        sys.stdout.write(records.format_record(record))


def run_sample(args) -> None:
    rng = random.Random(args.rand_seed)
    source = records.iter_records(args.files)
    if args.number is not None:
        _write(sample_count(source, args.number, rng))
    else:
        _write(sample_fraction(source, args.proportion, rng))


def run_shuffle(args) -> None:
    everything = list(records.iter_records(args.files))
    random.Random(args.rand_seed).shuffle(everything)
    _write(everything)


def _add_seed(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-s", "--rand-seed", type=int, default=DEFAULT_SEED,
                        help="random seed (default: %(default)s)")


def register(subparsers) -> None:
    sample = subparsers.add_parser("sample", help="sample records at random")
    _common.add_input_files(sample)
    amount = sample.add_mutually_exclusive_group(required=True)
    amount.add_argument("-n", "--number", type=_common.non_negative_int,
                        help="exact number of records")
    amount.add_argument("-p", "--proportion", type=_fraction,
                        help="probability of keeping each record (0..1)")
    _add_seed(sample)
    sample.set_defaults(func=run_sample)

    shuffle = subparsers.add_parser("shuffle", help="shuffle record order")
    _common.add_input_files(shuffle)
    _add_seed(shuffle)
    shuffle.set_defaults(func=run_shuffle)
