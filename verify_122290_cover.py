"""Verify a CRT cover for the interior of an even gap, defaulting to 122,290."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys


GAP = 122290
CERTIFICATE = Path(__file__).with_name("gap_122290_cover.json")


def small_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def verify_fermat_witness(witness):
    offset, value, number = witness
    try:
        from gmpy2 import powmod
    except ImportError:
        powmod = pow
    if powmod(value, number - 1, number) == 1:
        raise ValueError(f"Compositeness check failed at offset {offset}")


def verify_remaining_offsets(covered, gap, left, witnesses, *, jobs=1):
    offsets = set()
    fermat_witnesses = []
    for offset, kind, value in witnesses:
        if (type(offset) is not int or not 0 < offset < gap
                or covered[offset] or offset in offsets):
            raise ValueError(f"Unexpected interior offset {offset}")
        number = left + offset
        if type(value) is not int or not 1 < value < number:
            raise ValueError(f"Invalid compositeness witness at offset {offset}")
        if kind == "fermat":
            fermat_witnesses.append((offset, value, number))
        elif kind != "divisor" or number % value:
            raise ValueError(f"Compositeness check failed at offset {offset}")
        offsets.add(offset)

    if jobs > 1 and len(fermat_witnesses) > 1:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            for _ in executor.map(verify_fermat_witness, fermat_witnesses, chunksize=4):
                pass
    else:
        for witness in fermat_witnesses:
            verify_fermat_witness(witness)
    for offset in offsets:
        covered[offset] = 1


def verify_cover(path, expected_gap=GAP, *, row=None, interior=(), jobs=1):
    data = json.loads(path.read_text(encoding="ascii"))
    if data["gap"] != expected_gap or expected_gap < 2 or expected_gap % 2:
        raise ValueError("Wrong gap size")

    residue, modulus = 0, 1
    covered = bytearray(expected_gap + 1)
    previous = 1
    for prime, value in data["classes"]:
        if prime <= previous or not small_prime(prime):
            raise ValueError(f"Out-of-order or nonprime modulus {prime}")
        if not 0 < value < prime or (value + expected_gap) % prime == 0:
            raise ValueError(f"An endpoint is blocked modulo {prime}")
        first = (-value) % prime
        if first <= expected_gap:
            covered[first::prime] = b"\x01" * ((expected_gap - first) // prime + 1)
        residue += modulus * ((value - residue) * pow(modulus, -1, prime) % prime)
        modulus *= prime
        previous = prime

    if not 0 <= residue < modulus or modulus <= previous:
        raise ArithmeticError("Unexpected CRT representative or modulus")
    for prime, value in data["classes"]:
        if residue % prime != value:
            raise ArithmeticError(f"CRT check failed at {prime}")
    if interior:
        if type(row) is not int or row < 1:
            raise ValueError("Offset witnesses require a positive row multiplier")
        verify_remaining_offsets(
            covered, expected_gap, residue + row * modulus, interior, jobs=jobs
        )
    if covered[0] or covered[expected_gap] or covered.count(1) != expected_gap - 1:
        raise ArithmeticError("Interior coverage failed")

    print(f"Verified all {expected_gap - 1:,} interior offsets")
    print(f"{len(data['classes']):,} distinct prime moduli, largest {previous}")
    print(f"W has {len(str(modulus)):,} decimal digits")
    print(f"SHA-256 of decimal W {sha256(str(modulus).encode('ascii')).hexdigest()}")
    print(f"SHA-256 of decimal a {sha256(str(residue).encode('ascii')).hexdigest()}")
    return residue, modulus


def main():
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--gap", type=int, default=GAP, help="Expected even gap size")
    parser.add_argument("--candidate", type=Path, help="Check both endpoint integers from a JSON file")
    parser.add_argument("--row", type=int, help="Positive row multiplier for a + row*W")
    parser.add_argument("--prove", action="store_true", help="Prove both endpoints using python-flint")
    args = parser.parse_args()
    if (args.row is not None and args.row < 1) or ((args.prove or args.candidate) and args.row is None):
        parser.error("Use a positive --row with --prove or --candidate")

    candidate = json.loads(args.candidate.read_text(encoding="ascii")) if args.candidate else None
    residue, modulus = verify_cover(
        args.certificate, args.gap,
        row=args.row, interior=candidate.get("interior_composites", ()) if candidate else (),
    )
    if args.row is not None:
        left = residue + args.row * modulus
        right = left + args.gap
        if args.candidate:
            actual = (candidate["gap"], candidate["row"], candidate["left"], candidate["right"])
            if actual != (args.gap, args.row, str(left), str(right)):
                raise ArithmeticError("Candidate disagrees with the CRT cover")
            print(f"Matched both endpoints at row {args.row}")
        if args.prove:
            from flint import fmpz

            if not fmpz(left).is_prime() or not fmpz(right).is_prime():
                raise ArithmeticError("At least one endpoint is composite")
            print(f"Row {args.row} has two prime endpoints and gap {args.gap}")
        elif not args.candidate:
            print(f"Row {args.row} left endpoint")
            print(left)
            print(f"Row {args.row} right endpoint")
            print(right)


if __name__ == "__main__":
    main()
