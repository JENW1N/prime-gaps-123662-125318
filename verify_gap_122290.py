"""Verify a recorded gap candidate, optionally with two ECPP certificates."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys

from verify_122290_cover import GAP, verify_cover

try:
    from gmpy2 import invert, mpz
except ImportError:
    mpz = int

    def invert(value, modulus):
        return pow(value, -1, modulus)


HERE = Path(__file__).resolve().parent
RECORDS = {
    122290: (
        530520,
        "gap_122290_cover_b.json",
        "2d7a21344e3dd26dfb667033b31d7f40a1459b266b3c837ded335fd3e0a27d05",
        "b450c29e7ac6bbe29083f66901baeccdd30e24fb3c2eb5ac82afc02ad1982a47",
    ),
    123662: (
        1523845,
        "gap_123662_cover.json",
        "b3acacbe505500ae76899e82ef63372cf82ad50295bc9d594bc3e9d14bcaf6cd",
        "eee9ce65ccc0ffd6cbc4d81f9b63df0f20805466363dcb21fff433c3d2911eee",
    ),
    125318: (
        5046505,
        "gap_125318_cover.json",
        "efa03b02aeef1364bab3d967db59202cdeb1e8da4586a315520dca1be3a4d011",
        "b1f5bcd673d8612de5eaa89f624f9bac563e3544e910927ecb5461977402a156",
    ),
}
ROW = RECORDS[GAP][0]
ENDPOINT_HASHES = RECORDS[GAP][2:]
SMALL_LIMIT = 1 << 64
# Sorenson and Webster, Math. Comp. 86 (2017), 985–1003.
SMALL_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_small_prime(n):
    if not 2 <= n < SMALL_LIMIT:
        return False
    for base in SMALL_BASES:
        if n % base == 0:
            return n == base
    odd_part = n - 1
    twos = 0
    while odd_part % 2 == 0:
        odd_part //= 2
        twos += 1
    for base in SMALL_BASES:
        residue = pow(base, odd_part, n)
        if residue in (1, n - 1):
            continue
        for _ in range(twos - 1):
            residue = residue * residue % n
            if residue == n - 1:
                break
        else:
            return False
    return True


def add_points(left, right, a, n):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 + y2) % n == 0:
            return None
        if y1 != y2:
            raise ValueError("Incompatible reductions of points with the same x-coordinate")
        numerator, denominator = 3 * x1 * x1 + a, 2 * y1
    else:
        numerator, denominator = y2 - y1, x2 - x1
    try:
        slope = numerator * invert(denominator % n, n) % n
    except (ZeroDivisionError, ValueError) as exc:
        raise ValueError("A point-addition denominator is not a unit") from exc
    x3 = (slope * slope - x1 - x2) % n
    return x3, (slope * (x1 - x3) - y1) % n


def multiply_point(k, point, a, n):
    result = None
    for bit in bin(k)[2:]:
        result = add_points(result, result, a, n)
        if bit == "1":
            result = add_points(result, point, a, n)
    return result


def large_enough(q, n):
    # q > (n**(1/4) + 1)**2, checked using integer arithmetic.
    positive_part = q * q + 6 * q + 1 - n
    return positive_part > 0 and positive_part**2 > 16 * (q + 1)**2 * q


def verify_step(step):
    index, n, cofactor, a, point, q = step
    n = mpz(n)
    a = mpz(a) % n
    x, y = (mpz(coordinate) % n for coordinate in point)
    b = (y * y - x * (x * x + a)) % n
    if gcd((4 * a**3 + 27 * b**2) % n, n) != 1:
        raise ValueError(f"The curve has bad reduction at step {index}")
    reduced = multiply_point(cofactor, (x, y), a, n)
    if reduced is None or multiply_point(q, reduced, a, n) is not None:
        raise ValueError(f"The point does not certify step {index}")


def verify_ecpp(certificate, expected, *, partial=False, jobs=1):
    if type(certificate) is int:
        if certificate != expected or not is_small_prime(certificate):
            raise ValueError("Invalid single-prime certificate")
        return certificate, 0
    if not isinstance(certificate, list) or not certificate:
        raise ValueError("Expected a nonempty list of ECPP steps")

    steps = []
    child = expected
    for index, step in enumerate(certificate, start=1):
        if not isinstance(step, list) or len(step) != 5:
            raise ValueError(f"Malformed step {index}")
        n, trace, cofactor, a, point = step
        if (not isinstance(point, list) or len(point) != 2
                or any(type(value) is not int for value in (n, trace, cofactor, a, *point))):
            raise ValueError(f"Malformed coefficients at step {index}")
        if n != child or n <= 3 or gcd(n, 6) != 1 or trace**2 >= 4 * n:
            raise ValueError(f"Invalid modulus or trace at step {index}")
        order = n + 1 - trace
        if cofactor < 1 or order % cofactor:
            raise ValueError(f"Invalid cofactor at step {index}")
        child = order // cofactor
        if child < 2 or not large_enough(child, n):
            raise ValueError(f"Invalid child at step {index}")
        steps.append((index, n, cofactor, a, point, child))

    if not partial and not is_small_prime(child):
        raise ValueError("Certificate does not end at a proved prime below 2**64")
    if jobs == 1:
        for step in reversed(steps):
            verify_step(step)
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            for _ in executor.map(verify_step, steps, chunksize=1):
                pass
    return child, len(steps)


def read_ecpp(path):
    contents = path.read_text(encoding="ascii").strip()
    if contents.startswith("c = ") and contents.endswith(";"):
        contents = contents[4:-1]
    return json.loads(contents)


def verified_endpoints(gap, cover=None, candidate_path=None, *, jobs=1):
    """Reconstruct the endpoints and check the interior and endpoint hashes."""
    row, cover_file, left_hash, right_hash = RECORDS[gap]
    cover = cover or HERE / cover_file
    candidate_path = candidate_path or HERE / f"gap_{gap}_candidate.json"
    candidate = json.loads(candidate_path.read_text(encoding="ascii"))
    if candidate["gap"] != gap or candidate["row"] != row:
        raise ValueError("Unexpected gap size or multiplier")
    a, modulus = verify_cover(
        cover, gap, row=row, interior=candidate.get("interior_composites", ()), jobs=jobs
    )
    left = a + row * modulus
    endpoints = (left, left + gap)
    for side, endpoint, expected_hash in zip(("left", "right"), endpoints, (left_hash, right_hash)):
        decimal = str(endpoint)
        if decimal != candidate[side]:
            raise ValueError(f"Incorrect {side} endpoint in the candidate file")
        if sha256(decimal.encode("ascii")).hexdigest() != expected_hash:
            raise ValueError("Unexpected endpoint reconstructed from the cover")
    return endpoints


def main():
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gap", type=int, choices=RECORDS, default=GAP)
    parser.add_argument("--cover", type=Path, help="Override the recorded cover file")
    parser.add_argument("--candidate", type=Path, help="Override the recorded candidate file")
    parser.add_argument("--left", type=Path, help="Left endpoint ECPP certificate")
    parser.add_argument("--right", type=Path, help="Right endpoint ECPP certificate")
    parser.add_argument("--jobs", type=int, default=1, help="Workers for certificate checks")
    args = parser.parse_args()
    if (args.left is None) != (args.right is None):
        parser.error("Supply both ECPP certificates together")
    if args.jobs < 1:
        parser.error("--jobs must be positive")

    gap = args.gap
    endpoints = verified_endpoints(gap, args.cover, args.candidate, jobs=args.jobs)
    if args.left is None:
        try:
            from gmpy2 import is_probab_prime
        except ImportError:
            parser.error("Install gmpy2 to check the probable-prime endpoints")
        for side, endpoint in zip(("left", "right"), endpoints):
            if not is_probab_prime(endpoint, 32):
                raise ValueError(f"The {side} endpoint is composite")
            print(f"{side.capitalize()} endpoint passed GMP's probable-prime test", flush=True)
        print(f"Verified gap {gap:,} with probable-prime endpoints and a composite interior")
        return

    for endpoint, path in zip(endpoints, (args.left, args.right)):
        _, count = verify_ecpp(read_ecpp(path), endpoint, jobs=args.jobs)
        print(f"Proved a {len(str(endpoint)):,}-digit endpoint prime in {count} ECPP steps")
    print(f"Proved consecutive primes at distance {gap:,}")


if __name__ == "__main__":
    main()
