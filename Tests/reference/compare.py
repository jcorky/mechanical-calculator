#!/usr/bin/env python3
"""Checks regenerated reference cases against the committed ones.

Tolerances and provenance must match exactly; expressions too, except that numbers printed from
numpy may differ by 1e-13 relative. Expected values must agree to half of each case's own
tolerance: numpy, BLAS and libm builds differ in the last digits of long computations, which is
not a change in the case. Usage: compare.py <committed dir> <regenerated dir>
"""
import math
import pathlib
import re
import sys

NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def allowed(expected, tol):
    kind, _, size = tol.partition(":")
    size = float(size)
    budget = size * abs(expected) if kind == "rel" else size
    # Half the case's tolerance (long integrations amplify libm's last digits; the 200-orbit
    # leapfrog energy moves 1e-7 relative between builds against a 1e-6 tolerance).
    return max(budget / 2.0, 4e-16 * abs(expected))


def same_expression(a, b):
    """Equal text, with numbers printed from numpy allowed to differ in their last digits."""
    if a == b:
        return True
    if NUMBER.sub("#", a) != NUMBER.sub("#", b):
        return False
    for x, y in zip(NUMBER.findall(a), NUMBER.findall(b)):
        fx, fy = float(x), float(y)
        if fx != fy and abs(fx - fy) > 1e-13 * max(abs(fx), abs(fy)):
            return False
    return True


def number(text):
    try:
        return float(text)
    except ValueError:
        return None


def main():
    old_dir, new_dir = map(pathlib.Path, sys.argv[1:3])
    problems = 0
    for old in sorted(old_dir.glob("*.tsv")):
        new = new_dir / old.name
        a = old.read_text().splitlines()
        b = new.read_text().splitlines() if new.exists() else []
        if len(a) != len(b):
            print(f"{old.name}: {len(a)} lines committed, {len(b)} regenerated")
            problems += 1
            continue
        for n, (x, y) in enumerate(zip(a, b), 1):
            if x == y:
                continue
            cx, cy = x.split("\t"), y.split("\t")
            same_shape = len(cx) == len(cy) and len(cx) >= 3 and x[:1] != "#"
            if same_shape and same_expression(cx[0], cy[0]) and cx[2] == cy[2] and cx[3:] == cy[3:]:
                ex, ey = number(cx[1]), number(cy[1])
                if ex is not None and ey is not None and (ex == ey or (math.isfinite(ex) and abs(ex - ey) <= allowed(ex, cx[2]))):
                    continue
            print(f"{old.name}:{n}: changed\n  - {x[:200]}\n  + {y[:200]}")
            problems += 1
    for new in sorted(new_dir.glob("*.tsv")):
        if not (old_dir / new.name).exists():
            print(f"{new.name}: new file not committed")
            problems += 1
    print(f"{problems} difference(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
