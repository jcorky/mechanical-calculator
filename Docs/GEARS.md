# How the gears compute

The calculator has one arithmetic unit: a stepped-drum calculating machine in the line of
Leibniz's Stepped Reckoner (1694) and Thomas de Colmar's Arithmometer (1820). Every number the
calculator touches is an `mcalc::Real`, and every `+`, `−`, `×`, `÷` and `√` on a `Real` is one of the
machine operator's procedures below. Nothing else does arithmetic: the elementary functions, the
linear algebra and the 400-odd library functions are all sequences of these five.

## The machine

| Part | Size | Role |
|---|---|---|
| Setting register | 20 digits | one Leibniz drum per digit: a cylinder with nine teeth of stepped lengths. The set digit `d` slides a pinion along its drum to where it meets `d` teeth per revolution. |
| Accumulator | 40 numeral wheels | on the carriage. Drum `i` drives wheel `i + shift`. Each wheel has a tens-carry lever that advances the next wheel by one when it passes 9 → 0 (or 0 → 9 reversed). |
| Counter | 21 numeral wheels | on the carriage. The wheel at the carriage position moves one step per crank turn. |
| Carriage | shifts 0 … 20 | slides the accumulator and counter past the drums, one place per step. |
| Crank, reverser, bell | | one turn adds the setting into the accumulator at the carriage position (or, reversed, subtracts it). A borrow out of the top wheel is an overdraft, and rings the bell. |

The model is `SteppedDrumMachine` in `Source/MechanicalCalculator/Calc/GearMachine.h`. It is literal:
`Turn()` advances each wheel by its drum's digit, then trips the carry levers in order up the
register, exactly as the levers ripple on the real mechanism.

## The operator's procedures

A value is kept between operations the way an operator writes it down: a sign, 18 significant
digits and a power of ten. Before an operation they set the operands into the machine at carriage
positions that line up the decimal points. After it they read the result and round it to 18
digits by the usual rule: add 5 under the 19th digit, then read 18.

| Operation | Procedure | Turns |
|---|---|---|
| `a + b`, `a − b` | Set the operand with the larger exponent, shift it up by the exponent gap (at most 20) and turn once. Set the other at shift 0 and turn once, reversed for a difference of magnitudes, larger minus smaller. | 2 |
| `a × b` | Set the multiplicand. For each multiplier digit `k` (units first), shift to `k` and turn that many times. The counter then shows the multiplier and the accumulator the exact 36-digit product. | the multiplier's digit sum |
| `a ÷ b` | Load the dividend into the top of the accumulator and set the divisor. At each place from shift 20 down: turn reversed until the bell rings, one turn forward to restore, then shift down. The counter reads 19 quotient digits. | Σ (digit + 2) |
| `√a` | Toepler's method: load `a` (made even in its power of ten), then at each column pair subtract the odd numbers 20y + 1, 20y + 3, … until the bell. The count is the next root digit. | Σ (digit + 2) |
| a series | The terms are turned into the accumulator without clearing it (`RealAccumulator`), and the sum is read once at the end. A product can be turned in whole (`AddProduct`). | one per term |

When the exponents differ by more than 20 places, the smaller operand would land below the
accumulator's last wheel, so the sum is the larger operand. Division by zero rings the bell forever:
the result is `inf`, or `nan` for 0/0. These are the only special values.

## Fast path and proof

Turning thousands of wheels per keystroke is too slow for a game, so each procedure also has a
fast path. It computes the same register contents with integer arithmetic: the exact product, the
same digit-by-digit repeated subtraction, and the same Toepler digit pairs. `Tests/TestGear.cpp`
runs 3000 random operand pairs through both paths, including near-cancellations and runs of 9s and
0s that stress the carries, and requires identical digits for all five operations.

When the Unreal layer evaluates an expression, it attaches a `GearRecorder`. The first operations
then run on the literal machine, and every set, shift and turn is recorded with full register
snapshots. `MechanismPlayer` (`Keyboard/`) turns that record into motion: the crank, the drums, the
pinions sliding to their digits, each numeral wheel stepping through the stepped teeth, the carry
levers tripping in turn, the carriage gliding, and the bell. The rest of a long computation takes
the fast path, which gives the same result.

## Elementary functions and constants

`GearMath.h` builds `exp`, `log`, the trigonometric, hyperbolic and inverse functions, `erf`,
`erfc`, `lgamma` and `tgamma` from the five procedures:

- **exp:** x = k ln 10 + r, where 10^k is a carriage shift, and the Taylor series for r is summed in
  the accumulator.
- **log:** x = m · 10^e and m = 2^j y, then 2 atanh((y − 1)/(y + 1)).
- **sin and cos:** reduction by k · π/2, with π split in two (Cody and Waite) and both products
  turned in whole, followed by Taylor series.
- **atan:** two angle-halvings, then its series.
- **erf, erfc:** an all-positive series for erf, and a continued fraction for erfc.
- **lgamma, tgamma:** Stirling's series after a recurrence shift, and reflection for x < ½.

π, ln 2 and ln 10 are computed on the gears the first time they are used: Machin's formula and
atanh series, with each term formed by exact long division in a 72-digit fixed-point register. So
they are correctly rounded to 18 digits.

Against `<cmath>`, the worst relative error over the test grid is 1.9 × 10⁻¹⁵ (lgamma next to its
zero at 1). Elsewhere it is a few units in the 17th digit.

## Consequences you can see

- Decimal inputs are exact: `0.1 + 0.2 == 0.3`. A typed number is entered digit for digit, and a
  `double` literal in the library code is taken as its shortest decimal.
- Nothing overflows on paper: `1e300 * 1e300` is `1e600`.
- Results carry 18 significant digits rather than about 16 for binary doubles. They can differ from
  the JavaScript reference values in the last digits wherever a computation is ill-conditioned. The
  few reference cases where this exceeded the original tolerance are listed in the generators,
  with the reason and an independent high-precision check of which value is closer.
- Each result on the tape is followed by the work it took, for example `~ 1290 turns, 431 shifts`.
