// The number type of the gear-driven calculator. Every arithmetic operation on a Real is carried
// out by the stepped-drum machine (GearMachine.h): a 20-digit setting register with a Leibniz drum
// per digit, a 40-digit accumulator and a 20-digit counter on a shifting carriage.
//
//   a + b, a - b   the operator sets one operand, turns it into the accumulator at the carriage
//                  position that aligns the decimal points, sets the other and turns it in
//                  forward (add) or reversed (subtract)
//   a * b          multiplicand set; for each multiplier digit, that many turns, then a shift
//   a / b          dividend in the accumulator, divisor set; turn backwards until the bell rings
//                  (overdraft), one turn forward to restore, shift; the counter reads the quotient
//   sqrt(a)        Toepler's method: subtract the odd numbers 1, 3, 5, ... column pair by column pair
//
// A value is held the way the operator writes it down between operations: a sign, 18 significant
// decimal digits and a power of ten. Results are rounded to 18 digits by the operator's rule
// (add 5 under the 19th digit, then read 18). Infinity and NaN stand for the bell ringing
// without end (division by zero) and for an invalid procedure.
//
// To keep keystrokes interactive the machine's procedures have a fast path that computes the
// same register contents with integer arithmetic. The tests run both on the same operands and
// require identical digits; GearContext can route operations through the literal wheel-by-wheel
// model instead, which is what the Unreal layer records to animate the mechanism.
#pragma once

#include <cstdint>
#include <string>

namespace mcalc
{
class Real
{
public:
	static constexpr int Digits = 18;            // significant digits held between operations

	Real() = default;
	Real(int X) { SetInteger(static_cast<long long>(X)); }
	Real(long X) { SetInteger(static_cast<long long>(X)); }
	Real(long long X) { SetInteger(X); }
	Real(unsigned X) { SetUnsigned(static_cast<unsigned long long>(X)); }
	Real(unsigned long X) { SetUnsigned(static_cast<unsigned long long>(X)); }
	Real(unsigned long long X) { SetUnsigned(X); }
	// A double literal or input is taken as the shortest decimal that reads back as that double,
	// so 0.1 is entered on the keys as 0.1, not as the binary fraction nearest to it.
	Real(double X);

	static Real FromParts(bool Negative, uint64_t Mantissa, int Exponent);   // value = +-Mantissa * 10^(Exponent - 17)
	static Real Infinity(bool Negative = false);
	static Real NaN();
	static Real Epsilon();                       // 1e-17, the spacing of 18-digit values at 1
	static Real Parse(const char* Text, size_t* Consumed);   // decimal text, exactly (rounded to 18 digits)

	bool IsNaN() const { return Kind == KNaN; }
	bool IsInf() const { return Kind == KInf; }
	bool IsFinite() const { return Kind == KFinite; }
	bool IsZero() const { return Kind == KFinite && Mant == 0; }
	bool IsNegative() const { return Neg && !IsNaN() && !IsZero(); }
	bool SignBit() const { return Neg; }
	uint64_t Mantissa() const { return Mant; }
	int Exponent() const { return Exp; }        // floor(log10 |x|) for finite non-zero x

	Real operator-() const;
	Real operator+() const { return *this; }
	Real& operator+=(const Real& B) { *this = *this + B; return *this; }
	Real& operator-=(const Real& B) { *this = *this - B; return *this; }
	Real& operator*=(const Real& B) { *this = *this * B; return *this; }
	Real& operator/=(const Real& B) { *this = *this / B; return *this; }

	friend Real operator+(const Real& A, const Real& B) { return Add(A, B, false); }
	friend Real operator-(const Real& A, const Real& B) { return Add(A, B, true); }
	friend Real operator*(const Real& A, const Real& B) { return Mul(A, B); }
	friend Real operator/(const Real& A, const Real& B) { return Div(A, B); }
	friend bool operator==(const Real& A, const Real& B) { return Compare(A, B) == 0; }
	friend bool operator!=(const Real& A, const Real& B) { return Compare(A, B) != 0; }
	friend bool operator<(const Real& A, const Real& B) { return Compare(A, B) == -1; }
	friend bool operator>(const Real& A, const Real& B) { return Compare(A, B) == 1; }
	friend bool operator<=(const Real& A, const Real& B) { const int C = Compare(A, B); return C == -1 || C == 0; }
	friend bool operator>=(const Real& A, const Real& B) { const int C = Compare(A, B); return C == 1 || C == 0; }

	// Reading the registers out (not arithmetic): for display, loop counts and the engine.
	double ToDouble() const;
	long long ToLong() const;   // truncates toward zero; saturates; NaN -> 0
	int ToInt() const;
	// %g-style: Sig significant digits, exponent shown as e6 / e-7 when needed.
	std::string ToString(int Sig = Digits) const;

	static Real Sqrt(const Real& A);             // Toepler's method on the machine
	static Real Abs(const Real& A) { Real R = A; R.Neg = false; return R; }
	static Real Trunc(const Real& A);
	static Real Floor(const Real& A);
	static Real Ceil(const Real& A);
	static Real Round(const Real& A);             // half away from zero
	static Real Scale10(const Real& A, int Power);   // exact: moves the decimal point
	static int Compare(const Real& A, const Real& B);   // -1, 0, 1, or 2 if unordered (NaN)

private:
	enum KindType : uint8_t
	{
		KFinite,
		KInf,
		KNaN,
	};
	uint64_t Mant = 0;   // 0, or 18 digits: 1e17 <= Mant < 1e18
	int Exp = 0;
	bool Neg = false;
	KindType Kind = KFinite;

	void SetInteger(long long X);
	void SetUnsigned(unsigned long long X);
	static Real Add(const Real& A, const Real& B, bool Subtract);
	static Real Mul(const Real& A, const Real& B);
	static Real Div(const Real& A, const Real& B);
};

// A running sum kept in the accumulator without reading it off between terms: the operator turns
// each term in at its place and rounds to 18 digits only once, at the end. Series (exp, sin,
// atan, the constants) are summed this way, as they would be at the machine. The accumulator holds
// 40 digits; terms below its 40th place relative to the largest term are lost, as on the wheels.
class RealAccumulator
{
public:
	RealAccumulator();
	void Add(const Real& Term);
	void Sub(const Real& Term) { Add(-Term); }
	// Multiplies into the accumulator without clearing it: the exact 36-digit product is added.
	void AddProduct(const Real& A, const Real& B, bool Subtract = false);
	Real Result() const;

private:
	// Positive and negative parts, each value = Digits * 10^LowExp, kept as base-1e9 limbs.
	static constexpr int Limbs = 6;
	uint64_t Pos[Limbs];
	uint64_t NegPart[Limbs];
	long long LowExp = 0;
	bool Started = false;
	bool Special = false;
	Real SpecialValue;
	void Rebase(long long NewLowExp);
	// Adds Value * 10^ValueLowExp (Value up to 36 digits as base-1e18 limbs).
	void AddDigits(const uint64_t* Value18, long long ValueLowExp, bool Negative);
};

// Work done by the machine on this thread (each operation counted once, on either path).
struct GearStats
{
	unsigned long long Additions = 0;      // add and subtract procedures
	unsigned long long Multiplications = 0;
	unsigned long long Divisions = 0;
	unsigned long long Roots = 0;
	unsigned long long Turns = 0;          // crank revolutions
	unsigned long long Shifts = 0;         // carriage steps
	unsigned long long Carries = 0;        // tens-carries (literal path only; estimated otherwise)
	unsigned long long Total() const { return Additions + Multiplications + Divisions + Roots; }
};

class GearRecorder;   // GearMachine.h

namespace GearContext
{
GearStats& Stats();                // this thread's counters
void ResetStats();
// While a recorder is attached, operations run on the literal machine and are recorded until
// the recorder says it has enough; the rest take the fast path (identical results).
void Attach(GearRecorder* Recorder);
GearRecorder* Attached();
} // namespace GearContext
} // namespace mcalc
