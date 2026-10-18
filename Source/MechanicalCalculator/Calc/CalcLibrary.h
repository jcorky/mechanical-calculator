// Values, the function registry and the helpers every library file uses to declare functions.
// A value is a matrix (a scalar is 1 x 1). A function that returns several named numbers returns a
// "record": a row vector whose elements carry field names, so `h = hohmann(r1, r2)` prints
// `dv1 = ...  dv2 = ...` and `h.dv1` or `h(1)` reads one back.
#pragma once

#include "CalcMatrix.h"

#include <string>
#include <utility>
#include <vector>

namespace mcalc
{
struct Value
{
	Matrix M;
	std::vector<std::string> Fields;   // empty, or one name per element (row-major)

	Value() = default;
	explicit Value(Real X) : M(Matrix::Scalar(X)) {}
	explicit Value(Matrix InM) : M(std::move(InM)) {}

	bool IsScalar() const { return M.IsScalar(); }
	bool IsRecord() const { return !Fields.empty(); }
	Real Scalar() const { return M.ScalarValue(); }
	int FieldIndex(const std::string& Name) const;
};

struct Result
{
	bool Ok = true;
	Value Val;
	std::string Error;

	static Result Of(Real X) { Result R; R.Val = Value(X); return R; }
	static Result Of(const Matrix& M) { Result R; R.Val = Value(M); return R; }
	static Result Of(const Value& V) { Result R; R.Val = V; return R; }
	static Result Fail(const std::string& Message) { Result R; R.Ok = false; R.Error = Message; return R; }
	// Named outputs, in order: Result::Record({{"dv1", a}, {"dv2", b}}).
	static Result Record(const std::vector<std::pair<std::string, Real>>& Items);
};

enum class ParamKind
{
	Number,   // must be a scalar
	Matrix,   // any shape (vectors and scalars included)
};

struct Param
{
	std::string Name;
	std::string Unit;        // "m", "m/s", "rad", "-" ...; shown in help
	ParamKind Kind = ParamKind::Number;
	bool HasDefault = false;
	Real Default = 0.0;
};

// Arguments as the function sees them: every declared parameter is present (defaults filled in).
class Call
{
public:
	explicit Call(const std::vector<Value>& InArgs) : Args(InArgs) {}
	size_t Count() const { return Args.size(); }
	Real N(size_t I) const { return I < Args.size() ? Args[I].Scalar() : 0.0; }
	int Int(size_t I) const;
	const Matrix& M(size_t I) const { return Args[I].M; }
	const Value& V(size_t I) const { return Args[I]; }

private:
	const std::vector<Value>& Args;
};

using NativeFn = Result (*)(const Call& C);

enum AngleUse
{
	AngleNone = 0,
	AngleArgs = 1,     // the first argument is an angle (sin, cos, tan)
	AngleResult = 2,   // the result is an angle (asin, acos, atan, atan2)
};

struct FunctionInfo
{
	std::string Name;
	std::string Category;    // keyboard library page, e.g. "orbits"
	std::string Source;      // where in Where It Isn't the relation lives, e.g. "orbits/ §2 Hohmann"
	std::string Summary;     // one line
	std::string Returns;     // what comes back, with units
	std::vector<Param> Params;
	bool Variadic = false;   // extra trailing arguments allowed (all Matrix kind)
	int Angle = 0;           // AngleArgs / AngleResult: elementary trig follows the DEG/RAD key
	NativeFn Fn = nullptr;

	std::string Signature() const;
	size_t RequiredCount() const;
};

class Registry
{
public:
	// The registry with every library registered (built on first use).
	static const Registry& Default();

	void Add(const FunctionInfo& Info);
	const FunctionInfo* Find(const std::string& Name) const;
	const std::vector<FunctionInfo>& All() const { return Functions; }
	std::vector<std::string> Categories() const;
	std::vector<const FunctionInfo*> InCategory(const std::string& Category) const;
	std::vector<std::string> Completions(const std::string& Prefix) const;
	// Names registered more than once (the later definition wins); the tests require none.
	const std::vector<std::string>& Duplicates() const { return Duplicated; }

private:
	std::vector<FunctionInfo> Functions;
	std::vector<std::string> Duplicated;
};

// Declaration helpers for library files.
inline Param Num(const char* Name, const char* Unit) { Param P; P.Name = Name; P.Unit = Unit; return P; }
inline Param Num(const char* Name, const char* Unit, Real Default) { Param P; P.Name = Name; P.Unit = Unit; P.HasDefault = true; P.Default = Default; return P; }
inline Param Mat(const char* Name, const char* Unit) { Param P; P.Name = Name; P.Unit = Unit; P.Kind = ParamKind::Matrix; return P; }

void Def(Registry& R, const char* Name, const char* Category, const char* Source, const char* Summary, const char* Returns, std::vector<Param> Params, NativeFn Fn);

// Shared constants, as the operator enters them: exact decimal values (CODATA / IAU / ISA,
// matching Where It Isn't's modules), and pi from the gears.
namespace cst
{
inline const Real& Pi() { return gm::Pi(); }
inline Real Tau() { return Real(2) * gm::Pi(); }
inline Real G0() { return Real(9.80665); }              // standard gravity, m/s^2
inline Real C() { return Real(299792458); }             // speed of light, m/s
inline Real Boltzmann() { return Real(1.380649e-23); }  // J/K
inline Real T0Noise() { return Real(290); }             // reference noise temperature, K
inline Real MuEarth() { return Real(3.986004418e14); }  // m^3/s^2
inline Real REarth() { return Real(6378.137e3); }       // equatorial radius, m
inline Real GammaAir() { return Real(1.4); }
inline Real CpAir() { return Real(1005); }              // J/(kg K)
inline Real RAir() { return Real(287.05); }             // J/(kg K)
} // namespace cst

// One registration function per library file; Registry::Default() calls them all.
void RegisterCoreMath(Registry& R);       // CalcBuiltins.cpp: elementary functions, statistics, matrices
void RegisterControlLibrary(Registry& R); // LibControl.cpp: Riccati/Mobius (ch03), dynamics and Lyapunov (ch08), MPC and approach benches
void RegisterFlightLibrary(Registry& R);  // LibFlight.cpp: equation sheet, atmosphere, flight model, planform, game
void RegisterRadarLibrary(Registry& R);   // LibRadar.cpp: radar equation, radar basics, lens antennas
void RegisterOrbitsLibrary(Registry& R);  // LibOrbits.cpp: Kepler, Hohmann, Lambert, CW, INS (Allan, Schuler, coning)
void RegisterOpsLibrary(Registry& R);     // LibOps.cpp: logistics (ch13 and bench), graphs (ch09), SLAM bench
void RegisterAtlasLibrary(Registry& R);   // LibAtlas.cpp: chapters 1, 2, 4, 5, 6 (groups, topology, Hilbert spaces, measure, optimization)
void RegisterFieldsLibrary(Registry& R);  // LibFields.cpp: chapters 7, 10, 11, 12 (lattices and GNSS, category theory, procedural systems, expansion fields)
} // namespace mcalc
