#include "CalcLibrary.h"

#include <algorithm>
#include <cmath>

namespace mcalc
{
int Value::FieldIndex(const std::string& Name) const
{
	for (size_t i = 0; i < Fields.size(); ++i)
	{
		if (Fields[i] == Name)
		{
			return static_cast<int>(i);
		}
	}
	return -1;
}

Result Result::Record(const std::vector<std::pair<std::string, double>>& Items)
{
	Result R;
	R.Val.M = Matrix(1, static_cast<int>(Items.size()));
	for (size_t i = 0; i < Items.size(); ++i)
	{
		R.Val.M.V[i] = Items[i].second;
		R.Val.Fields.push_back(Items[i].first);
	}
	return R;
}

int Call::Int(size_t I) const
{
	const double X = N(I);
	if (!std::isfinite(X))
	{
		return 0;
	}
	return static_cast<int>(std::lround(X));
}

std::string FunctionInfo::Signature() const
{
	std::string S = Name + "(";
	bool InOptional = false;
	for (size_t i = 0; i < Params.size(); ++i)
	{
		if (Params[i].HasDefault && !InOptional)
		{
			S += "[";
			InOptional = true;
		}
		if (i > 0)
		{
			S += ", ";
		}
		S += Params[i].Name;
	}
	if (Variadic)
	{
		S += Params.empty() ? "..." : ", ...";
	}
	if (InOptional)
	{
		S += "]";
	}
	return S + ")";
}

size_t FunctionInfo::RequiredCount() const
{
	size_t N = 0;
	for (const Param& P : Params)
	{
		if (!P.HasDefault)
		{
			++N;
		}
	}
	return N;
}

void Registry::Add(const FunctionInfo& Info)
{
	for (FunctionInfo& Existing : Functions)
	{
		if (Existing.Name == Info.Name)
		{
			Duplicated.push_back(Info.Name + " (" + Existing.Category + " and " + Info.Category + ")");
			Existing = Info;
			return;
		}
	}
	Functions.push_back(Info);
}

const FunctionInfo* Registry::Find(const std::string& Name) const
{
	for (const FunctionInfo& F : Functions)
	{
		if (F.Name == Name)
		{
			return &F;
		}
	}
	return nullptr;
}

std::vector<std::string> Registry::Categories() const
{
	std::vector<std::string> Out;
	for (const FunctionInfo& F : Functions)
	{
		if (std::find(Out.begin(), Out.end(), F.Category) == Out.end())
		{
			Out.push_back(F.Category);
		}
	}
	return Out;
}

std::vector<const FunctionInfo*> Registry::InCategory(const std::string& Category) const
{
	std::vector<const FunctionInfo*> Out;
	for (const FunctionInfo& F : Functions)
	{
		if (F.Category == Category)
		{
			Out.push_back(&F);
		}
	}
	return Out;
}

std::vector<std::string> Registry::Completions(const std::string& Prefix) const
{
	std::vector<std::string> Out;
	for (const FunctionInfo& F : Functions)
	{
		if (F.Name.compare(0, Prefix.size(), Prefix) == 0)
		{
			Out.push_back(F.Name);
		}
	}
	std::sort(Out.begin(), Out.end());
	return Out;
}

void Def(Registry& R, const char* Name, const char* Category, const char* Source, const char* Summary, const char* Returns, std::vector<Param> Params, NativeFn Fn)
{
	FunctionInfo Info;
	Info.Name = Name;
	Info.Category = Category;
	Info.Source = Source;
	Info.Summary = Summary;
	Info.Returns = Returns;
	Info.Params = std::move(Params);
	Info.Fn = Fn;
	R.Add(Info);
}

const Registry& Registry::Default()
{
	static const Registry Instance = [] {
		Registry R;
		RegisterCoreMath(R);
		RegisterControlLibrary(R);
		RegisterFlightLibrary(R);
		RegisterRadarLibrary(R);
		RegisterOrbitsLibrary(R);
		RegisterOpsLibrary(R);
		RegisterAtlasLibrary(R);
		RegisterFieldsLibrary(R);
		return R;
	}();
	return Instance;
}
} // namespace mcalc
