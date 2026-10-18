using UnrealBuildTool;

public class MechanicalCalculatorTarget : TargetRules
{
	public MechanicalCalculatorTarget(TargetInfo Target) : base(Target)
	{
		Type = TargetType.Game;
		DefaultBuildSettings = BuildSettingsVersion.Latest;
		IncludeOrderVersion = EngineIncludeOrderVersion.Latest;
		ExtraModuleNames.Add("MechanicalCalculator");
	}
}
