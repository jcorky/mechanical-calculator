using UnrealBuildTool;

public class MechanicalCalculatorEditorTarget : TargetRules
{
	public MechanicalCalculatorEditorTarget(TargetInfo Target) : base(Target)
	{
		Type = TargetType.Editor;
		DefaultBuildSettings = BuildSettingsVersion.Latest;
		IncludeOrderVersion = EngineIncludeOrderVersion.Latest;
		ExtraModuleNames.Add("MechanicalCalculator");
	}
}
