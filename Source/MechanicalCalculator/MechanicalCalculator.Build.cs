using System.IO;
using UnrealBuildTool;

public class MechanicalCalculator : ModuleRules
{
	public MechanicalCalculator(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;

		PublicDependencyModuleNames.AddRange(new string[]
		{
			"Core",
			"CoreUObject",
			"Engine",
			"InputCore",
			"ProceduralMeshComponent",
		});

		if (Target.bBuildEditor)
		{
			// Building and saving the master materials from code (KeyboardMaterials.cpp).
			PrivateDependencyModuleNames.Add("AssetRegistry");
		}

		// The engine-independent code: the calculator (Calc), the meshes (Geometry) and the
		// keyboard's firmware, switch, sound and backlight models (Keyboard). Tests/ builds the
		// same files natively with CMake.
		PublicIncludePaths.AddRange(new string[]
		{
			Path.Combine(ModuleDirectory, "Calc"),
			Path.Combine(ModuleDirectory, "Geometry"),
			Path.Combine(ModuleDirectory, "Keyboard"),
		});
	}
}
