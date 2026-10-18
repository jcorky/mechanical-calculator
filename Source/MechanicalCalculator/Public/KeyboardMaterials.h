// Materials without .uasset files. The project ships code only, so the two master materials, one
// opaque and one translucent, both with BaseColor / Metallic / Roughness / Specular / Emissive /
// VertexEmissive / Opacity parameters, come from the first of:
//   1. saved assets in /Game/MechanicalCalculator/Materials/ (made by the console command
//      `mc.SaveMaterials` in the editor; needed for packaged builds);
//   2. in an editor build, masters built in memory from material expressions;
//   3. the engine's BasicShapeMaterial, tinted through its Color parameter (no metal or glow).
#pragma once

#include "CoreMinimal.h"

class UMaterialInstanceDynamic;
class UMaterialInterface;

struct FKeyboardSurface
{
	FLinearColor BaseColor = FLinearColor(0.8f, 0.8f, 0.8f);
	float Metallic = 0.0f;
	float Roughness = 0.5f;
	float Specular = 0.5f;
	FLinearColor Emissive = FLinearColor::Black;
	float VertexEmissive = 0.0f;   // multiplies the vertex colour into emission (the RGB backlight)
	float Opacity = 1.0f;
	bool bTranslucent = false;
};

namespace KeyboardMaterials
{
UMaterialInstanceDynamic* Make(UObject* Outer, const FKeyboardSurface& Surface);
void Apply(UMaterialInstanceDynamic* Material, const FKeyboardSurface& Surface);
// True when the masters support metal, glow and translucency (cases 1 and 2 above).
bool HasFullMasters();
} // namespace KeyboardMaterials
