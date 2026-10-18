// The physical layout: 104 keys, unique ids and engine keys, no two caps overlapping.
#include "KeyLayout.h"
#include "TestFramework.h"

#include <set>

TEST(layout_full_size_ansi)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	CHECK(L.Keys.size() == 104);
	std::set<std::string> Ids, Engine;
	for (const mkb::KeyDef& K : L.Keys)
	{
		CHECK(Ids.insert(K.Id).second);
		if (!K.EngineKey.empty())
		{
			CHECK(Engine.insert(K.EngineKey).second);
		}
		CHECK(K.X >= 0 && K.Y >= 0 && K.X + K.W <= L.Width + 1e-9 && K.Y + K.H <= L.Height + 1e-9);
		CHECK(K.ProfileRow >= 1 && K.ProfileRow <= 4);
		if (K.Action == mkb::KeyAction::Insert)
		{
			CHECK(!K.Insert.empty());
		}
	}
	for (size_t i = 0; i < L.Keys.size(); ++i)
	{
		for (size_t j = i + 1; j < L.Keys.size(); ++j)
		{
			const mkb::KeyDef &A = L.Keys[i], &B = L.Keys[j];
			const bool Apart = A.X + A.W <= B.X + 1e-9 || B.X + B.W <= A.X + 1e-9 || A.Y + A.H <= B.Y + 1e-9 || B.Y + B.H <= A.Y + 1e-9;
			if (!Apart)
			{
				std::printf("    overlap %s %s\n", A.Id.c_str(), B.Id.c_str());
			}
			CHECK(Apart);
		}
		// The display window sits clear of every key.
		const mkb::KeyDef& K = L.Keys[i];
		const mkb::Rect& D = L.Display;
		CHECK(K.X + K.W <= D.X || D.X + D.W <= K.X || K.Y + K.H <= D.Y || D.Y + D.H <= K.Y);
	}
	CHECK(L.Find("Space") && L.Find("Space")->NeedsStabilizer());
	CHECK(L.FindByEngineKey("Enter") == L.Find("Enter"));
	CHECK(L.Find("NumEnter")->H == 2);
}
