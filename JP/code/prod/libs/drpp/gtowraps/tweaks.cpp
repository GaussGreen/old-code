#include "tweaks.h"
#include "drexception.h"

Tweak toTweak (DRString n)
{
	DRString temp;

	upper_string name (n);
	for (int i = 0; i < NumOfTweaks; i++) {
			temp = TweakStrings[i];
			if (name == temp) return (Tweak) i;
	}
	throw DRException("Could Not Find Tweak. Check spelling.");
	return ALL;
}

DRString toString (Tweak tweak)
{
	DRString ans;
	ans = TweakStrings[(int) tweak];
	return ans;
}
