#include "ICMKernel\\random\icm_RandomRNG-Str.h"
#include "ICMKernel\util\icm_macro.h"


ICM_RandomRNG_Str::ICM_RandomRNG_Str() 
{ 
	Init(1000);
}

ICM_RandomRNG_Str::ICM_RandomRNG_Str(long InitialSeed)
{
	Init( InitialSeed);
}

void ICM_RandomRNG_Str::Init(long InitialSeed)
{
		itsInitialSeed = InitialSeed;
		itsCurrentSeed = InitialSeed;

}


ICM_RandomRNG_Str::ICM_RandomRNG_Str(const ICM_RandomRNG_Str& ref):ICM_RandomRan1(ref)
{
}

ARM_Object * ICM_RandomRNG_Str::Clone(void)
{
	return new ICM_RandomRNG_Str(*this);
}

ICM_RandomRNG_Str& ICM_RandomRNG_Str::operator=(const ICM_RandomRNG_Str& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomRNG_Str(); 
		new(this)ICM_RandomRNG_Str(ref);
	}
	return *this; 
}


ICM_RandomRNG_Str::~ICM_RandomRNG_Str() {}


void ICM_RandomRNG_Str::reset(){
	Init(1000);
}
