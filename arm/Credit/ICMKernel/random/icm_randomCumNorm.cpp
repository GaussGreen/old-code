#include "ICMKernel\\random\icm_randomCumNorm.h"
#include "ICMKernel\util\icm_macro.h"


ICM_RandomCumNorm::ICM_RandomCumNorm() 
{ 
	Init();
}

ICM_RandomCumNorm::ICM_RandomCumNorm(const ARM_Vector* pvUnifRand):ICM_RandomLaws(pvUnifRand)
{}

ICM_RandomCumNorm::ICM_RandomCumNorm(const ICM_RandomCumNorm& ref):ICM_RandomLaws(ref)
{}

// randomGen is not cloned !! 
ARM_Object * ICM_RandomCumNorm::Clone(void)
{
	return new ICM_RandomCumNorm(*this);
}

ICM_RandomCumNorm& ICM_RandomCumNorm::operator=(const ICM_RandomCumNorm& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomCumNorm(); 
		new(this)ICM_RandomCumNorm(ref);
	}
	return *this; 
}


ICM_RandomCumNorm::~ICM_RandomCumNorm() 
{
// nothing
}




