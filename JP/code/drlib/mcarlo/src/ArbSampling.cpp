//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids DR
//
//   Filename    : AdaptiveSampling.cpp
//
//   Description : 
//
//   Date        : Oct 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ArbSampling.hpp"

DRLIB_BEGIN_NAMESPACE



IObject* ArbSampling::defaultArbSampling()
{
	return new ArbSampling();
};

CClassConstSP const ArbSampling::TYPE =
	CClass::registerClassLoadMethod(
		"ArbSampling", 
		typeid(ArbSampling), 
		load);

void
ArbSampling::load(CClassSP& clazz)
{
    clazz->setPublic();
    
	REGISTER(ArbSampling, clazz);

	SUPERCLASS(CObject);

	IMPLEMENTS(IMFSampling);
	
	FIELD(points, "Points on the grid");

	FIELD(samples, "Number of samples");

	EMPTY_SHELL_METHOD(defaultArbSampling);
}

void 
ArbSampling::validatePop2Object()
{
	
};

bool ArbSamplingLoad()
{
	return ArbSampling::TYPE != NULL;
};


DRLIB_END_NAMESPACE