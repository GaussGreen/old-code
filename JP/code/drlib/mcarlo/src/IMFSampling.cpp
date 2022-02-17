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
#include "edginc/IMFSampling.hpp"

DRLIB_BEGIN_NAMESPACE


static void IMFSamplingLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IMFSampling, clazz);
    EXTENDS(IObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet, in particular allow creation of arrays
}


CClassConstSP const IMFSampling::TYPE = 
    CClass::registerInterfaceLoadMethod("IMFSampling", 
                                        typeid(IMFSampling), 
                                        IMFSamplingLoad);

bool IMFSamplingLoad()
{
	return IMFSampling::TYPE != NULL;
};

DRLIB_END_NAMESPACE