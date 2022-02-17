
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IQMCAssetRNG.hpp
//
//   Description : During diffusion assets need correlated numbers and possibly extra sources of randomness
//                 This interface encapsulates this functionality.
//
//
//----------------------------------------------------------------------------


#ifndef IQMCASSETRNG
#define IQMCASSETRNG

#include "edginc/IRNGGenerator.h"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class MCARLO_DLL IQMCAssetRNG
{
    public:
        virtual const double * getCorrelatedRandoms(size_t factor = 0) = 0; // get random numbers for diffusion that are correlated between assets
        virtual IUniformRNGGenSP  getExtra() = 0;      // get uniform numbers we can use as we pleased
        virtual ~IQMCAssetRNG() {}
};

DECLARE_REF_COUNT(IQMCAssetRNG);

DRLIB_END_NAMESPACE

#endif
