//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRNGManager.hpp
//
//   Description : random number manager mainly for MCPathConfigSRMGen
//
//
//----------------------------------------------------------------------------


#ifndef IRNG_MANAGER_HPP

#define IRNG_MANAGER_HPP

#include "edginc/refCountPtr.hpp"
#include "edginc/IRNGGenerator.h" //IUniformRNGGenSP
#include "edginc/IQMCAssetRNG.hpp"

DRLIB_BEGIN_NAMESPACE


class MCARLO_DLL IQMCRNGManager
{
public:

    virtual void seekToPath(size_t path) = 0; ///< prepare all needed RNGs for a path

    /* This is the interface we eventually want assets to use. It gives access to extra source of random numbers */
    virtual IQMCAssetRNGSP  getAssetRNG(size_t assedIdx) = 0; ///< get RNGs needed for the given asset

    /** Legacy interface to facilitate transition until we find how to handle IR/FX relation expressed in calls to fx->begin(...) that implies we currently need to pass rngMgr around */
    virtual const double *  getCorrelatedRandoms(size_t index, size_t factor = 0) = 0; ///< return pointer to correlated numbers for this index; currently index should be randomIndex stored inside asset.

    virtual IUniformRNGGenSP    getSharedGen() = 0; ///< returns extra generator of uniform numbers

    virtual ~IQMCRNGManager() {} // for SP

};
DECLARE_REF_COUNT(IQMCRNGManager);

DRLIB_END_NAMESPACE

#endif

