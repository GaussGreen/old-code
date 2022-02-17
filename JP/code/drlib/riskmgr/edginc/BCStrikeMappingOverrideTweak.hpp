//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingOverrideTweak.hpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//                 Override parameter
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_BCSTRIKEMAPPINGTWEAK_HPP
#define QLIB_BCSTRIKEMAPPINGTWEAK_HPP

#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

struct BCStrikeMappingOverrideTwk: public virtual IScalarTweak {
    virtual ~BCStrikeMappingOverrideTwk();

    /** Shifts the original StrikeMappingOverride according to the type
     * of tweak, and returns the new StrikeMappingOverride */
    virtual double applyShift(double& unadjustedValue) = 0;
};

DRLIB_END_NAMESPACE


#endif
