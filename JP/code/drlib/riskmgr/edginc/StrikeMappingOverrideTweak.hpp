//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : StrikeMappingOverrideTweak.hpp
//
//   Description : Tweaks the overridden strike mapping parameter in base correlation models
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_STRIKE_MAPPING_OVERRIDE_TWEAK_H
#define QLIB_STRIKE_MAPPING_OVERRIDE_TWEAK_H

#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

struct StrikeMappingOverrideTwk: public virtual IScalarTweak {
    virtual ~StrikeMappingOverrideTwk() {}
};

DRLIB_END_NAMESPACE

#endif //QLIB_STRIKE_MAPPING_OVERRIDE_TWEAK_H
