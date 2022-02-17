//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : StrikeMappingTweak.hpp
//
//   Description : Strike mapping tweak for base correlation
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_SRIKE_MAPPING_TWEAK_H
#define QLIB_SRIKE_MAPPING_TWEAK_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

struct StrikeMappingTwk: public virtual IScalarTweak {
    virtual ~StrikeMappingTwk() {}
};

DRLIB_END_NAMESPACE

#endif //QLIB_SRIKE_MAPPING_TWEAK_H
