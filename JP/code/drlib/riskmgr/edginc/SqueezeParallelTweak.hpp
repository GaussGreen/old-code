//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : SqueezeParallelTweak.hpp
//
//   Description : Multiplicative parallel tweak for squeeze functions
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_SQUEEZE_PARALLEL_TWEAK_H
#define QLIB_SQUEEZE_PARALLEL_TWEAK_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

struct SqueezeParallelTwk: public virtual IScalarTweak {
    virtual ~SqueezeParallelTwk() {}
};

DRLIB_END_NAMESPACE

#endif //QLIB_SQUEEZE_PARALLEL_TWEAK_H
