//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CompressionRatioTweak.hpp
//
//   Description : Tweak for compression ratio
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_COMPRESSION_RATIO_TWEAK_H
#define QLIB_COMPRESSION_RATIO_TWEAK_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

struct CompressionRatioTwk: public virtual IScalarTweak {
    virtual ~CompressionRatioTwk() {}
};

DRLIB_END_NAMESPACE

#endif //QLIB_COMPRESSION_RATIO_TWEAK_H
