//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : WeightedInstrumentTweak.hpp
//
//   Description : Tag class for shift in instrument weight 
//                 (WeightedInstrumentCollection)
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#ifndef DRLIB_WEIGHTED_INSTRUMENT_TWEAK_H
#define DRLIB_WEIGHTED_INSTRUMENT_TWEAK_H

#include "edginc/Void.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct RISKMGR_DLL WeightedInstrumentTweak: CObject {
    static CClassConstSP const TYPE;
    WeightedInstrumentTweak(IntArrayConstSP instruments = IntArrayConstSP());
    ~WeightedInstrumentTweak();

    IntArrayConstSP instruments;

    typedef Void Qualifier;
    //Is continuous, i.e. tweaks to it can be made arbitrarily small
    enum { discrete = 0};
};

DRLIB_END_NAMESPACE

#endif
