//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeterministicYieldCurve.hpp
//
//   Description : Basically a YieldCurve with no Vol
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_DETERMINISTICYIELDCURVE_HPP
#define EDR_DETERMINISTICYIELDCURVE_HPP
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/** Basically a traditional YieldCurve with Vol */
class MARKET_DLL IDeterministicYieldCurve: public virtual IYieldCurve{
public:
    static CClassConstSP const TYPE; // in YieldCurve.cpp

    virtual ~IDeterministicYieldCurve();

    // just a marker interface at the moment
private:
    static void load(CClassSP& clazz); // in YieldCurve.cpp
};

// support for wrapper class
typedef MarketWrapper<IDeterministicYieldCurve> IDeterministicYieldCurveWrapper;

DRLIB_END_NAMESPACE
#endif
