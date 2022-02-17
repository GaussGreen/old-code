//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YieldAsset.hpp
//
//   Description : YieldAsset interface
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_YIELDASSET_HPP
#define EDR_YIELDASSET_HPP
#include "edginc/GeneralAsset.hpp"

DRLIB_BEGIN_NAMESPACE

/** A YieldAsset covers IR type of assets eg 3M LIBOR */
class MARKET_DLL IYieldAsset: public virtual IGeneralAsset{
public:
    static CClassConstSP const TYPE; // in MarketFactor.cpp

    virtual ~IYieldAsset();

    //// perhaps a specialised getProcessedVol method?

private:
    static void load(CClassSP& clazz); // in MarketFactor.cpp
};

DRLIB_END_NAMESPACE
#endif
