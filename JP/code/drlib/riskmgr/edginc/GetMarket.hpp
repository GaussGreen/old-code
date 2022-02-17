//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GetMarket.hpp
//
//   Description : What an class must implement to get market data
//                 from the cache
//
//   Author      : Andrew J Swain
//
//   Date        : 24 September 2001
//
//
//----------------------------------------------------------------------------

#ifndef GETMARKET_HPP
#define GETMARKET_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE
class IModel;
class MarketData;

/** What an class must implement to get market data from the cache */
class RISKMGR_DLL IGetMarket: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    virtual ~IGetMarket();

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market) = 0;

protected:
    IGetMarket();
private:
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
