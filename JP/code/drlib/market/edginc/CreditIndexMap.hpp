//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CreditIndexMap.hpp
//
//   Description : Abstract base class for mapping single names to an index
//
//   Author      : Gordon C Stephens
//
//   Date        : November 2005
//
//----------------------------------------------------------------------------

#ifndef CREDIT_INDEX_MAP_HPP
#define CREDIT_INDEX_MAP_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CreditIndexBase);

/** A class implementing ICreditIndexMap has the ability to
    associate an index to a single name. */
class MARKET_DLL ICreditIndexMap : public virtual IGetMarket
{
    public:

        static CClassConstSP const TYPE;

        /** Return the credit index associated to the single name
            as defined by this map */
        virtual CreditIndexBaseConstSP getIndex(const string& singleName) const = 0;

        /** Return the single credit index defined by the map */
        // maps that have the ability to represent more than one index
        // should return an empty sp
        virtual CreditIndexBaseConstSP getIndex() const = 0;

    private:

        static void load(CClassSP& clazz);
};

typedef smartPtr<ICreditIndexMap>      ICreditIndexMapSP;
typedef smartConstPtr<ICreditIndexMap> ICreditIndexMapConstSP;
typedef MarketWrapper<ICreditIndexMap> ICreditIndexMapWrapper;

DRLIB_END_NAMESPACE

#endif

