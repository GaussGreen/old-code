//----------------------------------------------------------------------------
//
//   Group       : Quanitative Research
//
//   Filename    : CreditIndexPreferred.hpp
//
//   Description : Proxy for Credit Index representations
//
//   Author      : Gordon Stephens
//
//   Date        : November 2005
//
//----------------------------------------------------------------------------

#ifndef CREDIT_INDEX_PREFERRED
#define CREDIT_INDEX_PREFERRED

#include "edginc/CreditIndex.hpp"
#include "edginc/CreditIndexMap.hpp"

DRLIB_BEGIN_NAMESPACE

/** CreditIndexPreferred provides a proxy by name for a
    real CreditIndex.
    It is useful for representing rolling index definitions
    in a constant manner. */
class MARKET_DLL CreditIndexPreferred : public CreditIndexBase,
                             public virtual ICreditIndexMap
{
    public:
        static CClassConstSP const TYPE;

        //------------------------
        // CreditIndexBase methods
        //------------------------

        /** Return the index basis adjustment */
        virtual CreditIndexBasisConstSP getIndexBasis() const;

        /** populate from market cache */
        virtual void getMarket(const IModel* model,
                               const MarketData* market);

        /** the name of the market object */
        virtual string getName() const;

        //------------------------
        // ICreditIndexMap methods
        //------------------------

        /** Return the credit index associated to the single name
            as defined by this map */
        virtual CreditIndexBaseConstSP getIndex(const string& singleName) const;

        /** Return the single credit index defined by the map */
        virtual CreditIndexBaseConstSP getIndex() const;

    private:
        CreditIndexPreferred();

        static void load(CClassSP& clazz);
        static IObject* defaultConstructor();

        //fields
        string             name;             //name of this Credit Index
        CreditIndexWrapper creditIndexToUse; //the real index of interest
};

DRLIB_END_NAMESPACE

#endif
