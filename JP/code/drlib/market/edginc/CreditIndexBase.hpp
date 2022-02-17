//----------------------------------------------------------------------------
//
//   Group       : Quanitative Research
//
//   Filename    : CreditIndexBase.hpp
//
//   Description : Base class for Credit Index representations
//
//   Author      : Gordon Stephens
//
//   Date        : November 2005
//
//----------------------------------------------------------------------------

#ifndef CREDIT_INDEX_BASE_HPP
#define CREDIT_INDEX_BASE_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/CreditIndexBasis.hpp"
#include "edginc/AdjustedCDSPSwithTweaking.hpp"
#include "edginc/AdjustedCDSParSpreads.hpp"

DRLIB_BEGIN_NAMESPACE

class CreditIndexBase;
typedef smartPtr<CreditIndexBase> CreditIndexBaseSP;
typedef smartConstPtr<CreditIndexBase> CreditIndexBaseConstSP;

/** A credit index models the relationship between an index curve and the
    constituent single name curves.
    It should be able to represent and/or calculate the basis spreads
    that result from a mismatch between the single name spreads
    and the index spreads.
*/
class MARKET_DLL CreditIndexBase : public MarketObject
{
    public:
        static CClassConstSP const TYPE;

        /** Return the index basis adjustment */
        virtual CreditIndexBasisConstSP getIndexBasis() const = 0;

        /** Create an index basis adjusted version of a curve
            for convenience */
        //NB should not need to be overridden
        AdjustedCDSPSwithTweakingSP adjustCurve(BootstrappableCDSParSpreadsSP unadjustedCurve) const;
        AdjustedCDSParSpreadsSP     adjustCurveNtwk(BootstrappableCDSParSpreadsSP unadjustedCurve) const;

   protected:
        CreditIndexBase(const CClassConstSP& clazz);

    private:
        static void load(CClassSP& clazz);
};

typedef MarketWrapper<CreditIndexBase> CreditIndexBaseWrapper;

DRLIB_END_NAMESPACE

#endif
