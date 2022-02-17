//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Author      : Charles Morcom
//
//   Filename    : ICDS.hpp
//
//   Description : Interface to define a CDS (credit default swap). This is
//                 something that has a fee leg and a protection leg, so that
//                 it is possible for external instruments to calculate things
//                 such as forward spreads, and to make bootstrapping easier.
//
//                 ICDSConvention is something that is capable of generating
//                 an ICDS when you give it a start date, an end date, and a
//                 coupon.
//
//                 ICDSWithConvention is a CDS which is also capable of acting
//                 as a CDS "template" to create others like it. This is useful
//                 for bootstrapping from a "similar" set of benchmarks, or
//                 for defining a vol cube which is valid for a whole family
//                 of CDS.
//
//   Date        : 27 January 2006
//
//----------------------------------------------------------------------------

#ifndef QR_ICDS_HPP
#define QR_ICDS_HPP

#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/ICreditVanillaInstrument.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICDS)
FORWARD_DECLARE(ICDSConvention)
FORWARD_DECLARE_WRAPPER(YieldCurve)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
FORWARD_DECLARE(IDiscountCurveRisky)

/** A CDSConvention is a CDS "template" that lets you generate similar CDS whose
    only differences are their start and end dates and their fee rates. */
class MARKET_DLL ICDSConvention : virtual public IObject {
public:
    /**Generate a new CDS with specified start and end dates and fee whose other
       characteristics are defined by this. "startDate" should be the start of protection;
       "endDate" should be the end of protection.*/
    virtual ICDSSP generateCDS(
        const DateTime&                     startDate,
        const DateTime&                     endDate,
        double feeRate
        ) const = 0;

    static CClassConstSP const TYPE;
    virtual ~ICDSConvention();
private:
    static void load(CClassSP& clazz);
};

/**Interface to define the behaviour expected of something that wants
 * to call itself a credit default swap. Note that this is, implicitly, 
 * a CDS with a constant fee. The convention is that a positive notional implies that
 * you are long protection and short risk. Thus if notional>0, you'd expect that
 * the value given default is positive. This matches the convention in CredDefSwap. */
class MARKET_DLL ICDS : public virtual ICreditVanillaInstrument,
             public virtual IHasCreditFeeLeg,
             public virtual IHasCreditContingentLeg {
public:
    static CClassConstSP const TYPE;
    virtual ~ICDS();

    /**Return the wrapper referring to the yield curve*/
    virtual YieldCurveWrapper getYieldCurveWrapper() const = 0;

    /**Return the wrapper describing the credit curve*/
    virtual ICDSParSpreadsWrapper getParSpreadsWrapper() const = 0;

    /**This has been overridden to set getPV = getFeeLegPV + getContingentLegPV.
     * If you override this, you should still ensure that this is the case!*/
    virtual double getPV(
        const DateTime&              valuationDate,
        const DateTime&              settlementDate,
        const IDiscountCurveRisky&   crv,
        const IDecretionCurveConstSP prepay,
        IForwardRatePricerSP         model,
        IBadDayAdjusterConstSP       bda
        ) const;

    /**This has been overridden to set getPV = getFeeLegPV + getContingentLegPV.
     * If you override this, you should still ensure that this is the case!*/
    virtual double getPV(
        const DateTime&              valuationDate,
        const IDiscountCurveRisky&   crv,
        const IDecretionCurveConstSP prepay,
        IForwardRatePricerSP         model,
        IBadDayAdjusterConstSP       bda
        ) const;

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
