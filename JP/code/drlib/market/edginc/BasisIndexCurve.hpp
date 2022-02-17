//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BasisIndexCurve.hpp
//
//   Description : CP, Fed Funds, Prime, PSA, and other index curves
//
//   Author      : Afshin Bayrooti
//
//   Date        : 23 Jan 2006
//
//
//---------------------------------------------------------------------------


#ifndef QLIB_BASIS_INDEX_CURVE_HPP
#define QLIB_BASIS_INDEX_CURVE_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketFactor.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

class CValueDateCollector;
class WrapperNameCollector;
class CModel;

class MARKET_DLL IBasisIndexCurve : public virtual IObject
{
public:
    static CClassConstSP const TYPE;

    virtual ~IBasisIndexCurve();

    /** Returns market name of the object (same as MarketObject::getName())*/ 
    virtual string getName() const = 0;

    /**
     * Get curve base date.
     */
    virtual DateTime getSpotDate() const = 0;

    /**
     * Get the payment date correspond to a reset date by following the
     * basis swap Libor leg convention
     */
    virtual DateTime getRefPaymentDate(const DateTime& resetDate ) const = 0;

    /**
     * Returns spotDate + n + 1 where n = spotOffset
     */
    virtual DateTime getFirstCurveDate() const = 0;

    /**
    * Returns an processed vol - which combines the vol market data with the
    * instrument data in the volRequest.
    */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const = 0;

    /**
    * Returns the basis reference (underlying libor) curve
    */
    virtual YieldCurveWrapper  getRefCurve() const = 0;

    /**
     * Calculate a forward rate for one basis payment period, using the basis 
     * day count convention and simple compounding.
     */
    virtual double fwd(const DateTime& startDate) const = 0;

    /** Calculate a forward rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double fwd(const DateTime&           lodate, 
                       const DateTime&           hidate,
                       const DayCountConvention* dcc,
                       int                       basis) const = 0;

    /**
     * Return par swap rates used to build basis curve.
     */
    virtual DoubleArraySP getParSwapRates() const = 0;

    /**
     * Get par swap rate for requested maturity.
     */
    virtual double parSwapRate(const DateTime& maturity) const = 0;

    /**
     * Get par spread for requested maturity date.
     *
     * Par spread = par swap rate of basis curve - par swap rate of libor curve
     */
    virtual double parSpread(const DateTime& maturity) const = 0;

    /**
    * Get time 0 fwd spread for requested reset and maturity dates.
    *
    * fwd spread = fwd rate of basis curve - fwd rate of libor curve
    */
    virtual double fwdSpread(const DateTime& reset, const DateTime& maturity) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef MarketWrapper<IBasisIndexCurve> IBasisIndexCurveWrapper;
DECLARE(IBasisIndexCurve);




DRLIB_END_NAMESPACE
#endif  // QLIB_BASIS_INDEX_CURVE_HPP




