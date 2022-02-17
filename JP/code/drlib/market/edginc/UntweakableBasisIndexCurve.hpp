//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : UntweakableBasisIndexCurve.hpp
//
//   Description : Basis index curve that is already bootstrapped and ready
//                 for use by the diffusion model
//
//   Author      : 
//
//   Date        : 23rd Aug 2006
//
//
//---------------------------------------------------------------------------


#ifndef QLIB_UNTWEAKABLE_BASIS_INDEX_CURVE_HPP
#define QLIB_UNTWEAKABLE_BASIS_INDEX_CURVE_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/BasisIndexCurve.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/SPCalib.hpp"
/*
#include "edginc/YieldCurveAdapter.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/RootFinder.hpp"
*/


DRLIB_BEGIN_NAMESPACE


/**
* Untweakable Basis index curve 
*
* One can only get forward rates off such curves - do not try to get
* discount factors, futures or other such quantities as they are 
* inconsistent with the methodology.
*/
class MARKET_DLL UntweakableBasisIndexCurve : public MarketObject, 
    virtual public IBasisIndexCurve,
    virtual public IMarketFactor,
    virtual public IGetMarket
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    UntweakableBasisIndexCurve();  // for reflection

    // this constructor assumes reference curve is a BootstrappedYieldCurve
    UntweakableBasisIndexCurve(
        const DateTime&     valueDate, 
        const DateTime&     baseDate, 
        IYieldCurveConstSP  basisCurve,
        DoubleArraySP       parSwapRates
    );

    /**
    * Methods inherit from IBasisIndexCurve
    */
    virtual DateTime getSpotDate() const;
    virtual DateTime getFirstCurveDate() const;
    virtual DateTime getRefPaymentDate(const DateTime& resetDate ) const;

    /**
    * Returns an processed vol - which combines the vol market data with the
    * instrument data in the volRequest.
    */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;

    /**
    * Calculate a forward rate for one basis payment period, using the basis 
    * day count convention and simple compounding.
    */
    virtual double fwd(const DateTime& startDate) const;

    /** Calculate a forward rate between two dates */
    virtual double fwd(const DateTime&           lodate, 
        const DateTime&           hidate,
        const DayCountConvention* dcc,
        int                       basis) const;

    /**
    * Return par swap rates used to build basis curve.
    */
    virtual DoubleArraySP getParSwapRates() const;

    /**
    * Get par swap rate for requested maturity.
    */
    virtual double parSwapRate(const DateTime& maturity) const;

    /**
    * Get par spread for requested maturity date.
    *
    * Par spread = par swap rate of basis curve - par swap rate of libor curve
    */
    virtual double parSpread(const DateTime& maturity) const;

    /**
    * Get time 0 fwd spread for requested reset and maturity dates.
    *
    * fwd spread = fwd rate of basis curve - fwd rate of libor curve
    */
    virtual double fwdSpread(const DateTime& reset, const DateTime& maturity) const;

    // accessors for basis curve inputs
    virtual YieldCurveWrapper getRefCurve() const;
    bool getIsAdditive() const;
    DayCountConventionConstSP getBasisSwapLiborDcc() const;
    DayCountConventionConstSP getBasisSwapBasisDcc() const;
    MaturityPeriodConstSP getBasisSwapLiborIvl() const;
    MaturityPeriodConstSP getBasisSwapBasisIvl() const;
    BadDayConventionConstSP getBadDayConvention() const;
    const HolidayWrapper& getHolidays() const;
    /*
    YieldCurveWrapper  getDiscountCurve() const;
    ExpiryArrayConstSP getExpiries() const;
    DoubleArrayConstSP getSpreads() const;
    StringArrayConstSP getInstruments() const;
    */

    /** Optimized equality test for performance */
    bool equalTo(const IObject* obj) const;

    /** Optimized hashCode for performance */
    int hashCode() const;

    /** get market data name */
    string getName() const;

    virtual void getMarket(const IModel* model, const MarketData* market);

    bool isAdditiveBasis() const { return isAdditive; }

    /** overrides CObject version to allow for easy default */
    bool accept(ICollector* collector) const;

    static void acceptValueDateCollector(
        const UntweakableBasisIndexCurve* ic,
        CValueDateCollector*               collector);

    static void acceptWrapperNameCollector(
        const UntweakableBasisIndexCurve* ic,
        WrapperNameCollector*              collector);


private:
    /** throws exception if data is invalid or missing */
    void validate(const string& method) const;

    // inputs
    DateTime                  today;
    DateTime                  valueDate;
    bool                      isAdditive;
    YieldCurveWrapper         reference;    // reference yield curve
    IYieldCurveConstSP        basisCurve;   // basis zero curve
    DoubleArraySP             parSwapRates; // [optional]

    // for curve construction
    string                    name;
    int                       spotOffset;
    SPCalibWrapper            spVol;
    IRVolBaseWrapper          irVol;
    DayCountConventionConstSP basisSwapLiborDcc;
    DayCountConventionConstSP basisSwapBasisDcc;
    MaturityPeriodConstSP     basisSwapLiborIvl;
    MaturityPeriodConstSP     basisSwapBasisIvl;
    BadDayConventionConstSP   badDayConvention;
    HolidayWrapper            holidays;
    
    /*
    ExpiryArrayConstSP        expiries;
    DoubleArrayConstSP        spreads;
    StringArrayConstSP        instruments;
    YieldCurveWrapper         discount;
    DayCountConventionConstSP moneyMarketDcc;
    IZeroCurveFactoryConstSP  zeroCurveFactory;
    DayCountConventionConstSP refSwapFixedDcc;
    DayCountConventionConstSP refSwapFloatDcc;
    MaturityPeriodConstSP     refSwapFixedIvl;
    MaturityPeriodConstSP     refSwapFloatIvl;
    */
};

DECLARE(UntweakableBasisIndexCurve);
typedef MarketWrapper<UntweakableBasisIndexCurve> UntweakableBasisIndexCurveWrapper;


DRLIB_END_NAMESPACE
#endif  // QLIB_UNTWEAKABLE_BASIS_INDEX_CURVE_HPP
