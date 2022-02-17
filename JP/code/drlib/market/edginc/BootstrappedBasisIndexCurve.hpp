//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BootstrappedBasisIndexCurve.hpp
//
//   Description : Basis index curve (as per ALIB GtoZCBasis).
//
//   Author      : Richard Appleton
//
//   Date        : 21st Feb 2006
//
//
//---------------------------------------------------------------------------


#ifndef QLIB_BOOTSTRAPPED_BASIS_INDEX_CURVE_HPP
#define QLIB_BOOTSTRAPPED_BASIS_INDEX_CURVE_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/BasisIndexCurve.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/YieldCurveAdapter.hpp"
#include "edginc/SPCalib.hpp"



DRLIB_BEGIN_NAMESPACE

struct BootstrappedBasisIndexCurveData;


/**
 * Basis index curve (cf. ALIB GtoZCBasis).
 *
 * One can only get forward rates off such curves - do not try to get
 * discount factors, futures or other such quantities as they are 
 * inconsistent with the methodology.
 */
class MARKET_DLL BootstrappedBasisIndexCurve : public MarketObject, 
                                    virtual public IBasisIndexCurve,
                                    virtual public IMarketFactor,
                                    virtual public IGetMarket
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    BootstrappedBasisIndexCurve();  // for reflection

    // this constructor assumes reference curve is a BootstrappedYieldCurve
    BootstrappedBasisIndexCurve(
        const string&               name,
        const YieldCurveWrapper&    discount,
        const YieldCurveWrapper&    reference,
        const ExpiryArray&          expiries,
        const DoubleArray&          spreads,
        const StringArray&          instruments,
        const bool                  isAdditive,
        const DayCountConvention&   basisSwapLiborDcc,
        const DayCountConvention&   basisSwapBasisDcc,
        const MaturityPeriod&       basisSwapLiborIvl,
        const MaturityPeriod&       basisSwapBasisIvl,
        const SPCalibWrapper&       spVol);

    // this constructor makes no assumptions re type of input zero curve
    BootstrappedBasisIndexCurve(
        const DateTime&             today,
        const string&               name,
        const YieldCurveWrapper&    discount,
        const YieldCurveWrapper&    reference,
        const ExpiryArray&          expiries,
        const DoubleArray&          spreads,
        const StringArray&          instruments,
        const bool                  isAdditive,

        // following parameters are for curve construction
        const int                   spotOffset,
        const IZeroCurveFactory&    zeroCurveFactory,
        const DayCountConvention&   moneyMarketDcc,
        const DayCountConvention&   refSwapFixedDcc,
        const DayCountConvention&   refSwapFloatDcc,
        const MaturityPeriod&       refSwapFixedIvl,
        const MaturityPeriod&       refSwapFloatIvl,
        const DayCountConvention&   basisSwapLiborDcc,
        const DayCountConvention&   basisSwapBasisDcc,
        const MaturityPeriod&       basisSwapLiborIvl,
        const MaturityPeriod&       basisSwapBasisIvl,
        const BadDayConvention&     swapBadDayConvention,
        const IRVolBaseWrapper&     irVol,
        const SPCalibWrapper&       spVol,
        const HolidayWrapper&       holidays);


    /**
     * Get curve base date.
     */
    virtual DateTime getSpotDate() const;
    virtual DateTime getFirstCurveDate() const;
    virtual DateTime getRefPaymentDate(const DateTime& resetDate ) const;

    /**
    * Returns an processed vol - which combines the vol market data with the
    * instrument data in the volRequest.
    */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;

    // accessors for basis curve inputs
    YieldCurveWrapper  getDiscountCurve() const;
    virtual YieldCurveWrapper  getRefCurve() const;
    virtual ExpiryArrayConstSP getExpiries() const;
    DoubleArrayConstSP getSpreads() const;
    virtual StringArrayConstSP getInstruments() const;
    bool               getIsAdditive() const;
    DayCountConventionConstSP getBasisSwapLiborDcc() const;
    DayCountConventionConstSP getBasisSwapBasisDcc() const;
    MaturityPeriodConstSP getBasisSwapLiborIvl() const;
    MaturityPeriodConstSP getBasisSwapBasisIvl() const;
    BadDayConventionConstSP getBadDayConvention() const;
    const HolidayWrapper& getHolidays() const;

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

    /** Optimized equality test for performance */
    bool equalTo(const IObject* obj) const;

    /** Optimized hashCode for performance */
    int hashCode() const;

    /** get market data name */
    string getName() const;

    virtual void getMarket(const IModel* model, const MarketData* market);

    /** overrides CObject version to allow for easy default */
    bool accept(ICollector* collector) const;

    static void acceptValueDateCollector(
        const BootstrappedBasisIndexCurve* ic,
        CValueDateCollector*               collector);

    static void acceptWrapperNameCollector(
        const BootstrappedBasisIndexCurve* ic,
        WrapperNameCollector*              collector);

    /** Build basis index curve */
    BootstrappedBasisIndexCurveData bootstrap() const;

private:
    /** access cached data, bootstrapping curve if required */
    BootstrappedBasisIndexCurveData get() const;

    /** throws exception if data is invalid or missing */
    void validate(const string& method) const;

    // inputs
    DateTime                  today;
    ExpiryArrayConstSP        expiries;
    DoubleArrayConstSP        spreads;
    StringArrayConstSP        instruments;
    bool                      isAdditive;
    YieldCurveWrapper         discount;
    YieldCurveWrapper         reference;

    // for curve construction
    string                    name;
    int                       spotOffset;
    IZeroCurveFactoryConstSP  zeroCurveFactory;
    DayCountConventionConstSP moneyMarketDcc;
    DayCountConventionConstSP refSwapFixedDcc;
    DayCountConventionConstSP refSwapFloatDcc;
    MaturityPeriodConstSP     refSwapFixedIvl;
    MaturityPeriodConstSP     refSwapFloatIvl;
    DayCountConventionConstSP basisSwapLiborDcc;
    DayCountConventionConstSP basisSwapBasisDcc;
    MaturityPeriodConstSP     basisSwapLiborIvl;
    MaturityPeriodConstSP     basisSwapBasisIvl;
    BadDayConventionConstSP   badDayConvention;
    IRVolBaseWrapper          irVol;
    SPCalibWrapper            spVol;
    HolidayWrapper            holidays;
};

DECLARE(BootstrappedBasisIndexCurve);
typedef MarketWrapper<BootstrappedBasisIndexCurve> BootstrappedBasisIndexCurveWrapper;

/*
 * Private object used for boostrapping basis curve.
 */
class MARKET_DLL BasisCurveFuncObj : public Func1D::NoDeriv
{
public:
    BasisCurveFuncObj(
        const string&               ccy,
        const string&               name,
        const DateTime&             today,
        const int                   spotOffset,
        const IYieldCurve&          discount,
        const IYieldCurve&          reference,
        const ExpiryArray&          expiries,
        const DoubleArray&          spreads,
        const StringArray&          instruments,
        const bool                  isAdditive,
        const IZeroCurveFactory&    factory,
        const DayCountConvention&   moneyMarketDcc,
        const DayCountConvention&   refSwapFixedDcc,
        const DayCountConvention&   refSwapFloatDcc,
        const MaturityPeriod&       refSwapFixedIvl,
        const MaturityPeriod&       refSwapFloatIvl,
        const DayCountConvention&   basisSwapLiborDcc,
        const DayCountConvention&   basisSwapBasisDcc,
        const MaturityPeriod&       basisSwapLiborIvl,
        const MaturityPeriod&       basisSwapBasisIvl,
        const BadDayConvention&     swapBadDayConvention,
        const IRVolBaseWrapper&     irVol,
        const SPCalibWrapper&       spVol,
        const HolidayWrapper&       holidays);

    void calculate();

    double operator()(double  x) const;

    IYieldCurve*  getZeroCurve() const;
    DoubleArraySP getParSwapRates() const;

private:
    double calcBasisLegPV(
        const DateTime&  swapMatDate, 
        const DateTime&  swapBasisPeriod,
        const ZeroCurve& basisCurve) const;

    void computeBasisForwardRates(
        const ZeroCurve& basisCurve, 
        DoubleArray&     basisPV) const;

    double calcLiborLegPV(
        const DateTime& swapMatDate, 
        const DateTime& swapLiborPeriod) const;

    void initSwaps();
    void initLiborLeg();
    void initBasisLeg();

    int dateArrayFloorIndex(
        const DateTimeArray& y,
        const DateTime&      x) const;


    // inputs
    const DateTime              valueDate;
    const ExpiryArray&          expiries;
    const DoubleArray&          spreads;
    const StringArray&          instruments;
    const bool                  isAdditive;
    const YieldCurveAdapter     discount;
    const YieldCurveAdapter     reference;

    // for curve construction
    const string                ccy;
    const string                name;
    const int                   spotOffset;
    const IZeroCurveFactory&    factory;
    ZeroCurveBenchmarkArray     benchmarks;
    const DayCountConvention&   mmDcc;
    const DayCountConvention&   refSwapFixedDcc;
    const DayCountConvention&   refSwapFloatDcc;
    const DayCountConvention&   basisSwapLiborDcc;
    const DayCountConvention&   basisSwapBasisDcc;
    const MaturityPeriod&       refSwapFixedIvl;
    const MaturityPeriod&       refSwapFloatIvl;
    const MaturityPeriod&       basisSwapLiborIvl;
    const MaturityPeriod&       basisSwapBasisIvl;
    const BadDayConvention&     badDayConv;
    const IRVolBaseWrapper&     irVol;
    const SPCalibWrapper&       spVol;
    const HolidayWrapper&       holidays;

    // for calculation
    DateTimeArray                  basisResetDates;
    DateTimeArray                  liborResetDates;
    DateTimeArray                  adjBasisResetDates;
    DateTimeArray                  adjLiborResetDates;
    DateTimeArray                  adjSwapDates;
    DoubleArray                    discBasisResetDateZeros; // discounts at basis resets
    DoubleArray                    discLiborResetDateZeros; // discounts at libor resets
    DoubleArray                    discSwapDateZeros;       // discounts at swap dates
    DoubleArray                    guessSwapRates;          // initial guess of basis par rates
    DoubleArray                    liborPV;                 // pv of libor leg
    DoubleArray                    basisAnn;                // annuity facts at basis resets
    DoubleArray                    yearFracLibor;           // year fractions of basis intervals
    DoubleArray                    yearFracBasis;           // year fractions of libor intervals
    int                            n;                       // index of instrument currently being fitted [0 based]
};



// add-ins for bootstrapping basis index curve

class MARKET_DLL BootstrappedBasisIndexCurveAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

    static void      load(CClassSP& clazz);
    static IObject*  defaultBootstrappedBasisIndexCurveAddin();

    BootstrappedBasisIndexCurveAddin();

    IObjectSP createBasisIndexCurve();

    YieldCurveSP            refCurve;
    YieldCurveSP            discountCurve;
    ExpiryArray             expiries;
    DoubleArray             spreads;
    StringArray             instruments;
    string                  basisSwapType;

    // following parameters are for curve construction
    string                  name;
    DayCountConventionSP    basisSwapLiborDcc;
    MaturityPeriodSP        basisSwapLiborIvl;
    DayCountConventionSP    basisSwapBasisDcc;
    MaturityPeriodSP        basisSwapBasisIvl;

    SPCalibWrapper          spVol;
    MarketDataConstSP       marketData;
};


class BootstrappedBasisIndexCurveLegacyAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

    static void      load(CClassSP& clazz);
    static IObject*  defaultBootstrappedBasisIndexCurveLegacyAddin();

    BootstrappedBasisIndexCurveLegacyAddin();

    IObjectSP createBasisIndexCurve();

    YieldCurveSP            refCurve;
    YieldCurveSP            discountCurve;
    ExpiryArray             expiries;
    DoubleArray             spreads;
    StringArray             instruments;
    string                  basisSwapType;

    // following parameters are for curve construction
    string                  name;
    DayCountConventionSP    basisSwapLiborDcc;
    MaturityPeriodSP        basisSwapLiborIvl;
    DayCountConventionSP    basisSwapBasisDcc;
    MaturityPeriodSP        basisSwapBasisIvl;

    DateTime                today;
    int                     spotOffset;
    IZeroCurveFactorySP     zeroCurveFactory;
    HolidayWrapper          holidays;
    DayCountConventionSP    moneyMarketDcc;
    IRVolBaseWrapper        irVol;
    string                  badDayConvention;
    DayCountConventionSP    refSwapFixedDcc;
    MaturityPeriodSP        refSwapFixedIvl;
    DayCountConventionSP    refSwapFloatDcc;
    MaturityPeriodSP        refSwapFloatIvl;

    SPCalibWrapper          spVol;
    MarketDataConstSP       marketData;
};


DRLIB_END_NAMESPACE
#endif  // QLIB_BOOTSTRAPPED_BASIS_INDEX_CURVE_HPP
