//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : BootstrappedYieldCurve.hpp
//
//   Description : A classical yield curve of cash, futures & swap rates
//
//   Author      : Richard Appleton
//
//   Date        : 25th November 2005
//
//
//----------------------------------------------------------------------------

#ifndef BOOTSTRAPPED_YIELD_CURVE_HPP
#define BOOTSTRAPPED_YIELD_CURVE_HPP

#include "edginc/DeterministicYieldCurve.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/PublicObject.hpp"
#include "edginc/RateShift.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/YieldNameCollector.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/RiskyCurve.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/ZeroPair.hpp"
#include "edginc/Duration.hpp"
#include "edginc/YCWeightedShift.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/IRRatePointwise.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL BootstrappedYieldCurve : public YieldCurve,
                      public virtual IDeterministicYieldCurve,
                      public virtual IPrivateObject,
                      public virtual IGetMarket,
                      public virtual IRestorableWithRespectTo<RateParallel>,
                      public virtual IRestorableWithRespectTo<RatePointwise>,
                      public virtual IRestorableWithRespectTo<IRRatePointwise>,
                      public virtual Theta::Shift,
                      public virtual YCWeightedShift::IShift,
                      public virtual IRiskyCurve,
                      public virtual Duration::IParHandlerWithClosedForm,
                      public virtual Duration::IParHandlerWithoutClosedForm
{
public:
    static CClassConstSP const TYPE;

    BootstrappedYieldCurve(
        const string&                  ccy,
        const string&                  name,
        const DateTime&                today, 
        int                            spotOffset,
        const HolidayWrapper&          holidays,
        const ZeroCurveBenchmarkArray& benchmarks,
        const IZeroCurveFactory&       factory,
        const DayCountConvention&      moneyMarketDcc,
        const MaturityPeriod*          futuresMaturity,
        const IRVolBaseWrapper&        irVol,
        const BadDayConvention*        badDayConvention,
        const DayCountConvention&      fixedDcc,
        const MaturityPeriod&          fixedIvl,
        const DayCountConvention*      floatDcc,
        const MaturityPeriod*          floatIvl,
        const ExpiryArray*             fixDates,  // if not NULL float rate has fixings (floatRateFixed = true)
        const DoubleArray*             fixRates,  // if passed must be same size as fixDates
        const DayCountConvention*      basisDcc,
        const MaturityPeriod*          basisIvl,
        const CurrencyBasisWrapper&    ccyBasis,
        bool                           isIndexCurve = false);

    /** overrides CObject version to allow for easy default */
    virtual bool accept(ICollector* collector) const;
 
	/** return either the growth or discount curve  */
    ZeroCurveSP get(bool growthCurve = true) const;

    /** make a risky curve from a credit spread curve */
    IYieldCurveSP makeRiskyCurve(
        const CreditSpreadCurve& spreadCurve,
        const DateTime*          maturityDate=NULL) const;

    /** make a credit spread curve from default rates */
	CreditSpreadCurveSP makeCreditSpreadCurve(
		const string&        name,
	    const CashFlowArray& defaultRates,
		double               recovery) const;

    IYieldCurveSP createForwardCurve(const DateTime& forwardDate) const;

    BootstrappedYieldCurve* getScaledCurve(const double& scalingFactor) const;

    virtual ~BootstrappedYieldCurve();
    
    /** @return Yield curve's currency */
    virtual string getCcy() const;
    
    /** @return Yield curve's name - used to identify sensitivities */
    virtual string getName() const;

    /** @return Yield curve's spot date */
    virtual DateTime getSpotDate() const;
    
    /** @return Yield curve's today date */
    virtual DateTime getToday() const;

    /** Useful accessor methods */
    virtual ExpiryArrayConstSP getExpiries() const;
    virtual StringArrayConstSP getInstruments() const;

    /** @return curve factory */
    IZeroCurveFactorySP getCurveFactory() const;

    /** @return curve benchmarks */
    const ZeroCurveBenchmarkArray& getBenchmarks() const;

    /** @return money market day count convention */
    DayCountConventionConstSP getMoneyMarketDcc() const;

    /** @return fixed leg day count convention */
    DayCountConventionConstSP getFixedDcc() const;

    /** @return floating leg day count convention */
    DayCountConventionConstSP getFloatDcc() const;

    /** @return basis leg day count convention */
    DayCountConventionConstSP getBasisDcc() const;

    /** @return fixed leg day period */
    MaturityPeriodConstSP getFixedIvl() const;

    /** @return float leg day period */
    MaturityPeriodConstSP getFloatIvl() const;

    /** @return basis leg day period */
    MaturityPeriodConstSP getBasisIvl() const;

    /** @return spot offset */
    int getSpotOffset() const;

    /** @return bad day convention */
    BadDayConventionConstSP getBadDayConvention() const ;

    /** @return holidays */
    const HolidayWrapper& getHolidays() const;

    /** @return fixing dates */
    ExpiryArrayConstSP getFixDates() const;

    /** @return IR vol */
    const IRVolBaseWrapper& getIrVol() const;

    /** strip out the rates and dates */
    virtual CashFlowArraySP getRatesAndDates() const;

    /** Compute discount factor between two dates
     * @param lodate Lower date
     * @param hidate Upper date
     * @return Discount factor between lodate & hidate
     */
    virtual double pv(const DateTime& lodate, const DateTime& hidate) const;

    /** Compute discount factor between value date and a date
     * @param date To get discount factor for
     * @return Discount factor between value date & given date
     */
    virtual double pv(const DateTime& date) const;

    /** Interpolate zero coupon rate at a date
     * @param date Interpolate zero coupon rate here
     * @return Zero coupon rate at given date
     */
    virtual double zero(const DateTime& date) const;   
     
    /** Returns tradeDate + n business days where n = spot offset */
    virtual DateTime settles(const DateTime& tradeDate) const;

    /** Returns rateMaturity->toDate(rateStart) bad day adjusted
        either using bad day convention supplied (if non null)
        together with the yield curve holidays or the yield curve's
        bad day convention and holidays */
    virtual DateTime rateMaturity(
        const DateTime&         rateStart,
        const Expiry*           rateMaturity,
        const BadDayConvention* bdc) const; // optional

    /** Interpolate a forward rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double fwd(const DateTime&           lodate, 
                       const DateTime&           hidate,
                       const DayCountConvention* dcc,
                       int                       basis) const;
                       
    /** Interpolate a forward rate between two dates 
     * isCMS=false => simple rate; isCMS=true =>compomding rate
     */
    virtual double fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMaturity,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const;

    /** Interpolate a forward rate */
    virtual double fwd(const DateTime&           payDate,
                       const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const;

    /** Get par swap rate */
    double parSwapRate(const DateTime& maturity) const;

    /** Clears locally cached zero curve */
    virtual void fieldsUpdated(const CFieldArray& fields);

    virtual IPublicObject* toPublicObject() const;

    /** Override clone method to copy non registered fields over */
    virtual IObject* clone() const;

    /** Same as usual clone() but optionally doesn't copy the local ZC cache */
    BootstrappedYieldCurve* cloneYieldCurve(bool withoutCache) const;

    /** Returns name identifying yield curve for rho parallel */
    virtual string sensName(RateShift* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(RateShift* shift);

    /** Returns name identifying yield curve for additive or 
     * multiplicative weighted shift */
    virtual string sensName(YCWeightedShift* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(YCWeightedShift* shift);

    /** Returns name identifying yield curve for rho parallel */
    virtual string sensName(const RateParallel* shift) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<RateParallel>& shift);
    /** Restores the object to its original form */
    virtual void sensRestore(const PropertyTweak<RateParallel>& shift);

    /** Returns the name of the yield curve - used to determine whether 
        to tweak the object */
    virtual string sensName(const RatePointwise* shift) const;   
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this  yield curve */
    virtual ExpiryWindowArrayConstSP sensQualifiers(const RatePointwise* shift) const;
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual TweakOutcome sensShift(const PropertyTweak<RatePointwise>& shift);
    void sensRestore(const PropertyTweak<RatePointwise>& shift);

    virtual bool sensShift(Theta* shift);

    /** Returns name identifying yield curve for IR delta pointwise */
    virtual string sensName(const IRRatePointwise* shift) const;
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this  yield curve */
    virtual ExpiryWindowArrayConstSP sensQualifiers(const IRRatePointwise* shift) const;
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual TweakOutcome sensShift(const PropertyTweak<IRRatePointwise>& shift);
    void sensRestore(const PropertyTweak<IRRatePointwise>& shift);

    /** Returns true if this and the given yield curve are identical */
    bool zeroCurveEquals(const IYieldCurve* yc2) const;

    /** Optimized hashCode for performance */
    int zeroCurveHash() const;
    
    /** Returns a key used to optimise repeated calculations of discount
        factors (or forward rates). The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative 
        year fraction (Act/365F) betweeen the two dates. */
    virtual IKey* logOfDiscFactorKey() const;

    /** Records name of iso code against yc name in market data object. 
        It is invoked ONCE only
        - immediately after this object is placed in the cache. */
    void initialise(MarketData* market);

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** risky discount factor - the bond spreads are implicit in the yield, so this
        will do the same as pv */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional) const;

    /** grab the dates used in the zero curve */
    virtual DateTimeArray zeroDates() const;
    
    /** drive which style of zero curve is used. Default is discounting */
    virtual void setProjectionCurve(bool useProjectionCurve = true) const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;
    
    //-----------------------------
    //Duration::IParhandler methods
    //-----------------------------

    /** 
     * Return a means of tweaking the MarketObject in a pointwise manner
     * assumption is that the results of this is a VectorResult.....
     */
    virtual VectorRiskPropertySensitivityConstSP getPointwiseTweaker() const;

    /** 
     * Return benchmarks for par instruments
     * these will determine the durations to calculate unless specified
     */
    virtual ExpiryArrayConstSP getParBenchmarks() const;

    /**
     * Return a par instrument for the specified maturity
     */
    virtual InstrumentSP getParInstrument(const ExpirySP maturity) const;

    /**
     * Return closed form solution
     */
    virtual ExpiryResultArraySP getDuration(const Duration* duration) const;

    /**
     * Adjust benchmark array (for forward or risky curve).
     */
    void adjust(
        ZeroCurveBenchmarkArray& adjusted,
        const DateTime&          date,
        const IYieldCurve*       yc) const;

    // IR-QR HACK
    bool isIndexCurve() const;

protected:
    BootstrappedYieldCurve(CClassConstSP clazz = TYPE);
    BootstrappedYieldCurve(const BootstrappedYieldCurve &rhs);
    BootstrappedYieldCurve& operator=(const BootstrappedYieldCurve& rhs);

    static void acceptValueDateCollector(BootstrappedYieldCurve*       yieldCurve, 
                                         CValueDateCollector* collector);
    static void acceptWrapperNameCollector(BootstrappedYieldCurve*       yc,
                                           WrapperNameCollector* collector);

    static void acceptYieldNameCollector(BootstrappedYieldCurve*       yc,
                                         YieldNameCollector* collector);

    /** interface to accessing zero curve - handles caching etc.
        Note - needs modification to support multi-threading */
    class MARKET_DLL ZCAccessor{
    public:
        /** returns cached zero curve - builds on demand if not there */
        const ZeroCurve& get(const BootstrappedYieldCurve& yc) const;
        
        /** Clears locally cached zero curve */
        void reset() const;

        /** Constructor */
        ZCAccessor();
        
    private:
        mutable ZeroPairConstSP zPair;
    };

    friend class ZCAccessor;
    friend class CashSwapCurveGlobalCache;
    friend class IrConverter; // jmf: UGLY!

    // build zero curve for current yield curve
    virtual ZeroPairSP zeroCurve() const;    

    // enforce existence of derived expiry data
    void checkCache();
    // build cache of ordered expiries
    void expiryCache();

    //  fields ///
    string                    ccy;
    string                    name;
    DateTime                  today;
    int                       spotOffset;
    HolidayWrapper            hols;
    IZeroCurveFactorySP       zcMethod;
    ZeroCurveBenchmarkArray   benchmarks;
    DayCountConventionConstSP moneyMarketDayCount;
    IRVolBaseWrapper          irVol;
    BadDayConventionSP        badDayConvention;
    DayCountConventionConstSP fixedDcc;
    MaturityPeriodConstSP     fixedIvl;
    DayCountConventionConstSP floatDcc;
    MaturityPeriodConstSP     floatIvl;
    ExpiryArrayConstSP        fixDates;
    DoubleArrayConstSP        fixRates;
    DayCountConventionConstSP basisDcc;
    MaturityPeriodConstSP     basisIvl;
    CurrencyBasisWrapper      ccyBasis;

    // just for futures
    MaturityPeriodSP          futMaturity;

    // derived data
    DateTime                  valueDate; // $unregistered
    ExpiryArraySP             allExpiries;
    bool                      nosort;

    // how to access zero curve (hides caching)
    ZCAccessor                zc; // $unregistered
    mutable bool              useProjectionCurve;

    // IR-QR HACK - we need 2 curve handles, one for fund and one for disc
    // Only way we can do this is to build 2 curves and set this flag which only
    // we will use
    bool                isIndexCurve_;

protected:
    static const string DEFAULT_FUTURES_MATURITY;

private:
    static void     load(CClassSP& clazz);
    static IObject* defaultBootstrappedYieldCurve();

    /** get first fixing rate */
    double firstFixingRate() const;
};


typedef smartPtr<BootstrappedYieldCurve> BootstrappedYieldCurveSP;
typedef smartConstPtr<BootstrappedYieldCurve> BootstrappedYieldCurveConstSP;


/** Addin for bootstrapping yield curves using benchmark array */
class MARKET_DLL BootstrappedCurveAddin: public CObject, public IPublicObject
{
public:
    static CClassConstSP const TYPE;

    BootstrappedCurveAddin(
        const string&                  ccy,
        const string&                  name,
        const DateTime&                today, 
        int                            spotOffset,
        const HolidayWrapper&          hols,
        const ZeroCurveBenchmarkArray& benchmarks,
        const IZeroCurveFactory&       factory,
        const string&                  moneyMarketDcc,   // may be ""
        const MaturityPeriod*          futuresMaturity,
        const IRVolBaseWrapper&        irVol,
        const string&                  badDayConvention, // may be ""
        const string&                  fixedDcc,
        const MaturityPeriod*          fixedIvl,
        const string&                  floatDcc,         // may be ""
        const MaturityPeriod*          floatIvl,
        const ExpiryArray*             fixDates,
        const DoubleArray*             fixRates,
        const string&                  basisDcc,         // may be ""
        const MaturityPeriod*          basisIvl,
        const CurrencyBasisWrapper&    ccyBasis,
        bool                           isIndexCurve);

    IPrivateObject* toPrivateObject() const;

private:
    static void      load(CClassSP& clazz);
    static IObject*  defaultBootstrappedCurveAddin();

    BootstrappedCurveAddin();

    DayCountConventionSP getDcc(const string& dcc) const;
    BadDayConventionSP   getBdc(const string& bdc) const;

    string                  ccy;
    string                  name;
    DateTime                today; 
    int                     spotOffset;
    HolidayWrapper          hols;
    ZeroCurveBenchmarkArray benchmarks;
    IZeroCurveFactorySP     factory;
    string                  moneyMarketDcc;
    MaturityPeriodSP        futuresMaturity;
    IRVolBaseWrapper        irVol;
    string                  badDayConvention;
    string                  fixedDcc;
    MaturityPeriodSP        fixedIvl;
    string                  floatDcc;
    MaturityPeriodSP        floatIvl;
    ExpiryArraySP           fixDates;
    DoubleArraySP           fixRates;
    string                  basisDcc;
    MaturityPeriodSP        basisIvl;
    CurrencyBasisWrapper    ccyBasis;
    bool                    isIndexCurve;
};


DRLIB_END_NAMESPACE

#endif  // BOOTSTRAPPED_YIELD_CURVE_HPP

