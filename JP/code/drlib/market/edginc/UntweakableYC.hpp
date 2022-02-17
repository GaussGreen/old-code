#ifndef UNTWEAKABLEYC_HPP
#define UNTWEAKABLEYC_HPP

#include "edginc/DeterministicYieldCurve.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/RateShift.hpp"
#include "edginc/YieldNameCollector.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/IRRatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL UntweakableYC: public YieldCurve,
                     public virtual IDeterministicYieldCurve,
                     public virtual ITweakableWithRespectTo<RateParallel>,
                     public virtual ITweakableWithRespectTo<RatePointwise>,
                     public virtual ITweakableWithRespectTo<IRRatePointwise>,
                     public virtual Theta::Shift,
                     public virtual YCWeightedShift::IShift
{
    friend class AvoidCompilerWarning;
    friend class IrConverter;

protected:
    /************* exported fields *************/
    string                  ccy;
    string                  name;
    DateTime                today;
    int                     spotOffset;
    DayCountConventionSP    moneyMarketDayCount;
    DayCountConventionSP    swapDayCount;
    int                     swapFrequency;
    ZeroCurveConstSP        zeroCurve;       ///< default [discount]
    ZeroCurveConstSP        growthZeroCurve; ///< used if setProjectionCurve()
    DateTime                valueDate;
    IRVolBaseWrapper        irVol;

    /************* transient fields *************/
    mutable bool            useDiscCurve;    ///< defaults to true

    /************* end of fields *************/

public:
    static CClassConstSP const TYPE;

    UntweakableYC(
        const DateTime&  valueDate, 
        const DateTime&  baseDate, 
        const ZeroCurve& discounting, 
        const ZeroCurve* estimating = NULL);

    virtual ~UntweakableYC(void){}

    /** overrides CObject version to allow for easy default */
    bool accept(ICollector* collector) const;

    /** Records name of iso code against yc name in market data object.
        It is invoked ONCE only
        - immediately after this object is placed in the cache. */
    void initialise(MarketData* market);

    /** populate from market cache */
    void getMarket(const IModel* model, const MarketData* market);

    /** @return Yield curve's currency */
    virtual string getCcy(void) const;

    /** @return Yield curve's name - used to identify sensitivities */
    virtual string getName(void) const;

    /** @return Yield curve's spot date */
    virtual DateTime getSpotDate(void) const;

    /** @return Yield curve's today date */
    virtual DateTime getToday(void) const;

    /** Useful accessor methods */
    ExpiryArrayConstSP getExpiries() const;
    StringArrayConstSP getInstruments() const;

    /** Returns tradeDate + n business days where n = spot offset */
    virtual DateTime settles(const DateTime& tradeDate) const;

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

    virtual double parSwapRate(const DateTime& maturity) const;

    /** Convenience method to return the curve dates */
    virtual DateTimeArraySP getDates() const;

    /** Interpolate zero coupon rate at a date
     * @param date Interpolate zero coupon rate here
     * @return Zero coupon rate at given date
     */
    virtual double zero(const DateTime& date) const;

    /** Returns rateMaturity->toDate(rateStart)  */
    virtual DateTime rateMaturity(
        const DateTime&         rateStart,
        const Expiry*           rateMaturity,
        const BadDayConvention* bdc) const;

    /** Interpolate a forward rate between two dates
     * see CompoundBasis for basis values
     */
    virtual double fwd(const DateTime&           lodate,
                       const DateTime&           hidate,
                       const DayCountConvention* dcc,
                       int                       basis) const;

    virtual double fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const;

    virtual double fwd(const DateTime&           payDate,
                       const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const;

     //// switch to using growth curve
    virtual void setProjectionCurve(bool useEstimatingCurve = true) const;

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

    /** Returns the name of the yield curve - used to determine whether
    to tweak the object */
    virtual string sensName(const IRRatePointwise* shift) const;

    /** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
    virtual ExpiryWindowArrayConstSP sensQualifiers(const IRRatePointwise* shift) const;

    /** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
    virtual TweakOutcome sensShift(const PropertyTweak<IRRatePointwise>& shift);

    virtual bool sensShift(Theta* shift);

    /** Returns a key used to optimise repeated calculations of discount
        factors (or forward rates). The calc method for this key returns the
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates. */
    virtual IKey* logOfDiscFactorKey(void) const;

    /** grab the dates used in the zero curve */
    virtual DateTimeArray zeroDates(void) const;

    /** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
    CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;

    // these to move to RiskyCurve
    /** this function calculates the discount factor based on the assumption
        that the on which the recovery is based is provided externally.
        This allows to use different methodologies
        (PV, face value + accrued etc.) to be included easily */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional) const;

    /** make a risky curve from a credit spread curve */
    virtual IYieldCurveSP makeRiskyCurve(
        const CreditSpreadCurve& spreadCurve,
        const  DateTime*    maturityDate=NULL) const;

	virtual CreditSpreadCurveSP makeCreditSpreadCurve(
		const string&        name,
	    const CashFlowArray& defaultRates,
		double               recovery) const;

    /** grab the dates used in the zero curve and corresponding rates */
    virtual CashFlowArraySP getRatesAndDates() const;

    virtual IYieldCurveSP createForwardCurve(const DateTime& forwardDate) const;

    // checks for any negative forward rate implied by the zero curve
    bool validate();

private:
    UntweakableYC(void);
    UntweakableYC(const UntweakableYC &rhs);
    UntweakableYC& operator=(const UntweakableYC& rhs);

    static IObject* defaultConstructor(void);

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static void acceptValueDateCollector(const UntweakableYC*  yieldCurve,
                                         CValueDateCollector* collector);
    static void acceptWrapperNameCollector(const UntweakableYC*       yc,
                                           WrapperNameCollector* collector);
    static void acceptYieldNameCollector(const UntweakableYC*   yc,
                                         YieldNameCollector* collector);
};

DRLIB_END_NAMESPACE

#endif
