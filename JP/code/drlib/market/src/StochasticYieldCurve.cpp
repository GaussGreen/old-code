//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StochasticYieldCurve.hpp
//
//   Description : Basically a traditional YieldCurve with Vol
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/StochasticYieldCurve.hpp"
#include "edginc/DeterministicYieldCurve.hpp"
#include "edginc/IRVolBase.hpp"

DRLIB_BEGIN_NAMESPACE

////// IStochasticYieldCurve /////////
IStochasticYieldCurve::~IStochasticYieldCurve(){}
void IStochasticYieldCurve::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IStochasticYieldCurve, clazz);
    EXTENDS(IYieldCurve);
}

CClassConstSP const IStochasticYieldCurve::TYPE =
CClass::registerInterfaceLoadMethod(
    "IStochasticYieldCurve", typeid(IStochasticYieldCurve), load);

/** Obvious implementation */
class StochasticYieldCurve: public YieldCurve,
                            public virtual IStochasticYieldCurve{
    IDeterministicYieldCurveWrapper yc; 
    IRVolBaseWrapper                irVol;
public:
    static CClassConstSP const TYPE; 

    virtual ~StochasticYieldCurve(){};

    virtual void getMarket(const IModel* model, const MarketData* market){
        yc.getData(model, market);
        irVol.getData(model, market);
    }

    /** drive which style of zero curve is used. Default is discounting */
    virtual void setProjectionCurve(bool useEstimatingCurve) const{
        yc->setProjectionCurve(useEstimatingCurve); // pass down
    }

    /** overrides CObject version to allow for easy default */
    bool accept(ICollector* collector) const{
        if (!CClass::invokeAcceptMethod(this, collector)){
            // if no method registered try yc
            if (!CClass::invokeAcceptMethod(yc.get(), collector)){
                // then try vol
                return irVol->accept(collector);
           }
        }
        return false;
    }

    /** @return Yield curve's currency */
    virtual string getCcy() const{
        return yc->getCcy();
    }
    
    /** @return Yield curve's name - used to identify sensitivities */
    virtual string getName() const{
        return yc.getName();
    }

    /** @return Yield curve's spot date */
    virtual DateTime getSpotDate() const{
        return yc->getSpotDate();
    }

    /** @return Yield curve's today date */
    virtual DateTime getToday() const {
        return yc->getToday();
    }
    
    /** Useful accessor methods */
    virtual ExpiryArrayConstSP getExpiries() const {
        return yc->getExpiries();
    }

    virtual StringArrayConstSP getInstruments() const {
        return yc->getInstruments();
    }

    /** Passes down to non stochastic YC */
    virtual DateTime settles(const DateTime& tradeDate) const{
        return yc->settles(tradeDate);
    }

    /** Compute discount factor between two dates
     * @param lodate Lower date
     * @param hidate Upper date
     * @return Discount factor between lodate & hidate
     */
    virtual double pv(const DateTime& lodate, 
                      const DateTime& hidate) const{
        return yc->pv(lodate, hidate);
    }
    
    /** Compute discount factor between value date and a date
     * @param date To get discount factor for
     * @return Discount factor between value date & given date
     */
    virtual double pv(const DateTime& date) const{
        return yc->pv(date);
    }
        
    /** Calculates present value to baseDate of supplied cash flows. Cash
        flows on or before baseDate (eg value date) are ignored. No
        ordering of the cashflows is assumed */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const{
        return yc->pv(cashFlows, baseDate);
    }
        
    /** Interpolate zero coupon rate at a date
     * @param date Interpolate zero coupon rate here
     * @return Zero coupon rate at given date
     */
    virtual double zero(const DateTime& date) const{ 
        return yc->zero(date);
    }

    /** Passes down to non stochastic yc */
    virtual DateTime rateMaturity(
        const DateTime&         rateStart,
        const Expiry*           rateMaturity,
        const BadDayConvention* bdc) const{ // optional
        return yc->rateMaturity(rateStart, rateMaturity, bdc);
    }

    /** Interpolate a forward rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double fwd(const DateTime&           lodate, 
                       const DateTime&           hidate,
                       const DayCountConvention* dcc,
                       int                       basis) const{
        return yc->fwd(lodate, hidate, dcc, basis);
    }

    /** Calculates a par swap rate for the swap starting at startDate,
     * maturing at maturityDate, with fixed day count convention 
     * fixedDayCountConv and fixed payments occuring at time 
     * intervals defined by interval. In other words, the routine calculates 
     * the fixed rate such that the present value of the fixed and 
     * floating sides are equal. The floating side is assumed to be at par.
     */
    virtual double couponRate(
        const DateTime&           startDate,   //(I) Date instrument begins at
        const DateTime&           maturityDate,//(I) Date instrument matures at
        const MaturityPeriod&     interval,    // (I) Time between payments 
        bool                      stubAtEnd,
        const DayCountConvention* dcc) const{
        return yc->couponRate(startDate, maturityDate, 
                              interval, stubAtEnd, dcc);
    }

    /** Interpolate a futures rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double future(const DateTime&           lodate, 
                          const DateTime&           hidate,
                          const DayCountConvention* dcc,
                          int                       basis,
                          double                    irVol) const{
        return yc->future(lodate, hidate, dcc, basis, irVol);
    }

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IKey* logOfDiscFactorKey() const{
        return yc->logOfDiscFactorKey();
    }

    /** grab the dates used in the zero curve */
    virtual DateTimeArray zeroDates() const{
        return yc->zeroDates();
    }

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest) const{
        return irVol->getProcessedVol(volRequest, this);
    }

    // these to move to RiskyCurve
    /** this function calculates the discount factor based on the assumption
        that the on which the recovery is based is provided externally. 
        This allows to use different methodologies 
        (PV, face value + accrued etc.) to be included easily */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional) const{
        return pv(lodate, hidate) * cashFlow;
    }


    double fwd(const DateTime&           refixDate, 
                           const Expiry*             rateMat,
                           const BadDayConvention*   bdc, // optional
                           const DayCountConvention* dcc,
                           const bool                isCMS) const{
        throw ModelException("StochasticYieldCurve::fwd","Not implemented");
    }

    double fwd(const DateTime&           payDate,
                           const DateTime&           refixDate, 
                           const Expiry*             rateMat,
                           const BadDayConvention*   bdc, // optional
                           const DayCountConvention* dcc,
                           const bool                isCMS) const{
        throw ModelException("StochasticYieldCurve::fwd","Not implemented");
    }

    /** make a risky curve from a credit spread curve */
    virtual IYieldCurveSP makeRiskyCurve(
        const CreditSpreadCurve& spreadCurve,
        const  DateTime*    maturityDate=NULL) const
    {
        throw ModelException("StochasticYieldCurve::makeRiskyCurve","Not implemented");
    }


    virtual CreditSpreadCurveSP makeCreditSpreadCurve(
	    const string&        name,
	    const CashFlowArray& defaultRates,
		double               recovery) const
    {
        throw ModelException("StochasticYieldCurve::makeCreditSpreadCurve","Not implemented");
    }


    virtual CashFlowArraySP getRatesAndDates() const
    {
        return yc->getRatesAndDates();
    }


    virtual IYieldCurveSP createForwardCurve(const DateTime& forwardDate) const
    {
        throw ModelException("StochasticYieldCurve::createForwardCurve","Not implemented");
    }

private:
    StochasticYieldCurve(): YieldCurve(TYPE){}
    static IObject* defaultConstructor(){
        return new StochasticYieldCurve();
    }
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(StochasticYieldCurve, clazz);
        SUPERCLASS(YieldCurve);
        IMPLEMENTS(IStochasticYieldCurve);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(yc, "the regular yield curve");
        FIELD(irVol, "swaption vols for this yield curve");
    }
        
};

CClassConstSP const StochasticYieldCurve::TYPE =
CClass::registerClassLoadMethod(
    "StochasticYieldCurve", typeid(StochasticYieldCurve), load);

DRLIB_END_NAMESPACE
