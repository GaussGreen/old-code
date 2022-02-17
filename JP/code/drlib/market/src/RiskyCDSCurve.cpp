//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskyCDSCurve.hpp
//
//   Description : contains a yield curve and a spread curve.
//
//   Author      : André Segger
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RiskyCDSCurve.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"


DRLIB_BEGIN_NAMESPACE

class RiskyCDSCurve::LogOfDiscFactorKey: public YieldCurve::IKey
{
public:
    LogOfDiscFactorKey(YieldCurve::IKey* pRiskless, const CleanSpreadCurve& pCleanSpreadCurve) : 
      riskless(pRiskless), cleanSpreadCurve(pCleanSpreadCurve)
    {
    }

    virtual double calc(const DateTime& newLoDate, const DateTime& newHiDate)
    {
        return riskless->calc(newLoDate, newHiDate) 
             * cleanSpreadCurve.getDefaultPV(newLoDate, newHiDate);
    }

private:
    auto_ptr<YieldCurve::IKey> riskless;
    const CleanSpreadCurve&    cleanSpreadCurve;
};



RiskyCDSCurve::RiskyCDSCurve(): YieldCurve(TYPE){};

RiskyCDSCurve::RiskyCDSCurve(const string&         myname, 
                             const IYieldCurve&    ccyCurve, 
                             const ICDSParSpreads& spreadCurve, 
                             const DateTime&       today) :
    YieldCurve(TYPE), riskFreeCurve(&ccyCurve), name(myname)
{
    // initialise the derived (transient) fields
    // (want fwd rates)
    cleanSpreadCurve = CDSHelper::getCleanSpreadCurve(spreadCurve, today, true);
    recovery = spreadCurve.getRecovery();
}

RiskyCDSCurve::RiskyCDSCurve(const string&           myname, 
                             const IYieldCurve&      ccyCurve, 
                             const CleanSpreadCurve& spreadCurve, 
                             double                  recovery) :
    YieldCurve(TYPE), riskFreeCurve(&ccyCurve), name(myname), 
        cleanSpreadCurve(&spreadCurve), recovery(recovery)
{
}

string RiskyCDSCurve::getCcy() const
{
    return riskFreeCurve->getCcy();
}

string RiskyCDSCurve::getName() const
{
    return name;
}


DateTime RiskyCDSCurve::getSpotDate() const {
    return riskFreeCurve->getSpotDate();
}


DateTime RiskyCDSCurve::getToday() const {
    return riskFreeCurve->getToday();
}

/** Useful accessor methods */
ExpiryArrayConstSP RiskyCDSCurve::getExpiries() const
{
    return riskFreeCurve->getExpiries();
}

StringArrayConstSP RiskyCDSCurve::getInstruments() const
{
    return riskFreeCurve->getInstruments();
}

/** Returns tradeDate + n business days where n = spot offset */
DateTime RiskyCDSCurve::settles(const DateTime& tradeDate) const{
    return riskFreeCurve->settles(tradeDate);
}


double RiskyCDSCurve::pv(const DateTime& lodate, const DateTime& hidate) const {
    return (riskFreeCurve->pv(lodate, hidate) * cleanSpreadCurve->getDefaultPV(lodate, hidate));
}


double RiskyCDSCurve::pv(const DateTime& date) const {
    return riskFreeCurve->pv(date) * cleanSpreadCurve->getDefaultPV(date);
}


/** calculate the risky PV based on the assumption that the bond recovers 
    recovery * recoveryNotional in case of default - this is the proper 2-state discounting */
double RiskyCDSCurve::riskyPV(const DateTime& lodate,
                              const DateTime& hidate,
                              double          cashFlow,
                              double          recoveryNotional) const {
    double defaultProb = cleanSpreadCurve->getDefaultProb(lodate, hidate);
    double pv = (1.0 - defaultProb ) * riskFreeCurve->pv(lodate, hidate) * cashFlow +
                defaultProb * recovery * recoveryNotional;

    return pv;
}

/** this function calculates the discount factor based on the assumption that the
    on which the recovery is based is provided externally. This allows to use
    different methodologies (PV, face value + accrued etc.) to be included easily  -
    this function will use the externally given recovery rather than the underlying
    risky curve's recovery */
double RiskyCDSCurve::riskyPV(const DateTime& lodate,
                              const DateTime& hidate,
                              double          cashFlow,
                              double          recoveryNotional,
                              bool            useAssetRecovery,
                              double          assetRecovery) const {
    double pv = 0.0;
    if ( useAssetRecovery ) {
       double defaultProb = cleanSpreadCurve->getDefaultProb(lodate, hidate);
       pv = (1.0 - defaultProb ) * riskFreeCurve->pv(lodate, hidate) * cashFlow +
            defaultProb * assetRecovery * recoveryNotional;

    } else {
       pv =  riskyPV(lodate, hidate, cashFlow, recoveryNotional);
    }

    return pv;
}


/** grab the dates used in the zero curve */
DateTimeArray RiskyCDSCurve::zeroDates() const
{
    return DateTime::merge(riskFreeCurve->zeroDates(), cleanSpreadCurve->getDates());
}


double RiskyCDSCurve::zero(const DateTime& date) const
{
    throw ModelException("RiskyCDSCurve::zero", "This method has not been implemented yet!");
    return riskFreeCurve->zero(date);
}


/**
 * Returns a key used to optimise repeated calculations of forward rates.
 */
YieldCurve::IKey* RiskyCDSCurve::logOfDiscFactorKey() const
{
    return new LogOfDiscFactorKey(riskFreeCurve->logOfDiscFactorKey(), *cleanSpreadCurve);
}

/** Passes down to risk free YC */
DateTime RiskyCDSCurve::rateMaturity(
    const DateTime&         rateStart,
    const Expiry*           rateMaturity,
    const BadDayConvention* bdc) const { // optional
    return riskFreeCurve->rateMaturity(rateStart, rateMaturity, bdc);
}

double RiskyCDSCurve::fwd(const DateTime&           lodate, 
                          const DateTime&           hidate,
                          const DayCountConvention* dcc,
                          int                       basis) const{
    double disc = pv(lodate, hidate);
    return RateConversion::discountToRate(disc, lodate, hidate, dcc, basis);
}

double RiskyCDSCurve::fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const{
    throw ModelException("RiskyCDSCurve::fwd","Not implemented");
}

double RiskyCDSCurve::fwd(const DateTime&           payDate,
                       const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const{
    throw ModelException("RiskyCDSCurve::fwd","Not implemented");
}

/** make a risky curve from a credit spread curve */
IYieldCurveSP RiskyCDSCurve::makeRiskyCurve(
    const CreditSpreadCurve& spreadCurve,
    const  DateTime*    maturityDate) const
{
    throw ModelException("RiskyCDSCurve::makeRiskyCurve","Not implemented");
}


CreditSpreadCurveSP RiskyCDSCurve::makeCreditSpreadCurve(
	const string&        name,
	const CashFlowArray& defaultRates,
	double               recovery) const
{
    throw ModelException("RiskyCDSCurve::makeCreditSpreadCurve","Not implemented");
}


CashFlowArraySP RiskyCDSCurve::getRatesAndDates() const
{
    throw ModelException("RiskyCDSCurve::getRatesAndDates","Not implemented");
}


IYieldCurveSP RiskyCDSCurve::createForwardCurve(const DateTime& forwardDate) const
{
    throw ModelException("RiskyCDSCurve::createForwardCurve","Not implemented");
}


/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
IVolProcessed* RiskyCDSCurve::getProcessedVol(
    const CVolRequest* volRequest) const{
    throw ModelException("RiskyCDSCurve::getProcessedVol",
                         "Not yet implemented for RiskyCDSCurve");
}

class RiskyCurveHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        // clazz->setPublic(); - NOT  visible to EAS/spreadsheet
        REGISTER(RiskyCDSCurve, clazz);
        SUPERCLASS(YieldCurve);
        IMPLEMENTS(IRiskyCurve);
        EMPTY_SHELL_METHOD(defaultRiskyCurve);
        FIELD(name,      "risky curve name");
        FIELD(riskFreeCurve,    "risk free yield curve");
        FIELD(cleanSpreadCurve, "transient field - clean spread curve");
        FIELD_MAKE_TRANSIENT(cleanSpreadCurve);
        FIELD(recovery,  "transient field -recovery");
        FIELD_MAKE_TRANSIENT(recovery);
    }

    static IObject* defaultRiskyCurve() {
        return new RiskyCDSCurve();
    }
};

CClassConstSP const RiskyCDSCurve::TYPE = CClass::registerClassLoadMethod(
             "RiskyCDSCurve", typeid(RiskyCDSCurve), RiskyCurveHelper::load);


DRLIB_END_NAMESPACE
