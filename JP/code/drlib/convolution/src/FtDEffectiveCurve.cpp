//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : Risky discount curve to price FtD contingent legs
//
// Date        : November 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/RiskyLogOfDiscFactorKey.hpp"
#include "edginc/ConditionalFtDLossIntegrand.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/FtDEffectiveCurve.hpp"


DRLIB_BEGIN_NAMESPACE

//default destructor
FtDEffectiveCurve::~FtDEffectiveCurve()
{}


/** Constructor with a risk free discount curve */
FtDEffectiveCurve::FtDEffectiveCurve(
        const DateTime&             valueDate,
        const YieldCurveConstSP     dsc,
        DateTimeArraySP             timeline,
        const vector<ICondLossDistributionsGenKeyArraySP>& keysByDate,
        IConditionalDefaultsModelSP condDefaultsModel,
        DoubleArraySP               namesLoss,
        const double                ftdNotional) :
    CObject(TYPE),
    valueDate(valueDate),
    timeline(timeline),
    keysByDate(keysByDate),
    condDefaultsModel(condDefaultsModel),
    namesLoss(namesLoss),
    ftdNotional(ftdNotional)
{
    //build the discount curve
    discount = IDiscountCurveSP::constCast(dsc); //used as is

    validate();
}

IObject* FtDEffectiveCurve::clone() const {
    static const string method("FtDEffectiveCurve::clone");
    throw ModelException(method, "JLHP method not yet implemented");
}


void FtDEffectiveCurve::validate() {
    static const string method("FtDEffectiveCurve::validate");

    //basic validation
    if (timeline->empty()) {
        throw ModelException(method, "Timeline must be non empty");
    }
    if ((*timeline)[0] != valueDate) {
        throw ModelException(method, "First point of the timeline must be today");
    }
}


/**Returns the value at paymentDate (and conditional on no default before then)
 * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
 * in case of default between startDate and endDate, and zero otherwise.
 * The default payment is made with a delay of recoveryDelay calendar days.
 * Whether integration is done continuously or discretely and what,
 * if any approximations, are used is up to the curve implementation. */
double FtDEffectiveCurve::protectionPV(const DateTime& paymentDate,   // jlhp - not used?
                                       const DateTime& startDt,
                                       const DateTime& endDt,
                                       RecoveryType    recTyp,
                                       double          recoveryDelay) const //number of delay days
{
    static const string method("FtDEffectiveCurve::protectionPV");

    if (recTyp != IDiscountCurveRisky::RECOVER_1) {
        throw ModelException(method, "RecoveryTypes other than 'Recover_1' "
                             "are not supported.");
    }

    if (condDefaultsModel->isCondFunctionIntegrationTimeDependent()) {
        throw ModelException(method, "Cannot price FtDs with model '" +
                             condDefaultsModel->getClass()->getName() +
                             "' (the integration accross market factors may "
                             "be time dependant, and the closed form solution "
                             "is therefore not valid).");
    }

    // get the integration dates
    DateTimeArrayConstSP integrationDates =
        EffectiveCurve::buildIntegrationDates(
            startDt,
            endDt,
            recoveryDelay,
            discount,
            timeline);
    
    ConditionalFtDLossIntegrand condFtDLossFunction(
        keysByDate,
        timeline,
        integrationDates,
        namesLoss,
        condDefaultsModel->marketFactorDimension(),
        valueDate,
        recoveryDelay,
        discount);

    double pv = condDefaultsModel->integrateCondFunction(
        &condFtDLossFunction, 
        ICondLossDistributionsGenKeyArraySP(), // Null
        valueDate); // dummy date (we know integration is NOT time dependent)

    // Normalize the pv by the ftdNotional
    pv /= ftdNotional; 

    return pv;
}

 
//----------------------------
// IDiscountCurveRisky methods
//----------------------------
/**Returns the probability that there are no default events from d1 to d2,
    conditional that there are no default events up to d1*/
double FtDEffectiveCurve::survivalProb(const DateTime& d1,
                                       const DateTime& d2) const
{
    static const string method("FtDEffectiveCurve::survivalProb");
    throw ModelException(method, "Not yet implemented.");
}

/**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
double FtDEffectiveCurve::survivalProb(const DateTime& dt) const
{
    static const string method("FtDEffectiveCurve::survivalProb(dt)");
    throw ModelException(method, "Not yet implemented.");
}

/**Returns the market recovery rate on defaulted assets given default at
    time defaultDate. This allows the possibility of a term-structure of recovery rates. */
double FtDEffectiveCurve::getRecovery(const DateTime& defaultDate) const
{
    static const string method("FtDEffectiveCurve::getRecovery");
    throw ModelException(method, "Not yet implemented.");
}

/**Recovery-rate for immediate default.*/
double FtDEffectiveCurve::getRecovery() const
{
    static const string method("FtDEffectiveCurve::getRecovery");
    throw ModelException(method, "Not yet implemented.");
}


/**Returns the value at paymentDate (and conditional on no default before then)
    * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
    * in case of default between startDate and endDate, and zero otherwise.
    * The default payment is made on recoveryDate, no matter when the default happens.
    * Whether integration is done continuously or discretely and what,
    * if any approximations,
    * are used is up to the curve implementation. */
double FtDEffectiveCurve::protectionPV(const DateTime&     paymentDate,
                                       const DateTime&     startDt,
                                       const DateTime&     endDt,
                                       RecoveryType        recTyp,
                                       const DateTime&     recoveryDate) const
{
    static const string method("FtDEffectiveCurve::protectionPV");
    throw ModelException(method, "Not yet implemented.");
}


/**Returns the value at paymentDate (and conditional on no default before then)
    * of a sequence of payments, with simple linear accrued-interest
    * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
    * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for
    * protectionPV. This allows the curve to compute default-accrual PVs itself
    * in a way which reflects the type of curve (e.g. flat forwards, etc.)
    * The accrual periods are the intervals between consecutive payment dates. To have a first
    * accrual period, you should set the first payment to zero, and then the first date is
    * the beginning of the first accrual period.
    * Payments before paymentDate are ignored, except to the extent that they affect accrued
    * interest due.
    * Any default payments are made with a delay of recoveryDelay calendar days after default.
    * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
    * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
    * this method should return the same value as pv(payments, paymentDate), inherited
    * from IDiscountCurve. */
double FtDEffectiveCurve::annuityPV(const CashFlowArray&    payments,
                                    const DateTime&         paymentDate,
                                    RecoveryType            accruedRecTyp,
                                    double                  recoveryDelay,
                                    DateTime                accrueStartDate) const
{
    static const string method("FtDEffectiveCurve::annuityPV");
    throw ModelException(method, "Not yet implemented.");
}

/**Returns the value at paymentDate (and conditional on no default before then)
    * of a sequence of payments, with simple linear accrued-interest
    * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
    * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for
    * protectionPV. This allows the curve to compute default-accrual PVs itself
    * in a way which reflects the type of curve (e.g. flat forwards, etc.)
    * The accrual periods are the intervals between consecutive payment dates. To have a first
    * accrual period, you should set the first payment to zero, and then the first date is
    * the beginning of the first accrual period.
    * Payments before paymentDate are ignored, except to the extent that they affect accrued
    * interest due.
    * Any default payments are made on recoveryDate, no matter when the default happens.
    * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
    * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
    * this method should return the same value as pv(payments, paymentDate), inherited
    * from IDiscountCurve. */
double FtDEffectiveCurve::annuityPV(const CashFlowArray&    payments,
                                    const DateTime&         paymentDate,
                                    RecoveryType            accruedRecTyp,
                                    const DateTime&         recoveryDate,
                                    DateTime                accrueStartDate) const 
{
    static const string method("FtDEffectiveCurve::annuityPV.");
    throw ModelException(method, "Not yet implemented.");
}


//------------------------------------------------------
// Riskless PV methods
//------------------------------------------------------

/** Compute discount factor between two dates. This is the price
    * at time date1 of a zero-coupon (discount) bond that pays 1
    * at time date2 if it still exists then.
    * Note that, if this bond is risky, this method
    * in IDiscountCurveRisky means the PV of such a bond which knocks
    * out on default with no recovery at all, and the price is
    * contingent on no default before date1.
    * @param date1 payment settlement date (conditional on existence at date1)
    * @param date2 payment/zero-coupon bond maturity date
    * @return Discount factor between date1 & date2
    */
double FtDEffectiveCurve::risklessPV(const DateTime& date1,
                                  const DateTime& date2) const
{
    return discount->pv(date1, date2);
}

/** Compute price for settlement today of a zero-coupon bond
    * maturing at date. Note settlement is to TODAY, not to
    * the curve spot date. This is because some curves
    * may have ambiguous spot-dates - for example should a combined
    * credit and rates curve have spot date T+1 or T+2?
    * @param date To get discount factor/PV for
    * @return Discount factor between today & given date
    */
double FtDEffectiveCurve::risklessPV(const DateTime& date) const
{
    return discount->pv(date);
}

/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double FtDEffectiveCurve::risklessPV(const CashFlowArray& cashFlows,
                                  const DateTime&      baseDate) const
{
    return discount->pv(cashFlows, baseDate);
}

//-----------------------
// IDiscountCurve methods
//-----------------------

/** @return Discounting curve's currency - typically this is the
    ISO code eg "USD" (although it is up to the client - you may
    assume however that yield curves in the same currency return the
    same value in getCcy(). Note it is NOT the name of the yield
    curve (which might be eg "GBP-LIBOR"). */
string FtDEffectiveCurve::getCcy() const
{
    return discount->getCcy();
}

/** Compute discount factor between two dates. This is the price
    * at time date1 of a zero-coupon (discount) bond that pays 1
    * at time date2 if it still exists then.
    * Note that, if this bond is risky, this method
    * in IDiscountCurveRisky means the PV of such a bond which knocks
    * out on default with no recovery at all, and the price is
    * contingent on no default before date1.
    * @param date1 payment settlement date (conditional on existence at date1)
    * @param date2 payment/zero-coupon bond maturity date
    * @return Discount factor between date1 & date2
    */
double FtDEffectiveCurve::pv(const DateTime& date1,
                             const DateTime& date2) const
{
    static const string method("FtDEffectiveCurve::pv");
    throw ModelException(method, "Not yet implemented.");
}

/** Compute price for settlement today of a zero-coupon bond
    * maturing at date. Note settlement is to TODAY, not to
    * the curve spot date. This is because some curves
    * may have ambiguous spot-dates - for example should a combined
    * credit and rates curve have spot date T+1 or T+2?
    * @param date To get discount factor/PV for
    * @return Discount factor between today & given date
    */
double FtDEffectiveCurve::pv(const DateTime& date) const
{
    static const string method("FtDEffectiveCurve::pv(date)");
    throw ModelException(method, "Not yet implemented.");
}

/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double FtDEffectiveCurve::pv(const CashFlowArray& cashFlows,
                             const DateTime&      baseDate) const
{
    double value = .0;
    for (int i=0; i < cashFlows.size(); ++i)
    {
       const DateTime cfDate = cashFlows[i].date;
		if (cfDate.isGreater(baseDate))
                    value += cashFlows[i].amount * pv(baseDate, cfDate);
    }
    return value;
}

/** return the bootstrapped dates */
DateTimeArray FtDEffectiveCurve::zeroDates() const
{
    static const string method("FtDEffectiveCurve::zeroDates");
    throw ModelException(method, "Not yet implemented.");
//     DateTimeArraySP mergedDates = discount->zeroDates();
//     return *(mergedDates.get());
}

// accessor methods for logOfDiscFactorKey
IDiscountCurve::IKey* FtDEffectiveCurve::getDiscountKey() const {
    return discount->logOfDiscFactorKey();
}

DefaultRates::IKey* FtDEffectiveCurve::getRiskyKey() const {
    static const string method("FtDEffectiveCurve::getRiskyKey");
    throw ModelException(method, "Not yet implemented.");
}

/** Returns a key used to optimise repeated calculations of discount
    factors/forward rate. The calc method for this key returns the 
    natural logarithm of the discount factor (or equivalently the
    product of the forward rate (continuous, Act/365F) and the negative
    year fraction (Act/365F) betweeen the two dates.
    The default implementation has no performance improvements. */
//this is questionable !
IDiscountCurve::IKey* FtDEffectiveCurve::logOfDiscFactorKey() const
{
    static const string method("FtDEffectiveCurve::logOfDiscFactorKey");
    throw ModelException(method, "Not yet implemented.");
}



/** Returns a DefaultRates object which gives access to 
    useful functionality including "default rates", aka clean default
    spreads. The information needed to bootstrap the clean spreads is
    obtained from this object.
    The DefaultRate returned is immutable and therefore it will not be
    changed - which means that there is no need to clone it */
DefaultRatesSP FtDEffectiveCurve::getDefaultRates() const {
    static const string method("FtDEffectiveCurve::getDefaultRates");
    throw ModelException(method, "Not yet implemented");
}


//-----------------------------------------------
// IGetMarket methods - exposed by IDiscountCurve
//-----------------------------------------------

//// populate from market cache
void FtDEffectiveCurve::getMarket(const IModel* model, const MarketData* market)
{
    //nothing to do
}

//-------------------------------
// FtDEffectiveCurve private methods
//-------------------------------

void FtDEffectiveCurve::load(CClassSP& clazz) {
    REGISTER(FtDEffectiveCurve, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDiscountCurveRisky);

    FIELD(valueDate, "Reference date");
    FIELD(timeline, "The curve dates");

    // unregistered    FIELD(keysByDate, "Keys used to obtain the 'names' survival probabilities");
    FIELD(condDefaultsModel, "Model to use to integrate accross market factors");
    FIELD(discount, "may be risk-free or risky");
    FIELD(namesLoss, "Per name loss given default");
}

CClassConstSP const FtDEffectiveCurve::TYPE = CClass::registerClassLoadMethod(
    "FtDEffectiveCurve", typeid(FtDEffectiveCurve), load);

DRLIB_END_NAMESPACE
