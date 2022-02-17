//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : CounterPartyCredit.hpp
//
//   Author      : essentially from ccm2 library
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_COUNTERPARTYCREDIT_HPP
#define QLIB_COUNTERPARTYCREDIT_HPP

#include "edginc/CCMBetaSens.hpp"
#include "edginc/CCMAbsoluteBetaTweak.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE
FORWARD_DECLARE_WRAPPER(CreditAsset);
FORWARD_DECLARE(SingleCreditAsset);
FORWARD_DECLARE(EffectiveCurve);

/** Captures information needed for counter party credit calculation. Not clear
    at this stage where this file should live - possibly in market. Although
    it is a bit trade specific - not least because of the beta. 
    This class is pretty weak - needs more useful methods. Currently used
    as a data container. */
class MARKET_DLL CounterPartyCredit: /*     CSA Description      */
    public CObject,
    public virtual IGetMarket/*,
    virtual public CCMBetaSens::IShift,
    virtual public TweakableWith<CCMAbsoluteBetaTweak> */
{
public:
    static CClassConstSP const TYPE;

    virtual ~CounterPartyCredit();

    //// basic validation to do upon object construction
    virtual void validatePop2Object();

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the name of the CreditAsset which holds the market data for
        the counterparty. */
    string getName() const;

    /** Returns the credit asset for the counterparty */
    SingleCreditAssetConstSP getAsset() const;

    /** Returns the number of entries in the threshold schedule */
    int nbThreshold() const;

    /**
     * Returns the recovery:
     * - field "nameRecovery" if defaultParamOverride=TRUE
     * - market recovery otherwise
     *  */
    double getRecovery() const;

    /** Returns the 'beta' for this name (which is a view onto the correlation
        of this name with other names in the trade) */
    double getBeta() const;

    /**
     * Returns the default date:
     * - field "nameDefDate" if specified
     * - market default date (if any) otherwise
     *  */
    DateTime getDefaultDate() const;

    /**
     * Calculates the conuterparty charge with collateral agreement if
     * applicable. * If there is no CSA, then the creditCharge is just
     * (1-Rcpty)*(PriceUncond-PriceCond). * When there is a CSA, then the
     * max amt lost is the threshold, so the CSA * charge is seen as a CDS
     * ctg leg with notional = threshold. The threshold * is time-varying
     * and the dates in input as nbrs of years after ValueDate.  * An
     * additional term comes from the possible diffusion of the MTM value
     * * between the time the threshold is reached and the time the
     * payment is * settled: "a" * The confidence factor
     * cpty->confidenceLevel is the probability that the collateral call *
     * works (depends of the confidence in the bankrupcy laws of the
     * country) * The creditCharge returned is the min (creditCharge with
     * CSA, creditCharge without CSA) to ensure CSA * never increases
     * creditCharge.  * * Differences are to be expected with CSM as the
     * timeline is not exactly the * same *
     */
    double calculateCSACapPrice(
        bool                   isLong,
        double                 price,
        double                 priceCond,
        double                 priceUncond,
        const EffectiveCurveSP cptyCurve,
        const DateTime&        today,
        const DateTime&        lastDay) const;
        
private:
    CounterPartyCredit(const CounterPartyCredit& rhs);
    CounterPartyCredit& operator=(const CounterPartyCredit& rhs);
    CounterPartyCredit();

    static IObject* defaultConstructor();

    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    /// fields ////
    CreditAssetWrapper counterPartyName;       /* counterparty name */

    CDoubleSP beta;
    int cltrCallDays;               /* collateral call frequency in days */
    int curePeriod;                 /* cure period (days) */
    DateTime nameDefDate;           /* cpty def date */
    DoubleArray thresholdEnd;       /* end date for the threshold schedule
                                       given as a interval (in years) */
    DoubleArray thresholdAmt;        /* CSA threshold amt until the 
                                        ThresholdEnd */

    bool defaultParamOverride;

    double nameRecovery;             /* cpty recovery */
    double confidenceLevel;          /* confidence level of the CSA */
    double callFreqMultiplier;       /* collateral call frequency multiplier */

};

DECLARE(CounterPartyCredit);

DRLIB_END_NAMESPACE

#endif
