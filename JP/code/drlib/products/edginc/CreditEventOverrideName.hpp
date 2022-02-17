//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CreditEventOverrideName.hpp
//
//   Description : Simple class used by credit instruments (e.g., CIS) in order 
//                 to override (at instrument level) a name's default-related 
//                 parameters.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITEVENTOVERRIDENAME_HPP
#define QLIB_CREDITEVENTOVERRIDENAME_HPP

#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/ICDSParSpreads.hpp" // For ICDSParSpreadsWrapper
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IForwardRatePricer);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);

/** Simple class used by credit instruments (e.g., CIS) in order to
 * override (at instrument level) a name's default-related parameters */
class PRODUCTS_DLL CreditEventOverrideName: public CObject,
                                            public virtual ICreditEventOverrideName
{
public:
    static CClassConstSP const TYPE;
    virtual ~CreditEventOverrideName();

    /** Name of this market object */
    virtual string getName() const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Returns the fee leg loss cashflows corresponding to a defaulted name,
     * assuming that the default settles in one date. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const;

   /** Returns the contingent leg loss cashflows corresponding to a defaulted 
     * name, assuming that the default settles in one date. */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const;

    /* Returns the overall recovery rate accross all settlements, such that
     * total loss = overall defaulted notional * (1-overall recovery rate) */
    virtual double getOverallRecoveryRate(double recoveryRate) const;

    /** Returns the percentage of the notional that has been triggered 
     * for protection */
    virtual double getOverallDefaultedNotionalFraction() const;

     /** GetMarket implementation*/
    virtual void getMarket(const IModel* model, const MarketData* market);

private:
    // For reflection
    static void load (CClassSP& clazz);
    CreditEventOverrideName();
    static IObject* defaultCreditEventOverrideName();
    
    // Fields
    ICDSParSpreadsWrapper name;       // Name of the curve the override applies to
    DateTime  eventDeterminationDate; // When fees and accrual stop
    DateTime  defaultPayDate;         // When fees and (1-R) are paid
    CDoubleSP recovery;               // optional, hence CDoubleSP
};

DECLARE (CreditEventOverrideName);

DRLIB_END_NAMESPACE

#endif
