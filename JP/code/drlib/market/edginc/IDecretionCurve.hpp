//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : IDecretionCurve.hpp
//
//   Description : A decretion curve interface for ABS
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#ifndef IDECRETIONCURVE_HPP
#define IDECRETIONCURVE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(Settlement);

/** Interface for objects representing Asset Backed CDS DecretionCurve objects.
 ** The decretion can either be loss or prepay that happens to the underlying, 
 ** which affects the principal of the CDS. Prepay simply decreases the 
 ** notional. Protection buyer will get compensated in case of loss.  
 */

class MARKET_DLL IDecretionCurve : public virtual IObject {
public:
    static CClassConstSP const TYPE;

    virtual ~IDecretionCurve();
    
    /** pv, which actually is balance here, can be relative to a start date */
    virtual double pv(const DateTime& startDate,
                      const DateTime& endDate) const = 0;
    
    /** pv (balance) relative to initial balance */
    virtual double pv(const DateTime& endDate) const = 0;

    /** return decretion speed on a date */
    virtual double getDecretionSpeed(const DateTime& date) const = 0;

    /** return if balances are stepwise or continuous */ 
    virtual bool isStepBalances() const = 0;

    /** Returns the step dates */
    virtual DateTimeArraySP getStepDates() const = 0;

    /** Returns the step dates */
    /** A modification for ABCDS, to be used in the integration timeline */
    /** Subtly different from the above - see the code */
    virtual DateTimeArraySP getAlternativeStepDates() const = 0;

    /** Returns current factor */
    virtual double getFactor(const DateTime& date) const = 0;

    /** Returns number of delayed business days for prepay */
    virtual SettlementConstSP getSettlement() const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<IDecretionCurve> IDecretionCurveConstSP;
typedef smartPtr<IDecretionCurve>      IDecretionCurveSP;
typedef MarketWrapper<IDecretionCurve> IDecretionCurveWrapper;
typedef vector<IDecretionCurveConstSP>      IDecretionCurveConstSPArray;
 
DRLIB_END_NAMESPACE
#endif
