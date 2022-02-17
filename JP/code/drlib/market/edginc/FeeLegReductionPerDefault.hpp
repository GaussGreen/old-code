//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Class containing all required information regarding a   
//                 reduction in notional in the fee leg (typically for CDOs).
//                 The reduction can be due to a loss or recovered notional.
//                 The information captured here shall permit the computation
//                 of the rebate payments corresponding to the reductions.
//                 The loss is defined as: 
//                    defaultedNotional * (1 - recoveryRate)
//                 And the recovered notional as:
//                    totalNameNotional - loss
//
//                 NOTE: No assumption is made about the recoveryRate - as long
//                 as the loss and recovered notional are computed correctly
//                 using the formulas above, the recoveryRate may well lay
//                 outside [0,1].
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FEELEGREDUCTIONPERDEFAULT_HPP
#define QLIB_FEELEGREDUCTIONPERDEFAULT_HPP

#include "edginc/Atomic.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);

class MARKET_DLL FeeLegReductionPerDefault : virtual public VirtualDestructorBase {
public:
    FeeLegReductionPerDefault(const DateTime& determinationDate, 
                              const DateTime& effectiveDate,
                              const DateTime& calculationDate,
                              const double defaultedNotional,
                              const double totalNameNotional,
                              const double recoveryRate);

    FeeLegReductionPerDefault(const FeeLegReductionPerDefault& reduction, 
                              const bool isLoss,
                              const double amount);

    const double getLossAmount() const;
    const double getRecoveredAmount() const;
    const double getReductionAmount(const bool wantLoss) const;

    static bool lessThanForLossDates(const FeeLegReductionPerDefaultSP x, 
                                     const FeeLegReductionPerDefaultSP y);

    /** Returns a CashFlowArray with all the reduction CashFlows (loss or
        recovered notional) inside the argument array of reductions */
    static CashFlowArraySP getReductions(
        const FeeLegReductionPerDefaultArrayConstSP x,
        const bool wantLoss);


    // FIELDS - note they are public
    DateTime determinationDate;

    /** Date when the actual loss takes place - in principle identical to
     calculation date, but may be shifted (e.g., in the fee leg of tranches
     it may be aligned to an accrue start date boundary)*/
    DateTime effectiveDate;

    /** Date when the actual loss is determined - if different from
        determinationDate, a rebate MAY have to be paid if the
        assumption made on determinationDate was not accurate */
    DateTime calculationDate; // aka rebate date

    double defaultedNotional; // Notional triggered for protection
    double totalNameNotional;
    double recoveryRate;

private:
    FeeLegReductionPerDefault(const FeeLegReductionPerDefault&);  // don't use
    FeeLegReductionPerDefault& operator=(const FeeLegReductionPerDefault&); // don't use
};

DRLIB_END_NAMESPACE

#endif
