//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Class containing all required information regarding a   
//                 loss in the contingent leg (typically for CDOs).
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONTINGENTLEGLOSSPERDEFAULT_HPP
#define QLIB_CONTINGENTLEGLOSSPERDEFAULT_HPP

#include "edginc/Atomic.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);


class MARKET_DLL CtgLegLossPerDefault : virtual public VirtualDestructorBase {
public:
    CtgLegLossPerDefault(const DateTime& defaultDate, 
                         const CashFlow& loss);

    CtgLegLossPerDefault(const DateTime& defaultDate, 
                         const DateTime& lossDate,
                         const double    lossAmount);

    static bool lessThanForLossDates(CtgLegLossPerDefaultSP x, 
                                     CtgLegLossPerDefaultSP y);

    static CtgLegLossPerDefaultArraySP produceLossesPerDefault(
        const DateTime&       defaultDate, 
        const CashFlowArraySP losses);

    /** Returns a CashFlowArray with all the reduction CashFlows inside 
        the argument array of reductions */
    static CashFlowArraySP getReductions(
        const CtgLegLossPerDefaultArrayConstSP x);

    // FIELDS - note they are public
    DateTime defaultDate;
    CashFlow loss;

private:
    CtgLegLossPerDefault(const CtgLegLossPerDefault&);  // don't use
    CtgLegLossPerDefault& operator=(const CtgLegLossPerDefault&); // don't use
};

DRLIB_END_NAMESPACE

#endif
