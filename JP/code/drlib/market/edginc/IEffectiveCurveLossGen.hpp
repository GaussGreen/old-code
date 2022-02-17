//----------------------------------------------------------------------------
//                                                                           
// Group       : CH Quantitative Research                                    
//                                                                           
// Description : Interface defining what a generator for 'effective loss 
//               curve' based engine needs to be able to do.
//                                                                           
// Date        : July 2006                                                   
//                                                                           
//----------------------------------------------------------------------------

#ifndef QLIB_IEFFECTIVECURVELOSSGEN_HPP
#define QLIB_IEFFECTIVECURVELOSSGEN_HPP

#include "edginc/Control.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ICreditLossGen.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(IFixedTrancheLossCalculator);
FORWARD_DECLARE(CounterPartyCredit);
FORWARD_DECLARE(YieldCurve);
FORWARD_DECLARE(BetaSkewGridPoint);
FORWARD_DECLARE(ConvolutionProduct);

/** What a generator for an 'effective loss curve' based engine needs to be
    able to do */
class MARKET_DLL IEffectiveCurveLossGen : virtual public ICreditLossGen {

public:
    virtual ~IEffectiveCurveLossGen();

    /** Generate the timeline on which the effective curve will be calculated.
        The lastObservationDate argument must be instrument's last observation 
        date for both the fee and contingent legs */
    virtual DateTimeArraySP generateTimeline(
        const DateTime& lastObservationDate) const = 0;

    /** Returns a set of IFixedTrancheLossCalculator which are capable
        of returning expected tranche losses along the specified
        timeline. The cpty indicates whether losses conditional on the
        counterparty surviving are required (null if not). The control
        and results parameters are there for any output requests that
        should be calculated.  */
    virtual void createLossCalculators(
        const DateTimeArray&                timeline,           // (I)
        CounterPartyCreditConstSP           cpty,               // (I)
        const DateTime&                     maturity,           // (I)
        Control*                            control,            // (I)
        Results*                            results,            // (I)
        bool                                recoverNotional,    // (I)
        YieldCurveConstSP                   discount,           // (I)
        IFixedTrancheLossCalculatorConstSP& lossCalculator,               // (O)
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,  // (O)
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,          // (O)
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const = 0;// (O)

    /** Returns all the points on the skew surface to which this
        model for the supplied instrument is sensitive.
        A null return value is ok - it is interpreted as the greek being 
        NotApplicable */
    virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
        OutputNameConstSP outputName,
        const DateTime&   maturity,
        YieldCurveConstSP discount,
        bool              tweakAll,
        const int         numberOfNeighbours) const = 0;

    /** When to stop tweaking the loss config(s) associated to this loss gen */
    virtual DateTime maxLossConfigsTweakableDate(
        const DateTime& creditEndDate) const = 0;
    
    /** Recovery of notional from the top of the portfolio requires
        an additional call to the convolution.
        It only has effect if the upper strike is greater than the
        sum of (name notional * name recovery).
        Therefore we can avoid making this additional call with the
        models assistance if the product displays the correct
        characteristics. Typically this means the model must be of
        a fixed recovery type */
    virtual const bool modelRecoveredNotional() const = 0;

    /** Returns how the effectiveCurve should be interpolated. This should
        really disappear when we return a DefaultRate object say in the 
        calculateEffectiveCurve() method. */
    virtual const string& getLossInterpolation() const = 0;

protected:
    IEffectiveCurveLossGen();

private:    
    IEffectiveCurveLossGen(const IEffectiveCurveLossGen& rhs); // don't use
    IEffectiveCurveLossGen& operator=(const IEffectiveCurveLossGen& rhs); // don't use
};

DECLARE_REF_COUNT(IEffectiveCurveLossGen);

DRLIB_END_NAMESPACE

#endif
