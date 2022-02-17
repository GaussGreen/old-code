//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : This class just wraps the existing IConvolutionModel
//               and holds the extra data needed to implement the 
//               IEffectiveCurveLossGen interface. 
//
// Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHELOSSCURVEGEN_HPP
#define QLIB_TRANCHELOSSCURVEGEN_HPP

#include "edginc/SingleCreditAsset.hpp"
#include "edginc/IEffectiveCurveLossGen.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/ConvolutionModel.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(ConvolutionEngine);

class CONVOLUTION_DLL TrancheLossCurveGen : virtual public IEffectiveCurveLossGen 
{

public:
    ~TrancheLossCurveGen();

    TrancheLossCurveGen(IConvolutionModelSP convolutionModel,
                        ConvolutionEngineConstSP engineForInnerLossConfigs,
                        CreditTrancheLossConfigConstSP tranche);

    /** Generate the timeline on which the effective curve will be 
        calculated. No customization implemented here, just forward 
        to the convolution model */
    virtual DateTimeArraySP generateTimeline(
        const DateTime& lastObservationDate) const; 

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
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const;// (O)      

    /** Returns all the points on the skew surface to which this
        model for the supplied instrument is sensitive.
        A null return value is ok - it is interpreted as the greek being 
        NotApplicable
        The lastObservationDate argument must be instrument's last observation 
        date for both the fee and contingent legs */
    virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
        OutputNameConstSP outputName,
        const DateTime&   lastObservationDate,
        YieldCurveConstSP discount,
        bool              tweakAll,
        const int         numberOfNeighbours) const;

    /** When to stop tweaking the loss config(s) associated to this loss gen */
    virtual DateTime maxLossConfigsTweakableDate(
        const DateTime& creditEndDate) const;

    /** Recovery of notional from the top of the portfolio requires
        an additional call to the convolution.
        It only has effect if the upper strike is greater than the
        sum of (name notional * name recovery).
        Therefore we can avoid making this additional call with the
        models assistance if the product displays the correct
        characteristics. Typically this means the model must be of
        a fixed recovery type */
    virtual const bool modelRecoveredNotional() const;

    /** Returns how the effectiveCurve should be interpolated. This should
        really disappear when we return a DefaultRate object say in the 
        calculateEffectiveCurve() method. */
    virtual const string& getLossInterpolation() const;

private:
    IConvolutionModelSP convolutionModel;
    ConvolutionEngineConstSP engineForInnerLossConfigs;
    CreditTrancheLossConfigConstSP tranche;
};

DECLARE_REF_COUNT(TrancheLossCurveGen);

DRLIB_END_NAMESPACE

#endif
