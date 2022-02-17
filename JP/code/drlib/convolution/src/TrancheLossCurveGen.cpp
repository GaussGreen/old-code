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

#include "edginc/config.hpp"
#include "edginc/TrancheLossCurveGen.hpp"
#include "edginc/ConvolutionEngine.hpp"

DRLIB_BEGIN_NAMESPACE

TrancheLossCurveGen::~TrancheLossCurveGen()
{}

TrancheLossCurveGen::TrancheLossCurveGen(
        IConvolutionModelSP convolutionModel,
        ConvolutionEngineConstSP engineForInnerLossConfigs,
        CreditTrancheLossConfigConstSP tranche) :
    convolutionModel(convolutionModel), 
    engineForInnerLossConfigs(engineForInnerLossConfigs),
    tranche(tranche)
{}


/** Generate the timeline on which the effective curve will be 
    calculated. No customization implemented here, just forward 
    to the convolution model */
DateTimeArraySP TrancheLossCurveGen::generateTimeline(
    const DateTime& lastObservationDate) const 
{
    return convolutionModel->generateTimeline(tranche->getToday(),
                                              lastObservationDate);
}


/** Returns a set of IFixedTrancheLossCalculator which are capable
    of returning expected tranche losses along the specified
    timeline. The cpty indicates whether losses conditional on the
    counterparty surviving are required (null if not). The control
    and results parameters are there for any output requests that
    should be calculated.  */
void TrancheLossCurveGen::createLossCalculators(
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
    IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const // (O)
{
    return convolutionModel->createLossCalculators(timeline,
                                                   tranche,
                                                   cpty,
                                                   maturity,
                                                   control,
                                                   results,
                                                   recoverNotional,
                                                   discount,
                                                   lossCalculator,
                                                   recoveredNotionalCalculator,
                                                   conditionalLossCalc,
                                                   conditionalRecNtnlCalc);
}
        

/** Returns all the points on the skew surface to which this
    model for the supplied instrument is sensitive.
    A null return value is ok - it is interpreted as the greek being 
    NotApplicable
    The lastObservationDate argument must be instrument's last observation 
    date for both the fee and contingent legs */
BetaSkewGridPointArrayConstSP TrancheLossCurveGen::getSensitiveBetaSkewPoints(
    OutputNameConstSP outputName,
    const DateTime&   lastObservationDate,
    YieldCurveConstSP discount,
    bool              tweakAll,
    const int         numberOfNeighbours) const
{         
    return convolutionModel->getSensitiveBetaSkewPoints(outputName, 
                                                        lastObservationDate,
                                                        tranche, 
                                                        discount, 
                                                        tweakAll,
                                                        numberOfNeighbours);
}


/** When to stop tweaking the loss config(s) associated to this loss gen */
DateTime TrancheLossCurveGen::maxLossConfigsTweakableDate(
    const DateTime& creditEndDate) const 
{
    DateTime lastExpiry;
    DateTime maxDate(tranche->getToday());
    int numNames = tranche->numInnerLossConfigs();
    for (int i=0; i < numNames; ++i) {
        SingleCreditAssetConstSP myAsset(tranche->nameAsset(i));
        lastExpiry = myAsset->getParSpreadCurve()->stopTweaking(creditEndDate);
        if (lastExpiry.isGreater(maxDate)) {
            maxDate = lastExpiry;
        }            
    }
    return maxDate;
}


/** Recovery of notional from the top of the portfolio requires
    an additional call to the convolution.
    It only has effect if the upper strike is greater than the
    sum of (name notional * name recovery).
    Therefore we can avoid making this additional call with the
    models assistance if the product displays the correct
    characteristics. Typically this means the model must be of
    a fixed recovery type */
const bool TrancheLossCurveGen::modelRecoveredNotional() const {
    return convolutionModel->modelRecoveredNotional(
        tranche/*, engineForInnerLossConfigs*/); // do we need the 2nd parameter?
}


/** Returns how the effectiveCurve should be interpolated. This should
    really disappear when we return a DefaultRate object say in the 
    calculateEffectiveCurve() method. */
const string& TrancheLossCurveGen::getLossInterpolation() const {
    return convolutionModel->getLossInterpolation();
}

DRLIB_END_NAMESPACE
