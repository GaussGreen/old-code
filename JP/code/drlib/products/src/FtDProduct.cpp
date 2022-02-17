//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : View into instrument as required for FtDs
//
// Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/ConvolutionModel.hpp"
#include "edginc/NToDefaultLossConfig.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/IEffectiveCurveGen.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/NullDecretionCurve.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/FtDProduct.hpp"
#include "edginc/Actual365F.hpp"

DRLIB_BEGIN_NAMESPACE


/** Product assumptions:
     - No roll off dates - to be relaxed later
     - All underlyings to have the same notional
     "payAtTheEnd" FtDs are not supported */

/** Class to ease output requests */
class FtDProduct::FtDOutputs {
public:
    FtDOutputs() : ctgPv(0.0), feePv(0.0), accInt(0.0)
    {}

    // Returns the price as ctgPv - feePv - accInt
    const double price() const {
        return ctgPv - feePv - accInt;
    }

    // Fields - note they are public
    double ctgPv;
    double feePv;
    double accInt;
};



/** Constructor for a ModelProduct that prices using a one-factor 
    Gaussian copula */
FtDProduct::FtDProduct(ConvolutionEngineConstSP convolutionEngine,
                       const GeneralisedCDO*    inst) :
    convolutionEngine(convolutionEngine), inst(inst)
{}


NToDefaultLossConfigConstSP FtDProduct::getFtDLossCfg(const GeneralisedCDO* inst) {
    static const string method("FtDProduct::getFtDLossCfg");

    ICreditLossConfigConstSP lossCfg = inst->getLossConfig();
    if (NToDefaultLossConfig::TYPE->isInstance(lossCfg)) {
        NToDefaultLossConfigConstSP ftdLossCfg( 
            DYNAMIC_CONST_CAST(NToDefaultLossConfig, lossCfg.get()));

        if (ftdLossCfg->getDefaultNumber() == 1) {
            return ftdLossCfg;
        }
    }
    throw ModelException(method, "The loss config object in the instrument "
                         "is not an FtD (i.e., is not an NToDefaultLossConfig "
                         "with N=1)");
}


//++++++++++++  IGeneralisedConvolutionProduct methods
//
/** Invoked by model to compute the price.
    Prices the FTD product using a semi-analytical closed form model
    and stores the result in the "results" object */
void FtDProduct::price(const IConvolutionModel* convolutionModel,
                       Control*                 control, 
                       Results*                 results,
                       IForwardRatePricerSP     fwdRateModel) const
{
    static const string method("FtDProduct::price");

    try {
        BoolArraySP     payAsYouGoArray;
        IntArraySP      numDelayDaysArray;
        DateTimeArraySP startDates;
        DateTimeArraySP endDates;
        DateTimeArraySP paymentDates;
        if (inst->cLeg.get()) {
            inst->cLeg->getPaymentInformation(payAsYouGoArray, 
                                              numDelayDaysArray,
                                              startDates, 
                                              endDates, 
                                              paymentDates);
            
            for (int i=0; i < payAsYouGoArray->size(); ++i) {
                if ((*payAsYouGoArray)[i] == false) {
                    throw ModelException(method, "FtDs can not yet be priced as "
                                         "'pay at the end'. Set the contingent "
                                         "leg's 'pay as you go' field/s to true.");
                    // If we were to support it, would need to implement
                    // methods in the FtDEffectiveCurve which are not 
                    // currently implemented
                }
            }
        } // else will use null parameters


        NToDefaultLossConfigConstSP ftdLossCfg = FtDProduct::getFtDLossCfg(inst);

        const DateTime& defaultDate = ftdLossCfg->getDefaultDate();
        const DateTime& valueDate = inst->getValueDate();
        YieldCurveConstSP discount = getDiscount();

        double outstandingNotional = ftdLossCfg->notional();
        CashFlowArray cLegNotionalReductions; // Empty
        BoolArray cLegPayNotionalReductions;  // Empty
        IDiscountCurveRiskySP feeEffCurve;
        IDiscountCurveRiskySP ctgEffCurve;
        FtDOutputs myOutput;
        if (defaultDate.empty()) { // ie, this FtD has NOT defaulted
            ICreditLossGenSP lossGen = 
                convolutionEngine->effCurveGenerator(ftdLossCfg,
                                                     inst->cptyInfo,
                                                     inst->getRecoverNotional(),
                                                     IModelConfigMapperConstSP());
            IEffectiveCurveGenSP effCurveGen =
                DYNAMIC_POINTER_CAST<IEffectiveCurveGen>(lossGen);
            if (!effCurveGen) {
                throw ModelException(method,
                                     "Internal error: The loss generator is not of "
                                     "type IEffectiveCurveGen");
            }
            
            // This is where all the interesting (and slow) stuff takes place
            effCurveGen->getFeeAndCtgEffectiveCurve(feeEffCurve, 
                                                    ctgEffCurve, 
                                                    discount,
                                                    convolutionModel,
                                                    lastObservationDate());
            
            // Got all required inputs, so just go and price the instrument 
            if (inst->fLeg.get()) {
                // JLHP - get the prepay info from the NToDefaultLossConfigConstSP... 
                // (if we ever put it there)

                IDecretionCurveConstSP prepay(
                    new NullDecretionCurve(ftdLossCfg->getName()));

                // if (payAccruedFee) { JLHP - add to GeneralisedCDO or assume always true
                DateTime earliestAccrualStartOrRiskyDate;
                DateTime latestAccrualEndOrRiskyDate;
                if (inst->cLeg.get()) {
                    earliestAccrualStartOrRiskyDate = 
                        inst->cLeg->firstObservationStartDate();
                    latestAccrualEndOrRiskyDate = 
                        inst->cLeg->lastObservationEndDate();
                }
                else {
                    // There is no contingent leg - set some defaults
                    earliestAccrualStartOrRiskyDate = valueDate;
                    latestAccrualEndOrRiskyDate = inst->fLeg->getLastPayDate();
                }

                // JLHP - Add dcc to GeneralisedCDO or assume this default
                DayCountConventionSP act365F(new Actual365F());
                myOutput.feePv = 
                    inst->fLeg->getFeeLegPV(valueDate, 
                                            valueDate, // settlementDate ? JLHP
                                            earliestAccrualStartOrRiskyDate,
                                            latestAccrualEndOrRiskyDate,
                                            *discount,
                                            *feeEffCurve, 
                                            prepay, 
                                            false, // includeAccrued
                                            act365F, // accrualDCC, not used here
                                            fwdRateModel);

                myOutput.accInt = 
                    inst->fLeg->getFeeLegAI(valueDate,
                                            valueDate, // settlementDate ? JLHP
                                            earliestAccrualStartOrRiskyDate,
                                            latestAccrualEndOrRiskyDate,
                                            act365F, // accrualDCC
                                            *feeEffCurve,
                                            discount,
                                            prepay,
                                            fwdRateModel);
            }
               
            // cLeg priced further down below, regardless of default status
        }
        ///////////////////////////////////////////////////////////////////////e
        else { // The FtD has defaulted
            throw ModelException(method, "Defaulted FtDs can not yet be priced.");

            outstandingNotional = 0.0; // since the name has defaulted

            if (inst->fLeg.get()) {
                AccrualPeriodArrayConstSP accrualPeriods; // null since not used

                FeeLegReductionPerDefaultArraySP fLegPastReductions =
                    ftdLossCfg->historicFeeLegReductions(
                        inst->triggerDelay,
                        inst->defaultToCalculationDelay,
                        inst->temporaryLossAmount->doubleValue(),
                        inst->lastTriggerDate,
                        accrualPeriods,
                        IBadDayAdjusterConstSP::attachToRef(inst),
                        true); // recoverNotional (hardcoded in FtDs like in CDSs)

                if (!!fLegPastReductions) {
                    if (fLegPastReductions->size() != 1) {
                        throw ModelException(method, "Internal error: In an FtD "
                                             "all the fee leg notional is "
                                             "expected to go down at once.");
                    }
                }

                CashFlowArraySP feeLegReductions;

                // Since the name has defaulted, the ctgEffCurve is now really riskless:
                // there are no more expected losses
                feeEffCurve.reset(new EffectiveCurve(
                    valueDate,
                    discount,
                    DateTimeArray(1, valueDate), // No more losses after default
                    DoubleArray(1, 0.0),   // No more losses after default
                    convolutionModel->getLossInterpolation()));
                
                double dummyDouble;
                DoubleArray dummyDoubleArray;
                myOutput.accInt = 
                    inst->fLeg->price(valueDate,
                                      valueDate, // valDateCF
                                      feeEffCurve,
                                      inst->discountYieldCurveWrapper(),
                                      0, // lowStrike
                                      ftdLossCfg->notional(), // highStrike
                                      0, // outstandingNotional,
                                      *feeLegReductions,
                                      dummyDouble,  // pastTrancheLosses
                                      dummyDouble,  // riskyDurationTotal
                                      dummyDouble,  // riskyNotionalsMean
                                      false, // do not compute debug stuff
                                      dummyDoubleArray, // debugUnitPrice,   
                                      dummyDoubleArray, // debugUnitHistPrice
                                      CashFlowArraySP(), // no rebate payments
                                      payAsYouGoArray,   // JLHP - check how this is used in an FtD - only here for accrual, not for "non-defaulted" fee legs
                                      numDelayDaysArray, // 
                                      startDates,        // 
                                      endDates,          // 
                                      paymentDates,      //
                                      IBadDayAdjusterConstSP::attachToRef(inst),
                                      fwdRateModel);
            }
            if (inst->cLeg.get()) {
                CtgLegLossPerDefaultArraySP cLegLosses =
                    ftdLossCfg->historicContingentLegLosses(
                        inst->triggerDelay,
                        inst->defaultToCalculationDelay,
                        inst->lastTriggerDate,
                        IBadDayAdjusterConstSP::attachToRef(inst),
                        inst); // IProtectionProvider

                // JLHP this is duplicated from GeneralisedCDO - factor it out
                int numCtgPayments = cLegLosses->size();
                cLegNotionalReductions.resize(numCtgPayments);
                cLegPayNotionalReductions.resize(numCtgPayments);
                double cumulativeAmount = 0.0;
                for (int i=0; i < numCtgPayments; ++i) {
                    // Here we need the cumulative losses on the contingent leg, so
                    // add them up manually... while keeping track of which ones have
                    // to be paid
                    cumulativeAmount += (*cLegLosses)[i]->loss.amount;
                    cLegNotionalReductions[i] = 
                        CashFlow((*cLegLosses)[i]->loss.date, cumulativeAmount);
                    cLegPayNotionalReductions[i] =
                        inst->isDateCoveredForProtection((*cLegLosses)[i]->defaultDate);
                }

                // Since the name has defaulted, the ctgEffCurve is now really riskless:
                // there are no more expected losses
                ctgEffCurve.reset(new EffectiveCurve(
                    valueDate,
                    discount,
                    DateTimeArray(1, valueDate), // No more losses after default
                    DoubleArray(1, 0.0), // No more losses after default
                    convolutionModel->getLossInterpolation()));
            }
        }

        if (inst->cLeg.get()) {
            // Price the contingent leg with the parameters computed above
            DoubleArray dummyDebugPrices; //empty
            myOutput.ctgPv = 
                inst->cLeg->price(ftdLossCfg->notional(), // initial notional
                                  outstandingNotional,
                                  valueDate, // today
                                  valueDate, // valDateCF
                                  ctgEffCurve,
                                  cLegNotionalReductions,
                                  cLegPayNotionalReductions,
                                  false, // computeDebugPrices
                                  dummyDebugPrices,
                                  dummyDebugPrices,
                                  IBadDayAdjusterConstSP::attachToRef(inst));
        }

        const double pv = myOutput.price();
        results->storePrice(pv, discount->getCcy());

        storeOutputRequests(myOutput, control, results);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    return;
}


void FtDProduct::storeOutputRequests(FtDOutputs myOutput, 
                                     Control*   control, 
                                     Results*   results) const
{
    OutputRequest* request; // reused across different OutputRequests

    request = control->requestsOutput(OutputRequest::FEE_LEG_FV);
    if (request) {
        results->storeRequestResult(request, myOutput.feePv);
    }

    request = control->requestsOutput(OutputRequest::CONTINGENT_LEG_FV);
    if (request) {
        results->storeRequestResult(request, myOutput.ctgPv);
    }

    request = control->requestsOutput(OutputRequest::ACCRUED_INT_FV);
    if (request) {
        results->storeRequestResult(request, myOutput.accInt);
    }
    return;    
}


/** Returns the object that defines the losses for this product */
ICreditLossConfigConstSP FtDProduct::getLossConfig() const {
    return inst->getLossConfig();
}

/** Returns the last observation date for credit related data */
DateTime FtDProduct::lastObservationDate() const {
    return inst->lastObservationDate();
}

/** Returns the max of maxDate and the last date from when a yield
    curve is used (eg for discounting or rate estimation) by
    payoff method.  This method does not have to worry about eg what
    dates are used when the clean spread curves are built. Used to
    control when to stop tweaking */
DateTime FtDProduct::lastYCSensDate(const DateTime& maxDate) const {
    return inst->lastYCSensDate(maxDate);    
}

/** Returns the [riskless] curve for discounting (could remove this if
    we had better methods on ICDSParSpreads - only used for calculating
    durations) */
YieldCurveConstSP FtDProduct::getDiscount() const {
    return inst->discountYieldCurveWrapper().getSP();
}
//
//------------  IGeneralisedConvolutionProduct methods


DRLIB_END_NAMESPACE
