//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : OptOnConvBond.cpp
//
//   Description : Option on Convertible Bond
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 25, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/OptOnConvBond.hpp"

#include "edginc/B30E360.hpp"
#include "edginc/Results.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Addin.hpp"
#include "edginc/BondParams.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE

static CFieldConstSP assetField;
static CClassConstSP derivativeAssetClass;


// copy market data relevant to the instrument
void OptOnConvBond::GetMarket(const IModel* model, const CMarketDataSP market){

    market->GetReferenceDate(valueDate);
    CAsset::getAssetMarketData(model, market.get(), "N", "", cvbAsset);

    if (DerivativeAsset::TYPE->isInstance(cvbAsset.get())) {
        DerivativeAsset* asset = dynamic_cast<DerivativeAsset*>(cvbAsset.get());
        cvb = ConvBondSP::dynamicCast(asset->getInstrument());
    } else {
        throw ModelException("OptOnConvBond::GetMarket", "Asset must be an AssetCVB");
    }

    cvbAsset->getMarket(model, market.get());
    (const_cast<IModel*>(model))->getInstrumentAndModelMarket(market.get(), cvb.get());

    if (!!bondSettle) {
        bondSettle->getMarket(model, market.get());
    }
}

void OptOnConvBond::validatePop2Object()
{
    static const string method = "OptOnConvBond::validatePop2Object";
    try{
        // convert the day count convention strings to objects
        swapDCC = DayCountConventionSP(DayCountConventionFactory::make(swapDCCString));
        lockoutDCC = DayCountConventionSP(DayCountConventionFactory::make(lockOutParameters->lockoutDCCString));
   
        swapInterval = MaturityPeriodSP(new MaturityPeriod(swapIntervalString));

    }
    catch (exception& e) {
       throw ModelException(e, method);
    }
            
    return;
}

void OptOnConvBond::Validate() {
    static const string method = "OptOnConvBond::Validate";

    try {
        // only allow DeutcheBank style for now
        if (!(ocbType == "D" || ocbType == "F" || ocbType == "Y" || ocbType == "P" || ocbType == "C" || ocbType == "B" )) {
            throw ModelException(method, "Only ocbTypes D, F, C and Y are implemented");
        }

        if (ocbType == "D" || ocbType == "C" || ocbType == "Y") {
            if (exerSched->lastDate() > swapMaturity) {
                throw ModelException(method, "swap maturity must be on or after the last OptOnConvBond exercise date.");
            }

            if (swapMaturity > cvb->bond->getUnadjMaturityDate()) {
                throw ModelException(method, "swap maturity must be on or before bond maturity.");
            }
        }

        if ( ocbType == "Y" && !assetSwapStrike ) {
            throw ModelException(method, "Missing assetSwapStrike object!");
        }

        if ( ocbType == "Y" && !FixedYieldParameters::TYPE->isInstance(assetSwapStrike.get())) {
            throw ModelException(method, "assetSwapStrike object must be of type FixedYieldParameters for fixed yield asset swaps!");
        }

        if (cvb->DECS == true || cvb->PERCS == true) {
            throw ModelException(method, "underlying convertible cannot be a mandatory.");
        }

        if ( !!lockOutParameters ) {
            lockOutParameters->validate();
        }

        if ( unwindNewSwap == false ) {
            // check that a swap fixing schedule is provided
            if ( !swapFixings ) {
                throw ModelException(method, "A swap fixing schedule must be provided if unwindNewSwap is false.");
            }

            // check that the last payment date coincides with the swap maturity
            if ( (*swapFixings)[swapFixings->size()-1]->getPaymentDate() != swapMaturity ) {
                throw ModelException(method, "The last payment date of the swap fixings (" + 
                                    (*swapFixings)[swapFixings->size()-1]->getPaymentDate().toString()
                                    + ") must be equal to the swap maturity (" + swapMaturity.toString() + ").");
            }
        }

        if ( ocbType == "C" ) {
            int count;
            string interval;
            swapInterval->decompose(count, interval);

            CashFlowArray tmpSpreadFlows = SwapTool::cashflows(
                valueDate.rollDate(-1),
                swapMaturity,
                !shortFrontStub,
                -spread,
                count,           // interval = count periods
                interval,          // e.g. Y, M, W, D
                swapDCC.get());

            if (tmpSpreadFlows.size() > 0 && tmpSpreadFlows[0].date < valueDate && tmpSpreadFlows[0].date != swapLastFixDate) {
                throw ModelException(method, "Please set last swap fixing date (" + 
                    swapLastFixDate.toString() + ")");
            }
        }

        

        cvb->Validate();

        if (!!bondSettle) {
            cvb->bond->setSettlement(bondSettle);
        }
            
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Pull out the component assets & correlations from the market data */
void OptOnConvBond::getMarket(const IModel* model, const MarketData* market){
    try{

//        bond->getMarket(model, market);
    } 
    catch (exception& e){
        throw ModelException(e, "OptOnConvBond::getMarket", "Failed for Bond");
    }
}


// -- Stuff every product needs to do -- //

/** what's today ? */
DateTime OptOnConvBond::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime OptOnConvBond::endDate(const Sensitivity* sensControl) const {
   DateTime cvbEndDate = cvb->endDate(sensControl);
   DateTime ocbEndDate = (exerSched->lastDate() > swapMaturity)?exerSched->lastDate():swapMaturity;

   return ((cvbEndDate>ocbEndDate)?cvbEndDate:ocbEndDate);
}

bool OptOnConvBond::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "OptOnConvBond::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}

void OptOnConvBond::getStrike(DateTime          baseDate, 
                              const double&     bondFloor, 
                              bool              includeLockout, 
                              bool*             isExercisable, 
                              double*           strike) const {
    static const string method = "OptOnConvBond::getStrike";
    try {
        if (ocbType == "D" || ocbType == "C") {
            // see if we're passed swap maturity
            if (baseDate > swapMaturity) {
                *isExercisable = false;
                *strike = 0.;

            } else {
                double exerPct;
                DateTime IRSwapStart;       // Interest rate swap start date
             
                // see if it's exercisable today
                try { // need a try block as an error is thrown if it's not exercisable
                    exerPct = exerSched->interpolate(baseDate);
                    *isExercisable = true;
                    IRSwapStart = baseDate;

                } catch (exception& ) {
                    *isExercisable = false;
                    // we should still calculate the strike as if it were exercisable today for info purposes
                    // use the exer level on the first day
                    exerPct = exerSched->interpolate(exerSched->firstDate()); 

                    if (baseDate < exerSched->firstDate()) {
                        // only consider cashflow payments between the first exercise date and swap maturity
                        // but then PV these to baseDate
                        IRSwapStart = exerSched->firstDate();
                    } else {
                        IRSwapStart = baseDate;
                    }
                }    
                // The IR swap doesn't start until the CVB coupons start accruing (datedDate)
                if (IRSwapStart < cvb->bond->getAccrualStartDate()) {
                    IRSwapStart = cvb->bond->getAccrualStartDate();
                }

                double floatLeg = getFloatLegPV(IRSwapStart);
                double fixedLeg = cvb->bond->couponsPV(IRSwapStart, swapMaturity, cvb->discount.getSP());
                floatLeg *= cvb->discount->pv(baseDate, IRSwapStart);
                fixedLeg *= cvb->discount->pv(baseDate, IRSwapStart);
                
                fixedLeg += swapBackEndFee*cvb->discount->pv(baseDate, swapMaturity);

                *strike = exerPct*swapNotional + floatLeg + fixedLeg;

                // calculate break of funds for Credit Suisse style swaps
                if ( ocbType == "C" ) {
                    double currentRate         = 0.0;
                    double currentDayCountFrac = 0.0;
                    int count;
                    string interval;
                    swapInterval->decompose(count, interval);

                    CashFlowArray spreadFlows = SwapTool::cashflows(
                        baseDate.rollDate(-1),
                        swapMaturity,
                        !shortFrontStub,
                        -spread,
                        count,           // interval = count periods
                        interval,          // e.g. Y, M, W, D
                        swapDCC.get());

                    if ( baseDate < swapLastFixDate ) {
                        throw ModelException(method, "last swap fixing date is in the future");
                    } else if ( spreadFlows.size() > 0 && spreadFlows[0].date == baseDate) {
                        currentRate         = 0.0;
                        currentDayCountFrac = 0.0;
                    } else if ( spreadFlows.size() > 0 ){
                        currentRate         = cvb->discount->fwd(baseDate,
                                                                 spreadFlows[0].date,
                                                                 swapDCC.get(),
                                                                 CompoundBasis::SIMPLE);
                        currentDayCountFrac = swapDCC->years(baseDate, spreadFlows[0].date);
                    } else {
                        currentRate         = 0.0;
                        currentDayCountFrac = 0.0;
                    }


                    if (swapLastFixRate > currentRate) {
                        double breakOfFunds;
                        breakOfFunds = (swapLastFixRate - currentRate + swapBreakOfFundsRate) * 
                            currentDayCountFrac                                    *
                            swapNotional;
                        *strike += breakOfFunds;
                    }


                }
            }

        } else if (ocbType == "F") {
            double exerPct;
            // see if it's exercisable today
            try { // need a try block as an error is thrown if it's not exercisable
                exerPct = exerSched->interpolate(baseDate);
                *isExercisable = true;
            } catch (exception& ) {
                *isExercisable = false;
                // we should still calculate the strike as if it were exercisable today for info purposes
                // use the exer level on the last day
                exerPct = exerSched->interpolate(exerSched->lastDate()); 
            }    
 
            // *strike = exerPct*cvb->bond->getFaceValue();
            *strike = exerPct;
        } else if (ocbType == "Y") {
            if ( baseDate <= swapMaturity ) {
                *isExercisable = true;
                // *strike        = assetSwapStrike->getStrike(baseDate, cvb, cvb->discount.getSP());
            } else {
                *isExercisable = false;
                *strike        = 0.0;
            }
        } else if (ocbType == "P") {
            // see if it's exercisable today
            try { // need a try block as an error is thrown if it's not exercisable
                *isExercisable = true;
            } catch (exception& ) {
                *isExercisable = false;
            }

            // make the risky curve
            IYieldCurveSP tmpCSC = CreditSpreadCurve(spread).makeRiskyCurve(*cvb->discount.get());
            YieldCurveSP riskyCurve(dynamic_cast<YieldCurve*>(tmpCSC.get()));

            /** value of all coupons paid between fromDate and to Date on toDate */
            // *strike = cvb->bond->couponsPV(baseDate, swapMaturity, riskyCurve);
            *strike = cvb->bond->couponsPV(baseDate, cvb->bond->getMaturityDate(), riskyCurve) + 
                      cvb->bond->redemptionPV(baseDate, riskyCurve);
        } else if (ocbType == "B") {
            if (baseDate > swapMaturity) {
                *isExercisable = false;
                *strike = 0.;

            } else {
                double exerPct;
                // see if it's exercisable today
                try { // need a try block as an error is thrown if it's not exercisable
                    exerPct = exerSched->interpolate(baseDate);
                    *isExercisable = true;
                } catch (exception& ) {
                    *isExercisable = false;
                }    
                *strike = bondFloor;
            }
        }

        if (includeLockout == true && *isExercisable == true) {
            *strike += lockoutFee(baseDate);
        }

        if (ocbType != "F" && ocbType != "Y") {
           // Adjust so that strike is calculated with respect to the strike at the option's maturity instead
           // of the bond face value
           *strike *= strikeAtOptionMat / cvb->bond->getFaceValue();
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double OptOnConvBond::lockoutFee(DateTime  baseDate) const {

    static const string method = "OptOnConvBond::lockoutFee";

    try {
        double lockout = -1e99;

        if ( lockOutParameters->lockoutType == "N" || 
             baseDate.getDate() > lockOutParameters->lockoutDate.getDate()) {
            lockout = 0.;
        } else if (lockOutParameters->lockoutType == "S") {
            lockout = lockOutParameters->lockoutRate * lockOutParameters->lockoutNotional * 
                      lockoutDCC->years(baseDate, lockOutParameters->lockoutDate);
        } else if (lockOutParameters->lockoutType == "P") {
            lockout = lockOutParameters->lockoutRate * 
                      lockOutParameters->lockoutNotional * 
                      lockoutDCC->years(baseDate, lockOutParameters->lockoutDate) * 
                      cvb->discount->pv(baseDate, lockOutParameters->lockoutDate);
        } else if (lockOutParameters->lockoutType == "F") {
            lockout = lockOutParameters->lockoutRate * lockOutParameters->lockoutNotional;
        }
                
        return lockout;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double OptOnConvBond::getFloatLegPV(DateTime  baseDate) const {
    static const string method = "OptOnConvBond::getFloatLegPV";

    try {
        double floatLegPV = 0.; 

        if (baseDate < swapMaturity) {
            int count;
            string interval;
            swapInterval->decompose(count, interval);
            CashFlowArray spreadFlows = SwapTool::cashflows(
                baseDate,
                swapMaturity,
                !shortFrontStub,
                -spread,
                count,           // interval = count periods
                interval,          // e.g. Y, M, W, D
                swapDCC.get());

            int i;
            if ( ocbType == "C" ) {
                int count;
                string interval;
                swapInterval->decompose(count, interval);

                CashFlowArray tmpSpreadFlows = SwapTool::cashflows(
                    baseDate.rollDate(-1),
                    swapMaturity,
                    !shortFrontStub,
                    -spread,
                    count,           // interval = count periods
                    interval,          // e.g. Y, M, W, D
                    swapDCC.get());

                if ( tmpSpreadFlows.size() > 0 && tmpSpreadFlows[0].date > baseDate) {
                    double currentRate         = cvb->discount->fwd(baseDate,
                                                             spreadFlows[0].date,
                                                             swapDCC.get(),
                                                             CompoundBasis::SIMPLE);
                    double currentDayCountFrac = swapDCC->years(baseDate, spreadFlows[0].date);

                    spreadFlows[0].amount +=
                        (currentRate - swapLastFixRate) * currentDayCountFrac;

                }
            }
            
            for (i=0; i<spreadFlows.size(); i++) {
                floatLegPV += spreadFlows[i].amount * cvb->discount->pv(baseDate, spreadFlows[i].date);
            }
            
            floatLegPV -= 1.;
            
            floatLegPV *= swapNotional;
        } else {
            floatLegPV = 0.;
        }

        return floatLegPV;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool OptOnConvBond::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc. */
    return cvb->avoidVegaMatrix(model);
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP OptOnConvBond::getSensitiveStrikes(OutputNameConstSP outputName,
                                                 const IModel*      model)
{
    return cvb->getSensitiveStrikes(outputName, model);
}

bool OptOnConvBond::priceDeadInstrument(CControl* control, CResults* results) const
{   
    static const string method = "OptOnConvBond::priceDeadInstrument";

    DateTime matDate = exerSched->lastDate();
    DateTime cvbMatDate = cvb->bond->getMaturityDate();
    bool   successful = false;


    if (valueDate > matDate) {
        results->storePrice(0., cvb->discount->getCcy());
        successful = true;
    } else if (valueDate == matDate && valueDate == cvbMatDate) { // need to fix for struck
            
        double convertPrice, bondValue;
        double  convRatio, convCash, parity, callLevel, putLevel;
        bool    isConvertible, isCallable, isPuttable, isAdjustedForAccrued;

        cvb->getConversionInfo(valueDate, &isConvertible, &convRatio, &convCash);
        if (isConvertible == true) {
            parity = cvb->asset->fwdValue(valueDate) * convRatio + convCash;
        } else {
            parity = 0.;
        }

        bondValue = cvb->bond->getRedemption() + cvb->getAccruedAtDate(valueDate);

        OptOnConvBond* ncThis = const_cast<OptOnConvBond*>(this);
        try {
            ncThis->strikeAtOptionMat = cvb->bond->getNotional(exerSched->lastDate());
        } catch (exception&) {
            ncThis->strikeAtOptionMat = cvb->bond->getFaceValue();
        }

        bool isHardCall;
        cvb->getCallLevel(valueDate, 0.0, &isCallable, &isAdjustedForAccrued, &isHardCall, &callLevel);
        cvb->getPutLevel(valueDate, &isPuttable, &putLevel);

        if ( isPuttable ) {
            /* when we put on the maturity date we always get the coupon */
            putLevel += (cvb->getCouponAtMat)?0.0:cvb->getAccruedAtDate(valueDate);

            bondValue = Maths::max(bondValue, putLevel);
        }

        double issuerValue = (isCallable)?Maths::min(callLevel, bondValue):bondValue;
        convertPrice = Maths::max(parity, issuerValue);

        bool isExercisable;
        double strike;
        getStrike(valueDate, bondValue, true, &isExercisable, &strike);

        cvb->recordOutputRequests(control, results, convertPrice, bondValue, true);
        double ocbValue = Maths::max(convertPrice-strike, 0.);
        results->storePrice(ocbValue, cvb->discount->getCcy());
        recordOutputRequests(control, results, convertPrice, bondValue, ocbValue);

        successful = true;
    } else {
        // a little validation
  //      if (riskFreeCoupons == true) {
  //          throw ModelException(method, "riskFreeCoupons == true only implemented for closed-form DECS");
  //      }
    }

    return successful;
}

void OptOnConvBond::recordOutputRequests(Control* control, Results* results, const double& convertPrice, 
                                         const double& convertBondFloor, double theoOCBValue) const
{

    static const string method = "OptOnConvBond::recordOutputRequests";
    OutputRequest* request = NULL;

    if ( control->isPricing() ) {

        // DELAY_PRICE -- don't bother as convert settlement is ignored anyway

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       exerSched->lastDate(),
                                       valueDate,
                                       cvb->asset.get());

        if (control->requestsOutput(OutputRequest::OPTION_ON_CONVERTIBLE_BOND_FLOOR, request)) {
            bool isExercisable;
            double bondFloor;

            getStrike(valueDate, 
                      convertBondFloor,
                      false, 
                      &isExercisable, 
                      &bondFloor);
            results->storeRequestResult(request, bondFloor);
        }

        if (control->requestsOutput(OutputRequest::OPTION_ON_CONVERTIBLE_LOCKOUT_FEE, request)) {
            results->storeRequestResult(request, lockoutFee(valueDate));
        }
        
        if (control->requestsOutput(OutputRequest::OPTION_ON_CONVERTIBLE_STRIKE, request)) {
            bool isExercisable;
            double strike;

            getStrike(valueDate, 
                      convertBondFloor,
                      true, 
                      &isExercisable, 
                      &strike);
            results->storeRequestResult(request, strike);
        }

        if (control->requestsOutput(OutputRequest::CONVERTIBLE_PRICE_DIRTY, request)) {
            results->storeRequestResult(request, convertPrice);
        }
        if (control->requestsOutput(OutputRequest::CONVERTIBLE_PRICE_CLEAN, request)) {
            results->storeRequestResult(request, convertPrice - cvb->getAccruedAtDate(valueDate));
        }
        if (control->requestsOutput(OutputRequest::CLEAN_PRICE, request)) {
            results->storeRequestResult(request, convertPrice - convertBondFloor);
        }
        if (control->requestsOutput(OutputRequest::THEO_PRICE, request)) {
            results->storeRequestResult(request, theoOCBValue);
        }
    }
}

double OptOnConvBond::intrinsicValue(ResultsSP results) const
{
    bool isExercisable;
    double bondFloor;
    double mtmCVBPrice;

    // get the (dirty) MTM price
    if (DerivativeAsset::TYPE->isInstance(cvbAsset.get())) {
        const DerivativeAsset* asset = dynamic_cast<const   DerivativeAsset*>(cvbAsset.get());
        try {
            mtmCVBPrice = asset->getMTM();
        } catch (exception& e) {
            throw ModelException(e, "OptOnConvBond::GetMarket"  );
        }
    } else {
        throw ModelException("OptOnConvBond::GetMarket", "Asset must be an AssetCVB");
    }

    double strike = 0.0;
    if (ocbType == "B") {
        IObjectConstSP ocbStrike  = 
            results->retrieveRequestResult(OutputRequest::OPTION_ON_CONVERTIBLE_STRIKE);
        if (CDouble::TYPE->isInstance(ocbStrike)){
            strike = (dynamic_cast<const CDouble&>(*ocbStrike)).doubleValue();
        } else {
            throw ModelException("OptOnConvBond::intrinsicValue", 
                                 "Need to have OPTION_ON_CONVERTIBLE_STRIKE calculated before requesting intrinsic value");
        } 
    }

    getStrike(valueDate, 
              strike,
              false, 
              &isExercisable, 
              &bondFloor);

    double mtmPrice = mtmCVBPrice - bondFloor;
    mtmPrice = Maths::max(mtmPrice, 0.);

    return mtmPrice;
}

/** IScaleOutputs method */
void OptOnConvBond::scaleOutputs(CControlSP control, ResultsSP unscaledResults)
{
    try {
         // retrieve clean price and store it as fair value
         const string ccy = unscaledResults->getCcyName();
         double fv        = unscaledResults->retrievePrice();
         unscaledResults->storePrice(fv * 100, ccy);
    }
    catch (exception& ) {
        // do nothing
    }

     // scale everything by 1/faceValue
     double scalingFactor = 1.0 /  cvb->bond->getFaceValue();

     unscaledResults->scalePostProcess(control, scalingFactor);
}

/** IIntrinsicMTM method */
void OptOnConvBond::calculateIntrinsic(CControlSP control, ResultsSP results)
{
    double mtmPrice = intrinsicValue(results);
    results->storePrice(mtmPrice, cvb->discount->getCcy());
}

ObjectArraySP OptOnConvBond::calculateFloatingCashFlows(const DateTime& baseDate) const
{
    static const string method = "OptOnConvBond::calculateFloatingCashFlows";

    ObjectArraySP outputMatrix;

    try {
        if (baseDate < swapMaturity) {
            int count;
            string interval;
            swapInterval->decompose(count, interval);

            DateTimeArraySP cashFlowDates = DateTimeArraySP(SwapTool::paymentDates(
                baseDate,
                swapMaturity,
                count,           // interval = count periods
                interval,          // e.g. Y, M, W, D
                !shortFrontStub));

            outputMatrix = ObjectArraySP(new ObjectArray(6));
            DateTimeArraySP dates(new DateTimeArray(cashFlowDates->size()));
            DoubleArraySP   libor(new DoubleArray(cashFlowDates->size()));
            DoubleArraySP   liborPlusSpread(new DoubleArray(cashFlowDates->size()));
            DoubleArraySP   interest(new DoubleArray(cashFlowDates->size()));
            DoubleArraySP   discFact(new DoubleArray(cashFlowDates->size()));
            DoubleArraySP   pv(new DoubleArray(cashFlowDates->size()));

            double df, yearFrac, amount;

            for (int i=0; i<cashFlowDates->size(); i++) {
                df = cvb->discount->pv(baseDate, (*cashFlowDates)[i]);
                (*dates)[i] = (*cashFlowDates)[i];
                if ( i == 0 ) {
                    (*libor)[i] = cvb->discount->fwd(
                                    baseDate,
                                    (*cashFlowDates)[i],
                                    swapDCC.get(),
                                    CompoundBasis::SIMPLE);
                       yearFrac = swapDCC->years(baseDate, (*cashFlowDates)[i]);
                } else {
                    (*libor)[i] = cvb->discount->fwd(
                                   (*cashFlowDates)[i-1], 
                                   (*cashFlowDates)[i],
                                   swapDCC.get(),
                                   CompoundBasis::SIMPLE);
                       yearFrac = swapDCC->years((*cashFlowDates)[i-1], (*cashFlowDates)[i]);
                }

                (*liborPlusSpread)[i] = (*libor)[i] + spread;

                amount = (*liborPlusSpread)[i] * yearFrac * swapNotional;
                (*interest)[i] = amount;
                (*discFact)[i] = df;

                (*pv)[i] = amount * df;
            }

            (*outputMatrix)[0] = dates;
            (*outputMatrix)[1] = libor;
            (*outputMatrix)[2] = liborPlusSpread;
            (*outputMatrix)[3] = interest;
            (*outputMatrix)[4] = discFact;
            (*outputMatrix)[5] = pv;

        } else {
            throw ModelException(method, "Cannot calculate floating leg for expired swaps");
        }

        return outputMatrix;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

ObjectArraySP OptOnConvBond::calculateFixedCashFlows(const DateTime& baseDate) const
{
    static const string method = "OptOnConvBond::calculateFixedCashFlows";

    ObjectArraySP outputMatrix;

    try {
        if (baseDate < swapMaturity) {
            CashFlowArraySP coupons = cvb->getAssetSwapLeg(swapMaturity, true);

            outputMatrix = ObjectArraySP(new ObjectArray(5));
            DateTimeArraySP dates(new DateTimeArray(coupons->size()));
            DoubleArraySP   rates(new DoubleArray(coupons->size()));
            DoubleArraySP   interest(new DoubleArray(coupons->size()));
            DoubleArraySP   discFact(new DoubleArray(coupons->size()));
            DoubleArraySP   pv(new DoubleArray(coupons->size()));

            double df;

            for (int i=0; i<coupons->size(); i++) {
                (*dates)[i] = (*coupons)[i].date;
                df = cvb->discount->pv(baseDate, (*coupons)[i].date);
                (*interest)[i] = (*coupons)[i].amount;
                (*discFact)[i] = df;
                (*pv)[i] = (*interest)[i] * df;

                if (BondParams::TYPE->isInstance(cvb->bond.get())) {
                    BondParams* bond = dynamic_cast<BondParams*>(cvb->bond.get());
                    (*rates)[i] = bond->couponPct;
                } else {
                    (*rates)[i] = 0.0;
                }

            }

            (*outputMatrix)[0] = dates;
            (*outputMatrix)[1] = rates;
            (*outputMatrix)[2] = interest;
            (*outputMatrix)[3] = discFact;
            (*outputMatrix)[4] = pv;

        } else {
            throw ModelException(method, "Cannot calculate fixed floating leg for expired swaps");
        }

        return outputMatrix;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

ObjectArraySP OptOnConvBond::calculateFixedAccrued(const DateTime& baseDate) const
{
    static const string method = "OptOnConvBond::calculateFixedAccrued";

    DateTime startDate;
    DateTime endDate;
    double   couponRate;
    double   accruedInterest;

    try {
        ObjectArraySP outputMatrix;
        if (baseDate < swapMaturity) {
            cvb->bond->getCouponsDetails(baseDate, startDate, endDate,couponRate, accruedInterest);

            outputMatrix = ObjectArraySP(new ObjectArray(4));
            (*outputMatrix)[0] = DateTimeSP( new DateTime(startDate));
            (*outputMatrix)[1] = DateTimeSP( new DateTime(endDate));
            (*outputMatrix)[2] = CDoubleSP(CDouble::create(couponRate));
            (*outputMatrix)[3] = CDoubleSP(CDouble::create(accruedInterest));

        } else {
            throw ModelException(method, "Cannot calculate accrued interest for expired swaps");
        }

        return outputMatrix;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// for reflection
OptOnConvBond::OptOnConvBond(): CInstrument(TYPE) {
    // initialize optional parameters
    ocbType = "D";
    spread = 0.;
    swapIntervalString      = "3M";
    swapNotional            = 100.;
    swapDCCString           = "act/360";
    swapBackEndFee          = 0.;
    swapLastFixRate         = 0.;
    swapBreakOfFundsRate    = 0.;
    bondFloor               = 0.0;;
    strikeAtOptionMat       = 0.0;
    unwindNewSwap           = true;
    shortFrontStub          = true;
}

bool OptOnConvBond::recurse(const CFieldConstSP& field,
                            const CClassConstSP& targetClass) const
{
    // this gets called as part of the tweaking and allows us to specify 
    // whether the fields within this class should be tweaked or not. The
    // target class indicates what is being shifted
    if ( field == assetField &&  targetClass != derivativeAssetClass ) {
        return false;
    }

    return true;
}

/** Returns the name of the instrument's discount currency.
 * This model does NOT have a yield curve for the local currency,
 * it assumes it is the same as the underlying. We should really
 * add the domestic yield curve and return its name here - by now
 * just return an empty string */
string OptOnConvBond::discountYieldCurveName() const {
    return "";
}


class OptOnConvBondHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(OptOnConvBond, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(FD1F::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(IRiskyPricer);
        IMPLEMENTS(IIntrinsicMTM);
        IMPLEMENTS(IScaleOutputs);
        IMPLEMENTS(ObjectIteration::IOverride);
        EMPTY_SHELL_METHOD(defaultOptOnConvBond);
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(cvbAsset, "the convertible asset that the option is on");
        FIELD(ocbType, "D, C, P, or F. D is standard");
        FIELD(exerSched, "exercise schedule");
        FIELD(spread, "spread for ocb");
        FIELD_MAKE_OPTIONAL(spread);
        FIELD(swapMaturity, "maturity of swap");
        FIELD_MAKE_OPTIONAL(swapMaturity);
        FIELD(swapIntervalString, "swap floating leg payment interval.");
        FIELD_MAKE_OPTIONAL(swapIntervalString);
        FIELD(shortFrontStub,     "true if front stub of the floating leg is short");
        FIELD_MAKE_OPTIONAL(shortFrontStub);
        FIELD(swapNotional, "swap notional");
        FIELD_MAKE_OPTIONAL(swapNotional);
        FIELD(swapDCCString, "swap floating leg day count convention");
        FIELD_MAKE_OPTIONAL(swapDCCString);
        FIELD(swapBackEndFee, "any back end fee applying to swap as percent of cvb face");
        FIELD_MAKE_OPTIONAL(swapBackEndFee);
        FIELD(unwindNewSwap,  "true = new IR swap to unwind, false = settle accrued interest");
        FIELD_MAKE_OPTIONAL(unwindNewSwap);
        FIELD(swapFixings, "the fixings of the underlying IR swap");
        FIELD_MAKE_OPTIONAL(swapFixings);
        FIELD(swapLastFixRate, "last rate at which swap floating leg fixed. Type C only");
        FIELD_MAKE_OPTIONAL(swapLastFixRate);
        FIELD(swapLastFixDate, "last date swap floating leg fixed. Type C only");
        FIELD_MAKE_OPTIONAL(swapLastFixDate);
        FIELD(swapBreakOfFundsRate, "swap break of funds rate. Type C only");
        FIELD_MAKE_OPTIONAL(swapBreakOfFundsRate);
        FIELD(lockOutParameters, "lock out information");
        FIELD_MAKE_OPTIONAL(lockOutParameters);        
        FIELD(assetSwapStrike, "asset swap strike object (currently only applicable to type Y asset swaps");
        FIELD_MAKE_OPTIONAL(assetSwapStrike);
        // FIELD(fixedYieldParameters, "fixed yield parameters (only applicable to type Y asset swaps");
        // FIELD_MAKE_OPTIONAL(fixedYieldParameters);
        FIELD(bondSettle, "Override the settlement of the bond underlying");
        FIELD_MAKE_OPTIONAL(bondSettle);
        // transient fields
        FIELD(swapDCC, "derived");
        FIELD_MAKE_TRANSIENT(swapDCC);
        FIELD(swapInterval, "derived");
        FIELD_MAKE_TRANSIENT(swapInterval)
        FIELD(lockoutDCC, "derived");
        FIELD_MAKE_TRANSIENT(lockoutDCC)
        FIELD(bondFloor, "derived");
        FIELD_MAKE_TRANSIENT(bondFloor);
        FIELD(strikeAtOptionMat, "strike at maturity of the option on convert");	
        FIELD_MAKE_TRANSIENT(strikeAtOptionMat)
        FIELD(cvb, "derived");
        FIELD_MAKE_TRANSIENT(cvb);
        FIELD_MAKE_TWEAKABLE(cvb);

        assetField           = clazz->getDeclaredField("cvbAsset");
        derivativeAssetClass = CClass::forName("DerivativeAsset");
    }

    static IObject* defaultOptOnConvBond(){
        return new OptOnConvBond();
    }
};

CClassConstSP const OptOnConvBond::TYPE = CClass::registerClassLoadMethod(
    "OptOnConvBond", typeid(OptOnConvBond), OptOnConvBondHelper::load);

bool  OptOnConvBondLoad() {
    return (OptOnConvBond::TYPE != 0);
   }


//---------------------------------------------------------------------------
//  Addin functions:
//  ASSET_SWAP_GET_FIXED_LEG
//  ASSET_SWAP_GET_FLOATING_LEG
//---------------------------------------------------------------------------


class AssetSwapAddin: public CObject{
    static CClassConstSP const TYPE;

    // input parameters
    OptOnConvBondSP     assetSwap;
    CMarketDataSP       market;

    static IObjectSP assetSwapGetFloatingLeg(AssetSwapAddin* params) {
        static const string routine = "AssetSwapAddin::assetSwapGetFloatingLeg";
        try {
            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            OptOnConvBondSP swap(copy(params->assetSwap.get()));

            model.getInstrumentAndModelMarket(params->market.get(), swap.get());

            // roll value date forward by 0 days - this is to populate
            // any samples which should be set now, but are not
            // populated yet, for example when running overnight grids
            // for instruments which have a SOD sample
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(swap);

            ObjectArraySP floatingLeg = swap->calculateFloatingCashFlows(params->assetSwap->getValueDate());
            
            return floatingLeg;
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP assetSwapGetFixedLeg(AssetSwapAddin* params) {
        static const string routine = "AssetSwapAddin::assetSwapGetFixedLeg";
        try {
            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            OptOnConvBondSP swap(copy(params->assetSwap.get()));

            model.getInstrumentAndModelMarket(params->market.get(), swap.get());

            // roll value date forward by 0 days - this is to populate
            // any samples which should be set now, but are not
            // populated yet, for example when running overnight grids
            // for instruments which have a SOD sample
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(swap);

            ObjectArraySP fixedLeg = swap->calculateFixedCashFlows(params->assetSwap->getValueDate());
            
            return fixedLeg;
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP assetSwapGetFixedAccrued(AssetSwapAddin* params) {
        static const string routine = "AssetSwapAddin::assetSwapGetFixedAccrued";
        try {
            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            OptOnConvBondSP swap(copy(params->assetSwap.get()));

            model.getInstrumentAndModelMarket(params->market.get(), swap.get());

            // roll value date forward by 0 days - this is to populate
            // any samples which should be set now, but are not
            // populated yet, for example when running overnight grids
            // for instruments which have a SOD sample
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(swap);

            ObjectArraySP fixedLeg = swap->calculateFixedAccrued(params->assetSwap->getValueDate());
            
            return fixedLeg;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    AssetSwapAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(AssetSwapAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAssetSwapAddin);
        FIELD(assetSwap, "asset swap object");
        FIELD(market, "market object");
        Addin::registerClassObjectMethod("ASSET_SWAP_GET_FLOATING_LEG",
                                         Addin::CONV_BOND,
                                         "Returns the floating leg of an asset swap",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)assetSwapGetFloatingLeg);

        Addin::registerClassObjectMethod("ASSET_SWAP_GET_FIXED_LEG",
                                         Addin::CONV_BOND,
                                         "Returns the fixed leg of an asset swap",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)assetSwapGetFixedLeg);

        Addin::registerClassObjectMethod("ASSET_SWAP_GET_FIXED_ACCRUED",
                                         Addin::CONV_BOND,
                                         "Returns the accrued interest of the fixed leg of an asset swap",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)assetSwapGetFixedAccrued);

    }

    static IObject* defaultAssetSwapAddin(){
        return new AssetSwapAddin();
    }
    
};
   
CClassConstSP const AssetSwapAddin::TYPE = CClass::registerClassLoadMethod(
    "AssetSwapAddin", typeid(AssetSwapAddin), load);



/** calculate the asset swap strike on a given date */
double FixedYieldParameters::getStrike(const DateTime&    baseDate,
                                       ConvBondConstSP    cvb,
                                       YieldCurveConstSP  discountCurve) const
{
    return getStrike(baseDate);
}

double FixedYieldParameters::getStrike(const DateTime&    baseDate) const
{
    double yield = getYield();
    double notional = 0.0;

    if ( yieldType == 0 || yieldType == 1) {
        BondParamsSP fixedYieldBond(new BondParams());

        fixedYieldBond->faceValue           = faceValue;
        fixedYieldBond->redemptionPct       = 1.0;
        fixedYieldBond->couponPct           = couponRate;
        fixedYieldBond->frequency           = couponFrequency;
        fixedYieldBond->maturityDate        = endDate;
        fixedYieldBond->datedDate           = startDate;
        fixedYieldBond->dayCountConvString  = "B30E/360";

        // it's NOT an OID, as we do not know the yield yet
        fixedYieldBond->isAnOID             = false;
        fixedYieldBond->yieldForOID         = false;

        // we should handle bad day conventions in some way and take in a holiday
        fixedYieldBond->badDayConvString    = "N";
        fixedYieldBond->endOfMonthAdj       = false;
        fixedYieldBond->eomIgnoreLeapYear   = true;

        fixedYieldBond->oddLastShort        = true;

        fixedYieldBond->validatePop2Object();
        fixedYieldBond->initialize();
        fixedYieldBond->isValid = true;

        // convert from zero coupon annualised
        DayCountConventionSP inDCC(new B30E360());
        double yieldToUse;
        if ( yieldType == 1) {
            yieldToUse = RateConversion::rateConvert(startDate,
                                        endDate,
                                        yield,
                                        inDCC.get(),
                                        CompoundBasis::ANNUAL,
                                        inDCC.get(), // correct ???
                                        couponFrequency);
        } else {
            yieldToUse = yield;
        }
        notional = fixedYieldBond->priceFromYield(yieldToUse, true, baseDate);
    } else if ( yieldType == 2 ) {
		double redemptionValue = faceValue;

        DayCountConventionSP accrualDCC(new B30E360());
		double dayFrac = accrualDCC->years(baseDate,endDate);
        notional = (redemptionValue * ( 1.0 + couponRate * dayFrac)) / (1.0 + yield * dayFrac);
    }

    return notional;
}

bool FixedYieldParameters::isExercisable(const DateTime& baseDate) const
{
    if ( baseDate >= startDate && baseDate <= endDate) {
        return true;
    } else {
        return false;
    }
}


FixedYieldParameters::FixedYieldParameters() : AssetSwapStrike(TYPE), couponRate(0.0), couponFrequency(0), issuePrice(0.0),
                                               faceValue(100.0), yieldOverride(0.0)
{
}

void FixedYieldParameters::validatePop2Object()
{
    static const string method = "FixedYieldParameters::validatePop2Object";
    if ( yieldType < 0 || yieldType > 2 ) {
        throw ModelException(method, "Yield type must be 0 (standard ISMA), 1 (zero coupon bond annualised) or 2 (Simple Yield)");
    }

    if ( startDate >= endDate ) {
        throw ModelException(method, "Start date (" + startDate.toString() + ") must be before end date (" +
            endDate.toString() + ").");
    } 

    if ( !Maths::isPositive(yieldOverride) && issuePrice <= 0.0 ) {
        throw ModelException(method, "Issue Price(" + Format::toString(issuePrice) + ") must be positive");
    }
}

double FixedYieldParameters::getYield() const
{
    double yield = 0.0;
    if (Maths::isPositive(yieldOverride)) {
        yield = yieldOverride;
    } else {
        if ( yieldType == 0 || yieldType == 1) {
            BondParamsSP fixedYieldBond(new BondParams());

            fixedYieldBond->faceValue           = faceValue;
            fixedYieldBond->redemptionPct       = 1.0;
            fixedYieldBond->couponPct           = couponRate;
            fixedYieldBond->frequency           = couponFrequency;
            fixedYieldBond->maturityDate        = endDate;
            fixedYieldBond->datedDate           = startDate;
            fixedYieldBond->dayCountConvString  = "B30E/360";

            // it's NOT an OID, as we do not know the yield yet
            fixedYieldBond->isAnOID             = false;
            fixedYieldBond->yieldForOID         = false;

            // we should handle bad day conventions in some way and take in a holiday
            fixedYieldBond->badDayConvString    = "N";
            fixedYieldBond->endOfMonthAdj       = false;
            fixedYieldBond->eomIgnoreLeapYear   = true;

            fixedYieldBond->oddLastShort        = true;

            fixedYieldBond->validatePop2Object();
            fixedYieldBond->initialize();
            fixedYieldBond->isValid = true;

            // this yield is the YTM based on a 30/360 day count fraction and the user defined
            // coupon frequency
            yield = fixedYieldBond->yieldToMaturity(issuePrice, true, startDate);

            if ( yieldType == 1 ) {
                DayCountConventionSP inDCC(new B30E360());
                // zero coupon annualised
                yield = RateConversion::rateConvert(startDate,
                                                    endDate,
                                                    yield,
                                                    inDCC.get(),
                                                    couponFrequency,
                                                    inDCC.get(), // correct ???
                                                    CompoundBasis::ANNUAL);
            }
        } else  if ( yieldType == 2 ) {
		    double redemptionValue = 100;

            DayCountConventionSP accrualDCC(new B30E360());
		    double dayFrac = accrualDCC->years(startDate,endDate);
		    yield = (couponRate * 100 + (redemptionValue - issuePrice) / dayFrac) / issuePrice;
        }
    }
    return yield;
}

class FixedYieldParametersHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(FixedYieldParameters, clazz);
        SUPERCLASS(AssetSwapStrike);
        EMPTY_SHELL_METHOD(defaultFixedYield);

        FIELD(couponRate,      "the coupon rate");
        FIELD(couponFrequency, "the coupon frequency");
        FIELD(issuePrice,      "the issue price");
        FIELD(faceValue,       "face value at maturity");
        FIELD_MAKE_OPTIONAL(faceValue);
        FIELD(startDate,       "the bonds start date");
        FIELD(endDate,         "the bonds end date");
        FIELD(yieldType,       "type of yield");
        FIELD(hols,            "holidays");
        FIELD_MAKE_OPTIONAL(hols);
        FIELD(yieldOverride,   "prespecified yield will override the input data");
        FIELD_MAKE_OPTIONAL(yieldOverride);


        Addin::registerConstructor("FIXED_YIELD_PARAMETERS",
                                   Addin::CONV_BOND,
                                   "Creates a handle to a fixed yield parameters object",
                                   FixedYieldParameters::TYPE);
    }

    static IObject* defaultFixedYield(){
        return new FixedYieldParameters();
    }
};


CClassConstSP const FixedYieldParameters::TYPE = CClass::registerClassLoadMethod(
    "FixedYieldParameters", typeid(FixedYieldParameters), FixedYieldParametersHelper::load);


class FixedYieldAddin: public CObject{
    static CClassConstSP const TYPE;

    // input parameters
    FixedYieldParametersSP     yieldParams;
    DateTime                   valueDate;

    static double getYield(FixedYieldAddin* params) {
        static const string routine = "FixedYieldAddin::getYield";
        try {
            return params->yieldParams->getYield();
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static double getNotional(FixedYieldAddin* params) {
        static const string routine = "FixedYieldAddin::getNotional";
        double notional = 0.0;
        try {
            notional = params->yieldParams->getStrike(params->valueDate);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
        return notional;
    }


    /** for reflection */
    FixedYieldAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(FixedYieldAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFixedYieldAddin);
        FIELD(valueDate,   "value date");
        FIELD(yieldParams, "fixed yield parameters");
        Addin::registerClassDoubleMethod("GET_FIXED_YIELD",
                                         Addin::CONV_BOND,
                                         "Returns the fixed yield",
                                         TYPE,
                                         (Addin::DoubleMethod*)getYield);


        Addin::registerClassDoubleMethod("GET_FIXED_YIELD_STRIKE",
                                         Addin::CONV_BOND,
                                         "Returns the fixed yield",
                                         TYPE,
                                         (Addin::DoubleMethod*)getNotional);
    }

    static IObject* defaultFixedYieldAddin(){
        return new FixedYieldAddin();
    }
};
   
CClassConstSP const FixedYieldAddin::TYPE = CClass::registerClassLoadMethod(
    "FixedYieldAddin", typeid(FixedYieldAddin), load);

LockOutParameters::LockOutParameters() : CObject(TYPE) 
{
}

void LockOutParameters::validate()
{
    static const string method = "LockOutParameters::validate";
    if (!(lockoutType == "N" || lockoutType == "S" || lockoutType == "P" || lockoutType == "F")) {
        throw ModelException(method, "lockoutType must be N, S, P, or F");
    }
}


class LockOutParametersHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(LockOutParameters, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLockOutParameters);
        FIELD(lockoutType,        "the lockout type");
        FIELD(lockoutNotional,    "the lockout notional");
        FIELD(lockoutDate,        "the lockout date");
        FIELD(lockoutRate,        "the lockout rate");
        FIELD(lockoutDCCString,   "the lockout day count convention");
    }

    static IObject* defaultLockOutParameters(){
        return new LockOutParameters();
    }
};

CClassConstSP const LockOutParameters::TYPE = CClass::registerClassLoadMethod(
    "LockOutParameters", typeid(LockOutParameters), LockOutParametersHelper::load);


class GrossSpreadAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    BondSP              bond;
    YieldCurveConstSP   yieldCurve;
    CreditSpreadCurveSP creditSpreadCurve;
    DateTime            valueDate;
    DateTime            workoutDate;
    double              workoutLevel;
    double              bondPrice;
    int                 frequency;
    DayCountConventionSP floatingDCC;
    double              swapPrice;
    double              fundingSpread;
    CMarketDataSP       market;

    static double calculateAssetSwapSpread(GrossSpreadAddin* params) {
        static const string routine = "GrossSpreadAddin::calculateAssetSwapSpread";
        try {
            CMarketDataSP       market              = params->market;

            // initialise the bond
            CClosedFormLN model("VolSurface");
            params->bond->getMarket(&model, market.get());

	    if (market.get())
            {
                // May need to fetch holidays for Business/252 DCC
                market->fetchForNonMarketObjects(params->floatingDCC, 
						 IModelConstSP(new NonPricingModel()),
						 "Business252");
            } 

            double spread = ConvBond::calculateAssetSwapSpread(
                                    params->bond,
                                    params->yieldCurve,
                                    params->valueDate,
                                    params->workoutDate,
                                    params->workoutLevel,
                                    params->bondPrice,
                                    params->frequency,
                                    true,
                                    params->floatingDCC);

            return spread;
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    GrossSpreadAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GrossSpreadAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGrossSpreadAddin);
        FIELD(bond,                 "Bond object");
        FIELD(yieldCurve,           "Risk free yield curve");
        FIELD(creditSpreadCurve,    "Credit spread curve");
        FIELD(valueDate,     "Value date");
        FIELD(workoutDate,   "Workout date");
        FIELD(workoutLevel,  "Workout level");
        FIELD(bondPrice,     "Bond price");
        FIELD(frequency,     "Frequency of floating leg");
        FIELD(floatingDCC,   "Day count convention of floating leg");
        FIELD(swapPrice,     "Swapped price");
        FIELD(fundingSpread, "spread required to fund the cash payments");
        FIELD(market,               "Market Data Cache");

        Addin::registerClassDoubleMethod("GET_GROSS_SPREAD",
                                         Addin::CONV_BOND,
                                         "Calculates an asset swap gross spread",
                                         TYPE,
                                         (Addin::DoubleMethod*)calculateAssetSwapSpread);

    }

    static IObject* defaultGrossSpreadAddin(){
        return new GrossSpreadAddin();
    }
    
};
   
CClassConstSP const GrossSpreadAddin::TYPE = CClass::registerClassLoadMethod(
    "GrossSpreadAddin", typeid(GrossSpreadAddin), load);


class AssetSwapSpreadAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    BondSP              bond;
    YieldCurveConstSP   yieldCurve;
    CreditSpreadCurveSP creditSpreadCurve;
    DateTime            valueDate;
    DateTime            workoutDate;
    double              workoutLevel;
    double              bondPrice;
    int                 frequency;  
    DayCountConventionSP floatingDCC;
    bool                stubAtFront;
    CMarketDataSP       market;

    static double calculateAssetSwapSpread(AssetSwapSpreadAddin* params) {
        static const string routine = "AssetSwapSpreadAddin::calculateAssetSwapSpread";
        try {
            CMarketDataSP       market              = params->market;

            // initialise the bond
            CClosedFormLN model("VolSurface");
            params->bond->getMarket(&model, market.get());

	    if (market.get())
            {
                // May need to fetch holidays for Business/252 DCC
                market->fetchForNonMarketObjects(params->floatingDCC, 
						 IModelConstSP(new NonPricingModel()),
						 "Business252");
            } 

            double aswSpread = ConvBond::calculateAssetSwapSpread(params->bond,
                                                                  params->yieldCurve,
                                                                  params->valueDate,
                                                                  params->workoutDate,
                                                                  params->workoutLevel,
                                                                  params->bondPrice,
                                                                  params->frequency,
                                                                  params->stubAtFront,
                                                                  params->floatingDCC);
        
            return aswSpread;
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    static double calculateZSpread(AssetSwapSpreadAddin* params) {
        static const string routine = "AssetSwapSpreadAddin::calculateZSpread";
        try {
            CMarketDataSP       market              = params->market;

            // initialise the bond
            CClosedFormLN model("VolSurface");
            params->bond->getMarket(&model, market.get());

	    if (market.get())
            {
                // May need to fetch holidays for Business/252 DCC
                market->fetchForNonMarketObjects(params->floatingDCC, 
						 IModelConstSP(new NonPricingModel()),
						 "Business252");
            } 

            double zSpread = ConvBond::calculateZSpread(params->bond,
                                                        params->yieldCurve,
                                                        params->valueDate,
                                                        params->workoutDate,
                                                        params->workoutLevel,
                                                        params->bondPrice,
                                                        params->frequency,
                                                        params->stubAtFront,
                                                        params->floatingDCC);
        
            return zSpread;
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    AssetSwapSpreadAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(AssetSwapSpreadAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAssetSwapSpreadAddin);
        FIELD(bond,                 "Bond object");
        FIELD(yieldCurve,           "Risk free yield curve");
        FIELD(creditSpreadCurve,    "Credit spread curve");
        FIELD(valueDate,     "Value date");
        FIELD(workoutDate,   "Workout date");
        FIELD(workoutLevel,  "Workout level");
        FIELD(bondPrice,     "Bond price");
        FIELD(frequency,     "Frequency of floating leg");
        FIELD(floatingDCC,   "Day count convention of floating leg");
        FIELD(stubAtFront,   "true if short stub is at the front");
        FIELD(market,               "Market Data Cache");

        Addin::registerClassDoubleMethod("GET_ASSET_SWAP_SPREAD",
                                         Addin::CONV_BOND,
                                         "Calculates an asset swap spread",
                                         TYPE,
                                         (Addin::DoubleMethod*)calculateAssetSwapSpread);

        Addin::registerClassDoubleMethod("GET_Z_SPREAD",
                                         Addin::CONV_BOND,
                                         "Calculates an Z-Spread (spread over yields not over zero-rates)",
                                         TYPE,
                                         (Addin::DoubleMethod*)calculateZSpread);
    }

    static IObject* defaultAssetSwapSpreadAddin(){
        return new AssetSwapSpreadAddin();
    }
    
};
   
CClassConstSP const AssetSwapSpreadAddin::TYPE = CClass::registerClassLoadMethod(
    "AssetSwapSpreadAddin", typeid(AssetSwapSpreadAddin), load);

DRLIB_END_NAMESPACE

