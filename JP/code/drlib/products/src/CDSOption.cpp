//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSOption.cpp
//
//   Description : Option on a CDS
//
//   Author      : Charles Morcom
//
//   Date        : 12 December 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/CDSOption.hpp"
#include "edginc/IDiscountCurve.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/CDSVolRequestSimpleEuropean.hpp"
#include "edginc/CDSVolProcessedSimpleEuropean.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/ICDS.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CredDefSwap.hpp"

DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  CDSOptionLoad() {
    return (CDSOption::TYPE != 0);
   }

const string CDSOption::EXERCISE_TYPE_EUROPEAN = "EXERCISE_TYPE_EUROPEAN";
const string CDSOption::EXERCISE_TYPE_AMERICAN = "EXERCISE_TYPE_AMERICAN";
const string CDSOption::EXERCISE_TYPE_BERMUDAN = "EXERCISE_TYPE_BERMUDAN";
const string CDSOption::STRIKE_TYPE_SPREAD = "STRIKE_TYPE_SPREAD";
const string CDSOption::STRIKE_TYPE_CLEAN_PRICE = "STRIKE_TYPE_CLEAN_PRICE";
const string CDSOption::STRIKE_TYPE_DIRTY_PRICE = "STRIKE_TYPE_DIRTY_PRICE";

CDSOption::~CDSOption(){}

/* Default constructor */
CDSOption::CDSOption() : CInstrument(TYPE), isExercised(false) {
};

string CDSOption::discountYieldCurveName() const {
    return underlying->getYieldCurveWrapper().getName();
}

/** Do some asset specific validation */
void CDSOption::Validate() {
        
    static const string method("CDSOption::Validate");

    /*=========================================================================
     * CHECK EXERCISE TYPE AND SCHEDULE
     *=======================================================================*/
    /* Exercise type must be OK */
    if (exerciseType!=CDSOption::EXERCISE_TYPE_AMERICAN
        && exerciseType!=CDSOption::EXERCISE_TYPE_EUROPEAN
        && exerciseType!=CDSOption::EXERCISE_TYPE_BERMUDAN) {
            throw ModelException(method, "exerciseType must have one of the "
                "following three values: \"EXERCISE_TYPE_AMERICAN\", "
                "\"EXERCISE_TYPE_EUROPEAN\", or \"EXERCISE_TYPE_BERMUDAN\".");
    }
    /* Check that there is an exercise schedule and that it is not empty */
    if (!exerciseSchedule || exerciseSchedule->getDateArray().size()<=0) {
        throw ModelException(method, 
            "There must be a non-empty exercise schedule.");
    } 
    /* Check that all dates are increasing */
    int i=0;
    DateTime lastDt;
    for (i=0; i<exerciseSchedule->getDateArray().size(); i++) {
        DateTime dt = exerciseSchedule->getDateArray()[i];
        if (i>0) {
            if (dt<=lastDt) {
                throw ModelException(method,"Exercise schedule "
                    "dates must be in strictly increasing order.");
            }
            lastDt=dt;
        }
    }
    /* Check that interpolation type is OK */
    if (exerciseType==CDSOption::EXERCISE_TYPE_AMERICAN &&
        exerciseSchedule->getInterp()==Schedule::INTERP_NONE) {
        throw ModelException(method, 
            "Option with American exerciseType must have schedule "
            "interpolation type \"INTERP_STAIRS\" or \"INTERP_LINEAR\". "
            "\"INTERP_NONE\" is not meaningful.");
    }
    /*=========================================================================
     * FOR NOW - ONLY ALLOW EUROPEAN OPTIONS, BUT REMOVE THIS LATER!
     *=======================================================================*/
    if (exerciseSchedule->getDateArray().size()!=1 || exerciseType!=CDSOption::EXERCISE_TYPE_EUROPEAN) {
        throw ModelException(method, 
            "Only european options have been implemented so far: the type must be "
            "EXERCISE_TYPE_EUROPEAN, and there must be exactly one exercise date/strike "
            " in the schedule.");
    }

    /*=========================================================================
     * CHECK STRIKE TYPE
     *=======================================================================*/
        if (strikeType!=CDSOption::STRIKE_TYPE_SPREAD &&
            strikeType!=CDSOption::STRIKE_TYPE_CLEAN_PRICE &&
            strikeType!=CDSOption::STRIKE_TYPE_DIRTY_PRICE) {
                throw ModelException(method, 
                    "The strikeType must have one of three values: "
                    "\"STRIKE_TYPE_SPREAD\", \"STRIKE_TYPE_CLEAN_PRICE\", "
                    "or \"STRIKE_TYPE_DIRTY_PRICE\"");
        }

    /*=========================================================================
     * IF IT'S AN INDEX OPTION, IT MAY NOT KNOCK-OUT ON DEFAULT!
     *=======================================================================*/
    if (isIndexOption && koOnDefault) {
        throw ModelException(method, "A CDSOption on an index underlying may "
            "not knock out on default.");
    }

    /*=========================================================================
     * IF IT'S AN INDEX SPREAD-OPTION, IT MAY NOT BE PHYSICALLY SETTLED!
     * This is because of the BB CDSW unwind value convention
     *=======================================================================*/
    if (isIndexOption && strikeType==CDSOption::STRIKE_TYPE_SPREAD && 
        !isCashSettled) {
        throw ModelException(method, "A CDSOption on an index underlying "
            "with a spread strike may not be physically settled.");
    }

    /*========================================================================
     * The possible underlying CDS instruments have different conventions
     * so here we force short protection, long fees upon the underlying
     * Also assign the usingCDS flag which will be used during pricing
     *========================================================================*/
    // make sure that the notional of the underlying is -1.0
    // so that the underlying is long risk/short protection
    if (underlying->getNotional()!=-1.0) {
        underlying->setNotional(-1.0);
    }
    //validate the underlying
    underlying->Validate();
}


/** Get the asset , par curve and discount market data */
void CDSOption::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);

    /*=========================================================================
     * GET THE DISCOUNT CURVE 
     *=======================================================================*/
    discount = underlying->getYieldCurveWrapper();
    discount.getData(model, market);
    //    rfCrv = YieldCurveSP(discount.get());

    /*=========================================================================
     * GET THE CREDIT SPREAD CURVE
     *=======================================================================*/
   cdsParSpreads = underlying->getParSpreadsWrapper();
//     ICDSParSpreads::getMarketData(
//         model,
//         market.get(),
//         discount.getName(),
//         cdsParSpreads);
//     crv = ICDSParSpreadsSP(cdsParSpreads.get());
    cdsParSpreads.getData(model,market);
    
    /*=========================================================================
     * SORT OUT EXERCISE AND PREMIUM SETTLEMENT
     *=======================================================================*/
    if (!premiumSettlement) {
        // default is settle same day as exercise T + 0 - Mehdi
      //This relates to settling the cash on the CDS upon exercise.
      //It should be T + 3 but we are defaulting it to T + 0 as it may depend whether
      // the option is spread or price based. User can override this.
      premiumSettlement = SettlementSP(new RollingSettlement(0,
							     cdsParSpreads.get()->getHolidays()));
    }
    premiumSettlement->getMarket(model, market.get());
    
    // exerciseSettlement relates to settling the underlying upon exercise.
    // For a CDS, this means the start date of the CDS: T +1D
    // Note that only the flows that are on or after this date
    // will be taken into account when calc the forward annuity.
    if (!exerciseSettlement) {
        exerciseSettlement = SettlementSP(new RollingSettlement(1,
            HolidayConstSP(Holiday::noHolidays())));
    }
    exerciseSettlement->getMarket(model, market.get());

    /*=========================================================================
     * ENSURE MARKET INITIALIZATION FOR THE UNDERLYING
     *=======================================================================*/
    underlying->GetMarket(model,market);

    // call instrument specific getMarket routine if applicable
    if (IGetMarket::TYPE->isInstance(this)) {
        IGetMarket* imnt = dynamic_cast<IGetMarket*>(this);
        imnt->getMarket(model,market.get());
    }
}

DateTime CDSOption::getValueDate() const {
    return valueDate;
}

/** Find next exercise and strike after valueDate, or dateOverride if not null. Return
    true if found exercise, else false if expired.*/
bool CDSOption::findNextExercise(
                                 const DateTime*    dateOverride, 
                                 DateTime*          exDt, 
                                 double*            strk
                                 ) const {

    const static string method = "CDSOption::findNextExercise";

    DateTime dt = (dateOverride ? *dateOverride : valueDate);

    // find index of next date in schedule <= valueDate
    const DateTimeArray& dta = exerciseSchedule->getDateArray();
    const DoubleArray& vla = exerciseSchedule->getValueArray();
    if (dta.size()<=0) {
        // this has been checked in Validate(), but I'm being careful.
        throw ModelException(method, "Empty exercise schedule!");
    }
    int i = 0;
    // Slow! Ugly! Naughty Charles!
    while (i<dta.size() && dta[i]<dt) {
        i++;
    }
    if (i>=dta.size()) {
        // you are past the last exercise - return false to indicate expiry.
        return false;
    }
    
    if (this->exerciseType==CDSOption::EXERCISE_TYPE_EUROPEAN ||
        this->exerciseType==CDSOption::EXERCISE_TYPE_BERMUDAN ||
        i==0 /* If start of American exercise schedule still forward */) {

        // no interpolation - just find next date.
        *exDt = dta[i];
        *strk = vla[i];

    } else {

        // American - you must interpolate
        *exDt = dt; // it's exercisable now for sure
        if (exerciseSchedule->getInterp()==Schedule::INTERP_STAIRS) {
            // strike is same as last one
            *strk = vla[i];
        } else if (exerciseSchedule->getInterp()==Schedule::INTERP_LINEAR) {
            const DateTime& lastDate = dta[i-1];
            const DateTime& nextDate = dta[i];
            double dateIvl = lastDate.yearFrac(nextDate);
            if (dateIvl<=0) {
                // this has been checked, but just to be sure...
                throw ModelException(method,
                    "Exercise schedule with disordered or indistinct dates: not a good idea!");
            }
            double nextStrike = vla[i];
            double lastStrike = vla[i-1];
            *strk = lastStrike + 
                (nextStrike-lastStrike)*lastDate.yearFrac(valueDate)/dateIvl;
        } else {
            // this has been checked, but just to be sure...
            throw ModelException(method, 
                "Option with American exerciseType must have schedule "
                "interpolation type \"INTERP_STAIRS\" or \"INTERP_LINEAR\". "
                "\"INTERP_NONE\" is not meaningful.");
        }
    }
    return true;
}

/** Prices option in closed form using par-spread curve provided and and a
    forward price/spread distribution at exercise defined by the model. 
    Note that this prices to the next (or current)
    exercise date: if the option is American or Bermudan, later options will be ignored
    silently. Should this throw an exception?*/
void CDSOption::priceClosedForm(CResults* results, Control* control, const Model* model) const {
    static const string method = "CDSOption::priceClosedForm";

    /**************************************************************************
     * TODO: There seems to be no way of tracking index defaults at present in 
     * CDSIndexParSpreads. For now, just handle crv.defaulted(). Will need to
     * add code to handle defaults in indices and strike 
     * adjustments, etc. later, when default tracking has been implemented
     * Charles Morcom, February 7, 2006
     *************************************************************************/

    try {
        // retrieve the model to be used in the calculation
        // of fee leg forward rates
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>
            (const_cast<Model*>(model));
        if (!ihfrp)
        {
            throw ModelException(method,
                "Model must implement IHasForwardRateModel");
        }
        IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

        // retrieve the pre-payment curve
        IDecretionCurveConstSP prepay = cdsParSpreads->getPrepayCurve();

        OutputRequest* request = 0;

        /* FOR NOW, isIndexOption IS A PROVIDED FIELD: MAYBE IT SHOULD BE READ FROM THE CURVE? */
        const string& ccy = cdsParSpreads->getCcy();
        bool isDefaulted = cdsParSpreads->defaulted();
        bool strikeIsSpread = (strikeType==CDSOption::STRIKE_TYPE_SPREAD);

        /*=====================================================================
         * FIND THE RELEVANT EXERCISE DATE AND STRIKE
         *===================================================================*/
        DateTime nextEx;
        double nextStrk;
        if (isExercised) {
            nextEx = exerciseDate;
            DateTime dummyDate; // should be same as exerciseDate!
            findNextExercise(&exerciseDate, &dummyDate, &nextStrk);
        } else {
            // check schedule to find next exercise
            if (!findNextExercise(0, &nextEx, &nextStrk)) {
                // expired and not exercised
                results->storePrice(0.0, ccy);
                return;
            }
        }
        request = control->requestsOutput(OutputRequest::EXPECTED_STRIKE);
        if (request) {
            results->storeRequestResult(request, nextStrk);
        } 
        request = control->requestsOutput(OutputRequest::EXPECTED_TIME_TO_EXERCISE);
        if (request) {
            results->storeRequestResult(request, valueDate.yearFrac(nextEx));
        } 

        // Forward CDS will start on this date
        DateTime ulSettlementDate = exerciseSettlement->settles(nextEx); // T+1
        // Cash paid on exercise on this date
        DateTime payDate = premiumSettlement->settles(nextEx); // T+0
	
        /* Record payment date if cash-settled */
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && isCashSettled) {
            DateTimeArray dates(1, payDate);
            OutputRequestUtil::recordPaymentDates(control, results, &dates); 
        }
        double dfToExPmt = discount->pv(payDate); // risk-free discount factor from valueDate
        double dfToExercise = discount->pv(nextEx);
        double dfToULSettlement = discount->pv(ulSettlementDate);
        double value = 0.0;

        //DISCOUNT_FACTOR
        request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
        if (request) {
            results->storeRequestResult(request, dfToExercise);
        }

        /*=====================================================================
         * CREATE CONCRETE CDS THAT RESULTS FROM EXERCISE
         *===================================================================*/
        ICDSSP cds = underlying;
        IBadDayAdjuster* bdAdj = (IBadDayAdjuster*)dynamic_cast<IBadDayAdjuster*>(underlying.get());
        if (!bdAdj)
        {
            throw ModelException(method, "Underlying must implement IBadDayAdjuster");
        }
        IBadDayAdjusterConstSP bda(IBadDayAdjusterConstSP::attachToRef(bdAdj));
        
        ICreditFeeLegSP undFeeLeg = underlying->getFeeLeg();
        // fee leg must be of a fixed rate variety for this to work
        IFixedRateCreditFeeLeg* frcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(undFeeLeg.get());
        if (!frcfl)
        {
            throw ModelException(method, "Fee leg is not of a fixed fee variety "
                                         "(must implement IFixedRateCreditFeeLeg)");
        }

        /* The fee is always the next strike, unless the spread is fixed as in,
           index options and price-struck single-name options. */
        double fee;
        if (isIndexOption || !strikeIsSpread) {
            fee = frcfl->getRate();
        } else {
            fee = nextStrk;
            frcfl->setRate(fee);
        }

        /* If not a fixed underlying, create underlying with correct dates */
        /* Is there a better way to do this than with dynamic_cast<>? */
        /* TODO: This seems not to work. The generateCDS() method in CredDefSwap seems
         * to work fine (even for several tweaks), but invariably, there is a core dump later.
         * this feels like a memory deletion problem, perhaps, but my debugger will not work,
         * so I have been unable to check this. For now I have just commented out the case
         * where I regenerate the underlying for single-name spread options. For European
         * options, it is always possible to book the underlying correctly, so no functionality
         * is lost.
         * [Charles Morcom February 18, 2006] */
        if (ulMaturityDefinedByOption) {
            ICDSConvention* changeableCDS = dynamic_cast<ICDSConvention*>(cds.get());
            if (!changeableCDS) {
                throw ModelException(method,
                    "The underlying must be castable to an ICDSConvention.");
            }
            cds = changeableCDS->generateCDS(ulSettlementDate,
                ulExpiry->toDate(ulSettlementDate),
                fee);
        }
        /*------------------> BEGIN SECTION REMOVED BECAUSE OF UNSOLVED CORE-DUMPS */
        //} else if (!isIndexOption && strikeIsSpread) {
        //    /* A single-name CDS spread option is relative to an underlying
        //       issued on the settlement date */
        //    if (!cds->hasFiniteMaturity()) {
        //        throw ModelException(method,
        //            "The underlying must have a finite maturity!");
        //    }
        //    ICDSWithConvention* changeableCDS = dynamic_cast<ICDSWithConvention*>(cds.get());
        //    if (!changeableCDS) {
        //        throw ModelException(method,
        //            "The underlying must be castable to an ICDSWithConvention.");
        //    }
        //    cds = changeableCDS->generateCDS(ulSettlementDate,
        //        cds->getMaturity(),
        //        fee);
        //}
        /*------------------> END SECTION REMOVED BECAUSE OF UNSOLVED CORE-DUMPS*/

        /*=====================================================================
         * ADJUST STRIKE FOR ACCRUED INTEREST, IF NECESSARY
         *===================================================================*/
        /* If it's a dirty-price strike, you need to adjust it to a clean strike */
        //negative fees lead to -ve ai
        double accruedInterest = abs(cds->getAccruedInterest(ulSettlementDate, frModel));
        if (strikeType==CDSOption::STRIKE_TYPE_DIRTY_PRICE) {
            // note, notional is -1, so sign off AI is POSITIVE!
            nextStrk -= accruedInterest;
        }
        
        /*=====================================================================
         * CHECK IF EXERCISED ALREADY, AND COMPUTE VALUE
         *===================================================================*/
        if (isExercised) {
            if (isCashSettled) {
                // cash settlement
                double cashFlow;
                if (strikeIsSpread) {
                    cashFlow = max(0.0,(isCall ? exerciseValue : -exerciseValue));
                } else {
                    // if call: cost to exercise is -F, since F is paid by protection
                    // _buyer_ and call is right to sell protection. Strike is K=1-F
                    // so -F = K-1
                    cashFlow = max(0.0,(isCall 
                                        ? exerciseValue-(nextStrk+accruedInterest-1.0) 
                                        : -exerciseValue+(nextStrk+accruedInterest-1.0)) );
                }
                
                request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                if (request) {
                    OutputRequestUtil::KnownCashFlows kcf;
                    kcf.addKnownCashFlow(ccy, payDate, cashFlow*notional);
                    kcf.recordKnownCashFlows(control, results);
                }
                results->storePrice(cashFlow*notional*dfToExPmt, ccy);
                return;
            } else {
                // cash payment and
                // ul CDS value computed below
            }
        }
        
        if (isExercised && !isCashSettled && !strikeIsSpread) {
            // if physically settled, exercised, and price-struck, 
            // then there is an upfront cash-payment for the underlying
            // if call, buy risk for positive notional
            value += notional*(nextStrk + accruedInterest - 1.0);
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                OutputRequestUtil::KnownCashFlows kcf;
                kcf.addKnownCashFlow(ccy, payDate, notional*(nextStrk + accruedInterest -1.0));
                kcf.recordKnownCashFlows(control, results);
            }
        }

        /*=====================================================================
         * COMPUTE UL VALUE GIVEN DEFAULT FVd TO ulSettlementDate AND, HENCE,
         * THE OPTION VALUE GIVEN DEFAULT. YOU NEED 
         * THIS FOR DEFAULTED OPTIONS, AS WELL AS UNDEFAULTED OPTIONS
         *===================================================================*/
        double ulVGD;
        if (payAtDefault && !isDefaulted && 
            (ulMaturityDefinedByOption || (!isIndexOption && strikeIsSpread))) {
            /* Note that this breaks if the underlying is digital!*/
            ulVGD = -cdsParSpreads->protectionPV(
                valueDate, 
                valueDate, 
                ulSettlementDate, 
                IDiscountCurveRisky::RECOVER_1_MINUS_R,
                0);
            ulVGD  = ulVGD/((1-cdsParSpreads->survivalProb(valueDate, nextEx)));

        } else if (payAtDefault && !isDefaulted) {
            /* Default payment occurs whenever the default occurs */
            double sp = cdsParSpreads->survivalProb(valueDate, nextEx);
            ICreditContingentLegSP ctgLeg = cds->getContingentLeg();
            ulVGD = 
	      - ctgLeg->getContingentLegPV(nextEx, ulSettlementDate, *(cdsParSpreads.get()), bda)*sp
	      + ctgLeg->getContingentLegPV(valueDate, ulSettlementDate, *(cdsParSpreads.get()), bda);
            // what you have is the unconditional default value FVd to ulSettlementDate
            // if you divide this by the default probability, then you get the
            // FV conditional on default
            ulVGD  = ulVGD/(1-sp);

        } else {
            /* Default payments are as if exercise straight into defaulted underlying
               at the exercise date */
            ulVGD = cds->getPVGivenDefault(ulSettlementDate, ulSettlementDate, *(cdsParSpreads.get()), frModel);
        }
        double optionVGD;
        if (koOnDefault) {
            // no value at all if knock-out before default
            optionVGD = 0.0;
        } else if (strikeIsSpread) {
            // no fee to pay: just (1-R) payment
            optionVGD = max(0.0, (isCall ? 1.0 : -1.0) * ulVGD);
        } else {
            optionVGD = max(0.0, (isCall ? 1.0 : -1.0) * (ulVGD - (1 - nextStrk - accruedInterest)));
        }
        request = control->requestsOutput(OutputRequest::RECOVERY_VALUE);
        if (request) {
            results->storeRequestResult(request, notional*optionVGD);
        } 

        /*=====================================================================
         * SPECIAL HANDLING IF THERE HAVE BEEN DEFAULTS
         *===================================================================*/
        if (isDefaulted) {
                /*=============================================================
                 * INDEX DEFAULTS
		 * Disabled the seperation between indices and single names
		 * Indices dealt with as single names - Mehdi 08Oct06
                 *===========================================================*/	  
//             if (isIndexOption) {

//                 // TODO: implement default event handling for index options
//                 throw ModelException(method, 
//                     "Defaults for index options not handled yet!");
            
//             } else {
                /*=============================================================
                 * DEFAULTED SINGLE-NAME OPTION
                 *===========================================================*/
                if (isExercised) {
                    // Exercised with physical settlement already 
                    // note that the whole underlying is taken, since
                    // the option payment is already included in value
                    value += dfToExPmt*ulVGD;
                    results->storePrice(value*notional, ccy);
                    return;

                } else if (koOnDefault) {
                    // no value if KO and defaulted before exercised.
                    results->storePrice(0.0, ccy);
                    return;

                } else {
                    /* NO KO - PAY VALUE AT NEXT EXERCISE, ASSUMING
                       EXERCISE IS RATIONAL */
		  /* The below was not multiplied by the notional -
		     Fixed on 5Oct06 - Mehdi*/
                    results->storePrice(optionVGD*dfToExPmt*notional, ccy);
                    return;
                }
		// }
        }

        /*=====================================================================
         * IF EXERCISED AND PHYSICAL SETTLEMENT, THEN ADD THE VALUE OF THE 
         * UNDERLYING
         *===================================================================*/
        if (isExercised) {
            results->storePrice(dfToExPmt*cds->getPV(valueDate, ulSettlementDate, *(cdsParSpreads.get()), prepay, frModel, bda), ccy);
            return;
        }

        /*=====================================================================
         * IF HERE, NOT EXERCISED, NOT EXPIRED, AND NOT DEFAULTED, 
         * SO nextEx>=valueDate
         *===================================================================*/
        double survivalProbToEx = cdsParSpreads->survivalProb(valueDate, nextEx);
        // IMPLIED_DEFAULT_PROBABILITY
        request = control->requestsOutput(OutputRequest::IMPLIED_DEFAULT_PROBABILITY);
        if (request) {
            results->storeRequestResult(request, 1-survivalProbToEx);
        }

        /*=====================================================================
         * CALCULATE THE PROTECTION AND FEE FORWARDS
         *===================================================================*/
        /* These forwards are for unconditional settlement at ulSettlementDate, given
           no default by the exercise date; actual payments occur at payDate */
        ICreditContingentLegSP ctgLeg = cds->getContingentLeg();
        DateTime ulContingentStart = DateTime(ulSettlementDate.getDate(),0);
        double fwdProtVal = 
            ctgLeg->getContingentLegPV(ulContingentStart,
                                       ulContingentStart,
                                       *(cdsParSpreads.get()), bda);

        ICreditFeeLegSP feeLeg = cds->getFeeLeg();
        DateTime latestRiskyDate = ctgLeg->lastObservationEndDate();
        DateTime earliestRiskyDate = ctgLeg->firstObservationStartDate();

        //calculate partial accrual, and substract it from fwdFeeVal
        //valued as of nextEx date, conditional upon survival to ulSettlementDate
        //negative notional results in negative results
        double fwdFeeVal = abs(feeLeg->getFeeLegPV(
                                    nextEx, ulSettlementDate,
                                    earliestRiskyDate, latestRiskyDate,
                                    *(discount.getSP()), *(cdsParSpreads.get()), prepay,
                                    true, cds->getAccrualDcc(), false, frModel));

        fwdFeeVal -= accruedInterest;

        // fee leg must be of a fixed rate variety for this to work
        IFixedRateCreditFeeLeg* cdsfrcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(feeLeg.get());
        double fwdAnnuityVal = fwdFeeVal/cdsfrcfl->getRate();
        request = control->requestsOutput(OutputRequest::CDS_RISKY_DURATION);
        if (request) {
            results->storeRequestResult(request, fwdAnnuityVal);
        }
        if (fwdFeeVal<=0) {
            // note we are short protection with notional 1, so annuity val must be >0!
            throw ModelException(method,
                "Annuity for CDS option should be positive "
                "(i.e. call is right to buy risk), but is negative\n");
        }
        // clean forward
        double conditionalFwdPV = (fwdFeeVal + fwdProtVal) * survivalProbToEx;
        // need front protection value to compute unconditional settlement
        // adjustment for index options
        double frontProtVal = 0.0;
        if (isIndexOption) {
            frontProtVal = ulVGD*(1-survivalProbToEx);
            // remember that ulVGD = unco front protection value / default probability
        }
        /* This is unconditional if index swap, else conditional */
        double forwardSpread = -(fwdProtVal + frontProtVal/survivalProbToEx)/fwdAnnuityVal;
        /* This is the forward value of the CDS per notional, from the perspective
           of the protection SELLER. */
        request = control->requestsOutput(OutputRequest::FWD_UNDERLYING);
        if (request) {
            results->storeRequestResult(request, frontProtVal + conditionalFwdPV + accruedInterest);
        }
        /* CDS_CLEAN_PRICE is the forward "clean price".
           This is 1 - fee paid by protection buyer - accrued interest */
        request = control->requestsOutput(OutputRequest::CDS_CLEAN_PRICE);
        if (request) {
            results->storeRequestResult(request, 
                1 + frontProtVal + conditionalFwdPV);
        }

        /*=====================================================================
         * COMPUTE RELEVANT FORWARD FOR CLOSED-FORM OPTION PRICE
         *===================================================================*/
        double fwdThingForPayoffExpectation; // spread for spread option, else price.
        if (strikeIsSpread) {
            /*=================================================================
             * SPREAD-STRUCK OPTION
             *===============================================================*/
            
            if (isIndexOption) {

                /* Unpleasantly enough, getAccrueFee only exists in CDSParSpreads, not
                   ICDSParSpreads.*/
                //CDSParSpreads* ccrv = dynamic_cast<CDSParSpreads*>(crv.get());

                /* NEED TO ADJUST FORWARD TO ACCOUNT FOR THE CDSW FLAT_FORWARD
                UNWIND CONVENTION */
                /* Create a credit spread curve for a flat spread of the strike
                to the exercise date, with the same recovery as the original curve */
                ExpiryArraySP fcExpiries(new ExpiryArray());
                DoubleArraySP fcSpreads(new DoubleArray());
                DoubleArraySP fcUpfronts(new DoubleArray());
                ExpirySP nextExExp(new BenchmarkDate(cds->getMaturity()));
                fcExpiries->push_back(nextExExp);
                fcSpreads->push_back(nextStrk);
                fcUpfronts->push_back(0);
                CDSParSpreadsSP cdswCurve(new CDSParSpreads(
                    "CDSW_CURVE",
                    valueDate,
                    (cdsParSpreads->spotDate(valueDate)).daysDiff(valueDate),
                    cdsParSpreads->getRecovery(),
                    cdsParSpreads->getSwapFrequency(),
                    fcExpiries,
                    fcSpreads,
                    fcUpfronts,
                    cdsParSpreads->isFeeAccrued(),
                    cdsParSpreads->dayCountConv()->toString(),
                    cdsParSpreads->getBadDayConvention()->toString(),
                    cdsParSpreads->getHolidays(),
                    YieldCurveConstSP(discount.get()),
                    cdsParSpreads->getPrepayCurveObj()
               ));

                ICreditFeeLegSP cdsFeeLeg = cds->getFeeLeg();
                // fee leg must be of a fixed rate variety for this to work
                IFixedRateCreditFeeLeg* cdsfrcfl = dynamic_cast<IFixedRateCreditFeeLeg*>(cdsFeeLeg.get());
                if (!cdsfrcfl)
                {
                    throw ModelException(method, "Fee leg is not of a fixed fee variety "
                                                "(must implement IFixedRateCreditFeeLeg)");
                }

                DateTime lastRiskyDate = cdsFeeLeg->getLastPayDate();
                //negative fees result in negative annuity
                double strikeAnnuity = abs(cdsFeeLeg->getFeeLegPV(
                                            nextEx, ulSettlementDate,
                                            earliestRiskyDate, lastRiskyDate,
                                            *(discount.getSP()), *cdswCurve, prepay,
                                            true, cds->getAccrualDcc(), false, frModel)/cdsfrcfl->getRate());

		        strikeAnnuity -=accruedInterest / cdsfrcfl->getRate() ;
                forwardSpread += (nextStrk - fee)*(1.0 - strikeAnnuity/(fwdAnnuityVal*survivalProbToEx));

                /* add to output - bad name, but OK for testing */
                request = control->requestsOutput(OutputRequest::FORWARD_ADJUSTMENT);
                if (request) {
                    results->storeRequestResult(request, (nextStrk - fee)*
                                                (1.0 - strikeAnnuity/(fwdAnnuityVal*survivalProbToEx)));
                } 
            }
            fwdThingForPayoffExpectation = forwardSpread;

        } else {
            /*=================================================================
             * PRICE-STRUCK OPTION
             *===============================================================*/
            double adjustedFwdPV = conditionalFwdPV;
            if (isIndexOption) {
                // index swap: you need to adjust forward spread for defaults
                // the adjustedFwdPV should be the unconditional foward.
                adjustedFwdPV = adjustedFwdPV + frontProtVal;
            }
            
            // note p = 1+f not 1-f, since f is from perspective of
            // protection seller, not buyer!
            double forwardDirtyPrice = (1.0 + adjustedFwdPV);           
            fwdThingForPayoffExpectation = forwardDirtyPrice;
        }
        request = control->requestsOutput(OutputRequest::FORWARD_CDS_SPREAD);
        if (request) {
            results->storeRequestResult(request, forwardSpread);
        } 
        
        /*=====================================================================
         * COMPUTE FORWARD CLOSED-FORM PRICE OF OPTIONALITY
         *===================================================================*/
        /* YOU SHOULD PROBABLY CHECK THE underlying TO MAKE SURE THE VOLCUBE
           IS APPROPRIATE - OR SHOULD YOU JUST LET THIS GO AND INTERPOLATE ANYWAY? */
        /* How does it know what kind of vol cube to construct? You may have to
           include the model type in the vol request. Or has this been determined by
           the model type in the market data fetching? */
        CDSVolRequestSimpleEuropean volReq(
            nextStrk, fwdThingForPayoffExpectation, strikeIsSpread, 
            isCall, nextEx, cds->getMaturity(), 
            (dynamic_cast<const ClosedFormMultiQSmile*>(model) ? true : false));
        CDSVolProcessedSimpleEuropeanSP 
            procVolSP(dynamic_cast<CDSVolProcessedSimpleEuropean*>(cdsParSpreads->getProcessedVol(&volReq))); 
        double forwardPayoffExpectation = procVolSP->optionPrice();
        // IMPLIED VOL
        request = control->requestsOutput(OutputRequest::IND_VOL);
        if (request) {
            results->storeRequestResult(request, procVolSP->impliedVolatility());
        } 
        /* Use AVG_VOL field for ATM vol, for testing */
        request = control->requestsOutput(OutputRequest::AVG_VOL);
        if (request) {
            results->storeRequestResult(request, procVolSP->atmVolatility());
        } 
        /* Use AVG_FWD field for option moneyness, for testing */
        request = control->requestsOutput(OutputRequest::AVG_FWD);
        if (request) {
            results->storeRequestResult(request, procVolSP->moneyness());
        } 
        /* Delta */
		request = control->requestsOutput(OutputRequest::OPTION_DELTA);
        if (request) {
            results->storeRequestResult(request, 
                procVolSP->delta());
		}
		/* Gamma */
		request = control->requestsOutput(OutputRequest::OPTION_GAMMA);
        if (request) {
            results->storeRequestResult(request, 
                procVolSP->gamma());
		}
		/* Vega */
		request = control->requestsOutput(OutputRequest::OPTION_VEGA);
        if (request) {
            results->storeRequestResult(request, 
                procVolSP->vega());
		}
		
        
        /*=====================================================================
         * COMBINE WITH ANNUITY/2-STATE PRICING FOR FINAL ANSWER
         *===================================================================*/
        if (strikeIsSpread) {
            // pricing in annuity measure
            value = forwardPayoffExpectation * fwdAnnuityVal;
        } else {
            value = forwardPayoffExpectation;
        }
        request = control->requestsOutput(OutputRequest::OPTION_PRICE);
        if (request) {
            /* This is just the forward BS part of the option (n.b.
               includes annuity for spread options! */
            results->storeRequestResult(request, forwardPayoffExpectation);
        } 
        
        if (isIndexOption) {
            // just PV to today at risk-free rate - default probability
            // information was encoded in forward adjustment already for
            // price struck option. For spread struck option, the forward
            // annuity should also be multiplied by the survival probability
            value *= dfToExPmt * (strikeIsSpread ? survivalProbToEx : 1.0);
        } else {
            // 2-state pricing since payoff expectation was calculated
            // conditional on survival until exercise.
            value = dfToExPmt*(survivalProbToEx*value 
                + (1-survivalProbToEx)*optionVGD);
        }
        results->storePrice(value*notional, ccy);
        return;

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class CDSOptionHelper {
public:
    static IObject* defaultCDSOption();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CDSOption, clazz);
        SUPERCLASS(CInstrument);
		IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(ClosedFormBSImpliedSmile::IIntoProduct);
        IMPLEMENTS(ClosedFormMultiQSmile::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultCDSOption);

        FIELD(valueDate,                  "Valuation Date - default from market.");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(       premiumSettlement,          "Describes when the premium is payable. Default is rolling T+3");
        FIELD_MAKE_OPTIONAL(premiumSettlement);
        FIELD(notional,                   "Notional of the option.");
        FIELD(issueDate,                  "Start date of the option; if missing, default is valueDate.");
        FIELD_MAKE_OPTIONAL(issueDate);
        FIELD(koOnDefault,                "True if option knocks-out on underlying credit default.");
        FIELD(isCall,                     "True if call (right to buy risk/sell protection, else put (Right to sell risk/buy protection).");
        FIELD(isIndexOption,              "True if underlying name is an index rather than a single-name. Affects forward adjustment.");
        FIELD(exerciseType,               "Exercise type of option: european, american, or bermudan/multi-european. Must have one of the following values \"EXERCISE_TYPE_EUROPEAN\", \"EXERCISE_TYPE_BERMUDAN\", or \"EXERCISE_TYPE_AMERICAN\".");
        FIELD(       exerciseSchedule,           "Exercise schedule. Interpolation type matters only if the option is American, and then it must be INTERP_STAIRS or INTERP_LINEAR."); 
        FIELD(       exerciseSettlement,         "Describes underlying settlement after exercise date. Default is immediate settlement.");
        FIELD_MAKE_OPTIONAL(exerciseSettlement);
        FIELD(payAtDefault,               "If true, payments happen at moment of default. If not, then default payment waits until the next exercise.");
        FIELD(strikeType,                 "STRIKE_TYPE_SPREAD, STRIKE_TYPE_CLEAN_PRICE or STRIKE_TYPE_DIRTY_PRICE");
        FIELD(       underlying,                 "Underlying CDS that defines the CDS once the exercise date is known.");
        FIELD(ulMaturityDefinedByOption,  "If true, the underlying maturity is determined relative to the exercise date. If false, the underlying maturity is fixed in underlying.");
        FIELD(       ulExpiry,                   "Underlying maturity (used if ulMaturityDefinedByOption).");
        FIELD_MAKE_OPTIONAL(ulExpiry);
        FIELD(isCashSettled,              "True if cash settled, else physically settled.");
        FIELD(isExercised,                "True if the option has been exercised.");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(exerciseDate,               "Date that option was exercised.");
        FIELD_MAKE_OPTIONAL(exerciseDate);
        FIELD(exerciseValue,              "Value (forward) of whatever the strike is defined relative to (upfront-fee? par spread?) at exercise.");
        FIELD_MAKE_OPTIONAL(exerciseValue);
        FIELD(cdsParSpreads,              "Underlying of option");
        FIELD_MAKE_OPTIONAL(cdsParSpreads);
        FIELD(discount,                   "Discount curve");
        FIELD_MAKE_OPTIONAL(discount);
	FIELD(crv,"credit curve");
	FIELD_MAKE_TRANSIENT(crv);
	FIELD(rfCrv,"yield curve");
	FIELD_MAKE_TRANSIENT(rfCrv);
    }

};

IObject* CDSOptionHelper::defaultCDSOption() {
    return new CDSOption();
}

CClassConstSP const CDSOption::TYPE = 
    CClass::registerClassLoadMethod("CDSOption", typeid(CDSOption),CDSOptionHelper::load);



/*=============================================================================
 * Pricing Models
 *===========================================================================*/
class CDSOptionClosedFormBS : public ClosedFormBSImpliedSmile::IProduct {
private:
    const CDSOption* cdso; // a reference
    const ClosedFormBSImpliedSmile* model;

public:
    CDSOptionClosedFormBS(const CDSOption* cdso, const ClosedFormBSImpliedSmile* model): cdso(cdso), model(model){}
    void price(ClosedFormBSImpliedSmile* model,
               Control*         control, 
               CResults*        results) const {
        cdso->priceClosedForm(results, control, model);
    }
};
    
/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
ClosedFormBSImpliedSmile::IProduct* CDSOption::createProduct(ClosedFormBSImpliedSmile* model) const {
    return new CDSOptionClosedFormBS(this, model);
}

class CDSOptionClosedFormMultiQ : public ClosedFormMultiQSmile::IProduct{
private:
    const CDSOption* cdso; // a reference
    const ClosedFormMultiQSmile* model;

public:
    CDSOptionClosedFormMultiQ(const CDSOption* cdso, const ClosedFormMultiQSmile* model): cdso(cdso), model(model) {}
    void price(ClosedFormMultiQSmile* model,
               Control*         control, 
               CResults*        results) const {
        cdso->priceClosedForm(results, control, model);
    }
};
    
ClosedFormMultiQSmile::IProduct* CDSOption::createProduct(ClosedFormMultiQSmile* model) const{
    return new CDSOptionClosedFormMultiQ(this, model);
}

/*=============================================================================
 * Theta::Shift Interface
 *===========================================================================*/
bool CDSOption::sensShift(Theta* theta) {
    try {
        valueDate = theta->rollDate(valueDate);
    } catch (exception& e) {
        throw ModelException(e, "CDSOption::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}


DRLIB_END_NAMESPACE

