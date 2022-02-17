//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ConvBond.cpp
//
//   Description : Convertible Bond Instrument
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 27, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ConvBond.hpp"
#include "edginc/Class.hpp"
#include "edginc/Format.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Black.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/BondParams.hpp"
#include "edginc/CanBeRisky.hpp"
#include "edginc/Format.hpp"
#include "edginc/BondFloatNtl.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/RiskyCDSCurve.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/HaveEquity.hpp"
#include "edginc/AccrualCalendar.hpp"
#include "edginc/FD1FLV.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/VolSurface.hpp"

#include "edginc/IndexSpecEQ.hpp"

DRLIB_BEGIN_NAMESPACE

typedef struct _FloatingData
{
    YieldCurveConstSP       yieldCurve;
    DateTime                valueDate;
    DateTime                workoutDate;
    double                  bondDiff;
    int                     count;
    string                  interval;
    bool                    stubAtEnd;
    DayCountConventionSP    floatDCC;
    BondSP                  bond;
    double                  bondPrice;
    double                  faceValue;
} FloatingData;

static CFieldConstSP creditCurveField;

double spreadPV(double spread, void *liborStreamData)
{
    FloatingData* floaterStruct = (FloatingData*)liborStreamData;

    double          floatLegPV  = 0.;
    CashFlowArray   spreadFlows;

    if (floaterStruct->valueDate < floaterStruct->workoutDate) {
        spreadFlows = SwapTool::cashflows(
            floaterStruct->valueDate,
            floaterStruct->workoutDate,
            floaterStruct->stubAtEnd,
            spread,
            floaterStruct->count,           // interval = count periods
            floaterStruct->interval,          // e.g. Y, M, W, D
            floaterStruct->floatDCC.get());

            int i;
            for (i=0; i<spreadFlows.size()-1; i++) {
                floatLegPV += spreadFlows[i].amount * floaterStruct->yieldCurve->pv(floaterStruct->valueDate, spreadFlows[i].date);
            }

            if (spreadFlows.size() > 0) {
                floatLegPV += (spreadFlows[spreadFlows.size()-1].amount - 1.0) * floaterStruct->yieldCurve->pv(floaterStruct->valueDate, spreadFlows[spreadFlows.size()-1].date);
            }

    } else {
        floatLegPV = 0.;
    }

    return floaterStruct->faceValue*floatLegPV - floaterStruct->bondDiff;
}

double zSpreadBond(double spread, void *liborStreamData)
{
    FloatingData* floaterStruct = (FloatingData*)liborStreamData;

    BootstrappedYieldCurveSP yieldCurve(BootstrappedYieldCurveSP::dynamicCast((IObjectSP)floaterStruct->yieldCurve->clone()));

    yieldCurve->sensShift(PropertyTweak<RateParallel>(spread));

    YieldCurveSP yc = yieldCurve;

    double bondPV = floaterStruct->bond->presentValue(floaterStruct->valueDate, yc);

    yieldCurve->sensRestore(PropertyTweak<RateParallel>(spread));

    return floaterStruct->bondPrice - bondPV;
}

// copy market data relevant to the instrument
void ConvBond::GetMarket(const IModel* model, const CMarketDataSP market) {

    market->GetReferenceDate(valueDate);
    discount.getData(model, market);

    if (!creditSpreads.getName().empty()) {
        creditSpreads.getData(model, market);
        
#ifdef CDS_BACKWARD_COMPATIBILITY
        if (ICDSParSpreads::TYPE->isInstance(creditSpreads.get())) {
            ICDSParSpreads* parSpreads = 
                dynamic_cast<ICDSParSpreads*>(creditSpreads.get());
            ICDSParSpreadsSP cdsSpreads = 
                ICDSParSpreadsSP::attachToRef(parSpreads);
                
            BadDayConventionSP bdc(BadDayConventionFactory::make("None"));
            cdsSpreads->setBadDayConvention(bdc);
            cdsSpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));

        }
#endif
        
    }
    bond->getMarket(model, market.get());

    CAsset::getAssetMarketData(model, market.get(), ccyTreatment,
                               discount.getName(), asset);

    if (!!resetSchedule) {
        resetSchedule->setValueDate(valueDate);
    }

    if( !!contConversion )
        contConversion->getMarket(model, market.get(), this);
}

void ConvBond::validatePop2Object()
{
    static const string method = "ConvBond::validatePop2Object";

    try{
        if (payoffType == "D" || payoffType == "d") {
            DECS = true;
        }
        else if (payoffType == "P" || payoffType == "p") {
            PERCS = true;
        }

        if (dynamic_cast<BondFloatNtl*>(bond.get()))
            schedsArePcts = true;

    }
    catch (exception& e) {
       throw ModelException(e, method);
    }

    return;
}

void ConvBond::Validate() {
    static const string method = "ConvBond::Validate";

    try {
        if (!!softPutSchedule || !!softPutTriggerSchedule) {
            throw ModelException(method, "The soft put feature has been disabled - please contact Derivatives Research.");
        }

        if (!(ccyTreatment == "N" || ccyTreatment == "S" || ccyTreatment == "V" )) {
            throw ModelException(method, "only ccyTreatment types N, V and S currently supported");
        }
        if (!(callTreatment == "K" || callTreatment == "B")) {
            throw ModelException(method, "only callTreatment type K currently supported");
        }
        if (!(makeWholeType == "N" || makeWholeType == "C" || makeWholeType == "A")) {
            throw ModelException(method, "only makeWholeType N (None), C (PV of Coupons) or A (Amazon) are currently supported");
        }

        if (dividendPassThrough && dividendAdjusted ) {
             throw ModelException(method, "A bond can either be dividend pass through or dividend adjusted, but not both");
        }

        if (dividendPassThrough && !conversionRatios->isFlat()) {
             throw ModelException(method, "The conversion ratio schedule must be flat for dividend pass through bonds");
        }

        if (dividendAdjusted && !conversionRatios->isFlat()) {
             throw ModelException(method, "The conversion ratio schedule must be flat for dividend adjusted bonds");
        }

        if (alreadyCalled == true){

            if ( callTreatment != "B" ) {
                throw ModelException(method, "Already called bonds only supported in conjunction with"
                    "Black-on-call call treatment");
            }

            if ( callNotification <= 0) {
                throw ModelException(method, "Bond can't already have been called(" +
                    dateCallNotifSent.toString()  + ") as call notification period(" +
                    Format::toString(callNotification) + ") < 1.");
            }
            if (dateCallNotifSent > valueDate) {
                    throw ModelException(method, "Bond can't already have been called(" +
                        dateCallNotifSent.toString()  + ") after valuation date (" +
                        valueDate.toString() +")");
            }

            DateTime callMatDate = dateCallNotifSent.rollDate(callNotification);

            if ( callMatDate > bond->getMaturityDate() ) {
                throw ModelException(method, "Bond can't be called(" +
                    dateCallNotifSent.toString() + ") < " + Format::toString(callNotification) +
                    "(notification period) days before bond maturity(" + bond->getMaturityDate().toString() + ")");
            }

            if ( frenchExtendedConv == true ) {
                DateTime endDate = frenchExtConvInt->toDate(callMatDate);
                if (endDate < valueDate) {
                    throw ModelException(method, "Bond was called (" + dateCallNotifSent.toString() +
                                     ") and expired ( " + endDate.toString() + ") before valuation date (" +
                                     valueDate.toString() + ")");
                }
            }
        }

        /*
        if (frenchDivTreatment == true){
            throw ModelException(method, "French dividend treatment is not supported");
        }
        if (frenchExtendedConv == true){
            throw ModelException(method, "French extended conversion is not supported");
        }
        */

        if (!(!callSchedule) && callSchedule->length() > 0){
            if (callSchedule->isPositive() == false)
                throw ModelException(method, "Call schedule must be positive");
            if (callSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                callSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "Call schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && callSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "Call schedule times must be start of day");
            if (!Maths::isPositive(callOptimality)) {
                throw ModelException(method, "Call optimality (" + Format::toString(callOptimality) +
                                             ") must be positive");
            }
        }
        if (!(!softCallSchedule) && callTreatment != "K" && softCallSchedule->length() > 0){
            if (softCallSchedule->isPositive() == false)
                throw ModelException(method, "Soft-call schedule must be positive");
            if (softCallSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                softCallSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "Soft-call schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && softCallSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "Soft-call schedule times must be start of day");
        }
        if (!(!softCallTriggerSchedule) && softCallTriggerSchedule->length() > 0) {
            if (softCallTriggerSchedule->isPositive() == false)
                throw ModelException(method, "Soft-call trigger schedule must be positive");
            if (softCallTriggerSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                softCallTriggerSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "Soft-call trigger schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && softCallTriggerSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "Soft-call trigger schedule times must be start of day");
            if ( (!(!callSchedule) && callSchedule->length() > 0)
                 && softCallTriggerSchedule->firstDate() > callSchedule->firstDate()) {
                throw ModelException(method,
                    "first soft call date ("                             +
                    softCallTriggerSchedule->firstDate().toString()      +
                    ") must be before first regular call notification (" +
                    callSchedule->firstDate().toString()                 + ")");
            }
        }

        if (!(!softPutSchedule) && callTreatment != "K" && softPutSchedule->length() > 0){
            if (softPutSchedule->isPositive() == false)
                throw ModelException(method, "Soft-put schedule must be positive");
            if (softPutSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                softPutSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "Soft-put schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && softPutSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "Soft-put schedule times must be start of day");
        }
        if (!(!softPutTriggerSchedule) && softPutTriggerSchedule->length() > 0) {
            if (softPutTriggerSchedule->isPositive() == false)
                throw ModelException(method, "Soft-put trigger schedule must be positive");
            if (softPutTriggerSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                softPutTriggerSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "Soft-put trigger schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && softPutTriggerSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "Soft-put trigger schedule times must be start of day");
        }

        if (!(!conversionRatios) && conversionRatios->length() > 0) {
            if (conversionRatios->isNonNegative() == false)
                throw ModelException(method, "conversion schedule must be positive");
            if (conversionRatios->lastDate() > bond->getUnadjMaturityDate() &&
                conversionRatios->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "conversion schedule cannot end after bond maturity.\n"
                "Last conversion date:" + conversionRatios->lastDate().toString() +
                "\nBond Maturity:" + bond->getUnadjMaturityDate().toString());
            if (DEBUG_SOD_Only == true && conversionRatios->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "conversion schedule times must be start of day");
        }
        if (!(!convCashSchedule) && convCashSchedule->length() > 0) {
            if (convCashSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                convCashSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "conversion cash schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && convCashSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "conversion cash schedule times must be start of day");
        }
        if (!(!parityFloorSchedule) && parityFloorSchedule->length() > 0) {
            if (DECS == true || PERCS == true)
                throw ModelException(method, "parity floor schedule should not be specified when convert is DEC or PERC");
            if (parityFloorSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                parityFloorSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "parity floor schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && parityFloorSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "parity floor schedule times must be start of day");
        }

        if (!(!putSchedule) && putSchedule->length() > 0) {
            if (putSchedule->lastDate() > bond->getUnadjMaturityDate() &&
                putSchedule->lastDate() > bond->getMaturityDate())
                throw ModelException(method, "put schedule cannot end after bond maturity");
            if (DEBUG_SOD_Only == true && putSchedule->timesAreAll(DateTime::START_OF_DAY_TIME) == false)
                throw ModelException(method, "put schedule times must be start of day");
        }

        if( !!putToStock )
        {
            if ( !putSchedule || putSchedule->length() == 0)
                throw ModelException(method, "put schedule empty but is putToStock");
            if ( ccyTreatment == "S" || ccyTreatment == "P" )
                throw ModelException(method, "Does not allow currency feature if putToStock");
            if( putToStock->putToStockNotified && putToStock->putToStockPutDate < valueDate )
                throw ModelException(method, "Put date already passed for putToStock");
        }


        if (DEBUG_SOD_Only == true && bond->getMaturityDate().getTime() != DateTime::START_OF_DAY_TIME)
            throw ModelException(method, "bondMaturity must be start of day");


        if (getCouponAtMat == false) {
            // Test to see that the coupon does not go ex-div before maturity since you don't get it.
            // If you're supposed to get it early and then must pay it back upon conversion, set getCouponAtMat
            // to true and use negative conversion cash in the amount of the coupon after record date.
            CashFlowArraySP allCFs = bond->getCashFlows();
            DateTimeArraySP exDates = bond->getExCouponDates();
            if( (*allCFs)[allCFs->size()-1].date != (*exDates)[exDates->size()-1]  &&
                 !Maths::equals((*allCFs)[allCFs->size()-1].amount, bond->getRedemption())) {
                throw ModelException(method, "getCouponAtMat must be TRUE for bonds whose final coupon goes "
                    "ex-div before redemption is paid");
            }
        }


        if (!(softCallAdjustType == "N" || softCallAdjustType == "n" ||
              softCallAdjustType == "D" || softCallAdjustType == "d" ||
              softCallAdjustType == "R" || softCallAdjustType == "r" ||
              softCallAdjustType == "T" || softCallAdjustType == "t" ||
              softCallAdjustType == "E" || softCallAdjustType == "e" )) {
            throw ModelException(method, softCallAdjustType + " is not a supported softCallAdjustType");
        }

        if ( softCallAdjustType == "R" || softCallAdjustType == "r" ||
             softCallAdjustType == "T" || softCallAdjustType == "t" ||
             softCallAdjustType == "E" || softCallAdjustType == "e" ) {
            if ( softCallDaysN < 0 ) {
                throw ModelException(method, " softCallDaysN must be positive");
            }
            if ( softCallDaysM < 0 ) {
                throw ModelException(method, " softCallDaysM must be positive");
            }
            if ( softCallDaysM > softCallDaysN ) {
                throw ModelException(method, "softCallDaysM must not be greater than softCallDaysN");
            }
        }

        if (!!resetSchedule && resetSchedule->length() > 0 ) {
            if (triggerReset) {
                throw ModelException(method, "Trigger resets not allowed for resettable convertibles");
            }

            if ( callTreatment == "B" ) {
                throw ModelException(method, "Black on Call is not allowed for resettables");
            }

            if ( !conversionRatios->isFlat() ) {
                throw ModelException(method, "The conversion ratio schedule must be flat for resettable convertibles");
            }
        }

        // state dependent trigger validation
        if ( triggerReset == true ) {
            if ( resetType != "UP_RESET" &&
                 resetType != "DOWN_RESET") {
                throw ModelException(method, "resetType must be either UP_RESET or DOWN_RESET");
            }

            if ( resetType == "UP_RESET"                &&
                 resetTrigger <= asset->getSpot()       &&
                 valueDate >= resetObservationStartDate &&
                 valueDate >= resetObservationEndDate ) {
                throw ModelException(method, "Instrument has knocked into a non-resettable convertible.");
            }

            if ( resetType == "DOWN_RESET"              &&
                 resetTrigger >= asset->getSpot()       &&
                 valueDate >= resetObservationStartDate &&
                 valueDate >= resetObservationEndDate ) {
                throw ModelException(method, "Instrument has knocked into a non-resettable convertible.");
            }

            if (DECS == true || PERCS == true) {
                throw ModelException(method, "Reset trigger must be false for DECS and PERCS.");
            }
        }

        // -- Validate DECS and PERCS --
        if (DECS == true || PERCS == true) {
            // cannot be resettable
            if (!!resetSchedule && resetSchedule->length() > 0 ) {
                throw ModelException(method, "DECS and PERCS can not be resettable");
            }

            // get coupon at mat must be true
            if (getCouponAtMat == false) {
                throw ModelException(method, "getCouponAtMat must be true for mandatory convertibles");
            }
            // no conversion cash ???
            if (!(!convCashSchedule)) {
// Put this back in ???
//                if(convCashSchedule->isAlmostZero() == false) {
//                   throw ModelException(method, "conversion cash is not allowed with mandatory convertibles");
//                }
            }

            if (initialConvRatio <= 0) {
                throw ModelException(method, "initialConvRatio (" + Format::toString(initialConvRatio) + ") must be positive");
            }

            if (DECS == true) {
                if (initialPrice > 0 && fabs(initialPrice*initialConvRatio - bond->getFaceValue()) > .0005*bond->getFaceValue()) {
                    throw ModelException(method, "initialPrice*initialConvRatio must be equal to face value. Use initialPrice of " + Format::toString(bond->getFaceValue()/initialConvRatio));
                }

                if (minConvRatio <= 0) {
                   throw ModelException(method, "minConvRatio (" + Format::toString(minConvRatio) + ") must be positive");
                }

                if (minConvRatio > initialConvRatio) {
                    throw ModelException(method, "minConvRatio (" + Format::toString(minConvRatio) + ") must be less than initialConvRatio (" + Format::toString(initialConvRatio) + ")");
                }
            }
        }

		if ( CDSParSpreads::TYPE->isInstance(creditSpreads.get()) && riskFreeCoupons == true ) {
            throw ModelException(method, "Risk free coupons are not supported when pricing off a CDS spread curve");
        }

        if ( stockInEscrow && convertIntoIssuerShares ) {
           throw ModelException(method, "Stock can only be held in Escrow accounts for exchangeables");
        }

        if ( CDSParSpreads::TYPE->isInstance(creditSpreads.get()) && useJointDefaultCorrelation ) {
            throw ModelException(method, "Joint default prob can only be used with Z-spread curves");
        }

        if (useJointDefaultCorrelation && (jointDefaultCorrelation < 0.0 || jointDefaultCorrelation > 1.0 )) {
            throw ModelException(method, "Joint default prob must be between zero and one");
        }

        // addOnConvRatio only available for limitted type of conv bond
        if( hasAddOnConvRatios() )
        {
            if (triggerReset) {
                throw ModelException(method, "Does not allow AddOnConvRatio when convert trigger");
            }

            if (DECS == true || PERCS == true)
                throw ModelException(method, "Does not allow AddOnConvRatio when convert is DEC or PERC");

            if ( dividendAdjusted )
                throw ModelException(method, "Does not allow AddOnConvRatio when conversion is dividendAdjusted");

            if (dividendPassThrough )
                throw ModelException(method, "Does not allow AddOnConvRatio when conversion is dividendPassThrough");

            if ( ccyTreatment == "S" )
                throw ModelException(method, "Does not allow AddOnConvRatio when convert is currency struck");

            if( !!resetSchedule && resetSchedule->length() > 0 )
                throw ModelException(method, "Does not allow AddOnConvRatio when convert is resettable");
        }

        // contConversion only available for limitted type of conv bond
        if( hasContConversion() )
        {
            if (triggerReset) {
                throw ModelException(method, "Does not allow ContConversion when convert is trigger reset");
            }

            if (DECS == true || PERCS == true)
                throw ModelException(method, "Does not allow ContConversion when convert is DEC or PERC");

            if ( ccyTreatment == "S" )
                throw ModelException(method, "Does not allow ContConversion when convert is currency struck");

            if( !!resetSchedule && resetSchedule->length() > 0 )
                throw ModelException(method, "Does not allow ContConversion when convert is resettable");

            if ( hasAddOnConvRatios() )
                throw ModelException(method, "Does not allow ContConversion when has addOnConvRatios");

		}

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Pull out the component assets & correlations from the market data */
void ConvBond::getMarket(const IModel* model, const MarketData* market){
    try{

//        bond->getMarket(model, market);
    }
    catch (exception& e){
        throw ModelException(e, "ConvBond::getMarket", "Failed for Bond");
    }
}


// -- Stuff every product needs to do -- //

/** what's today ? */
DateTime ConvBond::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime ConvBond::endDate(const Sensitivity* sensControl) const {
    return bond->getMaturityDate(); //max((*cfl)[last].date, maturity);
}

/** Returns name identifying credit curve for adjusted credit rho parallel */
string ConvBond::sensName(AdjCreditSpreadRhoParallel* shift) const{
   static const string method = "ConvBond::sensName";
   try {    
        if (CreditSpreadCurve::TYPE->isInstance(creditSpreads.get())) {
            CreditSpreadCurveConstSP CScurve = 
                CreditSpreadCurveConstSP(dynamic_cast<const CreditSpreadCurve *>(creditSpreads.get()));
            return CScurve->getName();  
        } else {
            // don't report sensitivity if pricing off CDS par spread curve
            return "";
        }
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

// Set to risky growth
bool ConvBond::sensShift(AdjCreditSpreadRhoParallel* shift){
   static const string method = "ConvBond::sensShift";
   try {
       // cache riskyGrowth flag
       shift->setRiskyGrowth(riskyGrowth);
       riskyGrowth = true;
       return true; // need to shift the credit spread curve
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

/** Restores the object to its original form */
void ConvBond::sensRestore(AdjCreditSpreadRhoParallel* shift){
   static const string method = "ConvBond::sensShift";
   try {
       riskyGrowth = shift->getRiskyGrowth();
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

/** Returns name identifying credit curve for adjusted credit rho parallel */
string ConvBond::sensName(AdjCreditSpreadRhoPointwise* shift) const{
   static const string method = "ConvBond::sensName";
   try {    
        if (CreditSpreadCurve::TYPE->isInstance(creditSpreads.get())) {
            CreditSpreadCurveConstSP CScurve = 
                CreditSpreadCurveConstSP(dynamic_cast<const CreditSpreadCurve *>(creditSpreads.get()));
            return CScurve->getName();  
        } else {
            // don't report sensitivity if pricing off CDS par spread curve
            return "";
        }
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

ExpiryArrayConstSP ConvBond::sensExpiries(AdjCreditSpreadRhoPointwise* shift) const{
    return ExpiryArrayConstSP(); // return null SP - the credit spread curve will take care of this
}


// Set to risky growth
bool ConvBond::sensShift(AdjCreditSpreadRhoPointwise* shift){
   static const string method = "ConvBond::sensShift";
   try {
       // cache riskyGrowth flag
       shift->setRiskyGrowth(riskyGrowth);
       riskyGrowth = true;
       return true; // need to shift the credit spread curve
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

/** Restores the object to its original form */
void ConvBond::sensRestore(AdjCreditSpreadRhoPointwise* shift){
   static const string method = "ConvBond::sensShift";
   try {
       riskyGrowth = shift->getRiskyGrowth();
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

bool ConvBond::sensShift(Theta* shift) {
    try {
        DateTime oldValueDate = valueDate;
        valueDate = shift->rollDate(valueDate);

        if (!!resetSchedule) {
            double  newSpot      = (shift->useAssetFwds())?asset->fwdValue(valueDate):asset->getSpot();
            double  newFX        = 1.0;
            bool    isStruck     = false;

            if ( StruckEquity::TYPE->isInstance(getEquityAsset().get())) {
                StruckEquity* struckAsset = dynamic_cast<StruckEquity*>(getEquityAsset().get());
                newFX                     = (shift->useAssetFwds())?struckAsset->fxFwdValue(valueDate):struckAsset->getFXSpot();
                newSpot                  /= newFX;
                isStruck                  = true;
            }
            resetSchedule->rollDate(oldValueDate, valueDate, newSpot, isStruck, newFX);
        }

        if( hasContConversion() )
            contConversion->rollDate(shift, this);

        if( !!putToStock && putToStock->putToStockNotified )
            putToStock->putToStockSamples->roll(shift->getUtil(oldValueDate), 0, asset.get());

    }
    catch (exception& e) {
        throw ModelException(e, "ConvBond::sensShift(theta)");
    }
    return true; // our components have theta type sensitivity
}

// for reflection
ConvBond::ConvBond(): CInstrument(TYPE) {
    // initialize optional parameters
    ccyTreatment                = "N";
    redemptionAdjustForAccrued  = true;
    triggerAdjustForAccrued     = false;
    notResetCallTrigger         = false;
    callOptimality              = 1.;
    softCallAdjustType          = "D";
    softCallDaysM               = 20;
    softCallDaysN               = 30;
    callNotification            = 30;
    callTreatment               = "K";
    putAdjustForAccrued         = true;
    putNotification             = 0;
    convertIntoIssuerShares     = true;
    instSettle = InstrumentSettlementSP(new PhysicalSettlement());
    premiumSettle = InstrumentSettlementSP(   );
    getCouponAtMat              = true;
    convertOnCall               = true;
    currentConversionAvg        = 0;
    currentConversionWeight     = 0;
    isPreferred                 = false;
    riskyGrowth                 = false;
    makeWholeType               = "N";
    makeWholeAmount             = 0.;
    alreadyCalled               = false;

    payoffType                  = "S";
    DECS                        = false;
    PERCS                       = false;
    initialPrice                = -1.;
    convPrice                   = -1.;
    initialConvRatio            = -1.;
    minConvRatio                = -1.;

    decsHasCutoff               = false;
    decsCutoffLevel             = 0.0;
    cappedDECSTrigger           = 0.0;

    riskFreeCoupons             = false;
    continuousReset             = false;
    delayReset                  = 0;
    includePutSpread            = false;

    frenchDivTreatment          = false;
    frenchExtendedConv          = false;
    frenchExtConvInt = MaturityPeriodSP(new MaturityPeriod("0D"));

    schedsArePcts               = false;

    inDefault                   = false;
    recoveryPct                 = 0.;
    useAssetRecovery            = false;

    dividendPassThrough         = false;
    dividendAdjusted            = false;
    dividendPassThroughPct      = 1.0;

    triggerReset                = false;
    resetTrigger                = 0.0;
    resetConversionRatio        = 0.0;

    accelerateDECS              = false;
    accelerateDECSTriggerLevel  = 0.0;

    DEBUG_SOD_Only              = true;

    stockInEscrow               = false;

    useJointDefaultCorrelation  = false;
    jointDefaultCorrelation     = 0.0;
    stockCreditSpread           = 0.0;

    contingentConversion        = false;
    triggerActiveAtMat          = false;
    contingentConversionTrigger = 0.0;

    payAccruedUponDefault       = true;

    maxWithParity               = false;
    manualResetConvPricePct     = 1.0;

    dontScaleNonPrefPriceBy100              = false;
}

ConvBond::ConvBond(
             YieldCurveWrapper         discount,
             CreditCurveWrapper        creditSpreads,
             CAssetWrapper             asset,
             const DateTime&           valueDate,
             const string&             ccyTreatment,
             const double&             faceValue,
             MaturityPeriodSP          bondMaturity,
             HolidayWrapper            holidays,
             int                       bondFrequency,
             const string&             dayCountConvString,
             MaturityPeriodSP          callStartPeriod,
             const double&             callTriggerLevel,
             const double&             coupon,
             const double&             redemption,
             const double&             putLevel,
             MaturityPeriodSP          putStartPeriod,
             const double&             spotAtStart,
             const double&             conversionPremium,
             const DateTime&           knockInDate):
        CInstrument(TYPE), ccyTreatment(ccyTreatment),
        redemptionAdjustForAccrued(true), triggerAdjustForAccrued(true),
        notResetCallTrigger(false),
        callOptimality(1.), softCallAdjustType("D"), softCallDaysM(20),
        softCallDaysN(30), callNotification(30), callTreatment("K"),
        putAdjustForAccrued(true), putNotification(0),
        convertIntoIssuerShares(true),
        getCouponAtMat(true), convertOnCall(true), currentConversionAvg(0.0),
        currentConversionWeight(0.0), isPreferred(false), riskyGrowth(false),
        makeWholeType("N"), makeWholeAmount(0.), alreadyCalled(false),
        schedsArePcts(false), inDefault(false), recoveryPct(0.0), 
        useJointDefaultCorrelation(false),
        jointDefaultCorrelation(0.0), stockCreditSpread(0.0), payoffType("S"),
        DECS(false), PERCS(false), initialConvRatio(-1.), minConvRatio(-1.),
        initialPrice(-1.), convPrice(-1.), riskFreeCoupons(false), frenchDivTreatment(false),
        frenchExtendedConv(false), dividendPassThrough(false), dividendAdjusted(false), triggerReset(false),
        resetTrigger(0.0), resetConversionRatio(0.0), accelerateDECS(false), accelerateDECSTriggerLevel(0.0),
        dontScaleNonPrefPriceBy100(false), DEBUG_SOD_Only(true), stockInEscrow(false),  contingentConversion(false),
        triggerActiveAtMat(false), contingentConversionTrigger(0.0)
{
    instSettle               = InstrumentSettlementSP(new PhysicalSettlement());
    premiumSettle            = InstrumentSettlementSP(   );

    this->discount           = YieldCurveWrapper(copy(discount.get()));
    this->creditSpreads      = CreditCurveWrapper(copy(creditSpreads.get()));
    this->asset              = CAssetWrapper(copy(asset.get()));
    this->valueDate          = valueDate;
    this->ccyTreatment       = ccyTreatment;
    this->frenchExtConvInt   = MaturityPeriodSP(   );

    this->accelerateDECS              = false;
    this->accelerateDECSTriggerLevel  = 0.0;

    this->continuousReset             = false;
    this->delayReset                  = 0;
    this->includePutSpread            = false;

    this->maxWithParity               = false;
    this->manualResetConvPricePct     = 1.0;

    this->useAssetRecovery            = false;
    this->payAccruedUponDefault       = true;

    this->decsHasCutoff               = false;
    this->decsCutoffLevel             = 0.0;
    this->cappedDECSTrigger           = 0.0;

    this->dividendPassThroughPct      = 1.0;

    // create conversion ratio schedule
    DoubleArray   convSchedRates(0);
    DateTimeArray convSchedDates(0);
    double convRatio = (redemption * faceValue) / ( spotAtStart * (1 + conversionPremium));
    convSchedDates.push_back(knockInDate);
    convSchedRates.push_back(convRatio);
    convSchedDates.push_back(bondMaturity->toDate(knockInDate));
    convSchedRates.push_back(convRatio);
    conversionRatios = ScheduleSP( new Schedule(convSchedDates, convSchedRates, "L"));

    // create call schedule
    DoubleArray   callSchedRates(0);
    DateTimeArray callSchedDates(0);
    callSchedDates.push_back(callStartPeriod->toDate(knockInDate));
    callSchedRates.push_back(callTriggerLevel);
    callSchedDates.push_back(bondMaturity->toDate(knockInDate));
    callSchedRates.push_back(callTriggerLevel);

    callSchedule = ScheduleSP( new Schedule(callSchedDates, callSchedRates, "L"));

    // create put schedule
    if ( putLevel > 0.0 ) {
        DoubleArray   putSchedRates(0);
        DateTimeArray putSchedDates(0);
        putSchedDates.push_back(putStartPeriod->toDate(knockInDate));
        putSchedRates.push_back(putLevel);
        putSchedule = ScheduleSP( new Schedule(putSchedDates, putSchedRates, "N"));
    }

     // create a new underlying bond
    BondParams* tmpBond = new BondParams();

    tmpBond->faceValue          = faceValue;
    tmpBond->redemptionPct      = redemption;
    tmpBond->couponPct          = coupon;
    tmpBond->frequency          = bondFrequency;
    tmpBond->maturityDate       = bondMaturity->toDate(knockInDate);
    tmpBond->datedDate          = knockInDate;
    tmpBond->dayCountConvString = dayCountConvString;
    tmpBond->oddLastShort       = false;
    tmpBond->endOfMonthAdj      = true;
    tmpBond->eomIgnoreLeapYear  = false;
    tmpBond->badDayConvString   = "N";
    tmpBond->hols               = holidays;
    tmpBond->exDivDays          = 0;
    tmpBond->exDivRule          = "N";
    tmpBond->taxRate            = 0.0;

    bond = BondSP(tmpBond);
    bond->validatePop2Object();
}



ConvBond* ConvBond::createDECS(YieldCurveWrapper         discount,
                               CreditCurveWrapper        creditSpreads,
                               CAssetWrapper             asset,
                               const DateTime&           valueDate,
                               const string&             ccyTreatment,
                               const double&             faceValue,
                               MaturityPeriodSP          bondMaturity,
                               HolidayWrapper            holidays,
                               int                       bondFrequency,
                               const string&             dayCountConvString,
                               const double&             coupon,
                               const double&             redemption,
                               const double&             spotAtStart,
                               const double&             downsideProtection,
                               const double&             decsPremium,
                               const bool                divPassThrough,
                               const DateTime&           knockInDate)
{

    ConvBond* newDECS = new ConvBond();

    newDECS->instSettle               = InstrumentSettlementSP(new PhysicalSettlement());
    newDECS->premiumSettle            = InstrumentSettlementSP(   );

    newDECS->discount           = YieldCurveWrapper(copy(discount.get()));
    newDECS->creditSpreads      = CreditCurveWrapper(copy(creditSpreads.get()));
    newDECS->asset              = CAssetWrapper(copy(asset.get()));
    newDECS->valueDate          = valueDate;
    newDECS->ccyTreatment       = ccyTreatment;
    newDECS->frenchExtConvInt   = MaturityPeriodSP(   );

    // create conversion ratio schedule
    DoubleArray   convSchedRates(0);
    DateTimeArray convSchedDates(0);
    double convRatio = 1.0;
    convSchedDates.push_back(bondMaturity->toDate(knockInDate));
    convSchedRates.push_back(convRatio);
    newDECS->conversionRatios = ScheduleSP( new Schedule(convSchedDates, convSchedRates, "N"));

    // create DECS parameters
    newDECS->payoffType         = "D";
    newDECS->DECS               = true;
    newDECS->PERCS              = false;
    newDECS->initialPrice       = (1-downsideProtection)*spotAtStart;
    newDECS->initialConvRatio   = faceValue / newDECS->initialPrice;
    newDECS->convPrice          = spotAtStart * (1.0 + decsPremium);
    newDECS->minConvRatio       = faceValue / newDECS->convPrice;
    // newDECS->minConvRatio       = 1.0 / ( 1 + decsPremium);
    newDECS->riskFreeCoupons    = false;

    newDECS->accelerateDECS              = false;
    newDECS->accelerateDECSTriggerLevel  = 0.0;

    newDECS->dividendPassThrough         = divPassThrough;
    newDECS->dividendPassThroughPct      = 1.0;         
    newDECS->dividendAdjusted            = false;

    newDECS->decsHasCutoff               = false;
    newDECS->decsCutoffLevel             = 0.0;
    newDECS->cappedDECSTrigger           = 0.0;

     // create a new underlying bond
    BondParams* tmpBond = new BondParams();

    tmpBond->faceValue          = faceValue;
    tmpBond->redemptionPct      = redemption;
    tmpBond->couponPct          = coupon;
    tmpBond->frequency          = bondFrequency;
    tmpBond->maturityDate       = bondMaturity->toDate(knockInDate);
    tmpBond->datedDate          = knockInDate;
    tmpBond->dayCountConvString = dayCountConvString;
    tmpBond->oddLastShort       = false;
    tmpBond->endOfMonthAdj      = true;
    tmpBond->eomIgnoreLeapYear  = false;
    tmpBond->badDayConvString   = "N";
    tmpBond->hols               = holidays;
    tmpBond->exDivDays          = 0;
    tmpBond->exDivRule          = "N";
    tmpBond->taxRate            = 0.0;

    newDECS->bond = BondSP(tmpBond);
    newDECS->bond->validatePop2Object();

    newDECS->validatePop2Object();

    return newDECS;
}


void ConvBond::scaleOutputs(CControlSP control, ResultsSP unscaledResults)
{
    // only scale if it's not preferred
   if (!isPreferred) {

        // retrieve clean price and store it as fair value
       const string ccy = unscaledResults->getCcyName();
       OutputNameConstSP cleanPriceName(new OutputName(OutputRequest::CLEAN_PRICE_FOR_SCALING));
        double cleanPrice = unscaledResults->retrieveScalarGreek(
                Results::DEBUG_PACKET,
                    cleanPriceName);
       double scaleAmount = dontScaleNonPrefPriceBy100 ? 1.0 : 100.0;
       unscaledResults->storePrice(cleanPrice * scaleAmount, ccy);

       // scale the clean and dirty prices by scaleAmount
       OutputRequest* request = NULL;
       if ( control->requestsOutput(OutputRequest::DIRTY_PRICE, request) ) {
           unscaledResults->scale(Results::INSTRUMENT_PACKET,
                                  OutputRequest::DIRTY_PRICE,
                                  scaleAmount);
       }
       if ( control->requestsOutput(OutputRequest::CLEAN_PRICE, request) ) {
           unscaledResults->scale(Results::INSTRUMENT_PACKET,
                                  OutputRequest::CLEAN_PRICE,
                                  scaleAmount);
       }

       // scale everything by 1/faceValue
       double scalingFactor = 1.0 / bond->getFaceValue();

       unscaledResults->scalePostProcess(control, scalingFactor);
   }
}

class ConvBondHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ConvBond, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(FD1F::IIntoProduct);
        IMPLEMENTS(FD1FGeneric::IIntoProduct);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(IScaleOutputs);
        IMPLEMENTS(IRiskyPricer);
        IMPLEMENTS(IInstrumentAsAsset);
        IMPLEMENTS(AdjCreditSpreadRhoParallel::IRestorableShift);
        IMPLEMENTS(AdjCreditSpreadRhoPointwise::IRestorableShift);
        EMPTY_SHELL_METHOD(defaultConvBond);
        FIELD(discount, "identifies discount curve");
        FIELD(creditSpreads, "credit spread to discount curve");
        FIELD(asset, "the equity");
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(ccyTreatment, "N or V for normal, S for struck");
        FIELD_MAKE_OPTIONAL(ccyTreatment);

        FIELD(bond, "bond object");

        FIELD(putAdjustForAccrued, "receive accrued upon put");
        FIELD_MAKE_OPTIONAL(putAdjustForAccrued);
        FIELD(putNotification, "put notification period length. in business days");
        FIELD_MAKE_OPTIONAL(putNotification);
        FIELD(putToStock,                 "parameters for putToStock, ie. receive stock when put");
        FIELD_MAKE_OPTIONAL(putToStock);

        FIELD(redemptionAdjustForAccrued, "receive accrued upon call");
        FIELD_MAKE_OPTIONAL(redemptionAdjustForAccrued);
        FIELD(triggerAdjustForAccrued, "soft-call trigger level adjusted for accrued interest");
        FIELD_MAKE_OPTIONAL(triggerAdjustForAccrued);
        FIELD(notResetCallTrigger, "Call trigger do not reset");
        FIELD_MAKE_OPTIONAL(notResetCallTrigger);
        FIELD(callOptimality, "call optimality factor");
        FIELD_MAKE_OPTIONAL(callOptimality);
        FIELD(softCallAdjustType, "describes how to adjust barrier for soft-call conditions")
        FIELD_MAKE_OPTIONAL(softCallAdjustType);
        FIELD(softCallDaysM, "M for M out of N day above soft-call trigger condition")
        FIELD_MAKE_OPTIONAL(softCallDaysM);
        FIELD(softCallDaysN, "N for M out of N day above soft-call trigger condition")
        FIELD_MAKE_OPTIONAL(softCallDaysN);
        FIELD(callNotification, "call notification period");
        FIELD_MAKE_OPTIONAL(callNotification);
        FIELD(callTreatment, "K for knockout, B for BlackOnCall");
        FIELD_MAKE_OPTIONAL(callTreatment);

        FIELD(convertIntoIssuerShares, "conversion is into shares of the issuer");
        FIELD_MAKE_OPTIONAL(convertIntoIssuerShares);

        FIELD(instSettle, "settlement of the convertible");
        FIELD_MAKE_OPTIONAL(instSettle);
        FIELD(premiumSettle, "premium settlement -- only different from instSettle for synthetics");
        FIELD_MAKE_OPTIONAL(premiumSettle);
        FIELD(convertOnCall, "conversion only allowed when called");
        FIELD_MAKE_OPTIONAL(convertOnCall);
        FIELD(getCouponAtMat, "do you get the coupon if you convert at maturity");
        FIELD_MAKE_OPTIONAL(getCouponAtMat);
        FIELD(currentConversionAvg, "current average for conversion");
        FIELD_MAKE_OPTIONAL(currentConversionAvg);
        FIELD(currentConversionWeight, "current weight for conversion");
        FIELD_MAKE_OPTIONAL(currentConversionWeight);
        FIELD(isPreferred, "is the convertible a preferred?");
        FIELD_MAKE_OPTIONAL(isPreferred);
        FIELD(riskyGrowth, "use risky growth model");
        FIELD_MAKE_OPTIONAL(riskyGrowth);
        FIELD(makeWholeType, "make-whole type - N for None, C for PV of Coupons, A for Amazon");
        FIELD_MAKE_OPTIONAL(makeWholeType);
        FIELD(makeWholeDate, "date through which make-whole payment is made");
        FIELD_MAKE_OPTIONAL(makeWholeDate);
        FIELD(makeWholeAmount, "make-whole payment is this amount less coupons already paid");
        FIELD_MAKE_OPTIONAL(makeWholeAmount);
        FIELD(alreadyCalled, "has call notification already been sent");
        FIELD_MAKE_OPTIONAL(alreadyCalled);
        FIELD(dateCallNotifSent, "date call notification was sent");
        FIELD_MAKE_OPTIONAL(dateCallNotifSent);

        FIELD(softCallSchedule, "soft call schedule");
        FIELD_MAKE_OPTIONAL(softCallSchedule);
        FIELD(softCallTriggerSchedule, "soft call trigger schedule");
        FIELD_MAKE_OPTIONAL(softCallTriggerSchedule);
        FIELD(callSchedule, "hard call schedule");
        FIELD_MAKE_OPTIONAL(callSchedule);
        FIELD(putSchedule, "put schedule");
        FIELD_MAKE_OPTIONAL(putSchedule);
        FIELD(convCashSchedule, "conversion cash schedule");
        FIELD_MAKE_OPTIONAL(convCashSchedule);
        FIELD(conversionRatios, "conversion ratio schedule");
        FIELD(parityFloorSchedule, "parity floor schedule, allow to convert only when bond is at or below parity floor");
        FIELD_MAKE_OPTIONAL(parityFloorSchedule);

        FIELD(contConversion,	"periodic contingent conversion plus floor provision");
        FIELD_MAKE_OPTIONAL(contConversion);

        FIELD(addOnConvRatios,	"additional conversion ratio schedule and trigger info for embedded option");
        FIELD_MAKE_OPTIONAL(addOnConvRatios);

        FIELD(softPutSchedule, "soft call schedule");
        FIELD_MAKE_OPTIONAL(softPutSchedule);
        FIELD(softPutTriggerSchedule, "soft call trigger schedule");
        FIELD_MAKE_OPTIONAL(softPutTriggerSchedule);

        FIELD(payoffType, "S for straight, D for DECS, P for PERCS");
        FIELD_MAKE_OPTIONAL(payoffType);
        FIELD(DECS, "true for Dividend Enhanced Convertible Securities");
        FIELD_MAKE_TRANSIENT(DECS);
        FIELD(PERCS, "true for Preferred Equity Redemption Cumulative Stocks");
        FIELD_MAKE_TRANSIENT(PERCS);
        FIELD(initialPrice, "initial stock price for DECS and PERCS. Use -1 for face/initialConvRatio.");
        FIELD_MAKE_OPTIONAL(initialPrice);
        FIELD(initialConvRatio, "initial (maximum) conversion ratio for DECS and PERCS.");
        FIELD_MAKE_OPTIONAL(initialConvRatio);
        FIELD(convPrice, "conversion price for DECS. Use -1 for face/minConvRatio.");
        FIELD_MAKE_OPTIONAL(convPrice);
        FIELD(minConvRatio, "minimum conversion ratio for DECS.");
        FIELD_MAKE_OPTIONAL(minConvRatio);

        FIELD(decsHasCutoff, "whether a mandatory has an upper cutoff");
        FIELD_MAKE_OPTIONAL(decsHasCutoff);
        FIELD(decsCutoffLevel, "cutoff level for mandatory");
        FIELD_MAKE_OPTIONAL(decsCutoffLevel);
        FIELD(triggerStartDate, "observation start date for Knock out");
        FIELD_MAKE_OPTIONAL(triggerStartDate);
        FIELD(cappedDECSTrigger, "trigger for early exercise");
        FIELD_MAKE_OPTIONAL(cappedDECSTrigger);

        FIELD(riskFreeCoupons, "are the coupons guaranteed. For DECS only.");
        FIELD_MAKE_OPTIONAL(riskFreeCoupons);
        FIELD(continuousReset, "true if call spread can be exercised at any time (false: conversion ratio from conversion schedule applies before maturity");
        FIELD_MAKE_OPTIONAL(continuousReset);
        FIELD(delayReset, "delay between conversion ratio reset and exercise");
        FIELD_MAKE_OPTIONAL(delayReset);
        FIELD(includePutSpread, "price put spread if intrinsic is above option value");
        FIELD_MAKE_OPTIONAL(includePutSpread);

        FIELD(frenchDivTreatment, "true if don't get a dividend in the year you convert");
        FIELD_MAKE_OPTIONAL(frenchDivTreatment);
        FIELD(frenchExtendedConv, "true if you have extra time to decide whether to convert or take cash");
        FIELD_MAKE_OPTIONAL(frenchExtendedConv);
        FIELD(frenchExtConvInt, "length of French extended conversion period");
        FIELD_MAKE_OPTIONAL(frenchExtConvInt);

        FIELD(schedsArePcts, "call and put schedules are entered as percents of notional.");
        FIELD_MAKE_OPTIONAL(schedsArePcts);

        FIELD(inDefault, "is the bond issuer in default.");
        FIELD_MAKE_OPTIONAL(inDefault);
        FIELD(recoveryPct, "expected recovery upon default as percent of face value.");
        FIELD_MAKE_OPTIONAL(recoveryPct);
        FIELD(payAccruedUponDefault, "whether to settle accrued upon default");
        FIELD_MAKE_OPTIONAL(payAccruedUponDefault);
        FIELD(riskFreeCouponEndDate, "end date for risk free coupons.");
        FIELD_MAKE_OPTIONAL(riskFreeCouponEndDate);

        FIELD(dividendPassThrough, "true means stock dividends are paid as additional coupons?");
        FIELD_MAKE_OPTIONAL(dividendPassThrough);
        FIELD(dividendPassThroughPct, "Percentage of dividends to be passed on");
        FIELD_MAKE_OPTIONAL(dividendPassThroughPct);

        FIELD(dividendAdjusted, "true means conv ratio is increased to have constant parity when a dividend is paid");
        FIELD_MAKE_OPTIONAL(dividendAdjusted);

        FIELD(resetSchedule, "Reset schedule for resettable convertibles");
        FIELD_MAKE_OPTIONAL(resetSchedule);
        FIELD(maxWithParity, "whether the conversion can take place at parity");
        FIELD_MAKE_OPTIONAL(maxWithParity);
        FIELD(manualResetConvPricePct, "conversion price as percentage of spot for manual reset upon conversion");
        FIELD_MAKE_OPTIONAL(manualResetConvPricePct);

        FIELD(triggerReset,                          "true if there is a state dependent reset");
        FIELD_MAKE_OPTIONAL(triggerReset);
        FIELD(resetObservationStartDate,             "state dependent reset observation start date");
        FIELD_MAKE_OPTIONAL(resetObservationStartDate);
        FIELD(resetObservationEndDate,               "state dependent reset observation end date");
        FIELD_MAKE_OPTIONAL(resetObservationEndDate);
        FIELD(resetTrigger,                          "state dependent reset trigger level");
        FIELD_MAKE_OPTIONAL(resetTrigger);
        FIELD(resetConversionRatio,                  "new conversion ratio");
        FIELD_MAKE_OPTIONAL(resetConversionRatio);
        FIELD(resetType,                             "reset type (UP_RESET, DOWN_RESET)");
        FIELD_MAKE_OPTIONAL(resetType);

        FIELD(useAssetRecovery, "whether to use the asset recovery when pricing off CDS curves");
        FIELD_MAKE_OPTIONAL(useAssetRecovery);

        FIELD(accelerateDECS,                "callable DECS feature");
        FIELD_MAKE_OPTIONAL(accelerateDECS);
        FIELD(accelerateDECSTriggerLevel,    "trigger level for callable DECS");
        FIELD_MAKE_OPTIONAL(accelerateDECSTriggerLevel);
        FIELD(accelerateDECSStartDate,       "trigger observation start date");
        FIELD_MAKE_OPTIONAL(accelerateDECSStartDate);
        FIELD(accelerateDECSEndDate,         "trigger observation start date");
        FIELD_MAKE_OPTIONAL(accelerateDECSEndDate);

        FIELD(DEBUG_SOD_Only, "All critical points must fall at the start of day. Default is true.");
        FIELD_MAKE_OPTIONAL(DEBUG_SOD_Only);

        FIELD(stockInEscrow,             "whether stock is held in Escrow account");
        FIELD_MAKE_OPTIONAL(stockInEscrow);

        FIELD(useJointDefaultCorrelation,          "whether to use the default correlation adjustment");
        FIELD_MAKE_OPTIONAL(useJointDefaultCorrelation);
        FIELD(jointDefaultCorrelation,             "default correlation between underlying and issuer");
        FIELD_MAKE_OPTIONAL(jointDefaultCorrelation);
        FIELD(stockCreditSpread,                   "credit spread of the underlying stock");
        FIELD_MAKE_OPTIONAL(stockCreditSpread);

        FIELD(contingentConversion,             "true = contingent conversion trigger");
        FIELD_MAKE_OPTIONAL(contingentConversion);
        FIELD(triggerActiveAtMat,               "true = trigger is active at maturity");
        FIELD_MAKE_OPTIONAL(triggerActiveAtMat);
        FIELD(contingentConversionTrigger,      "Contingent conversion trigger level");
        FIELD_MAKE_OPTIONAL(contingentConversionTrigger);
        FIELD(dontScaleNonPrefPriceBy100, "true = price is not scaled for non preferred bonds");
        FIELD_MAKE_OPTIONAL(dontScaleNonPrefPriceBy100);

        // look up field for use on recurse
        creditCurveField = clazz->getDeclaredField("creditSpreads");
    }

    static IObject* defaultConvBond(){
        return new ConvBond();
    }
};

CClassConstSP const ConvBond::TYPE = CClass::registerClassLoadMethod(
    "ConvBond", typeid(ConvBond), ConvBondHelper::load);
bool   ConvBondLoad() {
    return ( ConvBond::TYPE != 0);
   }



double ConvBond::getStrike() const
{

    double strike;


    if (DECS == true || PERCS == true) {
        strike = bond->getFaceValue()/initialConvRatio;
    } else {
        bool   isConvertible;
        double convRatio;
        double convCash;
        DateTime lastConvDate;
        double putLevel;
        bool isPutable;
        double finalCashFlow;

        finalCashFlow = bond->getRedemption();

        // Use the face value instead if there is no redemption at maturity
        if (Maths::isZero(finalCashFlow)) {
            finalCashFlow = bond->getFaceValue();
        }
            
        finalCashFlow += getAccruedAtDate(bond->getMaturityDate()); // accrued will be coupon level if getCouponAtMat == false or 0 otherwise
        getPutLevel(bond->getMaturityDate(), &isPutable, &putLevel);
        if (isPutable == true) {
            putLevel += (putAdjustForAccrued || getCouponAtMat)?0.0:getAccruedAtDate(bond->getMaturityDate());
            finalCashFlow = Maths::max(finalCashFlow, putLevel);
        }

		bool isSpot = !hasAddOnConvRatios();
		double spotOrLevel = isSpot?0:finalCashFlow;

        getConversionInfo(bond->getMaturityDate(),
					  spotOrLevel,
                      &isConvertible,
                      &convRatio,
                      &convCash,
					  isSpot);

        if (Maths::isZero(convRatio)) {
            // find the last non-zero conv ratio date
            // not the most efficient algorithm but I don't think it will be used much
            lastConvDate = bond->getMaturityDate().rollDate(-1);
            while (lastConvDate >= valueDate) {
                 getConversionInfo(lastConvDate,
					  spotOrLevel,
                      &isConvertible,
                      &convRatio,
                      &convCash,
					  isSpot);
                 if (!Maths::isZero(convRatio)) {
                     break;
                 }
                 lastConvDate = lastConvDate.rollDate(-1);
            }

            if (Maths::isZero(convRatio)) {
                // the bond is no longer convertible and so the strike won't matter
                convRatio = 1.;
            }
        }

        // using redemption is a bit bogus is the lastConvDate is not maturity -- revisit if it becomes an issue
        strike = (finalCashFlow-convCash)/convRatio;
    }
    return strike;
}

bool ConvBond::hasAddOnConvRatios() const
{
	return (!!addOnConvRatios && addOnConvRatios->nbSchedule() );
}

void ConvBond::getConversionInfo(const DateTime& convDate,
                        bool    *isConvertible,
                        double  *ratio,
                        double  *cash,
                        double  *parityFloor) const
{
    static const string method = "ConvBond::getConversionRate";

	if( hasAddOnConvRatios() )
		throw ModelException(method, "Must pass in spot price when convert has AddOnConvRatio");

	getConversionInfo(convDate, 0 /* dummy spot*/, isConvertible, ratio, cash, false, parityFloor);
}

void ConvBond::getConversionInfo(const DateTime& convDate,
                        double  spotOrLevel,
                        bool    *isConvertible,
                        double  *ratio,
                        double  *cash,
                        bool    isSpot,
                        double  *parityFloor) const
{
    static const string method = "ConvBond::getConversionRate";

    if (!!conversionRatios) {
        if (conversionRatios->length() == 0) {
            throw ModelException(method, " Empty conversion ratio schedule - AS");
         //   *isConvertible = true;
         //   *ratio = conversionRatios->lastValue();
        } else if ( convDate >= conversionRatios->firstDate()) {
            if (!dividendAdjusted) {
                *isConvertible = conversionRatios->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                                          0.,
                                                                          bond->getMaturityDate(),
                                                                          convDate,
                                                                          *ratio);
            } else {
                // create the adjusted conversion schedule

                double      conversionRatio = conversionRatios->firstValue();
                double      spot            = asset->getSpot();
                DateTime    initialDate     = conversionRatios->firstDate();
                DateTime    issueDate       = bond->getAccrualStartDate();

                DividendListSP divSchedule = DividendCollector::divsBetweenDates(asset.get(),
                                                                                 issueDate,
                                                                                 issueDate,
                                                                                 bond->getMaturityDate());

                int divCount = 0;
                const DividendArray& divArray = divSchedule->getArray();
                while (divCount < divArray.size() && divArray[divCount].getPayDate() <= convDate) {
                    if ( divArray[divCount].getDivType() != Dividend::AMOUNT) {
                        throw ModelException(method,
                            "All dividends up to the maturity date must be discrete dividends in the dividends adjusted model");
                    }
                    conversionRatio /= (1. - divArray[divCount].getDivAmount()/spot);
                    ++divCount;
                }
                *isConvertible = true;
                *ratio         = conversionRatio;

            }
        } else {
            *isConvertible = false;
            *ratio = 0.0;
        }
    } else {
        throw ModelException(method, " Conversion ratio schedule is NULL");
    }


    if (!!convCashSchedule && convCashSchedule->length() > 0) {
        // don't care about return bool.
        convCashSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                 0.,
                                                 bond->getMaturityDate(),
                                                 convDate,
                                                 *cash);
    } else {
        // no cash upon conversion
        *cash = 0.;
    }

    if( parityFloor )
    {
        if( parityFloorSchedule.get() && parityFloorSchedule->length() > 0 &&
            parityFloorSchedule->interpolateCVBSchedule(
                Schedule::ADJUST_NONE,
                0.,
                bond->getMaturityDate(),
                convDate,
                *parityFloor) )
        {}
        else
        {
            // 100% by default
            *parityFloor = 1.;
        }
    }

    // process the addOnConversion if isConvertible from base conversion schedule
    if( hasAddOnConvRatios() && spotOrLevel > 0 && *isConvertible )
    {
        addOnConvRatios->adjustConvRatio(convDate, isSpot, spotOrLevel, ratio, cash);
    }
}

double ConvBond::getImpliedStrike(const DateTime& convDate,
                            double  level) const
{
    bool    isConvertible;
    double  ratio, cash;
    getConversionInfo(convDate, level, &isConvertible, &ratio, &cash, false/*is not spot*/);

    if( Maths::isZero(ratio) )
        ratio = 1;

    return (level - cash)/ratio;
}

double ConvBond::getMakeWholePayment( const DateTime & callDate ) const
{
    if( makeWholeType == "N" || callDate >= makeWholeDate )
        return 0.;

    if( makeWholeType == "C" )
        return bond->couponsPV( callDate, makeWholeDate, discount.getSP() );

    if( makeWholeType == "A" )
        return Maths::max( 0.,
            makeWholeAmount - bond->couponsPV( DateTime(0,5), callDate, YieldCurveConstSP() ) );

    return 0.;
}

void ConvBond::getCallLevel(const DateTime& callDate,
                            double   fwdVol,
                            bool    *isCallable,
                            bool    *isAdjustedForAccrued,
                            bool    *isHardCall,
                            double  *level) const
{
    static const string method = "ConvBond::getCallLevel";

    *isCallable = false;
    *isAdjustedForAccrued = false;
    *isHardCall = false;
    *level = 0.;

    DateTime callNotificationDate, adjFirstCallDate;

    if ( !!callSchedule && callSchedule->length() > 0) {
        if ( callTreatment == "B" ) {
            adjFirstCallDate = callSchedule->firstDate().rollDate(-callNotification);
        } else {
            adjFirstCallDate = callSchedule->firstDate();
        }
    }

    if ( !(!softCallTriggerSchedule) && softCallTriggerSchedule->length() > 0) {
        if ( callDate < softCallTriggerSchedule->firstDate() ) {
            // before softCalls
            *isCallable = false;
            *isAdjustedForAccrued = false;
            *level = 0.;
        } 
        else if ( !callSchedule ||
                  (callSchedule->length() == 0) ||
                  (callSchedule->length() > 0 && callDate < adjFirstCallDate)) {
            // between softCalls and hardCalls or maturity
            *isCallable = softCallTriggerSchedule->interpolateCVBSchedule(
                Schedule::ADJUST_NONE,
                bond->getRedemption(schedsArePcts), // not adjusted at mat so could pass in anything
                bond->getMaturityDate(),
                callDate,
                *level);

            if (*isCallable == true) {
                if (schedsArePcts == true) {
                    *level *= bond->getNotional(callDate);
                }
                
                double fwdRate;
                // default adjustment
                Actual365F fwdRateDCC;

                if (softCallAdjustType == "N" || softCallAdjustType == "n") {
                    // do nothing
                } else if (softCallAdjustType == "D" || softCallAdjustType == "d") {

                    /********************************************************************
                     * This one corresponds to a 20 out of 30 days soft call
                     * Soft call adjusted by :-
                     * ( 1.0 + adj ) * soft call
                     *
                     * where adj     = ( 3.75 + 0.95 * fwdRate / fwdVol ) * fwdVol / sqrt(tdpy)
                     * where fwdRate = risk free interest rate to maturity of soft call
                     * fwdVol        = volatility
                     * tdpy          = number of trading days per year which is 256.
                     *************************************************************************/
                    fwdRate = discount->fwd(callDate, callDate.rollDate(1), &fwdRateDCC, CompoundBasis::ANNUAL);

                    *level *= (1.0 + (3.75*fwdVol + 0.95*fwdRate)/16.);

                } else if (softCallAdjustType == "T" || softCallAdjustType == "t") {
                    double a, b, n, m, fwdVolDaily;
                    int daysTotal;

                    daysTotal = softCallDaysN;
                    n = (double) softCallDaysN;
                    m = (double) softCallDaysM;

                    fwdVolDaily = fwdVol / sqrt(253.0);
                    /*note that here the number of trading days is 253 */

                    if( daysTotal > 1 ){
                        a = 0.5972 + 1.5688/n;
                        b = 0.8311 - 3.7685/n ;

                        fwdRate = discount->fwd(callDate, callDate.rollDate(1), &fwdRateDCC, CompoundBasis::ANNUAL);

                        *level *=
                            exp( fwdRate * m/365.0 + fwdVolDaily * sqrt(n) * (a * (m/n) + b* (m*m*m)/(n*n*n)) );
                    }
                } else if (softCallAdjustType == "R" || softCallAdjustType == "r") {
                    double a, b, n, m, fwdVolDaily;
                    int daysTotal;

                    daysTotal = softCallDaysN;
                    n = (double) softCallDaysN;
                    m = (double) softCallDaysM;

                    fwdVolDaily = fwdVol / sqrt(253.0);
                    /* note that here the number of trading days is 253 */

                    if( daysTotal > 1 ) {
                        a = 9.11 * n /1000.0;
                        b = sqrt(2.0 * n / Maths::PI) - ( 0.656 + 1.01 * n / 100.0);

                        fwdRate = discount->fwd(callDate, callDate.rollDate(1), &fwdRateDCC, CompoundBasis::ANNUAL);

                        *level *=
                            exp( fwdRate * m/365.0 + fwdVolDaily *(a + b *(2.0 * (m-1.0)/(n-1.0) -1.0)));
                    }
                } else if (softCallAdjustType == "E" || softCallAdjustType == "e") {
                    double fwdVolDaily, m, n;

                    n = (double) softCallDaysN;
                    m = (double) softCallDaysM;

                    fwdVolDaily = fwdVol / sqrt(365.0);
                    /*note that here the number of trading days is 365, until
                      the effect of trading time is incorporated into the model*/

                    if( n == 30 && ( m == 20 || m == 30 )){
                        if ( m == 20 ){
                            *level *= exp(fwdVol * 0.2316);
                        }
                        else if ( m == 30 ){
                            *level *= exp(fwdVol * 0.3393);
                        }
                    }
                    else if (n > 1){
                        *level *= exp( fwdVol *( 0.050959 * (double) (m-1) / sqrt((double) n)
                                       + 0.0008774 * (double) n + 0.03055));
                    }
                } else {
                    // shouldn't be able to get here but no harm in validating again
                    throw ModelException(method, softCallAdjustType + " is not a supported softCallAdjustType");
                }

                double makeWholePayment = getMakeWholePayment( callDate );
                *level += makeWholePayment;

                if( triggerAdjustForAccrued == true
                    && ( makeWholeType != "C" || Maths::isZero( makeWholePayment ) ) )
                {
                    *isAdjustedForAccrued = true;
                    *level += getBondAccruedAtDate(callDate);
                }
                else
                    *isAdjustedForAccrued = false;

                // AM: Add on any conversion cash to the soft call level
                if (!!convCashSchedule && convCashSchedule->length() > 0)  {
                    double      convCash    = 0.0;
                    bool callDateInSchedule = convCashSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                 0.,
                                                 bond->getMaturityDate(),
                                                 callDate,
                                                 convCash);
                    if (callDateInSchedule) {
                       *level += convCash;
                    }
                }
            }
        }
    }

    if ( !(!callSchedule) && callSchedule->length() > 0 && callDate >= adjFirstCallDate ) {
        // need this to avoid writing over soft-call level
        // in hard call period

        // adjust call date for Black on call
        if ( callTreatment == "B" ) {
            callNotificationDate = callDate.rollDate(callNotification);
        } else {
            callNotificationDate = callDate;
        }

        *isCallable = callSchedule->interpolateCVBSchedule(
                                                    Schedule::ADJUST_NONE,
                                                    bond->getRedemption(schedsArePcts),
                                                    bond->getMaturityDate(),
                                                    callNotificationDate,
                                                    *level);

        if (*isCallable == true) {
            if (schedsArePcts == true) {
                *level *= bond->getNotional(callNotificationDate);
            }
            double makeWholePayment = getMakeWholePayment( callDate );
            *level += makeWholePayment;

            *level *= callOptimality;

            if( redemptionAdjustForAccrued == true
                && ( makeWholeType != "C" || Maths::isZero( makeWholePayment ) ) )
            {
                *isAdjustedForAccrued = true;
                *level += getBondAccruedAtDate(callNotificationDate);
            }
            else
                *isAdjustedForAccrued = false;

            *isHardCall = true;
        }
    }
}

void ConvBond::getCallRedemptionLevel(const DateTime& callDate,
                                      bool    *isCallable,
                                      bool    *isAdjustedForAccrued,
                                      double  *level) const
{
    static const string method = "ConvBond::getCallRedemptionLevel";

    *isCallable = false;
    *isAdjustedForAccrued = false;
    *level = 0.;

    if ( !(!softCallTriggerSchedule) && softCallTriggerSchedule->length() > 0) {
        if ( callDate < softCallTriggerSchedule->firstDate() ) {
            // before softCalls
            *isCallable = false;
            *isAdjustedForAccrued = false;
            *level = 0.;
        } else if ( !callSchedule ||
                   (callSchedule->length() == 0) ||
                   (callSchedule->length() > 0 && callDate < callSchedule->firstDate())) {
            // between softCalls and hardCalls or maturity
            *isCallable = softCallSchedule->interpolateCVBSchedule(
                                               Schedule::ADJUST_NONE,
                                               bond->getRedemption(schedsArePcts), // not adjusted at mat so could pass in anything
                                               bond->getMaturityDate(),
                                               callDate,
                                               *level);

            if (schedsArePcts == true) {
                *level *= bond->getNotional(callDate);
            }
        }
    }

    if ( !(!callSchedule) && callSchedule->length() > 0 && callDate >= callSchedule->firstDate()) {
        // need this to avoid writing over soft-call level
        // in hard call period
        *isCallable = callSchedule->interpolateCVBSchedule(
                                                    Schedule::ADJUST_NONE,
                                                    bond->getRedemption(schedsArePcts),
                                                    bond->getMaturityDate(),
                                                    callDate,
                                                    *level);

        if (*isCallable == true) {
            if (schedsArePcts == true) {
                *level *= bond->getNotional(callDate);
            }
        }
    }

    if (*isCallable == true) {
        if (redemptionAdjustForAccrued == true) {
            *isAdjustedForAccrued = true;
            *level += getBondAccruedAtDate(callDate);
        } else {
            *isAdjustedForAccrued = false;
        }
    }
}


void ConvBond::getPutLevel(const DateTime& putDate,
                           bool    *isPutable,
                           double  *level) const
{
    static const string method = "ConvBond::getPutLevel";

    if (!!putSchedule && putSchedule->length() > 0) {
        *isPutable = putSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                         bond->getRedemption(schedsArePcts),
                                                         bond->getMaturityDate(),
                                                         putDate,
                                                         *level);

        if (*isPutable == true) {
            if (schedsArePcts == true) {
                *level *= bond->getNotional(putDate);
            }

            if (putAdjustForAccrued == true) {
                *level += getBondAccruedAtDate(putDate);
            }
        }
    } else {
        *isPutable = false;
        *level     = 0.0;
    }
}

void ConvBond::getSoftPutLevel(const DateTime& putDate,
                               bool    *isPutable,
                               double  *trigger,
                               double  *level) const
{
    static const string method = "ConvBond::getSoftPutLevel";

    if (!!softPutSchedule && softPutSchedule->length() > 0) {
        *isPutable = softPutSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                         bond->getRedemption(schedsArePcts),
                                                         bond->getMaturityDate(),
                                                         putDate,
                                                         *level);

        if (*isPutable == true) {
            *isPutable = softPutTriggerSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                                        bond->getRedemption(schedsArePcts),
                                                                        bond->getMaturityDate(),
                                                                        putDate,
                                                                        *trigger);

            if (schedsArePcts == true) {
                *level *= bond->getNotional(putDate);
            }

            if (putAdjustForAccrued == true) {
                *level += getBondAccruedAtDate(putDate);
            }
        }
    } else {
        *isPutable = false;
        *level     = 0.0;
    }
}



/** return the next put date and level if possible. Returns canCalcYield = false if it is putable today */
void ConvBond::getFirstPutInfo(bool *canCalcYield, DateTime *firstPutDate, double *firstPutLevel) const
{
    static const string method = "ConvBond::getFirstPutDate";

    bool isPutable;

    *canCalcYield = false; // until we know it's true
    *firstPutDate = DateTime();
    *firstPutLevel = 0.;

    if (!putSchedule || putSchedule->length() == 0) {
        // keep what's set above
    } else {

        double tmpPutLevel;
        getPutLevel(valueDate, &isPutable, &tmpPutLevel);
        if (isPutable == true) {
            // we can't calc the yield put we know the date and level
            *canCalcYield = false;
            *firstPutDate = valueDate;
            *firstPutLevel = tmpPutLevel;
        } else if (putSchedule->getInterp() == Schedule::INTERP_NONE) {
            const DateTimeArray& putDates = putSchedule->getDateArray();
            int i;
            for (i=0; i<putDates.size(); i++) {
                if (putDates[i] > valueDate) {
                    *canCalcYield = true;
                    *firstPutDate = putDates[i];
                    getPutLevel(*firstPutDate, &isPutable, firstPutLevel);
                    break;
                }
            }
        } else { // we're stairs or linear and we're only putable in the future
            *canCalcYield = true;
            *firstPutDate = putSchedule->firstDate();
            getPutLevel(*firstPutDate, &isPutable, firstPutLevel);
        }
    }
}

double ConvBond::getAccruedAtDate(const DateTime& aiDate) const
{
    static const string method = "ConvBond::getAccruedAtDate";
    double ai;

    ai = getBondAccruedAtDate(aiDate);

    /*
    if (isPreferred == true) {
        ai = 0.;
    } else {
        ai = getBondAccruedAtDate(aiDate);
    }
    */
    return ai;
}

double ConvBond::getBondAccruedAtDate(const DateTime& aiDate) const
{
    static const string method = "ConvBond::getBondAccruedAtDate";
    double ai;

    if (getCouponAtMat == true) {
        ai = bond->getAccruedAtDate(aiDate);
    } else {
        // only compare the dates and not the times
        if (aiDate.getDate() == bond->getMaturityDate().getDate()) {
            ai = bond->getMaturityCoupon();
        } else {
            ai = bond->getAccruedAtDate(aiDate);
        }
    }
    return ai;
}

CashFlowArraySP ConvBond::getCashFlows(const DateTime& startDate) const
{
    static const string method = "ConvBond::getCashFlows";

    return bond->getCashFlows(startDate);
}

CashFlowArraySP ConvBond::getCoupons(const DateTime& startDate) const
{
    static const string method = "ConvBond::getCoupons";

    return bond->getCoupons(startDate);
}

/** Part of IInstrumentAsAsset interface.
    Returns the dates on which the instrument has to be held in order to
    hold the right to the corresponding coupon as returned by getCoupons() */
DateTimeArraySP ConvBond::getExCouponDates() const{
    return bond->getExCouponDates();
}

double ConvBond::yieldToMaturity(double price, bool clean) const
{
    return bond->yieldToMaturity(price, clean, valueDate);
}

double ConvBond::yieldToFirstPut(double price, bool clean) const
{
    static const string method = "ConvBond::yieldToFirstPut";

    try {
        bool canCalcYield;
        DateTime firstPutDate;
        double firstPutLevel;

        getFirstPutInfo(&canCalcYield, &firstPutDate, &firstPutLevel);

        if (canCalcYield == true) {
            return bond->yieldToFirstPut(price, clean, firstPutDate, firstPutLevel, valueDate);
        } else {
            throw ModelException(method, "can't calc yield to first put");
        }

    } catch (exception& e) {
        throw ModelException(e, method);
    }
    return true; // our components have theta type sensitivity

}

bool ConvBond::hasPut() const
{
    return (!(!putSchedule) && putSchedule->length() > 0);
}


void ConvBond::recordOutputRequests(Control* control, Results* results,
                                    double fairValue, double bondFloor,
                                    bool isAnOCB) const
{

    static const string method = "ConvBond::recordOutputRequests";
    OutputRequest* request = NULL;
    double cleanPrice;

    if ( control->isPricing() ) {

        // Bond value
        if ( control->requestsOutput(OutputRequest::NAKED_BOND_PRICE, request) ) {
             results->storeRequestResult(request, bondFloor);
        }

        // Accrued interest
        if ( control->requestsOutput(OutputRequest::ACCRUED_INTEREST, request) ) {
            if (isPreferred) {
                results->storeRequestResult(request, 0.0);
            } else {
                // accrued should be to the settlement date
                DateTime aiDate = bond->settles(valueDate);               
                results->storeRequestResult(request, getAccruedAtDate(aiDate));
            }
        }

        if ( control->requestsOutput(OutputRequest::PARITY, request) ) {
            bool isConvertible;
            double convRatio;
            DateTime nextConvDate;

           if ( (!!conversionRatios && conversionRatios->length() > 0) &&
                 conversionRatios->getNextDate(valueDate,
                                               bond->getMaturityDate(),
                                               nextConvDate)) {
                isConvertible = conversionRatios->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                                         0.,
                                                                         bond->getMaturityDate(),
                                                                         nextConvDate,
                                                                         convRatio);
                                // if resettable get the current ratio
                 if ( (!!resetSchedule && resetSchedule->length() > 0)) {
                    convRatio = resetSchedule->getCurrentConversionRatio(conversionRatios->firstValue(),
                                                             valueDate,
                                                             bond->getFaceValue());
                }

           } else if ( DECS == true ) {
               // no conversion schedule, so use minConvRatio instead
               // the parity calculation can be overriden if a conversion ratio schedule is entered
               // the conversion ratio schedule
               isConvertible = true;
               convRatio     = minConvRatio;
           } else {
               isConvertible = false;
               convRatio     = 0.0;
           }

           if(!!resetSchedule && resetSchedule->isFwdStart())
                results->storeRequestResult(request, convRatio);
           else 
                results->storeRequestResult(request, asset->getSpot() * convRatio);
        }

        // PARITY_CONVERSION_RATIO
        if ( control->requestsOutput(OutputRequest::PARITY_CONVERSION_RATIO, request) ) {
            double      convRatio   = 0.0;
            double      convCash    = 0.0;
            bool        isConvertible;
            DateTime    nextConvDate;

            if ( DECS == true ) {
               // no conversion schedule, so use initialConvRatio instead
               convRatio     = initialConvRatio;
            } else if (!!resetSchedule) {
                ResetScheduleSP adjustedResetSchedule(resetSchedule.clone());

                // convert the schedule to bond currency for currency struck bonds
                if ( StruckEquity::TYPE->isInstance(asset.get())) {
                    const StruckEquity* struckAsset = dynamic_cast<const StruckEquity*>(asset.get());

                    double       fxSpot      = struckAsset->getFXSpot();

                    // get all reset dates
                    DateTimeArray resetDates = adjustedResetSchedule->getDates();
                    DoubleArraySP fxFwds(new DoubleArray(resetDates.size()));
                    for(int i=0 ; i<fxFwds->size() ; ++i) {
                        if (resetDates[i] >= valueDate ) {
                            (*fxFwds)[i] = struckAsset->fxFwdValue(resetDates[i]);
                        } else {
                            (*fxFwds)[i] = 0.0;
                        }
                    }

                    // calculate initial conversion price
                    double initialCR        = conversionRatios->firstValue();
                    if (Maths::isZero(initialCR)) {
                        initialCR = DBL_EPSILON;
                    }
                    double initialConvPrice = bond->getFaceValue() / initialCR;

                    // convert the reset schedule into a currency struck schedule, if necessary
                    adjustedResetSchedule->preProcessSchedule(true, fxSpot, fxFwds, initialConvPrice);
                } else {
                    double       fxSpot      = 1.0;

                    // get all reset dates
                    DateTimeArray resetDates = adjustedResetSchedule->getDates();
                    DoubleArraySP fxFwds(new DoubleArray(resetDates.size()));
                    for(int i=0 ; i<fxFwds->size() ; ++i) {
                        (*fxFwds)[i] = 1.0;
                    }

                    // calculate initial conversion price
                    double initialCR        = conversionRatios->firstValue();
                    if (Maths::isZero(initialCR)) {
                        initialCR = DBL_EPSILON;
                    }
                    double initialConvPrice = bond->getFaceValue() / initialCR;

                    // convert the reset schedule into a currency struck schedule, if necessary
                    adjustedResetSchedule->preProcessSchedule(false, fxSpot, fxFwds, initialConvPrice);
                }

                convRatio = adjustedResetSchedule->getCurrentConversionRatio(
                                         conversionRatios->firstValue(),
                                         valueDate,
                                         bond->getFaceValue());
            } else if ( (!!conversionRatios && conversionRatios->length() > 0) &&
                 conversionRatios->getNextDate(valueDate,
                                               bond->getMaturityDate(),
                                               nextConvDate)) {
                double stockPrice = asset->fwdValue(nextConvDate);
                getConversionInfo(nextConvDate, stockPrice, &isConvertible, &convRatio, &convCash);
           } else {
               isConvertible = false;
               convRatio     = 0.0;
           }

           if(!!resetSchedule && resetSchedule->isFwdStart()){
               double refLvl = resetSchedule->getStartLevel(valueDate, 
                                                            asset->fwdValue(resetSchedule->getStartDate(valueDate)));
               results->storeRequestResult(request, convRatio/refLvl);
           }
           else 
               results->storeRequestResult(request, convRatio);

        }

        // THETA_ACCRUED_INTEREST
        double thetaAccInt = 0.0;
        if ( control->requestsOutput(OutputRequest::THETA_ACCRUED_INTEREST, request) ) {
            if (isPreferred) {
                thetaAccInt = 0.0;
            } 
            else {
                // determine next buisness day
                HolidaySP hols(Holiday::weekendsOnly());
                DateTime tomorrow = hols->addBusinessDays(getValueDate(), 1);
                double accruedToday      = getAccruedAtDate(getValueDate());
                double accruedTomorrow   = getAccruedAtDate(tomorrow);
                double couponsPV         = bond->couponsPV(getValueDate(),
                                                           tomorrow,
                                                           YieldCurveConstSP());
                thetaAccInt = accruedTomorrow - accruedToday + couponsPV;
            }
            results->storeRequestResult(request,thetaAccInt);
        }

        // IND_CDS_PAR_SPREAD (and deprecated CURRENT_SPREAD)

        if (control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD) ||
            control->requestsOutput(OutputRequest::CURRENT_SPREAD)) {
            IObjectSP result = creditSpreads->getCurrentSpreadOrUntweakable(
                                   getValueDate(), bond->getMaturityDate());
            request = control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD);
            if (request) {
                results->storeRequestResult(
                    request, result,
                    OutputNameConstSP(new OutputName(creditSpreads.getName())));
            }
            request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
            if (request) {
                results->storeRequestResult(request, result);
            }
        }

        // COUPON_DUE
		request = control->requestsOutput(OutputRequest::COUPON_DUE);
		if (request) {
		    AccrualCalendarArraySP coupons = bond->couponDue(valueDate);
		    if (coupons->empty()) {
		        results->storeNotApplicable(request);
		    } else {
		        results->storeRequestResult(request, coupons);
		    }
		}

		// ACCRUAL_CALENDAR
		request = control->requestsOutput(OutputRequest::ACCRUAL_CALENDAR);
		if (request) {
		    AccrualCalendarArraySP accData = bond->accrualCalendar(valueDate);
		    if (accData->empty()) {
		        results->storeNotApplicable(request);
		    } else {
		        results->storeRequestResult(request, accData);
		    }
        }


        if (!isAnOCB) {
            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                                             results,
                                             fairValue,
                                             valueDate,
                                             discount.get(),
                                             asset.get(),
                                             premiumSettle.get());

            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           bond->getMaturityDate(),
                                           valueDate,
                                           asset.get());

            // Option price
            if ( control->requestsOutput(OutputRequest::OPTION_PRICE, request) ) {
                results->storeRequestResult(request, fairValue - bondFloor);
            }

            // Dirty Price
            if ( control->requestsOutput(OutputRequest::DIRTY_PRICE, request) ) {
                results->storeRequestResult(request, fairValue);
            }

            // calculate clean price and output into DEBUG packet
            if (isPreferred) {
                cleanPrice = fairValue;
            } else {
                cleanPrice = fairValue - getAccruedAtDate(valueDate);
            }
		    OutputNameConstSP cleanPriceOutput(new OutputName(OutputRequest::CLEAN_PRICE_FOR_SCALING));
		    results->storeScalarGreek(cleanPrice, Results::DEBUG_PACKET, cleanPriceOutput);

            // Clean price
            if ( control->requestsOutput(OutputRequest::CLEAN_PRICE, request) ) {
               results->storeRequestResult(request, cleanPrice);
            }

            // STRIKE
            if ( control->requestsOutput(OutputRequest::STRIKE, request) ) {
               double strike = getStrike();

               // convert strike to underlying ccy for struck bonds

                if (StruckEquity::TYPE->isInstance(asset.getSP())) {
                    // cast to struck equity object
                    const IObject* obj           = dynamic_cast<const IObject*>(asset.get());
                    const StruckEquity* struckEq = dynamic_cast<const StruckEquity*>(obj);
                    if ( struckEq ) {
                        double spotFX = struckEq->getFX()->getSpot();
                        strike /= spotFX;
                    }
                }

                // this is a bit dodgy for closed form DECS
                results->storeRequestResult(request,strike);
            }

            // YTM
            if ( DECS == false && PERCS == false &&
                 control->requestsOutput(OutputRequest::YIELD_TO_MATURITY, request) ) {
                try { // don't fail if this doesn't work
                    double ytm = yieldToMaturity(fairValue, false);
                    results->storeRequestResult(request, ytm);
                } catch (exception& e) {
                    results->storeRequestResult(request, IObjectSP(new Untweakable(e)));
                    /*
                    results->storeNotApplicable(request);
                    ModelException myE = ModelException(e);
                    myE.addMsg(method + ": Failed to compute YIELD_TO_MATURITY. Continuing...\n");
                    myE.errorLog();
                    */
                }
            }

            // YTP
            if ( control->requestsOutput(OutputRequest::YIELD_TO_FIRST_PUT, request) ) {
                try { // don't fail if this doesn't work
                    bool     canCalcYield;
                    DateTime firstPutDate;
                    double   firstPutLevel;
                    getFirstPutInfo(&canCalcYield, &firstPutDate, &firstPutLevel);

                    if ( canCalcYield ) {
                        double ytp = yieldToFirstPut(fairValue, false);
                        results->storeRequestResult(request, ytp);
                    } else {
                        results->storeNotApplicable(request);
                    }
                } catch (exception& e) {
                    results->storeRequestResult(request, IObjectSP(new Untweakable(e)));
                    /*
                    results->storeNotApplicable(request);
                    ModelException myE = ModelException(e);
                    myE.addMsg(method + ": Failed to compute YIELD_TO_FIRST_PUT. Continuing...\n");
                    myE.errorLog();
                    */
                }
            }

            // KNOWN_CASHFLOWS
            if ( control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS) ) {
                CashFlowArraySP allCFs = bond->getKnownCashFlows();

                // exclude coupon payment if necessary
                if ( allCFs->size() > 0 ) {
                    if ( getCouponAtMat == false ) {
                        (*allCFs)[allCFs->size()-1].amount   -= ((*allCFs)[allCFs->size()-1].amount - bond->getRedemption());
                    }

                    // exclude redemption payment for mandatories
                    if ( isPreferred || DECS == true || PERCS == true ) {
                        (*allCFs)[allCFs->size()-1].amount   -= bond->getRedemption();
                    }


                    if (!Maths::isPositive((*allCFs)[allCFs->size()-1].amount)) {
                        allCFs->pop_back();
                    }
                }

                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        allCFs.get());
            }

            // PAYMENT_DATES
            if ( control->requestsOutput(OutputRequest::PAYMENT_DATES) ) {
                DateTimeArraySP dates = bond->getPaymentDates();

                if ( dates->size() > 0 ) {
                    // exclude redemption payment for mandatories
                    if ( isPreferred && getCouponAtMat == false ) {
                        dates->pop_back();
                    }
                }

                OutputRequestUtil::recordPaymentDates(control,results,dates.get());
            }

            if ( control->requestsOutput(OutputRequest::ASSET_SWAP_SPREAD, request) ) {

                try { // don't fail if this doesn't work
                    double bondPrice = bondFloor - getAccruedAtDate(valueDate);
                    double aswSpread = getAssetSwapSpread(bondPrice);
                    results->storeRequestResult(request, aswSpread);
                } catch (exception& e) {
                    results->storeRequestResult(request, IObjectSP(new Untweakable(e)));
                }
            }

            /*
            if ( control->requestsOutput(OutputRequest::Z_SPREAD) ) {
                double zSpread = getZSpread(bondFloor);
               results->storeRequestResult(request, aswSpread);
            }
            */
        }
        
        if ( hasContConversion() )
        {
            if( control->requestsOutput(OutputRequest::CONT_CONVERT_STATUS, request) )
            {
                bool flag = contConversion->isHistEnabled(getValueDate());
                results->storeRequestResult(request, flag);
            }
        }

    }
}

YieldCurveSP ConvBond::getRiskyCurve() const
{
    // create a risky curve
    IYieldCurveSP riskyCurve = creditSpreads.get()->makeRiskyCurve(*discount.get());

    if (!riskyCurve.get())
	{
        throw ModelException("ConvBond::getRiskyCurve", "failed to create risky curve");
    }

    return YieldCurveSP(dynamic_cast<YieldCurve*>(riskyCurve.get()));
}

YieldCurveSP ConvBond::getRiskyCurve(const DateTime* maturityDate) const
{
    // create a risky curve
    IYieldCurveSP riskyCurve = creditSpreads.get()->makeRiskyCurve(*discount.get(), maturityDate);

    if (!riskyCurve.get())
	{
        throw ModelException("ConvBond::getRiskyCurve", "failed to create risky curve");
    }

    return YieldCurveSP(dynamic_cast<YieldCurve*>(riskyCurve.get()));
}

YieldCurveSP ConvBond::getRiskyCurve(const DateTime* maturityDate, bool useJointDdefProb) const
{
    // create a risky curve
    IYieldCurveSP riskyCurve;
    if ( !useJointDefaultCorrelation ) {
        riskyCurve = creditSpreads.get()->makeRiskyCurve(*discount.get(), maturityDate);
    } else {
        // adjust the credit curve for the joint default probability
        CreditSpreadCurveSP csCurve(CreditSpreadCurveSP::dynamicCast((IObjectSP)creditSpreads->clone()));
        csCurve->addSpread(stockCreditSpread*(1.- jointDefaultCorrelation));
        riskyCurve = csCurve->makeRiskyCurve(*discount.get(), maturityDate);
    }

    if (!riskyCurve.get())
	{
        throw ModelException("ConvBond::getRiskyCurve", "failed to create risky curve");
    }

    return YieldCurveSP(dynamic_cast<YieldCurve*>(riskyCurve.get()));
}



/* now value the bond if it has been called. Use the same european approximation for what is
   actually an american that is used in black on call itself. most of the routines are reusable,
   but in this case we don't need to value calls at each node of a slice and we have to remember
   to include puts in the valuation. also unlike in black on call performance is not really an
   issue, so value a call for each day until call maturity. remember that call maturity may equal
   bond maturity, so have to cope with that below */
double ConvBond::alreadyCalledBlackApprox(Control* control, Results* results) const
{
    static const string method = "ConvBond::alreadyCalledBlackApprox";

    YieldCurveConstSP    optDiscZero;
    DateTime             callMatDate, endDate, callValuationDate;
    double               makeWholePayment, couponsThroughToday, couponsBeforeCall, makeWholeBeforeCall;
    double               couponsAfterCall, makeWholeAfterCall, putLevel, strike, convRatio, convCash;
    double               optDiscFactor, callValue, fwdPrice, variance, indicativeVol;
    double               riskyDiscFactorAfter, riskyDiscFactorBefore;
    bool                 isPuttable, isConvertible;
    double               cleanPrice;

    double alreadyCalledPrice = 0.0;

    YieldCurveConstSP riskyCurve = getRiskyCurve();

    // setup option valuation discount zero curve. set to risk free zero curve
    // (equity's growth curve) if issuer issues shares on conversion and risky
    // curve otherwise. NOTE - risky zero always used for cash discounting
    if ( convertIntoIssuerShares == true ) {
        optDiscZero = discount.getSP();
    } else {
        optDiscZero = riskyCurve;
    }

    callMatDate = dateCallNotifSent.rollDate(callNotification);

    // french conversion needs fixing
    if ( frenchExtendedConv == true ) {
        // throw ModelException(method, "AS: french extended conversion not implemented yet!");
        endDate = frenchExtConvInt->toDate(callMatDate);
    } else {
        endDate = callMatDate;
    }

    // get accrued interest at maturity of call (not call date) if required.
    // NOTE - if maturity of call is bond maturity this routine
    // takes account of get coupon at maturity so always
    // call it in that case

    bool    isCallable, isAdjustedForAccrued;
    double  callLevel;

    getCallRedemptionLevel(callMatDate,&isCallable,&isAdjustedForAccrued,&callLevel);

    // no adjusting for accrued during the Make-Whole period for eCouponDate style
    if ( makeWholeType == "C" && callMatDate <= makeWholeDate) {
        if ( endDate == bond->getMaturityDate() || isAdjustedForAccrued ) {
            callLevel -= getAccruedAtDate(endDate);
        }
    }

    // no adjusting for accrued during the Make-Whole period for eCouponDate style
    if (  (makeWholeType == "N" || makeWholeType == "A" || (makeWholeType == "C" && callMatDate > makeWholeDate)) &&
           !isAdjustedForAccrued && callMatDate == bond->getMaturityDate() ) {
        if ( endDate == bond->getMaturityDate() || isAdjustedForAccrued ) {
            callLevel += getAccruedAtDate(endDate);
        }
    }


    /* get the value of the make-whole payment */
    if (makeWholeType == "C") {
        // Get the PV of the coupons from today until the end of the make-whole period.
        makeWholePayment = bond->couponsPV(getValueDate(),
                                           makeWholeDate,
                                           discount.getSP());
    } else if (makeWholeType == "A") {
        makeWholePayment = 0;
        if ( valueDate < makeWholeDate ) {
            DateTime anAncientDate(0,5);
            couponsThroughToday = bond->couponsPV(anAncientDate,
                                                  valueDate,
                                                  YieldCurveConstSP());

            // the make-whole payment can never be negative.
            makeWholePayment = Maths::max(makeWholeAmount - couponsThroughToday,0.0);
        }
    } else { // no make-whole
        makeWholePayment = 0.;
    }

    // get value of each call. overall call worth MAX(calls)
    for ( callValuationDate = valueDate ; callValuationDate <= endDate ;
          callValuationDate = callValuationDate.rollDate(1)) {

         // AS: This is still wrong for French extended conversion
         optDiscFactor          = optDiscZero->pv(valueDate, callValuationDate);
         riskyDiscFactorBefore  = riskyCurve->pv(valueDate, callValuationDate);
         if ( callValuationDate > callMatDate ) {
             riskyDiscFactorAfter   = 1.0;
         } else {
             riskyDiscFactorAfter   = riskyCurve->pv(callValuationDate,callMatDate);
         }


        // get value today of coupons paid after today, but on or before call date (couponsBeforeCall).
        // NOTE - if callMatDate is bond maturity this takes account of get coupon at mat.
        //
        // NOTE - coupons always discounted at risky rate.
        //
        // NOTE - if valuing on bond maturity and get coupon at maturity is TRUE always get
        //        last coupon. note that in this case callAccrued above is zero so that last
        //        coupon does not affect conversion decision.
        if ( valueDate == bond->getMaturityDate() && getCouponAtMat ) {
            couponsBeforeCall = bond->getMaturityCoupon();
        } else {
            if (  makeWholeType == "N" ||
                 (makeWholeType == "C" && callValuationDate >= makeWholeDate) ||
                  makeWholeType == "A") {
                 DateTime dateToUse = (callValuationDate<callMatDate)?callValuationDate:callMatDate;
                 couponsBeforeCall = bond->couponsPV(valueDate, dateToUse, discount.getSP());

                if (makeWholeType == "C") {
                     makeWholeBeforeCall = bond->couponsPV(valueDate, makeWholeDate, discount.getSP());
                    couponsBeforeCall -= makeWholeBeforeCall;
                }
            }

            if ( makeWholeType == "N" ||
                (makeWholeType == "C" && makeWholeDate <= callMatDate && makeWholeDate <= callValuationDate) ||
                 makeWholeType == "A") {
                // get value on call date of coupons paid after call date,
                // but on or before call maturity (couponsAfterCall).
                DateTime dateToUse = (callValuationDate<callMatDate)?callValuationDate:callMatDate;
                couponsAfterCall = bond->couponsPV(callValuationDate, dateToUse, discount.getSP());

                if (makeWholeType == "C" &&
                    makeWholeDate > callValuationDate) {
                    makeWholeAfterCall = bond->couponsPV(callValuationDate, makeWholeDate, discount.getSP());
                    couponsAfterCall -= makeWholeAfterCall;
                }
            }
        }

        // get put level at call date if any puts. takes account of accrued interest on
        // call date (not call maturity)if required. can call this routine unlike calls
        // equivalent as accrued is paid same day as put

        getPutLevel(callValuationDate, &isPuttable, &putLevel);

        // get strike of call on call date. strike actually used below is in stock units
        //
        // NOTE - call level and accrued always discounted at
        //        risky rate to account for default risk.
        strike =
            callLevel*riskyDiscFactorAfter + couponsAfterCall;

        // get conversion ratio and cash
        bool isSpot = !hasAddOnConvRatios();
        double spotOrLevel = isSpot?0:Maths::max(putLevel, strike);
        getConversionInfo(callValuationDate, spotOrLevel, &isConvertible, &convRatio, &convCash, isSpot);

        // must take puts and conversion cash into account.
        // remember that both of those are paid at call date
        strike = Maths::max(putLevel, strike) - convCash;

        // handle zero conversion ratio
        if (Maths::isZero(convRatio)) {
            convRatio = 1.0;
            callValue = Maths::max(-strike, 0.0) * optDiscFactor;
        } else {
            // AS: this is truly awful - I should get all forwards in one go!
            fwdPrice = asset->fwdValue(callValuationDate);

            CVolRequestSP volRequest( new LinearStrikeVolRequest(
                                            strike/convRatio,
                                            valueDate,
                                            callValuationDate,
                                            false));

            LinearStrikeVolRequest* lsVolRequest = dynamic_cast<LinearStrikeVolRequest*>(volRequest.get());
            CVolProcessedBSSP volBS(asset->getProcessedVol(lsVolRequest));
            variance = volBS->CalcVar(valueDate, callValuationDate);

            if ( callValuationDate == endDate && valueDate < callValuationDate ) {
                indicativeVol = volBS->CalcVol(valueDate, callValuationDate);
            } else {
                // aproximate instantaneous vol
                DateTime nextDate = valueDate.rollDate(1);
                indicativeVol = volBS->CalcVol(valueDate, nextDate);
            }

            // call Black to value option
            callValue = Black::price(true /*isCall*/,
                                     fwdPrice,
                                     strike/convRatio,
                                     optDiscFactor,
                                     variance);
        }


        // NOTE - as above when strike is calculated, cash is discounted here at risky rate
        // to account for default risk --> at high stock levels when conversion is most likely
        // the option disc factor will be dominant and at low stock levels when conversion is
        // unlikely the risky disc factor is dominant which is required behaviour.
        callValue =
            convRatio * callValue + couponsBeforeCall +
            (strike + convCash)*riskyDiscFactorBefore;
        // returned price does not include the makeWolePayment when already called. This is
        // because the payment is made when the call notification is sent.

        alreadyCalledPrice = Maths::max(alreadyCalledPrice, callValue);
    }

    results->storePrice(alreadyCalledPrice, discount->getCcy());

    if ( control->isPricing() ) {
        OutputRequest* request = NULL;

        // Indicative vol
        if ( control->requestsOutput(OutputRequest::IND_VOL, request) ) {
             results->storeRequestResult(request, indicativeVol);
        }

        // Bond value
        if ( control->requestsOutput(OutputRequest::NAKED_BOND_PRICE, request) ) {
             results->storeRequestResult(request, alreadyCalledPrice);
        }

        // Accrued interest
        double accruedInterest = 0.0;
        if ( control->requestsOutput(OutputRequest::ACCRUED_INTEREST, request) ) {
            if (isPreferred) {
                accruedInterest = 0.0;
            } else {
                if ( valueDate < callMatDate ) {
                    accruedInterest = getAccruedAtDate(valueDate);
                } else {
                    accruedInterest = getAccruedAtDate(callMatDate);
                }
            }
            results->storeRequestResult(request, accruedInterest);
        }

        // Option price
        if ( control->requestsOutput(OutputRequest::OPTION_PRICE, request) ) {
            results->storeRequestResult(request, 0.0);
        }

        // Dirty Price
        if ( control->requestsOutput(OutputRequest::DIRTY_PRICE, request) ) {
            results->storeRequestResult(request, alreadyCalledPrice);
        }

        if (isPreferred) {
            cleanPrice = alreadyCalledPrice;
        } else {
            cleanPrice = alreadyCalledPrice - accruedInterest;
        }
        OutputNameConstSP cleanPriceOutput(new OutputName(OutputRequest::CLEAN_PRICE_FOR_SCALING));
	    results->storeScalarGreek(cleanPrice, Results::DEBUG_PACKET, cleanPriceOutput);

        // CLEAN_PRICE
        if ( control->requestsOutput(OutputRequest::CLEAN_PRICE, request) ) {
           results->storeRequestResult(request, cleanPrice);
        }

        // STRIKE
        if ( control->requestsOutput(OutputRequest::STRIKE, request) ) {
           double strike = getStrike();

           // convert strike to underlying ccy for struck bonds

            if (StruckEquity::TYPE->isInstance(asset.getSP())) {
                // cast to struck equity object
                const IObject* obj           = dynamic_cast<const IObject*>(asset.get());
                const StruckEquity* struckEq = dynamic_cast<const StruckEquity*>(obj);
                if ( struckEq ) {
                    double spotFX = struckEq->getFX()->getSpot();
                    strike /= spotFX;
                }
            }
            results->storeRequestResult(request,strike);
        }


        if ( control->requestsOutput(OutputRequest::PARITY, request) ) {
            bool isConvertible;
            double convRatio;
            DateTime nextConvDate;

           if ( (!!conversionRatios && conversionRatios->length() > 0) &&
                 conversionRatios->getNextDate(valueDate,
                                               bond->getMaturityDate(),
                                               nextConvDate)) {
                isConvertible = conversionRatios->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                                         0.,
                                                                         bond->getMaturityDate(),
                                                                         nextConvDate,
                                                                         convRatio);
                // if resettable get the current ratio
                if ( (!!resetSchedule && resetSchedule->length() > 0)) {
                    convRatio = resetSchedule->getCurrentConversionRatio(conversionRatios->firstValue(),
                                                             valueDate,
                                                             bond->getFaceValue());
                }
                                              
           } else if ( DECS == true ) {
               // no conversion schedule, so use minConvRatio instead
               // the parity calculation can be overriden if a conversion ratio schedule is entered
               // the conversion ratio schedule
               isConvertible = true;
               convRatio     = minConvRatio;
           } else {
               isConvertible = false;
               convRatio     = 0.0;
           }

            results->storeRequestResult(request, asset->getSpot() * convRatio);
        }

        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         alreadyCalledPrice,
                                         valueDate,
                                         discount.get(),
                                         asset.get(),
                                         premiumSettle.get());

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       bond->getMaturityDate(),
                                       valueDate,
                                       asset.get());


        // YTM
        if ( DECS == false && PERCS == false &&
             control->requestsOutput(OutputRequest::YIELD_TO_MATURITY, request) ) {
            try { // don't fail if this doesn't work
                double ytm = yieldToMaturity(alreadyCalledPrice, false);
                results->storeRequestResult(request, ytm);
            } catch (exception& e) {
                results->storeRequestResult(request, IObjectSP(new Untweakable(e)));
            }
        }

        // YTP
        if ( control->requestsOutput(OutputRequest::YIELD_TO_FIRST_PUT, request) ) {
            try { // don't fail if this doesn't work
                bool     canCalcYield;
                DateTime firstPutDate;
                double   firstPutLevel;
                getFirstPutInfo(&canCalcYield, &firstPutDate, &firstPutLevel);

                if ( canCalcYield ) {
                    double ytp = yieldToFirstPut(alreadyCalledPrice, false);
                    results->storeRequestResult(request, ytp);
                } else {
                    results->storeNotApplicable(request);
                }
            } catch (exception& e) {
                results->storeRequestResult(request, IObjectSP(new Untweakable(e)));
            }
        }

        // current spread
        if ( control->requestsOutput(OutputRequest::CURRENT_SPREAD, request) ) {
            double currentSpread = 0.0;
            try {
                currentSpread = creditSpreads->getCurrentSpread(getValueDate(),
                                                                bond->getMaturityDate());
                results->storeRequestResult(request,currentSpread);
            } catch (exception& e) {
                UntweakableSP untweakableSpread(new Untweakable(e));
                results->storeRequestResult(request,untweakableSpread);
            }
        }

        // THETA_ACCRUED_INTEREST
        double thetaAccInt = 0.0;
        if ( control->requestsOutput(OutputRequest::THETA_ACCRUED_INTEREST, request) ) {
            // determine next buisness day
            if (isPreferred) {
                thetaAccInt = 0.0;
            } 
            else {
                HolidaySP hols(Holiday::weekendsOnly());
                DateTime tomorrow = hols->addBusinessDays(getValueDate(), 1);
                double accruedToday = 0.0;
                if ( getValueDate() < callMatDate ) {
                    accruedToday      = getAccruedAtDate(getValueDate());
                } else {
                    accruedToday      = getAccruedAtDate(callMatDate);
                }
                double accruedTomorrow   = getAccruedAtDate(tomorrow);
                double couponsPV         = bond->couponsPV(getValueDate(),
                                                           tomorrow,
                                                           YieldCurveConstSP());
                thetaAccInt = accruedTomorrow - accruedToday + couponsPV;
            }

            results->storeRequestResult(request,thetaAccInt);
        }
    }

    return alreadyCalledPrice;

    // still to calculate:
}


SampleListSP ConvBond::getPutToStockSample(const DateTime putDate) const
{
    // get the holiday for the stock
    const IHaveEquity *haveEq = dynamic_cast<const IHaveEquity *>(asset.get());
    if( !haveEq )
        throw ModelException("ConvBond::getPutToStockSample", "Convert asset must support IHaveEquity (ie. EquityBase/ProtEquity/StruckEquity)");
    HolidayConstSP hols = haveEq->getEquity()->getMarketHolidays();

    int nbSampleDays = putToStock->putToStockSampleDays;
    DateTimeArray dates(nbSampleDays);
    DoubleArray   values(nbSampleDays);
    DoubleArray   weights(nbSampleDays);

    DateTime dateLoc = putDate;
    for(int i=(nbSampleDays-1); i>=0; i--)
    {
        // advance by 1 business day. use EOD
        dateLoc = DateTime(hols->addBusinessDays(dateLoc, -1).getDate(), DateTime::END_OF_DAY_TIME);
        dates[i] = dateLoc;
        values[i] = 0;
        weights[i] = 1.0/nbSampleDays;
    }
    
    // in case first date is on value date
    if( dates[0].getDate() == valueDate.getDate() )
        values[0] = asset->getSpot();

    return SampleListSP(new SampleList(dates, values, weights));
}

// calculate expected put pay off without discounting
double ConvBond::calcPutToStockFactor(const DateTime putDate, const SampleList *samples, const CVolRequest *volRequest) const
{
    static const string method = "ConvBond::calcPutToStockFactor";
    
    try {
        // just in case, during theta shift, put date becomes history
        double fwdX = (putDate>valueDate)?asset->fwdValue(putDate):asset->getSpot();
        
        // obtain BS vol for put level adjustment
        CVolProcessedSP  procVol(asset->getProcessedVol(volRequest));
        CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(procVol);
        
        // fwd, variance and covariance
        double fwdZ = samples->expectedAverage(asset.get(), valueDate);
        double varZ = 0, covarXZ = 0;
        
        if( samples->numFutureDates(valueDate, 0) > 0 )
        {
            DateTime startDate = samples->getFirstDate();
            if( startDate < valueDate ) startDate = valueDate;

            varZ = samples->averageVariance(volBS.get(), startDate, true); // usePastWeight
            
            DateTimeArray dtArray(1);
            DoubleArray dbArray(1);
            dtArray[0] = putDate;
            dbArray[0] = 1.0;
            SampleList closingSample(dtArray, dbArray, dbArray);
            covarXZ = samples->averageCovariance(&closingSample, volBS.get(), startDate);
        }

        return fwdX * exp(varZ - covarXZ) / fwdZ / putToStock->putToStockDiscount;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


double ConvBond::alreadyPutToStockValue(Control* control, Results* results) const
{
    static const string method = "ConvBond::alreadyPutToStockValue";
    
    try {
        YieldCurveConstSP riskyCurve = getRiskyCurve();
        YieldCurveConstSP optDiscZero = convertIntoIssuerShares ? discount.getSP() : riskyCurve;

        // first need to find the put put level
        const DateTime &putDate = putToStock->putToStockPutDate;
        bool putBool;
        double putLevel;
        getPutLevel(putDate, &putBool, &putLevel);
        if( !putBool )
            throw ModelException(method, "Not puttable on " + putDate.toString());

        LinearStrikeVolRequestSP volRequest(
                new LinearStrikeVolRequest(getStrike(),
                                           valueDate,
                                           bond->getMaturityDate(),
                                           false));
        double factor = calcPutToStockFactor(putDate, putToStock->putToStockSamples.get(), volRequest.get());
        double bondFlr = putLevel * riskyCurve->pv(valueDate, putDate);
        double putValue = putLevel * factor * optDiscZero->pv(valueDate, putDate);

        // output
        results->storePrice(putValue, discount->getCcy());
        recordOutputRequests(control, results, putValue, bondFlr, false);

        return putValue;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Part of IInstrumentAsAsset interface.
    Returns the yield curve used for discounting */
YieldCurveConstSP ConvBond::getDiscount() const{
    return discount.getSP();
}

/** Returns the bonds maturity date */
DateTime ConvBond::getBondMaturityDate() const{
    return bond->getMaturityDate();
}


bool ConvBond::priceDeadInstrument(CControl* control, CResults* results) const
{
    static const string method = "ConvBond::priceDeadInstrument";

    DateTime matDate = bond->getMaturityDate();
    double convertPrice;
    bool   successful = false;
    double bondValue = 0.0;

    bool    isConvertible, isCallable;
    double  convRatio, convCash;

    if ( !!putToStock && putToStock->putToStockNotified )
    {
        alreadyPutToStockValue(control, results);
        successful = true;
    } 
    else if ( alreadyCalled ) {
        alreadyCalledBlackApprox(control, results);
        successful = true;
    } else if ( valueDate > matDate ) {
        results->storePrice(0., discount->getCcy());
        recordOutputRequests(control, results, 0.0, 0.0, false);
        successful = true;
    } else if (inDefault == true) {
        double bondValue = recoveryPct*bond->getFaceValue();
        results->storePrice(bondValue, discount->getCcy());
        recordOutputRequests(control, results, bondValue, bondValue, false);
        successful = true;
    } else if (valueDate == matDate) { // need to fix for struck

        double stockPrice = asset->fwdValue(valueDate);

        getConversionInfo(valueDate, stockPrice, &isConvertible, &convRatio, &convCash);

        if (DECS == true) {
            //double face = bond->getFaceValue();

            /*
            if (stockPrice*initialConvRatio < face) {
                convertPrice = stockPrice*initialConvRatio;
            } else if (stockPrice*minConvRatio < face) {
                convertPrice = face;
            } else {
                convertPrice = stockPrice*minConvRatio;
            }
            */

            if (stockPrice < initialPrice) {
                convertPrice = stockPrice * initialConvRatio;
            } else if (stockPrice < convPrice){
                convertPrice = initialPrice * initialConvRatio;
            } else {
                convertPrice = initialPrice * initialConvRatio + (stockPrice - convPrice)*minConvRatio;
            }

            convertPrice += convCash;

        } else if (PERCS == true) {
            double face = bond->getFaceValue();

            if (stockPrice*initialConvRatio < face) {
                convertPrice = stockPrice*initialConvRatio;
            } else {
                convertPrice = face;
            }

            convertPrice += convCash;
        } else {
            double  parity, callLevel, putLevel;
            bool    isPuttable, isAdjustedForAccrued;

            getConversionInfo(valueDate, stockPrice, &isConvertible, &convRatio, &convCash);
            if (isConvertible == true) {
                parity = asset->fwdValue(valueDate) * convRatio + convCash;
            } else {
                parity = 0.;
            }

            bondValue = bond->getRedemption() + getAccruedAtDate(valueDate);

            /* get the forward vol - currently assumes log-normal vol */
            LinearStrikeVolRequestSP volRequest(
                new LinearStrikeVolRequest(getStrike(),
                                           valueDate,
                                           bond->getMaturityDate(),
                                           false));

            // interpolate the vol using our LN request
            CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));

            // calculate the instantaneous(-ish) vol
            double fwdVol;
            try {
                fwdVol = volBS->CalcVol(valueDate, DateTime(valueDate.getDate()+1,
                                                            valueDate.getDate()));
            }
            catch (exception& ) {
                fwdVol = 0.0;
            }

            bool isHardCall;
            getCallLevel(valueDate, fwdVol, &isCallable, &isAdjustedForAccrued, &isHardCall, &callLevel);
            getPutLevel(valueDate, &isPuttable, &putLevel);

            // callLevel += (getCouponAtMat)?0.0:getAccruedAtDate(valueDate);

            if ( isPuttable ) {
                /* when we put on the maturity date we always get the coupon */
                putLevel += (putAdjustForAccrued || getCouponAtMat)?0.0:getAccruedAtDate(valueDate);

                bondValue = Maths::max(bondValue, putLevel);
            }

            double issuerValue = (isCallable)?Maths::min(callLevel, bondValue):bondValue;
            convertPrice = Maths::max(Maths::max(parity,putLevel), issuerValue);
        }

        recordOutputRequests(control, results, convertPrice, bondValue, false);
        results->storePrice(convertPrice, discount->getCcy());

        successful = true;
    } else if (canPriceClosedForm()) {
        priceClosedForm(control, results, convertPrice, bondValue);
        recordOutputRequests(control, results, convertPrice, bondValue, false);

        successful = true;
    } else{
        // a little validation
        if (riskFreeCoupons == true) {
            throw ModelException(method, "riskFreeCoupons == true only implemented for closed-form DECS");
        }

        // Should see if we're above the call level here and value if we are.
        // It's a little bit tough though as the soft-call level depends on the vol.
    }

    return successful;
}


bool ConvBond::canPriceClosedForm() const
{
    static const string method = "ConvBond::canPriceClosedForm";

    try {
        // 1. it must be a DECS
        if (DECS == false) {
            return false;
        }

        if (decsHasCutoff == true) {
            return false;
        }

        // it cannot be a resettable
        if ( !!resetSchedule && resetSchedule->length() > 0 ) {
            return false;
        }

        if ( accelerateDECS == true ) {
            return false;
        }

        // 2. it must be convertible only at maturity or in the past
        const DateTimeArray& convDates = conversionRatios->getDateArray();

        if (convDates.getLength() == 0) {
            // not sure if this could happen but it would be ok
        }
        if (convDates.getLength() == 1) {
            // give three days leeway in case maturity was adjusted for bad days
            if (convDates[0].getDate() < bond->getMaturityDate().getDate()-3) {
                return false;
            }
        } else {
            // check that the interp type is none
            if (conversionRatios->getInterp() != Schedule::INTERP_NONE) {
                return false;
            }
            // make sure the penultimate date is in the past
            if (convDates[convDates.getLength()-2] >= valueDate) {
                return false;
            }
            // give three days leeway in case maturity was adjusted for bad days
            if (convDates[convDates.getLength()-1].getDate() < bond->getMaturityDate().getDate()-3) {
                return false;
            }
        }

        // 3. no calls
        if (!(!callSchedule) && callSchedule->length() > 0) {
            return false;
        }
        if (!(!softCallSchedule) && softCallSchedule->length() > 0) {
            return false;
        }

        if (!(!softCallTriggerSchedule) && softCallTriggerSchedule->length() > 0) {
            return false;
        }

        // 4. no puts
        if (!(!putSchedule) && putSchedule->length() > 0) {
            return false;
        }

        // // 5. can't be struck -- this could be allowed but would have to fix vol interp in priceClosedForm
        // //                       and in getSensitiveStrikes
        // if (ccyTreatment != "N" && ccyTreatment != "V") {
        //     return false;
        // }

        // 4. no E2C
        if (FirmAsset::TYPE->isInstance(asset.get())) {
            return false;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

    return true;
}


void ConvBond::priceClosedForm(CControl* control, CResults* results,
                                  double& price, double& bondValue) const
{
    static const string method = "ConvBond::priceClosedForm";

    try {
        int i;
        DateTime matDate = bond->getMaturityDate();

        // price live instrument
        double fwdPrice = asset->fwdValue(matDate);

        // choose how to interpolate the vol - go for traditional route for now
        LinearStrikeVolRequestSP volRequestLow(new LinearStrikeVolRequest(
                                                   initialPrice,
                                                   valueDate,
                                                   matDate,
                                                   false));
        // interpolate the vol using our LN request
        CVolProcessedBSSP volBSLow(asset->getProcessedVol(volRequestLow.get()));
        // calculate the variance
        double varianceLow = volBSLow->CalcVar(valueDate, matDate);

        // choose how to interpolate the vol - go for traditional route for now
        LinearStrikeVolRequestSP volRequestHigh(new LinearStrikeVolRequest(
                                                   convPrice,
                                                   valueDate,
                                                   matDate,
                                                   false));
        // interpolate the vol using our LN request
        CVolProcessedBSSP volBSHigh(asset->getProcessedVol(volRequestHigh.get()));
        // calculate the variance
        double varianceHigh = volBSHigh->CalcVar(valueDate, matDate);


        // create a risky curve
        YieldCurveConstSP riskyCurve = getRiskyCurve();

        // now add the coupon pv.
        YieldCurveConstSP couponYC;
        if (riskFreeCoupons == false) {
            couponYC = riskyCurve;
        } else {
            couponYC = YieldCurveConstSP::dynamicCast((IObjectConstSP)discount.getSP());
        }

        price = initialConvRatio*fwdPrice;

        price -= initialConvRatio*Black::price(true, fwdPrice,
            initialPrice, 1.0, varianceLow);

        price += minConvRatio*Black::price(true, fwdPrice,
            convPrice, 1.0, varianceHigh);

        if (RiskyCDSCurve::TYPE->isInstance(riskyCurve)) {
            // 2-state discounting required

            CashFlowArraySP bondCashFlows = bond->getCashFlows(valueDate);
            CashFlowArraySP dividendCashFlows(new CashFlowArray(0));

            // add the dividend schedule
            if ( dividendPassThrough ) {
                DividendListSP divSchedule = DividendCollector::divsBetweenDates(asset.get(),
                                                                                 valueDate,
                                                                                 valueDate,
                                                                                 matDate);

                const DividendArray& divArray = divSchedule->getArray();
                for (i=0 ; i<divArray.size() ; ++i) {
                    CashFlow cfl(divArray[i].getPayDate(),
                                 divArray[i].getDivAmount() * initialConvRatio * dividendPassThroughPct);
                    dividendCashFlows->push_back(cfl);
                }
            }

            CashFlowArraySP cashFlowStream = CashFlow::merge(bondCashFlows, dividendCashFlows);

            bondValue = 0;
            double face;
            double accrued;
            int i;
            bondValue = 0.0;
            double cashFlow;
            i         = (*cashFlowStream).size() - 1;

            // discount the risky coupons - no risk free coupons allowed yet
            // since we're looking at preferred dividends, we discount based on a recovery of 0.0
            while ( i >= 0) {
                face    = bond->getFaceValue();
                accrued = getAccruedAtDate((i==0)?valueDate:(*cashFlowStream)[i-1].date);
                cashFlow =  (*cashFlowStream)[i].amount + bondValue;
                if ( i == (*cashFlowStream).size() - 1) {
                    cashFlow -= bond->getRedemption();
                }
                bondValue = riskyCurve->riskyPV((i==0)?valueDate:(*cashFlowStream)[i-1].date,
                                            (*cashFlowStream)[i].date,
                                            cashFlow,
                                            cashFlow,
                                            useAssetRecovery,
                                            recoveryPct);

                if ( !convertIntoIssuerShares ) {
                    double amount = (*cashFlowStream)[i].amount -
                                    (( i == (*cashFlowStream).size() - 1)?bond->getRedemption():0.0);

                    price = riskyCurve->riskyPV((i==0)?valueDate:(*cashFlowStream)[i-1].date,
                                                (*cashFlowStream)[i].date,
                                                amount + price,
                                                bond->getFaceValue(),
                                                useAssetRecovery,
                                                recoveryPct);
                }
                --i;
            }

            if ( convertIntoIssuerShares ) {
                double df = discount->pv(valueDate, matDate);
                price *= df;
                price += bondValue;
            }

        } else {
            double df;
            if ( convertIntoIssuerShares ) {
                df = discount->pv(valueDate, matDate);
                price *= df;
            } else {
                price = riskyCurve->riskyPV(valueDate,
                                            matDate,
                                            price,
                                            bond->getFaceValue(),
                                            useAssetRecovery,
                                            recoveryPct);
            }

            CashFlowArraySP myCFs = bond->getCashFlows(valueDate);
            bondValue = 0;
            double pvCoupon;
            int i;
            for (i=0; i < (*myCFs).size(); i++) {
                if ((*myCFs)[i].date <= riskFreeCouponEndDate) {
                    pvCoupon = (*myCFs)[i].amount * discount->pv(valueDate, (*myCFs)[i].date);
                    bondValue += pvCoupon;
                } else {
                    pvCoupon = (*myCFs)[i].amount * couponYC->pv(valueDate, (*myCFs)[i].date);
                    bondValue += pvCoupon;
                }
            }

            // add the dividend schedule
            if ( dividendPassThrough ) {
                DividendListSP divSchedule = DividendCollector::divsBetweenDates(asset.get(),
                                                                                 valueDate,
                                                                                 valueDate,
                                                                                 matDate);

                const DividendArray& divArray = divSchedule->getArray();
                for (i=0 ; i<divArray.size() ; ++i) {
                    if (divArray[i].getPayDate() <= riskFreeCouponEndDate) {
                        pvCoupon = divArray[i].getDivAmount() * initialConvRatio * dividendPassThroughPct * discount->pv(valueDate, divArray[i].getPayDate());
                        bondValue += pvCoupon;
                    } else {
                        pvCoupon = divArray[i].getDivAmount() * initialConvRatio * dividendPassThroughPct * couponYC->pv(valueDate, divArray[i].getPayDate());
                        bondValue += pvCoupon;
                    }
                }
            }

            // subtract redemption if there is any.
            bondValue -= bond->getRedemption()*couponYC->pv(valueDate, matDate);

            price += bondValue;
        }

        results->storePrice(price, discount->getCcy());

        OutputRequest* request = NULL;
		// Duration
		if ( control->requestsOutput(OutputRequest::BOND_DURATION, request) ) {
            try {
			    results->storeRequestResult(request, bond->duration(valueDate,couponYC));
            } catch (exception& e) {
                UntweakableSP untweakableSpread(new Untweakable(e));
                results->storeRequestResult(request,untweakableSpread);
            }
		}
		// Convexity
		if ( control->requestsOutput(OutputRequest::BOND_CONVEXITY, request) ) {
            try {
			    results->storeRequestResult(request, bond->convexity(valueDate,couponYC));
            } catch (exception& e) {
                UntweakableSP untweakableSpread(new Untweakable(e));
                results->storeRequestResult(request,untweakableSpread);
            }
		}
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

    return;
}

void ConvBond::setRisky(bool flag)
{
    riskyGrowth = flag;
}

/** static function to calculate asset swap spreads */
double ConvBond::calculateAssetSwapSpread(BondSP              bond,
                                          YieldCurveConstSP   yieldCurve,
                                          DateTime            valueDate,
                                          DateTime            workoutDate,
                                          double              workoutLevel,
                                          double              bondPrice,
                                          int                 frequency,
                                          bool                stubAtFront,
					  DayCountConventionSP swapDCC)
{
    static const string routine = "ConvBond::calculateAssetSwapSpread";
    try {
        // create fixed leg of asset swap bond
        BondSP bondToPut(bond->getBondToPut(workoutDate, workoutLevel));
        // calculate risk-free pv of the coupons
        double bondPV = bondToPut->presentValue(valueDate, yieldCurve);
        // adjust bond PV for accrued interest
        bondPV -= bondToPut->getAccruedAtDate(valueDate);

        double               spread       = 0.0;
        int                  count        = 12 / frequency;
        string               interval     = "M";

        double priceTolPct = 1.e-10;
        FloatingData floaterStruct;

        floaterStruct.yieldCurve    = yieldCurve;
        floaterStruct.valueDate     = valueDate;
        floaterStruct.workoutDate   = workoutDate;
        floaterStruct.bondDiff      = bondPV - bondPrice;
        floaterStruct.count         = count;
        floaterStruct.interval      = interval;
        floaterStruct.floatDCC      = swapDCC;
        floaterStruct.bond          = bondToPut;
        floaterStruct.stubAtEnd     = !stubAtFront;
        floaterStruct.bondPrice     = bondPrice + bondToPut->getAccruedAtDate(valueDate);
        floaterStruct.faceValue     = bondToPut->getFaceValue();

        // calculate gross spread
        spread = zbrentUseful(
                 &spreadPV,                     /* (I) The function to find the root of */
                 &floaterStruct,                /* (I) Parameter block */
                0.00,                           /* (I) Lowvalue for x */ // be careful that yield/freq != -1
                4.,                             /* (I) High value for x */
                bondPV*priceTolPct);            /* (I) Tolerance */

        return spread;
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

double ConvBond::calculateZSpread(BondSP                bond,
                                  YieldCurveConstSP     yieldCurve,
                                  DateTime              valueDate,
                                  DateTime              workoutDate,
                                  double                workoutLevel,
                                  double                bondPrice,
                                  int                   frequency,
                                  bool                  stubAtFront,
				  DayCountConventionSP  swapDCC)
{
    static const string routine = "ConvBond::calculateZSpread";
    try {
        // create fixed leg of asset swap bond
        BondSP bondToPut(bond->getBondToPut(workoutDate, workoutLevel));
        // calculate risk-free pv of the coupons
        double bondPV = bondToPut->presentValue(valueDate, yieldCurve);
        // adjust bond PV for accrued interest
        bondPV -= bondToPut->getAccruedAtDate(valueDate);

        double               zSpread      = 0.0;
        int                  count        = 12 / frequency;
        string               interval     = "M";

        double priceTolPct = 1.e-10;
        FloatingData floaterStruct;

        floaterStruct.yieldCurve    = yieldCurve;
        floaterStruct.valueDate     = valueDate;
        floaterStruct.workoutDate   = workoutDate;
        floaterStruct.bondDiff      = bondPV - bondPrice;
        floaterStruct.count         = count;
        floaterStruct.interval      = interval;
        floaterStruct.floatDCC      = swapDCC;
        floaterStruct.stubAtEnd     = !stubAtFront;
        floaterStruct.bond          = bondToPut;
        floaterStruct.bondPrice     = bondPrice + bondToPut->getAccruedAtDate(valueDate);

        // calculate upper bond initial guess
        double  initialGuess = 4.0;
        double  result;
        bool    success = false;
        while (!success) {
            try {
                result = zSpreadBond(initialGuess, &floaterStruct);
                success = true;
            } catch (exception&) {
                initialGuess /= 2.0;
            }

        }

        // calculate Z-Spread
        zSpread = zbrentUseful(
                  &zSpreadBond,                  /* (I) The function to find the root of */
                  &floaterStruct,                /* (I) Parameter block */
                  0.00,                          /* (I) Lowvalue for x */ // be careful that yield/freq != -1
                  initialGuess,                  /* (I) High value for x */
                  bondPV*priceTolPct);           /* (I) Tolerance */

        return zSpread;
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


/** function to calculate asset swap spread for a convertible bond */
double ConvBond::getAssetSwapSpread(const double& bondPrice) const
{
    DateTime workoutDate;
    double   workoutLevel;
    bool     found = false;

    // determine workout date and level
    if (!!putSchedule && putSchedule->length() > 0) {
        found = putSchedule->getNextDate(valueDate.rollDate(1),
                                         bond->getMaturityDate(),
                                         workoutDate);

        if (found) {
            bool    isPuttable;
            getPutLevel( workoutDate,
                        &isPuttable,
                        &workoutLevel);

            if (!isPuttable) {
                // this should never happen
                throw ModelException("ConvBond::getAssetSwapSpread",
                                     "Internal error");
            }
        }
    }

    // if there is no more put date left, use maturity
    if (!found) {
        workoutDate     = bond->getMaturityDate();
        workoutLevel    = bond->getRedemption();
    }
    DayCountConventionSP swapDCC(DayCountConventionFactory::make("Actual/360"));
    double aswSpread = ConvBond::calculateAssetSwapSpread(
                                    bond,
                                    discount.getSP(),
                                    valueDate,
                                    workoutDate,
                                    workoutLevel,
                                    bondPrice,
                                    CompoundBasis::QUARTERLY,
                                    true /* stub at front */,
                                    swapDCC);

    return aswSpread;

}

/** function to calculate asset swap spread for a convertible bond */
double ConvBond::getZSpread(const double& bondPrice) const
{
    DateTime workoutDate;
    double   workoutLevel;
    bool     found;

    // determine workout date and level
    if (!!putSchedule && putSchedule->length() > 0) {
        found = putSchedule->getNextDate(valueDate.rollDate(1),
                                         bond->getMaturityDate(),
                                         workoutDate);

        if (found) {
            bool    isPuttable;
            getPutLevel( workoutDate,
                        &isPuttable,
                        &workoutLevel);

            if (!isPuttable) {
                // this should never happen
                throw ModelException("ConvBond::getAssetSwapSpread",
                                     "Internal error");
            }
        }
    }

    // if there is no more put date left, use maturity
    if (!found) {
        workoutDate     = bond->getMaturityDate();
        workoutLevel    = bond->getRedemption();
    }
    DayCountConventionSP swapDCC(DayCountConventionFactory::make("Actual/360"));
    double zSpread = ConvBond::calculateZSpread(
                                  bond,
                                  discount.getSP(),
                                  valueDate,
                                  workoutDate,
                                  workoutLevel,
                                  bondPrice,
                                  CompoundBasis::QUARTERLY,
                                  true /* stub at front */,
                                  swapDCC);

    return zSpread;
}


/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool ConvBond::avoidVegaMatrix(const IModel* model) {
    // AVOID for local vol, otherwise do it
    return (FD1FLV::TYPE->isInstance(model));
}

/** returns all strikes on the vol surface to which
    this instrument is sensitive */
DoubleArraySP ConvBond::getSensitiveStrikes(OutputNameConstSP outputName,
                                            const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("CVanilla::getSensitiveStrikes",
                             "VEGA_MATRIX is not valid for this instrument");
    }

    if (canPriceClosedForm() == false) {

        LinearStrikeVolRequestSP volRequest(
            new LinearStrikeVolRequest(getStrike(),
                                       valueDate,
                                       bond->getMaturityDate(),
                                       false));

        SensitiveStrikeDescriptor sensStrikeDesc;
        sensStrikeDesc.forwardOnly = false;

        asset->getSensitiveStrikes(volRequest.get(), outputName,
                                   sensStrikeDesc, sensStrikes);
    } else {
        // This would be wrong for ccy struck DECS but they presently aren't allowed to be priced closed form
        sensStrikes->resize(2);
        (*sensStrikes)[0] = initialPrice;
        (*sensStrikes)[1] = convPrice;
    }

    return sensStrikes;
}

/** returns the coupon stream of the bond from which the fixed leg of the asset swap bond floor is calculated */
CashFlowArraySP ConvBond::getAssetSwapLeg(const DateTime& swapMaturity, bool payAccruedAtMat) const
{
    CashFlowArraySP allCFs = bond->getCoupons(getValueDate());

    while ( allCFs->size() > 0 && (*allCFs)[allCFs->size()-1].date > swapMaturity)
        allCFs->pop_back();

    // remove
    if (payAccruedAtMat) {
        double accruedAtMat = getAccruedAtDate(swapMaturity);

        if ((*allCFs)[allCFs->size()-1].date == swapMaturity) {
            (*allCFs)[allCFs->size()-1].amount += accruedAtMat;
        } else {
            CashFlow matCashFlow(swapMaturity,accruedAtMat);
            allCFs->push_back(matCashFlow);
        }
    }

    return allCFs;
}

/** returns the coupon stream of the bond from which the fixed leg of the asset swap bond floor is calculated */
DateTimeArraySP ConvBond::getPaymentDates() const
{
    CashFlowArraySP cashFlows = bond->getCashFlows(valueDate);

    DateTimeArraySP paymentDates(new DateTimeArray(cashFlows->size()));

    for (int i=0; i<paymentDates->size(); ++i) {
        (*paymentDates)[i] = (*cashFlows)[i].date;
    }

    return paymentDates;

}

/** Part of IInstrumentAsAsset interface.
    Returns a date after which the instrument can no longer be used as an asset  */
DateTime ConvBond::maturityDate() const{
    return getBondMaturityDate();
}

/** Part of IInstrumentAsAsset interface.
    Returns the 'coupons' or payments that this instrument will make during
    its lifetime. This can include historic payments. */
CashFlowArraySP ConvBond::getCoupons() const{
    DateTime firstDate = getValueDate();
    DateTimeArraySP exDates(getExCouponDates());
    if (!exDates->empty()) {
        firstDate = exDates->front().rollDate(-1);
    }
    return getCoupons(firstDate);
}

/** Part of IInstrumentAsAsset interface.
    Returns the accured interest (if any) to date */
double ConvBond::getAccrued() const{
    return getAccruedAtDate(getValueDate());
}

CAssetSP ConvBond::getEquityAsset() {
    static const string method = "ConvBond::getEquityAsset";
    try {
        if (FirmAsset::TYPE->isInstance(asset.get()) ) {
            //CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
            FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(asset.get());            
            return firmAsset->getEquityAsset().getSP();
        } else {
            return asset.getSP();
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


/** Returns the name of the instrument's discount currency. */
string ConvBond::discountYieldCurveName() const {
    return discount.getName();
}




// ----------------------------------------------------------------------------------
// -- ADDINS -- ADDINS -- ADDINS -- ADDINS -- ADDINS -- ADDINS -- ADDINS -- ADDINS --
// ----------------------------------------------------------------------------------
class CVBTest: public CObject{
    typedef array<CVBTest> CVBTestArray;
public:
    static CClassConstSP const TYPE;
    DoubleArray doubArray;
    DateTimeArray dtArray;

    /** for reflection */
    CVBTest():  CObject(TYPE){}
    CVBTest(int i):  CObject(TYPE) {
        doubArray.resize(i);
        dtArray.resize(i);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CVBTest, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCVBTest);
        FIELD(dtArray, "date time array")
        FIELD(doubArray, "double array");
        clazz->setPrivate(); // hide this class
    }

    static IObject* defaultCVBTest(){
        return new CVBTest();
    }

};

typedef smartConstPtr<CVBTest> CVBTestConstSP;
typedef smartPtr<CVBTest> CVBTestSP;

CClassConstSP const CVBTest::TYPE = CClass::registerClassLoadMethod(
    "CVBTest", typeid(CVBTest), load);

class CVBTest2: public CObject{
    typedef array<CVBTest2> CVBTest2Array;
public:
    static CClassConstSP const TYPE;
    DoubleArray doubArray;
    DateTimeArray dtArray1;
    DateTimeArray dtArray2;

    /** for reflection */
    CVBTest2():  CObject(TYPE){}
    CVBTest2(int i):  CObject(TYPE) {
        doubArray.resize(i);
        dtArray1.resize(i);
        dtArray2.resize(i);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CVBTest2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCVBTest2);
        FIELD(dtArray1, "date time array")
        FIELD(doubArray, "double array");
        FIELD(dtArray2, "date time array")
        clazz->setPrivate(); // hide this class
    }

    static IObject* defaultCVBTest2(){
        return new CVBTest2();
    }

};

typedef smartConstPtr<CVBTest2> CVBTest2ConstSP;
typedef smartPtr<CVBTest2> CVBTest2SP;

CClassConstSP const CVBTest2::TYPE = CClass::registerClassLoadMethod(
    "CVBTest2", typeid(CVBTest2), load);


class ConvBondAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes two parameters - the yield curve and dates to
        get pv factor between
    */
    ConvBondSP convBond;
    DateTime startDate;
    CMarketDataSP   market;


    static IObjectSP addinGetCashFlows(ConvBondAddin* params){
        static const string routine = "ConvBondAddin::addinGetCashFlows";
        try {
            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            BondSP bond(copy(params->convBond->bond.get()));

            bond->getMarket(&model, params->market.get());

            /* roll value date forward by 0 days - this is to populate
             any samples which should be set now, but are not
             populated yet, for example when running overnight grids
             for instruments which have a SOD sample */
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(bond);

            CashFlowArraySP flows(bond->getCashFlows(params->startDate));
            CVBTestSP output(new CVBTest(flows->size()));

            for (int i = 0; i < flows->size(); i++) {
                output->dtArray[i] = (*flows)[i].date;
                output->doubArray[i] = (*flows)[i].amount;
            }

            return IObjectSP(output);
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP addinGetCoupons(ConvBondAddin* params){
        static const string routine = "ConvBondAddin::addinGetCoupons";
        try {
            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            BondSP bond(copy(params->convBond->bond.get()));

            bond->getMarket(&model, params->market.get());

            /* roll value date forward by 0 days - this is to populate
             any samples which should be set now, but are not
             populated yet, for example when running overnight grids
             for instruments which have a SOD sample */
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(bond);

            CashFlowArraySP flows(bond->getCoupons(params->startDate));
            CVBTestSP output(new CVBTest(flows->size()));

            for (int i = 0; i < flows->size(); i++) {
                output->dtArray[i] = (*flows)[i].date;
                output->doubArray[i] = (*flows)[i].amount;
            }


            return IObjectSP(output);
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP addinGetCashFlowsWithExDiv(ConvBondAddin* params){
        static const string routine = "ConvBondAddin::addinGetCashFlowsWithExDiv";
        try {

            CClosedFormLN model("VolSurface");

            // work on a copy of the instrument since we're possibly amending data
            // when doing the 0 day theta shift
            BondSP bond(copy(params->convBond->bond.get()));

            bond->getMarket(&model, params->market.get());

            /* roll value date forward by 0 days - this is to populate
             any samples which should be set now, but are not
             populated yet, for example when running overnight grids
             for instruments which have a SOD sample */
            ThetaSP thetaShift(new Theta(0, HolidaySP(Holiday::noHolidays())));
            thetaShift->applyScenario(bond);

            CashFlowArraySP flows(bond->getCashFlows());
            DateTimeArraySP exDates(bond->getExCouponDates());
            int i;
            // count the flows after the start date
            int nFlows;
            if (params->startDate.empty()) {
                nFlows = flows->size();
            } else {
                if (params->startDate > bond->getMaturityDate()) {
                    throw ModelException(routine, "startDate cannot be after maturityDate");
                }

                nFlows = 0;
                while(nFlows < flows->size()
                    && (*exDates)[flows->size() - nFlows - 1] > params->startDate) {
                    nFlows++;
                }
                if (nFlows < flows->size()) {
                    nFlows++; // show one historic
                }
            }

            // possibly an extra one for redemption
            int numCashFlows = ((*exDates)[flows->size()-1] == (*flows)[flows->size()-1].date)?nFlows:nFlows+1;

            CVBTest2SP output(new CVBTest2(numCashFlows));


            for (i=nFlows-2; i>=0; i--) {
                output->dtArray1[i]  = (*flows)[flows->size() - nFlows + i].date;
                output->doubArray[i] = (*flows)[flows->size() - nFlows + i].amount;
                output->dtArray2[i]  = (*exDates)[flows->size() - nFlows + i];
            }
            if ( numCashFlows == nFlows ) {
                output->dtArray1[nFlows - 1]  = (*flows)[flows->size()-1].date;
                output->doubArray[nFlows - 1] = (*flows)[flows->size()-1].amount;
                output->dtArray2[nFlows - 1]      = (*flows)[flows->size()-1].date;
            } else {
                output->dtArray1[nFlows - 1]  = (*flows)[flows->size()-1].date;
                output->doubArray[nFlows - 1] = (*flows)[flows->size()-1].amount - bond->getRedemption();
                output->dtArray2[nFlows - 1]  = (*exDates)[flows->size()-1];
                output->dtArray1[nFlows]      = (*flows)[flows->size()-1].date;
                output->doubArray[nFlows]     = bond->getRedemption();
                output->dtArray2[nFlows]      = (*flows)[flows->size()-1].date;
            }

            return IObjectSP(output);
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }
    /** for reflection */
    ConvBondAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ConvBondAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConvBondAddin);
        FIELD(convBond, "convertible bond object");
        FIELD(startDate, "return coupons that fall after this date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(market, "market object");
        Addin::registerClassObjectMethod("CVB_GET_CASHFLOWS",
                                         Addin::CONV_BOND,
                                         "Returns the cash flows",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)addinGetCashFlows);
        Addin::registerClassObjectMethod("CVB_GET_COUPONS",
                                         Addin::CONV_BOND,
                                         "Returns the coupons",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)addinGetCoupons);
        Addin::registerClassObjectMethod("CVB_GET_CASHFLOWS_WITH_EX_DIV",
                                         Addin::CONV_BOND,
                                         "Returns the cash flows and ex div dates",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)addinGetCashFlowsWithExDiv);
    }

    static IObject* defaultConvBondAddin(){
        return new ConvBondAddin();
    }

};

CClassConstSP const ConvBondAddin::TYPE = CClass::registerClassLoadMethod(
    "ConvBondAddin", typeid(ConvBondAddin), load);


class DECSVolAddin: public CObject{
    static CClassConstSP const TYPE;

    // input parameters
    ConvBondSP      convBond;
    CMarketDataSP   market;
    double          lowStrikeVol;
    double          highStrikeVol;

    static IObjectSP updateVol(DECSVolAddin* params) {
        static const string routine = "DECSVolAddin::updateVol";
        try {
            CClosedFormLN model("VolSurface");
            ConvBondSP convBond(copy(params->convBond.get()));
            
            // call model->getInstrumentAndModelMarket to initiate market data selection
            model.getInstrumentAndModelMarket(params->market.get(), convBond.get());

            if ( convBond->DECS == false ) {
                throw ModelException("DECSVolAddin::updateVol", "This addin function should only be used for DECS/MEDS instruments");
            }

            if ( convBond->ccyTreatment != "V" &&
                convBond->ccyTreatment != "N" ) {
                throw ModelException("DECSVolAddin::updateVol", "This addin function cannot currently be used for currency struck instruments");
            }


            // get vol surface name
            OutputNameArrayConstSP names(
                RiskProperty<VolParallel>().subjectNames(convBond));

            if (names->size() != 1) {
                throw ModelException("DECSVolAddin::updateVol", "Eror: Found multiple vol surfaces in instrument");
            }
            string surfaceName = (*names)[0]->toString();

            // get the existing vol surface
            MarketObjectSP  mktObject = params->market->GetData(surfaceName, VolSurface::TYPE);

            VolSurfaceSP    volSurf(VolSurfaceSP::dynamicCast(mktObject));
            ExpiryArrayConstSP expiries = volSurf->getExpiries();

            // set up data
            DoubleArray  strikes(2);
			strikes[0] = convBond->initialPrice;
			strikes[1] = convBond->convPrice;
            CDoubleMatrixSP vols(new DoubleMatrix(2,expiries->size()));
            for (int i=0 ; i<expiries->size() ; ++i) {
                (*vols)[0][i] = params->lowStrikeVol;
                (*vols)[1][i] = params->highStrikeVol;
            }

            volSurf = VolSurfaceSP(volSurf->update(strikes, vols, true));

            params->market->AddData(volSurf);

            return params->market;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    DECSVolAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DECSVolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDECSVolAddin);
        FIELD(convBond, "convertible bond object");
        FIELD(market, "market object");
        FIELD(lowStrikeVol,      "low strike vol");
        FIELD(highStrikeVol,     "high strike vol");
        Addin::registerClassObjectMethod("CVB_OVERRIDE_DECS_VOL",
                                         Addin::CONV_BOND,
                                         "creates a two strike vol surface and replaces the existing one in the market cache",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)updateVol);
    }

    static IObject* defaultDECSVolAddin(){
        return new DECSVolAddin();
    }

};

CClassConstSP const DECSVolAddin::TYPE = CClass::registerClassLoadMethod(
    "DECSVolAddin", typeid(DECSVolAddin), load);



//////////////////////////////////////////////////////////////
//															//
//			PutToStock              						//
//															//
//////////////////////////////////////////////////////////////

PutToStock::PutToStock() : CObject(TYPE), putToStockNotified(false), 
    putToStockSampleDays(0), putToStockDiscount(1.0)
{}

void PutToStock::validatePop2Object()
{
    static const string method = "PutToStock::validatePop2Object";

    try{
            if ( putToStockSampleDays < 1 )         
                throw ModelException(method, "Must be putToStockSampleDays >= 1");
            if ( !Maths::isPositive(putToStockDiscount) )
                throw ModelException(method, "putToStockDiscount must be > 0");

            if ( putToStockNotified )
            {
                if( !eqSampleDates || !eqSampleValues )
                    throw ModelException(method, "Empty sample dates and/or values when is putToStockNotified");
                if( eqSampleDates->size() != putToStockSampleDays || 
                    eqSampleValues->size() != putToStockSampleDays )
                    throw ModelException(method, "Sample dates and/or values size must equal putToStockSampleDays");

                // need to first create a sample list
                DoubleArray weights(putToStockSampleDays);
                for(int i=0; i<putToStockSampleDays; i++) weights[i] = 1.0/putToStockSampleDays;
                putToStockSamples = SampleListSP(new SampleList(*eqSampleDates, *eqSampleValues, weights));

                if ( !putToStockSamples || putToStockSamples->numDates(0) == 0 )
                    throw ModelException(method, "Empty samples when is putToStockNotified");
                if ( putToStockSamples->numDates(0) != putToStockSampleDays )
                    throw ModelException(method, "Num of sample elements != putToStockSampleDays");
                if ( !putToStockSamples->weightsSumToOne() )
                    throw ModelException(method, "Sample weights does not sum to 1");
                if ( putToStockPutDate.getDate() <= putToStockSamples->getLastDate().getDate() )
                    throw ModelException(method, "Put date must be after sample days for putToStock");
            }
    }
    catch (exception& e) {
       throw ModelException(e, method);
    }

    return;
}

/** Invoked when Class is 'loaded' */
void PutToStock::load(CClassSP& clazz)
{
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PutToStock, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPutToStock);
        FIELD(putToStockNotified,    "if is already put to stock");
        FIELD(putToStockSampleDays,  "put to stock sample period length. in business days");
        FIELD(putToStockDiscount,    "put to stock discount. default to 1.");
        FIELD(putToStockPutDate,     "put date if is already put to stock");
        FIELD_MAKE_OPTIONAL(putToStockPutDate);
        FIELD(eqSampleDates,                "stock price observation dates if already put notified");
        FIELD_MAKE_OPTIONAL(eqSampleDates);
        FIELD(eqSampleValues,               "stock prices if already put notified");
        FIELD_MAKE_OPTIONAL(eqSampleValues);
        // transient
        FIELD(putToStockSamples,            "");
        FIELD_MAKE_TRANSIENT(putToStockSamples);
}

IObject* PutToStock::defaultPutToStock()
{
        return new PutToStock();
}

CClassConstSP const PutToStock::TYPE = CClass::registerClassLoadMethod(
    "PutToStock", typeid(PutToStock), load);



//////////////////////////////////////////////////////////////
//															//
//			Addon Conversion Schedule						//
//															//
//////////////////////////////////////////////////////////////


AddOnConvRatio::AddOnConvRatio() : CObject(TYPE), interfaceEZ(false)
{}

void AddOnConvRatio::validatePop2Object()
{
    static const string method = "AddOnConvRatio::validatePop2Object";

    if( nbSchedule() == 0 )
            throw ModelException(method, "Min allowed nb of trigger/schedule is 1");
    if( nbSchedule() > 3 )
            throw ModelException(method, "Max allowed nb of trigger/schedule is 3");

    // if simple interface, shares input must match triggers
    if( interfaceEZ && (!shares || shares->size() != nbSchedule() ))
        throw ModelException(method, "Nb of shares must equal nb of triggers");

    // check that trigger is strictly increasing and positive
    for(int i=0; i<nbSchedule(); i++)
    {
        if( (*triggers)[i] <= 0 )
            throw ModelException(method, "Trigger level must be positive");
        if( i>0 && (*triggers)[i] <= (*triggers)[i-1] )
            throw ModelException(method, "Trigger level must be strictly increasing");

        ScheduleSP &sch = (i==0?schedule1:(i==1?schedule2:schedule3));
        // if simple interface, need to create schedules
        if( interfaceEZ )
        {
            DateTimeArray dates(2);
            DoubleArray values(2);
            dates[0] = startDate;
            dates[1] = endDate;
            values[0] = (*shares)[i];
            values[1] = (*shares)[i];
            sch = ScheduleSP(new Schedule(dates, values, Schedule::INTERP_LINEAR));
        }
        if( !sch )
            throw ModelException(method, "No schedule for trigger " + Format::toString(i+1));
    }
}

int AddOnConvRatio::nbSchedule() const
{
    return (!triggers)?0:triggers->size();
}

double AddOnConvRatio::getNextTrigger(double spot, double spotMax)
{
    if( spot > spotMax )
        throw ModelException("AddOnConvRatio::getNextTrigger", "Internal error. spotMax < spot");

    for(int i=0; i<nbSchedule(); i++)
        if( spot < (*triggers)[i] ) return (*triggers)[i];

    return spotMax;
}

void AddOnConvRatio::adjustConvRatio(const DateTime& convDate,
    bool isSpot, double spotOrLevel, double *ratio, double *cash) const
{
    double dRatio;
    DateTime dummyDate;

    // add incremental conv ratios. stop if either not triggered or not convertible
    int i=0;
    while(  i<nbSchedule() &&
            spotOrLevel > (isSpot?(*triggers)[i]:((*triggers)[i]*(*ratio)+(*cash))) )
	{
		const Schedule *sch = (i==0?schedule1:(i==1?schedule2:schedule3)).get();
		// stop if date outside range, or is not convertible
        if ( convDate < sch->firstDate() ||
			!sch->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                    0.,
                                                    dummyDate,
                                                    convDate,
                                                    dRatio) )
			break;

        *ratio += dRatio;
        *cash -= dRatio * (*triggers)[i];
        i++;
	}
	if( *ratio < 0 )
		throw ModelException("AddOnConvRatio::adjustConvRatio", "AddOnConvRatio input error. Conversion ratio negative at " + convDate.toString() );
}

void AddOnConvRatio::insertCritDates(DateTimeArray &critDates) const
{
    for(int i=0; i<nbSchedule(); i++)
    {
		const Schedule *sch = (i==0?schedule1:(i==1?schedule2:schedule3)).get();
        const DateTimeArray& conversionDates = sch->getDateArray();
        for (int j=0; j<conversionDates.size(); j++)
            critDates.push_back(conversionDates[j]);
    }
}

class AddOnConvRatioHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AddOnConvRatio, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAddOnConvRatio);
        FIELD(triggers,     "triggers for the schedules");
        FIELD(interfaceEZ,  "0 if use schedule1/2/3 as input (default), 1 if use simple interface");
        FIELD_MAKE_OPTIONAL(interfaceEZ);
        FIELD(startDate,    "start date to possibly add shares, if simple interface");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(endDate,      "end date to possibly add shares, if simple interface");
        FIELD_MAKE_OPTIONAL(endDate);
        FIELD(shares,       "additional shares corresponding to triggers, if simple interface");
        FIELD_MAKE_OPTIONAL(shares);
        FIELD(schedule1,    "1st schedule of additional conv ratios");
        FIELD_MAKE_OPTIONAL(schedule1);
        FIELD(schedule2,    "2nd schedule of additional conv ratios");
        FIELD_MAKE_OPTIONAL(schedule2);
        FIELD(schedule3,    "3rd schedule of additional conv ratios");
        FIELD_MAKE_OPTIONAL(schedule3);
    }

    static IObject* defaultAddOnConvRatio(){
        return new AddOnConvRatio();
    }
};

CClassConstSP const AddOnConvRatio::TYPE = CClass::registerClassLoadMethod(
    "AddOnConvRatio", typeid(AddOnConvRatio), AddOnConvRatioHelper::load);



//////////////////////////////////////////////////////////////////////
//                                                                  //
//          Contingent Conversion : implementation                  //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

bool ConvBond::hasContConversion() const
{
    return (!!contConversion && !contConversion->isTrivial() );
}

ContConversion::ContConversion() : CObject(TYPE), badDayConvString("F"),
frequency(4), numPeriod(0), periodDaysM(20), periodDaysN(30),
periodTrigger(-1), periodTrigAdjType("E"), postPeriodTrigger(-1), postPeriodEnabled(false), periodTrigActiveAtMat(false),
floorDaysM(20), floorDaysN(5), floorTrigger(-1), upEqBarrier4Flr(-1), dnEqBarrier4Flr(-1),
curPeriodEnabled(false), floorEnabled(false), nHistPeriod(0), ignoreFlrHist(true) 
{}

void ContConversion::validatePop2Object()
{
    static const string method = "ContConversion::validatePop2Object";

    try {

    badDayConv = BadDayConventionSP(BadDayConventionFactory::make(badDayConvString));

    if( firstPeriodStart.getTime() != DateTime::START_OF_DAY_TIME )
        throw ModelException(method, "First period start date must be start of day");

    if( hasPeriod() )
    {
        if( periodDaysM == 0 || periodDaysN == 0 )
            throw ModelException(method, "periodDaysM and periodDaysN can not be zero");
        if( periodDaysM > periodDaysN )
            throw ModelException(method, "periodDaysM must be <= periodDaysN");
        if( frequency != 1 && frequency != 2 && frequency != 4 && frequency != 12 )
            throw ModelException(method, "frequency must be 1, 2, 4 or 12");
        if( periodTrigAdjType != "N" && periodTrigAdjType != "E" )
            throw ModelException(method, "periodTrigAdjType type " + periodTrigAdjType + " invalid. Must be N or E");
    }

    if( hasFloor() )
    {
        if( !ignoreFlrHist )
            throw ModelException(method, "ignoreFlrHist = FALSE not supported");
        if( Maths::isPositive(floorTrigger - 1.0) )
            throw ModelException(method, "floorTrigger must be <= 1.");
        if( floorDaysM == 0 || floorDaysN == 0 )
            throw ModelException(method, "floorDaysM and floorDaysN can not be zero");
        if( upEqBarrier4Flr < dnEqBarrier4Flr )
            throw ModelException(method, "Floor up barrier < down barrier");
    }

    // allocate memory for samples to avoid always have to validate NULLness
    if( !eqSampleDates ) eqSampleDates = DateTimeArraySP(new DateTimeArray());
    if( !cvSampleDates ) cvSampleDates = DateTimeArraySP(new DateTimeArray());
    if( !eqSampleValues ) eqSampleValues = DoubleArraySP(new DoubleArray());
    if( !cvSampleValues ) cvSampleValues = DoubleArraySP(new DoubleArray());

    // some validation require market info and is left for getMarket()

    }
    catch (exception& e) {
       throw ModelException(e, "ConvBond::validatePop2Object");
    }
}


/** pull out holiday from the market data. also some one time processing */
void ContConversion::getMarket(const IModel* model, const MarketData* market, const ConvBond *cvb)
{
    static const string method = "ContConversion::getMarket";
    try{

        // get the holiday for the stock
        const IHaveEquity *haveEq = dynamic_cast<const IHaveEquity *>(cvb->asset.get());
        if( !haveEq )
            throw ModelException(method, "Convert asset must support IHaveEquity (ie. EquityBase/ProtEquity/StruckEquity)");
        hols = haveEq->getEquity()->getMarketHolidays();

        // set up end dates of the sample dates for each period. do trigger level setup
        // observation period ends at EOD
        DateTime periodStart = firstPeriodStart;
        periodStarts = DateTimeArraySP(new DateTimeArray(numPeriod+1));
        periodTriggers = DoubleArraySP(new DoubleArray(numPeriod+1));
        for(int i=0; i<=numPeriod; i++)
        {
            (*periodStarts)[i] = badDayConv->adjust(periodStart, hols.get());
            (*periodTriggers)[i] = (i<numPeriod?periodTrigger:postPeriodTrigger) * getConversionPrice(periodStart, cvb);
            periodStart = MaturityPeriod::toDate(12/frequency, "M", periodStart);       
        }

        // checking historical equity samples
        valueDate = cvb->getValueDate();
        validateSample(eqSampleDates.get(), eqSampleValues.get(), cvb->asset->getSpot());
        validateSample(cvSampleDates.get(), cvSampleValues.get());
 
        // postPeriod trigger
        if( valueDate < (*periodStarts)[numPeriod] && postPeriodEnabled )
            throw ModelException(method, "PostPeriod convertibility is true but postPeriod starts after valueDate");

        bondMatDate = cvb->getBondMaturityDate();

        // further historical related processing
        preprocess(cvb);
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void ContConversion::validateSample(const DateTimeArray *dates, DoubleArray *samples, double spot)
{
    string method = "ContConversion::validateSample(" + Format::toString(Maths::isZero(spot)?"cvb":"eq") + "Sample)";

    int dateSz = (!dates)?0:dates->size();
    int sampleSz = (!samples)?0:samples->size();

    if( dateSz != sampleSz )
        throw ModelException(method, "Dimension of date and sample disagree");

    if( dateSz == 0 ) return; // we know both samples are empty

    for(int i=0; i<dates->size(); i++)
    {
        // make sure it's ordered, not on the same date, but of same time of the day
        if( i>0 )
        {
            if( (*dates)[i].getDate() <= (*dates)[i-1].getDate() )
                throw ModelException(method, "Sample dates out of order");
            if( (*dates)[i].getTime() != (*dates)[i-1].getTime() )
                throw ModelException(method, "Sample dates must be of the same time of day");
        }
        
        if( !Maths::isPositive((*samples)[i]) && (*dates)[i] <= valueDate )
        {
            if( Maths::isZero((*samples)[i]) && (*dates)[i] == valueDate && !Maths::isZero(spot) )
                (*samples)[i] = spot;
            else
                throw ModelException(method, "Sample is not positive at " + (*dates)[i].toString());
        }       
    }
}

// one time set up of the sampling dates and other important dates
void ContConversion::preprocess(const ConvBond *cvb)
{
    static const string method = "ContConversion::preprocess";

    try {
        int i;

        // find out # of historical periods, including the period that has just started
        nHistPeriod = 0;
        for(i=0; i<=numPeriod; i++)
        {
            if( valueDate >= (*periodStarts)[i] ) nHistPeriod++;
        }

        // for 1st partial historical period which has full historical sampling, find out if is not eligible to convert in this period
        curPeriodEnabled = false;
        if( nHistPeriod > 0 && nHistPeriod <= numPeriod )
        {
            i = nHistPeriod-1;
            int m = getHistTrigCount(eqSampleDates.get(), eqSampleValues.get(), (*periodStarts)[i], periodDaysN, (*periodTriggers)[i]);
            if( m >= periodDaysM ) curPeriodEnabled = true;
        }

        // postperiod enable
        if( valueDate >= (*periodStarts)[numPeriod] && !postPeriodEnabled )
        {
            i=0; 
            while( i<eqSampleDates->size() && (*eqSampleDates)[i] < (*periodStarts)[numPeriod] )
                i++;
            while( i<eqSampleDates->size() && (*eqSampleDates)[i] <= valueDate && !postPeriodEnabled )
            {
                // sanity check
                if( i>0 && hols->businessDaysDiff((*eqSampleDates)[i-1], (*eqSampleDates)[i]) != 1 )
                    throw ModelException(method, "eqSample must be daily for postPeriod");

                double trigger = postPeriodTrigger * getConversionPrice((*eqSampleDates)[i], cvb);
                postPeriodEnabled = !Maths::isNegative((*eqSampleValues)[i] - trigger);
                i++;
            }
        }

        // floor
        floorEnabled = false;
        bool isConvertible;
        double convRatio, convCash;
        cvb->getConversionInfo(valueDate, &isConvertible, &convRatio, &convCash);
        if( !ignoreFlrHist && hasFloor() && isConvertible && valueDate >= firstPeriodStart )
        {
            DateTime end = valueDate;
            double spot = getHistSum(eqSampleDates.get(), eqSampleValues.get(), end, 1);
            double convPrice = getConversionPrice(end, cvb);
            if( spot <= upEqBarrier4Flr*convPrice && spot >= dnEqBarrier4Flr*convPrice )
            {
                for(i=0; i<floorDaysN; i++)
                {
                    double eqSum = getHistSum(eqSampleDates.get(), eqSampleValues.get(), end, floorDaysM);
                    double cvSum = getHistSum(cvSampleDates.get(), cvSampleValues.get(), end, floorDaysM);
                    if( cvSum < floorTrigger * (eqSum * convRatio + convCash * floorDaysM) ) 
                    {
                        floorEnabled = true;
                        break;
                    }
                    end = hols->addBusinessDays(end, -1);
                }
            }
        } // end of floor
        
        
    } catch (exception &e) {
        throw ModelException(e);
    }
}


void ContConversion::setVol(CVolProcessedBSConstSP vols)
{
    this->vols = vols;
}

// insert period dates and floor start/end dates
void ContConversion::insertCritDates(DateTimeArray &critDates) const
{
    for(int i=nHistPeriod; i<=numPeriod; i++)
    {
        critDates.push_back((*periodStarts)[i]);
    }
}

bool ContConversion::hasPeriod() const
{
    return ( numPeriod && Maths::isPositive(periodTrigger) );
}

bool ContConversion::hasPostPeriod() const
{
    return Maths::isPositive(postPeriodTrigger);
}

bool ContConversion::hasFloor() const
{
    return Maths::isPositive(floorTrigger) && Maths::isPositive(upEqBarrier4Flr);
}

bool ContConversion::isTrivial() const
{
    return (( !hasPeriod() && !hasPostPeriod() && !hasFloor() ) || isHistEnabled(bondMatDate) );
}

bool ContConversion::isHistEnabled(const DateTime &date) const
{
    return  date < firstPeriodStart ||
            ( curPeriodEnabled && date >= (*periodStarts)[nHistPeriod-1] &&
              date < (*periodStarts)[Maths::min(numPeriod, nHistPeriod)] ) ||
            ( date >= (*periodStarts)[numPeriod] && postPeriodEnabled ) ||
            ( date >= firstPeriodStart && date.getDate() == valueDate.getDate() && floorEnabled );
}

bool ContConversion::rollDate(Theta* shift, const ConvBond *cvb)
{
    static const string method = "ContConversion::sensShift";
    try {

        // calculate the new date
        DateTime newDate = shift->rollDate(valueDate);
        
        // shift historical equity samples. 
        // need to insert points if not already in there!
        vector<DateTime>::iterator iter1(eqSampleDates->begin());
        vector<double>::iterator iter2(eqSampleValues->begin());
        DateTime tmpDate = valueDate;
        int obsTime = eqSampleDates->size()?(*eqSampleDates)[0].getTime():DateTime::END_OF_DAY_TIME;
        while( tmpDate <= newDate && !postPeriodEnabled )
        {
            while(iter1 != eqSampleDates->end() && (*iter1).getDate() < tmpDate.getDate() )
            {
                ++iter1; ++iter2;
            }

            if( iter1 == eqSampleDates->end() || tmpDate.getDate() < (*iter1).getDate() )
            {
				iter1 = eqSampleDates->insert(iter1, DateTime(tmpDate.getDate(), obsTime)); ++iter1;
				iter2 = eqSampleValues->insert(iter2, 0); ++iter2;
            }
	
            tmpDate = hols->addBusinessDays(tmpDate, 1); 
        }
    
        const CAsset *asset = cvb->asset.get();
        const DateTimeArray &dates = *eqSampleDates;
        DoubleArray &values = *eqSampleValues;
        bool useSpot = shift->useAssetFwds()? false : true;
        for(int i=0; i<dates.size(); i++)
        {
            if ( ( dates[i].isGreater(valueDate) && !dates[i].isGreater(newDate)) ||
                ( dates[i].equals(valueDate)    && Maths::isZero(values[i]))    ) 
            {
                values[i] = useSpot?asset->getSpot():asset->fwdValue(dates[i]);
            } // end of if
        } 

        // should also shift cv sample, but cv sample is not enabled for now. so no work

        // roll today 
        valueDate = newDate;

        // further historical related processing using the new valuedate
        // make sure it's not the first 
        preprocess(cvb);

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    

    return false; // no other component needs to be tweaked
}

// get conversion ratio and devide by 
double ContConversion::getConversionPrice(const DateTime &date, const ConvBond *cvb) const
{
    bool isConvertible;
    double convRatio, convCash;
    cvb->getConversionInfo(date, &isConvertible, &convRatio, &convCash);

    // the trade is not well specified if it's not convertible when this function is called
    if( !isConvertible || Maths::isZero(convRatio) )
        throw ModelException("ContConversion::getConversionPrice", "Contingent conversion require "
        "conversion for date " + date.toString() + ", but is not convertible according to conv schedule");
    
    double par = cvb->bond->getFaceValue();
    return (par - convCash)/convRatio;
}

// we do m out of n vol adjustment to trigger. borrow code from softCall adj in convert
double ContConversion::getAdjTrigger(int m, int n,  // m out of past n days
                                double trigger, // trigger level in $ terms
                                double fwdRate, // annual fwd rate
                                double vol)     // % vol for n days
{
    if( m>n || n==0 )
        throw ModelException("ContConversion::getTriggerAdj", "Must have n>m and n>0");

    if( n == 30 && m == 20 )
        trigger *= exp(vol * 0.2316);
    else if (m > 1)
        trigger *= exp(vol * (0.050959 * (double) (m-1) / sqrt((double) n) + 0.0008774 * (double) n + 0.03055));

    return trigger;
}

// obtain quotes with datetime.date() >= start and datetime.date() < end. make sure datetime <= valueDate
DoubleArray ContConversion::getHistQuotes(const DateTimeArray *dates, const DoubleArray *samples, DateTime end, int nDays) const
{
    if( !dates || !samples )
        throw ModelException("ContConversion::getHistQuotes", "No date or sample provided");

    DoubleArray quotes;
    DateTime start = hols->addBusinessDays(end, -nDays);
    if( end > valueDate ) end = valueDate;
    nDays = hols->businessDaysDiff(start, end);
    if( nDays <=0 ) return quotes; // return empty array
    quotes.resize(nDays);

    int i=0, j=0;
    while( i<dates->size() && (*dates)[i].getDate() < start.getDate() ) i++;
    while( i<dates->size() && (*dates)[i].getDate() < end.getDate() )
    {
        if( j < nDays ) quotes[j] = (*samples)[i];
        i++; j++;
    }

    // sanity check, make sure there is no extra days or missing days
    if( nDays != j )
        throw ModelException("ContConversion::getHistQuotes", "There is missing or extra days in the sample between " + start.toString() + " and " + end.toString());

    return quotes;
}

// obtain historical trigger count with sample > trigger
int ContConversion::getHistTrigCount(const DateTimeArray *dates, const DoubleArray *samples, const DateTime &end, int nDays, double trigger, int *nHistDays) const
{
    DoubleArray quotes = getHistQuotes(dates, samples, end, nDays);

    int count=0;
    for(int i=0; i<quotes.size(); i++)
    {
        if( quotes[i] >= trigger) count++;
    }
    if( nHistDays ) *nHistDays = quotes.size();
    return count;
}

double ContConversion::getHistSum(const DateTimeArray *dates, const DoubleArray *samples, const DateTime &end, int nDays, int *nHistDays) const
{
    DoubleArray quotes = getHistQuotes(dates, samples, end, nDays);

    double sum = 0;
    for(int i=0; i<quotes.size(); i++)
        sum += quotes[i];

    if( nHistDays ) *nHistDays = quotes.size();
    return sum;
}

// adjust price depending on probability of being convertible or not
// don't call this routine if already postPeriodEnabled
void ContConversion::adjustPrices(
        const DateTime &date,
        int nStockStep, 
        const double *s,            // stock grid
        double *priceNCV,           // price if not convertible
        double *priceCV,            // price if convertible
        double callLevel,
        const ConvBond *cvb) const
{
    static const string method = "ContConversion::adjustPrices";
    try {

        if( postPeriodEnabled )
            throw ModelException(method, "Internal error. Can not call if already postPeriodEnabled");

        bool isConvertible;
        double convRatio, convCash;

        // this is to cure FD evolution problem when useFwdGrid=true. Due to spline interp,
        // if there is a kink in priceCV/NCV, even if priceNCV < priceCV, but interpolated value
        // may cause priceNCV to be slightly above priceCV around kink region.
        for(int n=0; n<=nStockStep; n++)
        {
            if( priceNCV[n] > priceCV[n] )
            {
                priceNCV[n] = 0.5*(priceCV[n] + priceNCV[n]);
                priceCV[n] = priceNCV[n];
            }
        }

        // if inside postPeriod, see if knock in conversion
        if( date > (*periodStarts)[numPeriod] )
        {
            if( hasPostPeriod() )
            {
                double kiTrigger = postPeriodTrigger * getConversionPrice(date, cvb);
                for(int n=nStockStep; n>=0; n--)
                {
                    if( s[n] < kiTrigger ) break;
                    priceNCV[n] = priceCV[n];
                }
            }
        }
        else if( date >= (*periodStarts)[0] ) // include postPeriodStart since may need intrinsic calculation for priceCV
        {
            // if period start, merge priceNCV and priceCV
            for(int i=nHistPeriod; i<=numPeriod; i++)
            {
                if( (*periodStarts)[i] == date )
                {
                    // here need to get the conv ratio and cash from the prev business day end of day! 
                    // it's used to get intrinsic value for priceCV
                    DateTime prevEOD = DateTime(hols->addBusinessDays(date, -1).getDate(), DateTime::END_OF_DAY_TIME);
                    cvb->getConversionInfo(prevEOD, &isConvertible, &convRatio, &convCash);

                    double trigger = (*periodTriggers)[i];

                    // for pure future period, do adjustment depending on if  
                    // the trigger sampling is partially historical or purely future
                    // first future period may be partial historical sample
                    int m = 0, n = 0;
                    if( i == nHistPeriod )
                        m = getHistTrigCount(eqSampleDates.get(), eqSampleValues.get(), (*periodStarts)[i], periodDaysN, (*periodTriggers)[i], &n);
                    m = periodDaysM - m;
                    n = periodDaysN - n;

                    if ( m<= 0 ) // definitely convertible
                    {
                        for(int j=0; j<=nStockStep; j++)
                            priceNCV[j] = priceCV[j];
                    } 
                    else if ( m>n ) // definitely not convertible
                    {
                        for(int j=0; j<=nStockStep; j++)
                        {
                            priceCV[j] = Maths::max(priceNCV[j], convRatio * s[j] + convCash);
                            if( !Maths::isNegative(callLevel) )
                                priceCV[j] = Maths::min(priceCV[j], callLevel);
                        }
                    } 
                    else // trigger is needed. may require adjustment if not the last period
                    {
                        if( periodTrigAdjType != "N" && i!=numPeriod)
                        {
                            Actual365F fwdRateDCC;
                            double fwdRate = cvb->getDiscount()->fwd(date, date.rollDate(1), &fwdRateDCC, CompoundBasis::ANNUAL);

                            DateTime start = hols->addBusinessDays((*periodStarts)[i], -periodDaysN);
                            double vol = vols->CalcVar(start>valueDate?start:valueDate, (*periodStarts)[i]);
                            vol = sqrt(vol);

                            trigger = getAdjTrigger(m, n, trigger, fwdRate, vol);

                        } // end of if(periodTrigAdjType)

                        for(int j=0; j<=nStockStep; j++)
                        {
                            double price = (s[j]>=trigger)?priceCV[j]:priceNCV[j];
                            priceNCV[j] = price;
                            priceCV[j] = Maths::max(price, convRatio * s[j] + convCash);
                            if( !Maths::isNegative(callLevel) )
                                priceCV[j] = Maths::min(priceCV[j], callLevel);
                        }
                    }

                    break;
                } // end of if(date is period start)
            } // end of date loop
        }

        // simple implementation for floor. set priceNCV to be at least trigger level of parity
        cvb->getConversionInfo(date, &isConvertible, &convRatio, &convCash);
        if( hasFloor() && isConvertible && date >= firstPeriodStart )
        {
            // need to adjust the trigger by historical sampling
            double offset = 0;
            if( !ignoreFlrHist && valueDate > hols->addBusinessDays(date, -floorDaysM) )
            {
                int n;
                double eqSum = getHistSum(eqSampleDates.get(), eqSampleValues.get(), date, floorDaysM, &n);
                double cvSum = getHistSum(cvSampleDates.get(), cvSampleValues.get(), date, floorDaysM);
                offset = (floorTrigger * (convRatio * eqSum + convCash * n) - cvSum) / (floorDaysM - n);
            }

            double convPrice = getConversionPrice(date, cvb);
            double upBarr = upEqBarrier4Flr * convPrice;
            double dnBarr = dnEqBarrier4Flr * convPrice;
            for(int n=0; n<=nStockStep; n++)
            {
                if( s[n] > upBarr ) break;
                if( s[n] >= dnBarr )
                {
                    double floor = floorTrigger * (convRatio * s[n] + convCash) + offset;
                    if( priceNCV[n] < floor ) priceNCV[n] = floor;
                    if( !Maths::isNegative(callLevel) )
                        priceNCV[n] = Maths::min(priceNCV[n], callLevel);
                }
            }
        }


    } catch (exception &e) {
        throw ModelException(e, method);
    }
}

class ContConversionHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ContConversion, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultContConversion);
        FIELD(badDayConvString,  "bad day addjustment N, P, F, or M")
        FIELD_MAKE_OPTIONAL(badDayConvString);
        FIELD(firstPeriodStart,  "Start date of 1st exchange period");
        FIELD(frequency,         "Exchange period frequency per year");
        FIELD(numPeriod,         "Num of exchange period");
        FIELD(periodDaysM,       "Num of days stock price needs to reach trigger to be convertible next period");
        FIELD(periodDaysN,       "Num of days prior to exchange period to sample stock price");
        FIELD(periodTrigger,     "Stock price need to reach trigger * conversion price. Use (-1) to turn off");
        FIELD(periodTrigAdjType, "Trigger adjustment type. N(one) or E(mpirical)");
        FIELD_MAKE_OPTIONAL(periodTrigAdjType);
        FIELD(postPeriodTrigger, "Stock price need to reach trigger * conversion price to knock in conversion thereafter, after end of the periods. Use (-1) to turn off");
        FIELD(postPeriodEnabled, "Stock price already reached trigger during postPeriod");
        FIELD_MAKE_OPTIONAL(postPeriodEnabled);
        FIELD(periodTrigActiveAtMat, "If trigger condition apply for conversion at maturity. Default is false");
        FIELD_MAKE_OPTIONAL(periodTrigActiveAtMat);
        FIELD(floorDaysM,        "Num of days to average convert price and stock price (parity)");
        FIELD_MAKE_OPTIONAL(floorDaysM);
        FIELD(floorDaysN,        "Length of conversion window after floor condition is met");
        FIELD_MAKE_OPTIONAL(floorDaysN);
        FIELD(floorTrigger,      "Avg convert price needs to be below trigger * avg stock price (parity) to be convertible the next day. Use (-1) to turn off");
        FIELD(upEqBarrier4Flr,   "Stock price must be below upper barrier * conversion price for floor privision to be considerred");
        FIELD(dnEqBarrier4Flr,   "Stock price must be below lower barrier * conversion price for floor privision to be considerred");
        FIELD(eqSampleDates,            "Historical stock price monitor dates. For periods and for floor");
        FIELD_MAKE_OPTIONAL(eqSampleDates);
        FIELD(eqSampleValues,           "Historical stock prices. For periods and for floor");
        FIELD_MAKE_OPTIONAL(eqSampleValues);
        FIELD(cvSampleDates,            "Historical convert price minotor dates. For floor");
        FIELD_MAKE_OPTIONAL(cvSampleDates);
        FIELD(cvSampleValues,           "Historical convert prices. For floor");
        FIELD_MAKE_OPTIONAL(cvSampleValues);
        // internal fields
        FIELD(valueDate,         "");
        FIELD_MAKE_TRANSIENT(valueDate);
        FIELD(bondMatDate,       "");
        FIELD_MAKE_TRANSIENT(bondMatDate);
        FIELD(curPeriodEnabled,  "");
        FIELD_MAKE_TRANSIENT(curPeriodEnabled);
        FIELD(floorEnabled,      "");
        FIELD_MAKE_TRANSIENT(floorEnabled);
        FIELD(nHistPeriod,       "");
        FIELD_MAKE_TRANSIENT(nHistPeriod);
        FIELD(periodStarts,             "");
        FIELD_MAKE_TRANSIENT(periodStarts);
        FIELD(periodTriggers,           "");
        FIELD_MAKE_TRANSIENT(periodTriggers);
        FIELD(hols,                     "");
        FIELD_MAKE_TRANSIENT(hols);
        FIELD(badDayConv,               "");
        FIELD_MAKE_TRANSIENT(badDayConv);
    }

    static IObject* defaultContConversion(){
        return new ContConversion();
    }
};

CClassConstSP const ContConversion::TYPE = CClass::registerClassLoadMethod(
    "ContConversion", typeid(ContConversion), ContConversionHelper::load);


DRLIB_END_NAMESPACE
