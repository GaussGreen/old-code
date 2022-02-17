//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Tracer.cpp
//
//   Description : Trigger Activated Convertible Security
//
//   Author      : André Segger
//
//   Date        : September 27, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ConvBond.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Tracer.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/DblBarrier.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

// copy market data relevant to the instrument
void Tracer::GetMarket(const IModel* model, 
                       const CMarketDataSP market){
    static const string method = "Tracer::GetMarket";
    try {
        this->market = CMarketDataSP::attachToRef(market.get());
        if (TracerModel::TYPE->isInstance(model)) {
            const TracerModel* umbrella=dynamic_cast<const TracerModel*>(model);

            market->GetReferenceDate(valueDate);
            discount.getData(umbrella->cvbModel.get(), market);
            creditSpreads.getData(umbrella->cvbModel.get(), market);
            frontBond->getMarket(umbrella->cvbModel.get(), market.get());
            
            CAsset::getAssetMarketData(umbrella->cvbModel.get(), market.get(), ccyTreatment, 
                                       discount.getName(), asset);
            
            CAsset::getAssetMarketData(umbrella->knockInModel.get(), market.get(), ccyTreatment, 
                                       discount.getName(), asset);
        }
        else {
            throw ModelException(method, "Model to price tracer must be a TracerModel");
        }

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void Tracer::validatePop2Object()
{
    try{
        // to do
    }
    catch (exception& e) {
       throw ModelException(e, "Tracer::validatePop2Object");
    }
            
    return;
}

void Tracer::Validate() {
    static const string method = "Tracer::Validate";

    try {
        // to do
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// -- Stuff every product needs to do -- //

/** what's today ? */
DateTime Tracer::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime Tracer::endDate(const Sensitivity* sensControl) const {

    // front bond parameters
    DateTime frontEndMaturity = frontBond->getMaturityDate();
    return bondMaturity->toDate(frontEndMaturity);
}

bool Tracer::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "Tracer::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}

// for reflection
Tracer::Tracer(): CInstrument(TYPE), trigger(0.0), conversionPremium(0.0), callTriggerLevel(0.0), coupon(0.0),
                  redemption(0.0), putLevel(0.0), bondFrequency(2), downsideProtection(0.0), decsPremium(0.0), 
                  decsCoupon(0.0), decsBondFrequency(2), decsDivPassThrough(false), backBondCapped(false), backBondCap(0.0),
                  knockedIn(false) {
    //cvbModel = IModelSP(new ClosedForm());
    // no code
}

class TracerHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Tracer, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(TracerModel::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultTracer);
        FIELD(discount,               "identifies discount curve");
        FIELD(creditSpreads,          "credit spread to discount curve");
        FIELD(asset,                  "the equity");
        FIELD(valueDate,              "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(ccyTreatment,           "N for normal, S for struck");
        FIELD_MAKE_OPTIONAL(ccyTreatment);
        FIELD(frontBond,                     "front bond object");
        FIELD(bondMaturity,                  "Maturity of front bond");
        FIELD(callStartPeriod,               "start of call schedule");
        FIELD(callTriggerLevel,       "call trigger level");
        FIELD(coupon ,                "back bond coupon" );
        FIELD(redemption,             "back bond redemption");
        FIELD(putLevel,               "put Level");
        FIELD(putStartPeriod,                "start of put schedule");
        FIELD(cbDates,                "Convertible valuation dates for rebate interpolation" );
        FIELD(trigger,                "Knock-in trigger" );
        FIELD(conversionPremium,      "Conversion Premium" );
        FIELD(dayCountConvString,     "back bond day count convention");
        FIELD(bondFrequency,          "back bond coupon frequency");
        FIELD(holidays,               "back bond holidays");
        FIELD(adjustRebateForAccrued, "adjust rebate for accrued interest");
        FIELD(backBondCapped,         "true if the back bond valuation is floored");
        FIELD(backBondCap,            "floor for back bond valuation");

        FIELD(canConvIntoMandatory,   "true if issuer can convert into a DECS");
        FIELD_MAKE_OPTIONAL(canConvIntoMandatory);
        FIELD(downsideProtection,     "Downside protection for DECS");
        FIELD_MAKE_OPTIONAL(downsideProtection);
        FIELD(decsPremium,            "Premium for DECS");
        FIELD_MAKE_OPTIONAL(decsPremium);
        FIELD(decsMaturity,                  "Maturity of DECS");
        FIELD_MAKE_OPTIONAL(decsMaturity);
        FIELD(decsHoliday,                   "Holiday for DECS");
        FIELD_MAKE_OPTIONAL(decsHoliday);
        FIELD(decsCoupon,                    "Coupon rate for DECS");
        FIELD_MAKE_OPTIONAL(decsCoupon);
        FIELD(decsDCC,                       "Day count convention for DECS");
        FIELD_MAKE_OPTIONAL(decsDCC);
        FIELD(decsBondFrequency,             "Coupon frequency for DECS");
        FIELD_MAKE_OPTIONAL(decsBondFrequency);
        FIELD(decsDivPassThrough,             "Coupon frequency for DECS");
        FIELD_MAKE_OPTIONAL(decsDivPassThrough);

        FIELD(knockedIn,              "whether the convertible has knocked in");
        FIELD_MAKE_OPTIONAL(knockedIn);
        FIELD(knockInDate,            "the date that convertible has knocked in");
        FIELD_MAKE_OPTIONAL(knockInDate);

        //FIELD(cvbModel,                      "model used to price back bond");
        //FIELD(knockInModel,                  "model used to price knock in");

        // transient fields
        FIELD(market,                        "Pointer to market data");
        FIELD_MAKE_TRANSIENT(market);
    }

    static IObject* defaultTracer(){
        return new Tracer();
    }
};

CClassConstSP const Tracer::TYPE = CClass::registerClassLoadMethod(
    "Tracer", typeid(Tracer), TracerHelper::load);
bool  TracerLoad() {
    return (Tracer::TYPE != 0);
}


   

void Tracer::recordOutputRequests(Control* control, Results* results, 
                                  double fairValue) const
{

    static const string method = "Tracer::recordOutputRequests";
    OutputRequest* request = NULL;

    if ( control->isPricing() ) {

        if (control->requestsOutput(OutputRequest::NAKED_BOND_PRICE, request)) {

            // create a risky curve
            IYieldCurveSP tmpCSC = creditSpreads.get()->makeRiskyCurve(*discount.get());
            YieldCurveSP risky(dynamic_cast<YieldCurve*>(tmpCSC.get()));

            double bondFloor = frontBond->presentValue(valueDate, risky);

            results->storeRequestResult(request, bondFloor);
        }
        // Indicative vol
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) {
            /*
            double             indVol     = 0.0;
            CVolRequestConstSP volRequest   = GetLNRequest();
            CVolRequestSP      ncVolRequest =CVolRequestSP::constCast(volRequest);

            LinearStrikeVolRequest* lsVolRequest = dynamic_cast<LinearStrikeVolRequest*>(ncVolRequest.get());
            // interpolate the vol using our LN request
            CVolProcessedBSSP volBS(cvb->asset->getProcessedVol(lsVolRequest));
            // calculate the indicative vol
            try {
                indVol = volBS->CalcVol(cvb->valueDate, cvb->bond->getMaturityDate());
            }
            catch (exception& ) {
                indVol = 0.0;
            }
            results->storeRequestResult(request, indVol);
            */
        }
    }
}

bool Tracer::priceDeadInstrument(CControl* control, CResults* results) const
{   
    static const string method = "Tracer::priceDeadInstrument";

    DateTime matDate = frontBond->getMaturityDate();
    bool   successful = false;

    if (valueDate > matDate) {
        results->storePrice(0., discount->getCcy());
        successful = true;
    };

    return successful;
}


/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool Tracer::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc. */
    return false;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP Tracer::getSensitiveStrikes(OutputNameConstSP outputName,
                                          const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("CVanilla::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    return sensStrikes;
}



/** Returns the name of the instrument's discount currency */
string Tracer::discountYieldCurveName() const {
    return discount.getName();
}


/////////////////////////////////////////////////////////
//           tree1f and fd product
/////////////////////////////////////////////////////////

/** private class */
class TracerProd: virtual public TracerModel::IProduct {
    
public:
    
    TracerProd(const Tracer* instr): tracer(instr) {

        // make the risky curve
     
        // create a risky curve
        IYieldCurveSP tmpCSC = tracer->creditSpreads.get()->makeRiskyCurve(*tracer->discount.get());
        risky = YieldCurveSP(dynamic_cast<YieldCurve*>(tmpCSC.get()));
    }

    virtual ~TracerProd();

    void price(TracerModel*  model,
               Control*      control, 
               CResults*     results) const;
    
private:
    const Tracer* tracer;
    YieldCurveSP risky;
};


TracerProd::~TracerProd() {
}

void TracerProd::price(TracerModel*  model,
                       Control*      control, 
                       CResults*     results) const
{
    double fairValue        = 0.0;
    OutputRequest* request  = NULL;

    if (tracer->knockedIn) {
        DateTime cbValueDate = tracer->knockInDate;
        double faceValue     = tracer->frontBond->getFaceValue();
        
        CAssetWrapper     assetToUse(copy(tracer->asset.get()));
        YieldCurveWrapper    ycToUse(copy(tracer->discount.get()));

        //-------------------------------------------------------------
        //  Value the convertible bond
        //-------------------------------------------------------------
        ConvBondSP cvbToPrice(new ConvBond(ycToUse, tracer->creditSpreads, assetToUse,
                                           tracer->valueDate, tracer->ccyTreatment, faceValue, 
                                           tracer->bondMaturity, tracer->holidays, tracer->bondFrequency, 
                                           tracer->dayCountConvString, tracer->callStartPeriod, 
                                           tracer->callTriggerLevel, tracer->coupon, tracer->redemption, 
                                           tracer->putLevel, tracer->putStartPeriod, tracer->trigger, 
                                           tracer->conversionPremium, cbValueDate));

        IModelSP cvbModel(model->cvbModel.clone());

        cvbToPrice->validatePop2Object();
        cvbToPrice->bond->getMarket(cvbModel.get(), tracer->market.get());

        // allow the model to collect factors from the instrument
        cvbModel->getMarket(tracer->market.get(), IInstrumentCollection::singleton(cvbToPrice));

        // create TweakGroup (holds what we want to tweak)
        TweakGroup   tweakGroup(cvbToPrice, cvbModel);
        SensMgr      sensMgr(&tweakGroup);

        CControlSP cvbControl(Control::makeFromFlags("", 0.0));
        ResultsSP results(tweakGroup.getModel()->Run(tweakGroup.getInstrument(), cvbControl.get()));

        fairValue = results->retrievePrice();
    } else {
        //---------------------------------------------------------
        // set up the convertible schedule and price all converts  
        //---------------------------------------------------------
        DateTimeArray   cvbPricingDates(0);
        DoubleArray     cvbPrices(0);
        DoubleArray     mandatoryPrices(0);

        cvbPricingDates.push_back(tracer->valueDate);

        int idx = 0;

        while ( idx < tracer->cbDates.size() && tracer->cbDates[idx] <= tracer->valueDate ) 
            ++idx;

        while ( idx < tracer->cbDates.size() && 
                tracer->cbDates[idx] <= tracer->frontBond->getMaturityDate()) {
            cvbPricingDates.push_back(tracer->cbDates[idx]);
            ++idx;
        }

        if ( cvbPricingDates[cvbPricingDates.size()-1] < tracer->frontBond->getMaturityDate() ) {
            cvbPricingDates.push_back(tracer->frontBond->getMaturityDate());
        }

        CControlSP cvbControl(Control::makeFromFlags("", 0.0));

        int i;
        for (i=0; i<cvbPricingDates.size() ; ++i) {
            DateTime cbValueDate = cvbPricingDates[i];
            double faceValue     = tracer->frontBond->getFaceValue();

            // create the forward yield curve
            IYieldCurveSP forwardCurve = tracer->discount.get()->createForwardCurve(cbValueDate);

            YieldCurveWrapper fwdCurveWrapper(dynamic_cast<YieldCurve*>(forwardCurve.get()));

            CAssetWrapper assetToUse(copy(tracer->asset.get()));

            /*-------------------------------------------------------------
              Value the convertible bond
              -------------------------------------------------------------*/
            ConvBondSP cvbToPrice(new ConvBond(fwdCurveWrapper, tracer->creditSpreads, assetToUse,
                                               tracer->valueDate, tracer->ccyTreatment, faceValue, 
                                               tracer->bondMaturity, tracer->holidays, tracer->bondFrequency, 
                                               tracer->dayCountConvString, tracer->callStartPeriod, 
                                               tracer->callTriggerLevel, tracer->coupon, tracer->redemption, 
                                               tracer->putLevel, tracer->putStartPeriod, tracer->trigger, 
                                               tracer->conversionPremium, cbValueDate));

            IModelSP cvbModel(model->cvbModel.clone());

            cvbToPrice->validatePop2Object();
            cvbToPrice->bond->getMarket(cvbModel.get(), tracer->market.get());

            // allow the model to collect factors from the instrument
            cvbModel->getMarket(tracer->market.get(), IInstrumentCollection::singleton(cvbToPrice));

            // create TweakGroup (holds what we want to tweak)
            TweakGroup tweakGroup(cvbToPrice, cvbModel);
            TweakGroupSP tweakGroupSP(TweakGroupSP::attachToRef(&tweakGroup));
            // set the new equity spot
            SpotLevel   spotScenario(tracer->trigger);
            spotScenario.findAndShift(tweakGroupSP, OutputNameConstSP());

            // roll to the new value date
            HolidaySP hols(Holiday::weekendsOnly());
            int busDayShift = hols->businessDaysDiff(tracer->valueDate,
                                                     cbValueDate);
            Theta       thetaShift(busDayShift, hols);
            thetaShift.applyScenario(tweakGroupSP);

            if (false) {
                string cvbFile = "c:\\cvb" + Format::toString(i) + ".xml";
                XMLWriter xml(cvbFile.c_str());
                cvbToPrice->write("OBJECT", &xml);
            }

            ResultsSP results(tweakGroup.getModel()->Run(tweakGroup.getInstrument(), cvbControl.get()));

            fairValue = results->retrievePrice();

            fairValue /= tracer->frontBond->getRedemption();

            cvbPrices.push_back(fairValue);

            /*-------------------------------------------------------------
              Value the mandatory 
              -------------------------------------------------------------*/
            if ( tracer->canConvIntoMandatory ) {

                ConvBondSP DECSToPrice(ConvBond::createDECS(fwdCurveWrapper, tracer->creditSpreads, assetToUse,
                                                            tracer->valueDate, tracer->ccyTreatment, faceValue,
                                                            tracer->decsMaturity, tracer->decsHoliday, tracer->decsBondFrequency, 
                                                            tracer->decsDCC, tracer->decsCoupon, 100.0, tracer->trigger,
                                                            tracer->downsideProtection, tracer->decsPremium, 
                                                            tracer->decsDivPassThrough, cbValueDate));

                IModelSP cvbModel(model->cvbModel.clone());

                DECSToPrice->validatePop2Object();
                DECSToPrice->bond->getMarket(cvbModel.get(), tracer->market.get());

                // allow the model to collect factors from the instrument
                cvbModel->getMarket(tracer->market.get(), IInstrumentCollection::singleton(DECSToPrice));

                // create TweakGroup (holds what we want to tweak)
                TweakGroup DECStweakGroup(DECSToPrice, cvbModel);
                TweakGroupSP decsTweakGroupSP(
                    TweakGroupSP::attachToRef(&DECStweakGroup));
                // set the new equity spot
                SpotLevel   DECSspotScenario(tracer->trigger);
                DECSspotScenario.findAndShift(decsTweakGroupSP, 
                                              OutputNameConstSP());

                // roll to the new value date
                thetaShift.applyScenario(decsTweakGroupSP);

                if (false) {
                    string decsFile = "c:\\decs" + Format::toString(i) + ".xml";
                    XMLWriter decs_xml(decsFile.c_str());
                    DECSToPrice->write("OBJECT", &decs_xml);
                }

                ResultsSP DECSresults(tweakGroup.getModel()->Run(DECStweakGroup.getInstrument(), cvbControl.get()));

                mandatoryPrices.push_back(DECSresults->retrievePrice() / tracer->frontBond->getRedemption() );
            }
        }

        // have to get all coupons for front bond and price 

        // create exercise schedule
        DoubleArray   exerSchedRates(0);
        DateTimeArray exerSchedDates(0);
        exerSchedDates.push_back(tracer->frontBond->getMaturityDate());
        exerSchedRates.push_back(tracer->trigger);
        ScheduleSP exerSchedule( new Schedule(exerSchedDates, exerSchedRates, "L"));

        // create upper barrier
        DoubleArray   upperBarrierRates(0);
        DateTimeArray upperBarrierDates(0);
        upperBarrierDates.push_back(tracer->valueDate);
        upperBarrierDates.push_back(tracer->frontBond->getMaturityDate());
        upperBarrierRates.push_back(tracer->trigger);
        upperBarrierRates.push_back(tracer->trigger);

        ScheduleSP upBarrierSchedule( new Schedule(upperBarrierDates, upperBarrierRates, "L"));

        // get all cash flows from the front bond
        CashFlowArraySP frontCashFlows = tracer->frontBond->getCashFlows(tracer->valueDate);

        // create schedule of back bond prices
        ScheduleSP backBondSchedule( new Schedule(cvbPricingDates, cvbPrices, "L"));

        // create schedule of back bond prices for the DECS
        ScheduleSP backBondDECSSchedule;
        if ( tracer->canConvIntoMandatory ) {
            backBondDECSSchedule = ScheduleSP(new Schedule(cvbPricingDates, mandatoryPrices, "L"));
        }

        // create upper rebate schedule
        DateTimeArray rebateSched1(0);
        DoubleArray   rebateSched2(0);


        double bondValue;
        double redemption = tracer->frontBond->getRedemption();

        rebateSched1.push_back(tracer->valueDate); // always add value date - final cash flow is maturity
        if ( tracer->canConvIntoMandatory ) {
            bondValue = Maths::min(backBondSchedule->interpolate(tracer->valueDate), backBondDECSSchedule->interpolate(tracer->valueDate));
        } else {
            bondValue = backBondSchedule->interpolate(tracer->valueDate);
        }
        if (tracer->backBondCapped) {
           bondValue = Maths::min(bondValue, tracer->backBondCap / tracer->frontBond->getRedemption());
        }
        if (tracer->adjustRebateForAccrued) {
            bondValue += tracer->frontBond->getAccruedAtDate(rebateSched1[rebateSched1.size()-1])/redemption;
        }
        rebateSched2.push_back(bondValue);

        for (i=0;i<frontCashFlows->size();++i) {
            if ( (*frontCashFlows)[i].date.rollDate(-1) > rebateSched1[rebateSched1.size()-1] ) {
                rebateSched1.push_back((*frontCashFlows)[i].date.rollDate(-1));
                // determin rebate at date
                if ( tracer->canConvIntoMandatory ) {
                    bondValue = Maths::min(backBondSchedule->interpolate(rebateSched1[rebateSched1.size()-1]), 
                                           backBondDECSSchedule->interpolate(rebateSched1[rebateSched1.size()-1]));
                } else {
                    bondValue = backBondSchedule->interpolate(rebateSched1[rebateSched1.size()-1]);
                }
                if (tracer->backBondCapped) {
                   bondValue = Maths::min(bondValue, tracer->backBondCap / tracer->frontBond->getRedemption());
                }
                if (tracer->adjustRebateForAccrued) {
                    bondValue += tracer->frontBond->getAccruedAtDate(rebateSched1[rebateSched1.size()-1])/redemption;
                }
                rebateSched2.push_back(bondValue);
            }
            if ( (*frontCashFlows)[i].date > rebateSched1[rebateSched1.size()-1] ) {
                rebateSched1.push_back((*frontCashFlows)[i].date);
                // determin rebate at date
                if ( tracer->canConvIntoMandatory ) {
                    bondValue = Maths::min(backBondSchedule->interpolate(rebateSched1[rebateSched1.size()-1]), 
                                           backBondDECSSchedule->interpolate(rebateSched1[rebateSched1.size()-1]));
                } else {
                    bondValue = backBondSchedule->interpolate(rebateSched1[rebateSched1.size()-1]);
                }
                if (tracer->backBondCapped) {
                   bondValue = Maths::min(bondValue, tracer->backBondCap / tracer->frontBond->getRedemption());
                }
                if (tracer->adjustRebateForAccrued) {
                    bondValue += tracer->frontBond->getAccruedAtDate(rebateSched1[rebateSched1.size()-1])/redemption;
                }
                rebateSched2.push_back(bondValue);
            }
        }

        ScheduleSP upRebateSchedule( new Schedule(rebateSched1, rebateSched2, "L"));

        // create the risky discount curve
        IYieldCurveSP riskyCurve = tracer->creditSpreads.get()->makeRiskyCurve(*tracer->discount.get());
        YieldCurveWrapper riskyCurveWrapper(dynamic_cast<YieldCurve*>(riskyCurve.get()));

        DoubleArraySP barrierPrices(new DoubleArray(0));

        // calculate probability weighted value of coupons
        for (i=0;i<frontCashFlows->size();++i) {
            // get cash flow expiry date
            const DateTime& cashFlowDate = (*frontCashFlows)[i].date;
            if (cashFlowDate > tracer->valueDate ) {
                // create exercise schedule
                DoubleArray   exerSchedRates(0);
                DateTimeArray exerSchedDates(0);
                /*
                if ( i == 0 && tracer->valueDate < cashFlowDate) {
                    exerSchedDates.push_back(tracer->valueDate);
                    exerSchedRates.push_back(0.0);
                }
                */
                exerSchedDates.push_back(cashFlowDate);
                exerSchedRates.push_back(0.0);
                ScheduleSP exerSchedule( new Schedule(exerSchedDates, exerSchedRates, "N"));

                InstrumentSettlementSP instSettle(new PhysicalSettlement());
                CDblBarrierSP knockOut(new CDblBarrier(tracer->valueDate,
                                       tracer->valueDate,         // start date
                                       false,                     // not fwd starting
                                       true,                      // one contract
                                       1.0,                       // notional
                                       0.0,                       // initialSpot
                                       instSettle,                // instSettle
                                       InstrumentSettlementSP(   ), // premiumSettle
                                       riskyCurveWrapper,         // use the risky curve for discounting
                                       tracer->asset,             // underlying asset
                                       tracer->ccyTreatment,      // ccy treatment
                                       "BINARY",                  // payoff mode
                                       exerSchedule,              // exercise schedule
                                       false,                     // no early exercise
                                       upBarrierSchedule,         // upper barrier levels
                                       "KO",                      // upper barrier type is KO
                                       ScheduleSP(   ),             // upper rebate schedule
                                       ScheduleSP(   ),             // lower Barrier
                                       "NA",                      // no lower barrier
                                       ScheduleSP(   ),             // no lower rebate
                                       true,                      // intraday monitoring
                                       true,                      // rebate at hit
                                       true,                      // don't scale rebate
                                       "ONCE_TOUCH",              // barrier depdence
                                       "BOTH"));                   // monitoring depdence

                IModelSP koModel(model->knockInModel.clone());

                /*
                knockOut->validatePop2Object();
                knockOut->GetMarket(koModel.get(), tracer->market);
                */

                // allow the model to collect factors from the instrument
                koModel->getMarket(tracer->market.get(), IInstrumentCollection::singleton(knockOut));

                ResultsSP results(koModel->Run(knockOut.get(), cvbControl.get()));
            
                double fairValue = results->retrievePrice();

                if (i == frontCashFlows->size()-1) {
                    barrierPrices->push_back(fairValue * ((*frontCashFlows)[i].amount - tracer->frontBond->getRedemption()));
                } else {
                    barrierPrices->push_back(fairValue * (*frontCashFlows)[i].amount);
                }
            } else {
                barrierPrices->push_back(0.0);
            }
        }

        // calculate value of redemption
        if ( frontCashFlows->size() > 0 ) {
            int redIdx = frontCashFlows->size()-1;
            // get cash flow expiry date
            const DateTime& redemptionDate = (*frontCashFlows)[redIdx].date;
            if ( redemptionDate > tracer->valueDate ) {
                // create exercise schedule
                DoubleArray   exerSchedRates(0);
                DateTimeArray exerSchedDates(0);
                exerSchedDates.push_back(redemptionDate);
                exerSchedRates.push_back(0.0);
                ScheduleSP exerSchedule( new Schedule(exerSchedDates, exerSchedRates, "N"));

                InstrumentSettlementSP instSettle(new PhysicalSettlement());
                CDblBarrierSP knockOut(new CDblBarrier(tracer->valueDate,
                                       tracer->valueDate,         // start date
                                       false,                     // not fwd starting
                                       true,                      // one contract
                                       1.0,                       // notional
                                       0.0,                       // initialSpot
                                       instSettle,                // instSettle
                                       InstrumentSettlementSP(   ), // premiumSettle
                                       riskyCurveWrapper,         // use the risky curve for discounting
                                       tracer->asset,             // underlying asset
                                       tracer->ccyTreatment,      // ccy treatment
                                       "BINARY",                  // payoff mode
                                       exerSchedule,              // exercise schedule
                                       false,                     // no early exercise
                                       upBarrierSchedule,         // upper barrier levels
                                       "KO",                      // upper barrier type is KO
                                       upRebateSchedule,          // upper rebate schedule
                                       ScheduleSP(   ),             // lower Barrier
                                       "NA",                      // no lower barrier
                                       ScheduleSP(   ),             // no lower rebate
                                       true,                      // intraday monitoring
                                       false,                     // rebate at hit
                                       true,                      // don't scale rebate
                                       "ONCE_TOUCH",              // barrier depdence
                                       "BOTH"));                  // monitoring depdence

                IModelSP koModel(model->knockInModel.clone());
                // tracer->knockInModel.clone());

                /*
                knockOut->validatePop2Object();
                knockOut->GetMarket(koModel.get(), tracer->market);
                */

                // allow the model to collect factors from the instrument
                koModel->getMarket(tracer->market.get(), IInstrumentCollection::singleton(knockOut));

                ResultsSP results(koModel->Run(knockOut.get(), cvbControl.get()));
            
                double fairValue = results->retrievePrice();

                barrierPrices->push_back(fairValue * tracer->frontBond->getRedemption());
            }
        }

        // calculate the price
        fairValue = 0.0;
        for (i=0;i<barrierPrices->size();++i) {
            fairValue += (*barrierPrices)[i];
        }

        // Bond value
        if ( control->requestsOutput(OutputRequest::TRACER_CONVERT_PRICE, request) ) {
            CashFlowArraySP convertPrices(new CashFlowArray(cvbPrices.size()));
            for (int i = 0; i < cvbPrices.size(); i++)
            {
                (*convertPrices)[i].date   = cvbPricingDates[i];
                (*convertPrices)[i].amount = cvbPrices[i];
            }
            results->storeRequestResult(request, convertPrices);
        }

        // Accrued interest
        if ( control->requestsOutput(OutputRequest::TRACER_DECS_PRICE, request) ) {
            if ( mandatoryPrices.size() > 0 ) {
                CashFlowArraySP DECSPrices(new CashFlowArray(mandatoryPrices.size()));
                for (int i = 0; i < mandatoryPrices.size(); i++)
                {
                    (*DECSPrices)[i].date   = cvbPricingDates[i];
                    (*DECSPrices)[i].amount = mandatoryPrices[i];
                }
                results->storeRequestResult(request, DECSPrices);
            }
        }
    }

    results->storePrice(fairValue, tracer->discount->getCcy());

    // calculate all output requests
    // take care of additional outputs
    if ( control && control->isPricing() ) {
        DateTime matDate = tracer->frontBond->getMaturityDate();

        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         fairValue,
                                         tracer->valueDate,
                                         tracer->discount.get(),
                                         tracer->asset.get(),
                                         0);
        // IND_VOL
        if (control->requestsOutput(OutputRequest::IND_VOL, request))  {
            // to do 
        }
        
        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       matDate,
                                       tracer->valueDate,
                                       tracer->asset.get());


        tracer->recordOutputRequests(control, results, fairValue);
    }
}

/** Implementation of CFDGridPass::IntoProduct interface */
TracerModel::IProduct* Tracer::createProduct(const TracerModel* model) const {
    return new TracerProd(this);
}

/////////////////////////////////////////////////////////////////////////////////
// Tracer model
/////////////////////////////////////////////////////////////////////////////////

void TracerModel::validatePop2Object(){
    // to test the model per instrument to price
}

// registration, invoked when class is 'loaded'
void TracerModel::load(CClassSP& clazz){
    clazz->setPublic(); 
    REGISTER(TracerModel, clazz);
    SUPERCLASS(CModel);
    EMPTY_SHELL_METHOD(defaultTracerModel);
    FIELD(cvbModel,"Model to price Convertible Bond");
    FIELD(knockInModel,"Model to price Knock-In Barrier");
}

// for TracerModel::IIntoProduct  
void TracerModel::loadIntoProduct(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(TracerModel::IIntoProduct, clazz);
    EXTENDS(Model::IModelIntoProduct);
}

IObject* TracerModel::defaultTracerModel(){
    return new TracerModel();
}

// constructor 
TracerModel::TracerModel():CModel(TYPE) {}

void TracerModel::Price(CInstrument* instrument,
                        CControl*    control,
                        CResults*    results) {
    static const string method = "Tracer::Price";
    IProduct* product = 0;
    try {
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException("Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support TracerModel::IntoProduct");
        }
        if (instrument->priceDeadInstrument(control, results)) {
            return; // done for a dead instrument
        }
        // cast to TracerModel::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        
        // create and price the product 
        product = intoProd.createProduct(this);
        product->price(this, control, results);
        delete product;
    } 
    catch (exception& e) {
        delete product;
        throw ModelException(e, method);
    }
}

MarketObjectSP TracerModel::GetMarket(const MarketData*    market,
                                      const string&        name,
                                      const CClassConstSP& type) const {
    static const string method = "TracerModel::GetMarket";
    try {
        return cvbModel->GetMarket(market, name, type);     
        
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

IModel::WantsRiskMapping TracerModel::wantsRiskMapping() const {
    return (cvbModel->wantsRiskMapping() == riskMappingDisallowed ||
            knockInModel->wantsRiskMapping() == riskMappingDisallowed) ?
                riskMappingDisallowed :
           (cvbModel->wantsRiskMapping() == riskMappingAllowed ||
            knockInModel->wantsRiskMapping() == riskMappingAllowed) ?
                riskMappingAllowed :
           riskMappingIrrelevant;
}

// Type registration
CClassConstSP const TracerModel::TYPE = 
    CClass::registerClassLoadMethod("TracerModel", typeid(TracerModel),
    TracerModel::load);

CClassConstSP const TracerModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("TracerModel::IIntoProduct",
                                    typeid(TracerModel::IIntoProduct), 0);

DRLIB_END_NAMESPACE
