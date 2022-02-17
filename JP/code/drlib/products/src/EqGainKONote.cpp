//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EqGainKONote.cpp
//
//   Description   Equity Gain KO Note.
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/EqGainKONote.hpp"
#include "edginc/Average.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/KOStubRule.hpp"
#include "edginc/KOFixedLeg.hpp"
#include "edginc/KOLiborLeg.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/CallableEquityKOSwap.hpp" // only for Barrier Events!

DRLIB_BEGIN_NAMESPACE


// AvgOutPerf Class, extenal data Class implementation //
void CEqGainKONote::AvgOutPerf::validatePop2Object() {
    static const string method("CEqGainKONote::AvgOutPerf::validatePop2Object");
    try {
        if (genPerfType=="F"||genPerfType=="C"||genPerfType=="P"||genPerfType=="S"){
            if (strikesPct.size() != 1){
                throw ModelException(method, "Number of strikesPct array should be 1 for C,S,P,F case.");
            }
        }
        else if(genPerfType=="CS"||genPerfType=="PS")
        {
            if (strikesPct.size() != 2){
                throw ModelException(method, "Number of strikesPct array should be 2 for CS,PS case.");
            }
        }
        else if (genPerfType=="BDF"){
            if (strikesPct.size() != 3){
                throw ModelException(method, "Number of strikesPct array should be 3 for BDF case.");
            }
        }
        else{
            throw ModelException(method, "Not valid perfType");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// AverageOut Class, internal class implementation //
CEqGainKONote::AverageOut::AverageOut(AvgOutPerfSP avgOutPerf){
    static const string method("CEqGainKONote::AvgOutPerf::setup");
    try {

        // copy data from external data class
        genPerfType = avgOutPerf->genPerfType;
        strikesPct =  avgOutPerf->strikesPct;
        participation = avgOutPerf->participation;

        int nSize = avgOutPerf->avgSamples.size();
        DateTimeArray avgDates(nSize);
        DoubleArray past(nSize);
        DoubleArray weights(nSize);
        for (int i = 0; i < nSize; i++){
            avgDates[i] = avgOutPerf->avgSamples[i].date;
            past[i] = avgOutPerf->avgSamples[i].amount;
            weights[i] = 1.0 / nSize;
        }
        avgOut = SampleListSP(new SampleList(avgDates,
                                             past,
                                             weights));
        // initilize internal data
        isCalls = BoolArray(0);
        strikes = DoubleArray(0);
        longShort = DoubleArray(0);

        // set up the above three array
        if(genPerfType == "F"){
            // Call - Put
            makePayoff(true,strikesPct[0],true);
            makePayoff(false,strikesPct[0],false);
        }
        else if(genPerfType=="C")
            makePayoff(true,strikesPct[0],true);
        else if(genPerfType=="P")
            makePayoff(false,strikesPct[0],true);
        else if(genPerfType=="S"){
            makePayoff(true,strikesPct[0],true);
            makePayoff(false,strikesPct[0],true);
        }
        else if(genPerfType=="CS"){
            makePayoff(true,strikesPct[0],true);
            makePayoff(true,strikesPct[1],false);
        }
        else if(genPerfType=="PS"){
            makePayoff(false,strikesPct[0],false);
            makePayoff(false,strikesPct[1],true);
        }
        else if(genPerfType=="BDF"){
            //long Put, Fwd (C-S), Short Call
            makePayoff(false,1.0+strikesPct[0],true);   // to use the same interface as general performance
            makePayoff(true,strikesPct[1],true);
            makePayoff(false,strikesPct[1],false);
            makePayoff(true,1.0+strikesPct[2],false);   // to use the same interface as general performance
        }
        else
            throw ModelException(method, "Not valid perfType");
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// helpers
class CEqGainKONoteHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(CEqGainKONote, clazz);
        SUPERCLASS(CDblBarrier);
        EMPTY_SHELL_METHOD(defaultEqGainKONote);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        IMPLEMENTS(KnownCashflows::IEventHandler);

        FIELD(SpotFixing,         "Initial Spot Price to decide number of contract for last put/call option.");

        FIELD(fixedLeg, "cash flow of fixed payment");
        FIELD(liborLeg, "payment stream for Float Leg");

        FIELD(ParticipationSchedule, "Participation of Equity Return and payment date");
        FIELD(CapStrikeArray, "Strike Level for Cap Level");
        FIELD_MAKE_OPTIONAL(CapStrikeArray);

        FIELD(CallSchedule, "Call Schedule with strike level.");
        FIELD_MAKE_OPTIONAL(CallSchedule);
        FIELD(isCallable, "Flag for Callable.  Default = false");
        FIELD_MAKE_OPTIONAL(isCallable);

        FIELD(AddStrikeArray1, "Additional Strike Level 1");
        FIELD_MAKE_OPTIONAL(AddStrikeArray1);
        FIELD(AddStrikeArray2, "Additional Strike Level 2");
        FIELD_MAKE_OPTIONAL(AddStrikeArray2);
        FIELD(PartStrike1, "Participation of Equity Return for Additional Strike 1");
        FIELD_MAKE_OPTIONAL(PartStrike1);
        FIELD(PartStrike2, "Participation of Equity Return for Additional Strike 2");
        FIELD_MAKE_OPTIONAL(PartStrike2);
        FIELD(PayType1, "Payoff type (C, P, F) for Additional Strike 1");
        FIELD_MAKE_OPTIONAL(PayType1);
        FIELD(PayType2, "Payoff type (C, P, F) for Additional Strike 2");
        FIELD_MAKE_OPTIONAL(PayType2);

        FIELD(pastValues, "past equity level (Do not Use!)");
        FIELD_MAKE_OPTIONAL(pastValues);
        FIELD(samples, "past equity level (Use This!)");
        FIELD_MAKE_OPTIONAL(samples);

        FIELD(isCalled, "if it's already called, true.");
        FIELD_MAKE_OPTIONAL(isCalled);
        FIELD(callDate, "if it's already called, needs the date.");
        FIELD_MAKE_OPTIONAL(callDate);

        FIELD(segDensity, "optional segment density to set more grid before next CritDates.  no use = 0 (default)");
        FIELD_MAKE_OPTIONAL(segDensity)

        FIELD(SwitchPlainSwap, "Price with PlainSwap or Not.  Default = True (with PlainSwap).");
        FIELD_MAKE_OPTIONAL(SwitchPlainSwap);
        //  When SwitchPlainSwap = false, There is no KNOWN_CASH_FLOW info in this product.

        FIELD(RemoveSettlementOption, "default = false. If final put/fwd not needed, make it true");
        FIELD_MAKE_OPTIONAL(RemoveSettlementOption);

        FIELD(hasAvgOut,"has Average Option?");
        FIELD_MAKE_OPTIONAL(hasAvgOut);
        FIELD(avgOutPerf, "You can add average-out option to payoff, after all early teminate dates.");
        FIELD_MAKE_OPTIONAL(avgOutPerf);

        // used internally
        FIELD(scalingStrike, "Scales strike for forward starting options");
        FIELD_MAKE_TRANSIENT(scalingStrike);
    }

static IObject* defaultEqGainKONote(){
        return new CEqGainKONote();
    }
};


CClassConstSP const CEqGainKONote::TYPE = CClass::registerClassLoadMethod(
    "EqGainKONote", typeid(CEqGainKONote), CEqGainKONoteHelper::load);
bool  CEqGainKONoteLoad() {
    return (CEqGainKONote::TYPE != 0);
}


CClassConstSP const CEqGainKONote::AvgOutPerf::TYPE = CClass::registerClassLoadMethod(
    "EqGainKONote::AvgOutPerf", typeid(CEqGainKONote::AvgOutPerf), CEqGainKONote::AvgOutPerf::load);

// constructor
CEqGainKONote::CEqGainKONote(): CDblBarrier(TYPE)
{
    isCallable = false;     //Flag for Callable.  Default = false
    PayType1 = "N";         // N explains "Same as PayMode"
    PayType2 = "N";         // N explains "Same as PayMode"
    Schedule pastValues (DateTimeArray(0),DoubleArray(0),"N");
//    CashFlowArray samples(CashFlowArray(0));
    isCalled = false;
    SwitchPlainSwap = true;
    RemoveSettlementOption = false;
    segDensity = 0;
    hasAvgOut = false;

    scalingStrike = 1.0;
};


void CEqGainKONote::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{

    int i, j;
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        int avSize = 0;
        DateTime avgMatDate = DateTime(0,0);
        if(hasAvgOut){
            avSize = avgOutPerf->avgSamples.size();
            avgMatDate = avgOutPerf->avgSamples[avSize-1].date;
        }

        OutputRequest* request = NULL;

        DateTime matDate= exerciseSchedule->lastDate();
        // DELAY_PRICE
         InstrumentUtil::delayPriceHelper(control,
                                          results,
                                          fairValue,
                                          valueDate,
                                          discount.get(),
                                          asset.get(),
                                          premiumSettle.get());
        // IND_VOL
        InstrumentUtil::recordIndicativeVol(control,results,indVol);

        // FWD_AT_MAT
        try{
            InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       matDate,
                                       valueDate,
                                       asset.get());
        }
        catch(exception&)
        {// continue if fwd failed - this hapens now for quanto asset with CEVj vol
        }

        //calculate Swap related Value and Store Values
        //libor
        double legResult;
        CashFlowArrayConstSP CflL = liborLeg->getCashFlowArray(valueDate, discount.get());
        CashFlowArraySP cfl (new CashFlowArray(0));
        if (CflL->size()>0)
       {
            for (i=0; i<liborLeg->getSize(); i++)
                    cfl->push_back((*CflL)[i]);
		    results->storeGreek(cfl, Results::DEBUG_PACKET, OutputNameSP(new OutputName("liborLeg")));
        }
        legResult = liborLeg->getPV(valueDate, discount.get());
        results->storeScalarGreek(legResult, Results::DEBUG_PACKET, OutputNameSP(new OutputName("liborLegValue")));
        // add to output
        request = control->requestsOutput(OutputRequest::LIBOR_LEG_VALUE);
        if (request) {
            results->storeRequestResult(request, legResult);
        }

        //fixed cash flow
        CashFlowArrayConstSP CflF = fixedLeg->getCashFlowArray();
        CashFlowArraySP cff (new CashFlowArray(0));
        if (CflF->size()>0)
        {
            for (j=0; j<fixedLeg->getSize(); j++)
                    cff->push_back((*CflF)[j]);
            results->storeGreek(cff, Results::DEBUG_PACKET, OutputNameSP(new OutputName("fixedLeg")));
        }
        legResult = fixedLeg->getPV(valueDate, discount.get());
		results->storeScalarGreek(legResult, Results::DEBUG_PACKET, OutputNameSP(new OutputName("fixedLegValue")));
        // add to output
        request = control->requestsOutput(OutputRequest::FIXED_LEG_VALUE);
        if (request) {
            results->storeRequestResult(request, legResult);
        }

        // Add PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArray date = ParticipationSchedule->getDates();

            if (CflL->size()>0)
            {
                for (i=0; i<liborLeg->getSize(); i++) {
                    date.push_back((*CflL)[i].date);
                }
            }
            if (CflF->size()>0)
            {
                for (j=0; j<fixedLeg->getSize(); j++) {
                    date.push_back((*CflF)[j].date);
                }
            }

            if (hasAvgOut)
                date.push_back(instSettle->settles(avgMatDate, asset.get()));

            OutputRequestUtil::recordPaymentDates(control,results,&date);
        }
        // Add KNOW_CASHFLOWS
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArraySP knownCfl = getKnownCashFlow();
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    knownCfl.get());   
        }

        // Add BARRIER_LEVEL
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        bool lowerOut = LowerBarBreached  && LowerBarType == "KO";
        bool upperOut = UpperBarBreached  && UpperBarType == "KO";
        bool isOut = ((lowerOut || upperOut) && BarrierDependence != "TWO_TOUCH")
                     || (lowerOut && upperOut);
        isOut = isCalled ? true : isOut;        //Adding Callable Info
        if (request && !fwdStarting && !isOut) {
            // report barrier levels over a date range
            DateTime upperDate = BarrierLevel::barrierWindow(valueDate);

            BarrierLevelArraySP levels(new BarrierLevelArray(0));

            bool isContUp = false;
            bool isContDn = false;
            if (IntraDayMonitor)
                isContUp = isContDn = true;
            else if (MonitoringDependence == "UPPER")
                isContDn = true;
            else if (MonitoringDependence == "LOWER")
                isContUp = true;

            if (LowerBarType != "NA") {
                // use economic barrier (if it exists)
                Schedule* s = lowerEcoBarrier.get() ? lowerEcoBarrier.get(): 
                    LowerBarrier.get();
                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(false,(*subset)[i].date,(*subset)[i].amount,isContDn);
                    levels->push_back(bl);
                }
            }

            if (UpperBarType != "NA") {
                // use economic barrier (if it exists)
                Schedule* s = upperEcoBarrier.get() ? upperEcoBarrier.get(): 
                    UpperBarrier.get();
                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(true,(*subset)[i].date,(*subset)[i].amount,isContUp);
                    levels->push_back(bl);
                }
            }

            if (!levels->empty()) {
                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       levels.get());
            }
        }            

    }
}

CashFlowArraySP CEqGainKONote::getKnownCashFlow() const{
    static const string method = "CEqGainKONote::getKnownCashFlow";
    try {     
        // these need to be in increasing date order to aggregate
        // hence the merge rather than push_back which handles this

        int i, j;

        int avSize = 0;
        DateTime avgMatDate = DateTime(0,0);
        if(hasAvgOut){
            avSize = avgOutPerf->avgSamples.size();
            avgMatDate = avgOutPerf->avgSamples[avSize-1].date;
        }

        //fixed cash flow
        CashFlowArrayConstSP CflF = fixedLeg->getCashFlowArray();
        CashFlowArrayConstSP CflL = liborLeg->getCashFlowArray(valueDate, discount.get());

        CashFlowArraySP equityCFL(new CashFlowArray(0));;                              
        // if there is known equity payoff, need to handle it.
        DateTimeArray payDates = ParticipationSchedule->getDates();
        DoubleArray part = ParticipationSchedule->getValues();
        DoubleArray pastLvl = DoubleArray(0);
        DoubleArray strikeArray = exerciseSchedule->getValues();
        DateTimeArray exeDates = exerciseSchedule->getDates();
        if (samples.size()==0) {
            for (i=0; i<exeDates.size(); i++)
                pastLvl.push_back(-1.0); //negative means not observed yet.
        }
        else {
            for (i=0; i<samples.size(); i++)   
                pastLvl.push_back(samples[i].amount);
        }
        bool isCap = CapStrikeArray.size() > 0 ? true : false ;
        bool isAdd1= PartStrike1.size() > 0 ? true : false ;
        bool isAdd2= PartStrike2.size() > 0 ? true : false ;
        bool isRebatePaidBefore = false;

        
        // for (i=0; i<payDates.size(); i++)
        double refLevel = SpotFixing<0? initialSpot : SpotFixing;
        for (i=0; i<pastLvl.size(); i++)
        {
            double strike = strikeArray[i];
            double strikeC = isCap ? CapStrikeArray[i] : strikeArray[i];    //strikeArray is used just for dummy
            double strike1 = isAdd1 ? AddStrikeArray1[i] : strikeArray[i];  //strikeArray is used just for dummy
            double strike2 = isAdd2 ? AddStrikeArray2[i] : strikeArray[i];  //strikeArray is used just for dummy
            
            if (exeDates[i] < valueDate)
            {// negative means not yet known.
                double eqGain  = part[i] * GetIntrinsic(pastLvl[i],strike ,isCall, PayoffMode != "FORWARD");
                if (isCap)
                        eqGain -= part[i] * GetIntrinsic(pastLvl[i],strikeC,isCall, PayoffMode != "FORWARD");
                if (isAdd1)
                        eqGain += PartStrike1[i] * GetIntrinsic(pastLvl[i],strike1, PayType1 == "C", PayType1 != "F");
                if (isAdd2)
                        eqGain += PartStrike2[i] * GetIntrinsic(pastLvl[i],strike2, PayType2 == "C", PayType2 != "F");

                eqGain *= notional/refLevel;

                CashFlow equityPay(payDates[i],eqGain);
                equityCFL->push_back(equityPay);
            }

            if ((UpperBarBreached || LowerBarBreached) && !isRebatePaidBefore)
            {// counting Rebate.
                if (payDates[i] >= UpperBarBreachDate || payDates[i] >= UpperBarBreachDate )
                {// Rebate payment
                    if (UpperBarBreached && !!UpperRebate){
                        if (UpperRebate->length()>0){
                            CashFlow equityPay(payDates[i],UpperRebate->interpolate(UpperBarBreachDate));
                            equityCFL->push_back(equityPay);
                        }
                    }
                    else if(!!LowerRebate) {
                        if (LowerRebate->length()>0){
                            CashFlow equityPay(payDates[i],LowerRebate->interpolate(LowerBarBreachDate));
                            equityCFL->push_back(equityPay);
                        }
                    }
                    isRebatePaidBefore = true;      // no more rebate in future payment.
                }
            }

        }

        // average payoff
        if(hasAvgOut){
            if(avgMatDate < valueDate){
                CEqGainKONote::AverageOutSP aOP = CEqGainKONote::AverageOutSP(
                                                                new CEqGainKONote::AverageOut(avgOutPerf)); 
                double avgPrice = avgValue(aOP, refLevel, refLevel, notional, valueDate);            
                DateTime avgPayDate  = instSettle->settles(avgMatDate, asset.get());
                equityCFL->push_back(CashFlow(avgPayDate, avgPrice));     
            }
            CashFlow::aggregate(*equityCFL);
        }
        
        // Pick up only known Fixed Leg.
        CashFlowArraySP knownCflF (new CashFlowArray(0));
        if (CflF->size()>0 && SwitchPlainSwap)
        {
            for (j=0; j<fixedLeg->getSize(); j++)
            {
                if(valueDate > fixedLeg->AccrueStartDates[j])
                    knownCflF->push_back((*CflF)[j]);
            }
        }
        // Pick up only known Libor Leg.
        CashFlowArraySP knownCflL (new CashFlowArray(0));
        if (CflL->size()>0 && SwitchPlainSwap)
        {
            for (i=0; i<liborLeg->getSize(); i++) {
                if(valueDate > liborLeg->AccrualDates[i])
                    knownCflL->push_back((*CflL)[i]);
            }
        }

        // Called case
        CashFlowArraySP knownCallK (new CashFlowArray(0));
        if (isCalled)
        {
            DoubleArray callStrikes = CallSchedule->getValues();
            for (i=0; i<payDates.size(); i++) {
                if (callDate < payDates[i]){
                    knownCallK->push_back(CashFlow(payDates[i], notional * callStrikes[i]));
                    break;
                }
            }
        }

        // now glue them all together
        CashFlowArraySP merge1(CashFlow::merge(knownCflF, knownCflL));
        CashFlowArraySP merge2(CashFlow::merge(merge1, equityCFL));
        CashFlowArraySP knownCfl(CashFlow::merge(merge2, knownCallK));

        return knownCfl;
    }
    catch (exception& e) 
	{
		throw ModelException(e, method);
	}
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CEqGainKONote::priceDeadInstrument(CControl* control, CResults* results) const
{
    DateTime  exerDate;

    static string method = "CEqGainKONote::priceDeadInstrument";
    bool isSameTime = IntraDayMonitor? false : true;
    int i;

    // check the past values
    int pastSize = samples.size();
    DateTimeArray pastDates = CashFlow::dates(samples);
    DoubleArray closings;
    DateTimeArray barDates;
    bool hasEcoUpBar = upperEcoBarrier.get() && !upperEcoBarrier->getDates().empty();
    bool hasEcoLoBar = lowerEcoBarrier.get() && !lowerEcoBarrier->getDates().empty();
    Schedule* ecoBar;

    for (i=0; i<pastSize; i++)
        closings.push_back(samples[i].amount);
    if (pastSize > 0){
        if (!UpperBarBreached && UpperBarrier->length()>0){
            ecoBar = hasEcoUpBar ? upperEcoBarrier.get() : UpperBarrier.get();
            barDates = ecoBar->getDates();
            for (i=0; i<pastSize && pastDates[i] < valueDate; i++){
                if (Neighbour(pastDates[i], barDates, 0, barDates.size()-1, 0)>=0){
                    if (closings[i] > ecoBar->interpolate(pastDates[i])*(1.0-FP_MIN)){
                        UpperBarBreached = true;
                        UpperBarBreachDate  = pastDates[i];
                    }
                }
            }
        }
        if (!LowerBarBreached && LowerBarrier->length()>0){
            ecoBar =hasEcoLoBar ? lowerEcoBarrier.get() : LowerBarrier.get();
            barDates = ecoBar->getDates();
            for (i=0; i<pastSize && pastDates[i] < valueDate; i++){
                if (Neighbour(pastDates[i], barDates, 0, barDates.size()-1, 0)>=0){
                    if (closings[i] < ecoBar->interpolate(pastDates[i])*(1.0+FP_MIN)){
                        LowerBarBreached = true;
                        LowerBarBreachDate  = pastDates[i];
                    }
                }
            }
        }
    }

    if (!fwdStarting)
    {
        double s = asset->getSpot();

        // check breach of barriers on value date for KO
        if ((UpperBarType == "KO") && (UpperBarBreached == false)
            && UpperBarrier->getDates()[0] <= valueDate
            && UpperBarrier->lastDate() >= valueDate) // barrier started and not finished.
        {
            if(UpperBarrier->getInterp() == "N")
            {
                for(i = 0; i < UpperBarrier->length(); i++)
                {
                    if (UpperBarrier->getDates()[i].equals(valueDate,isSameTime))
                    {
                        ecoBar = hasEcoUpBar ? upperEcoBarrier.get() : UpperBarrier.get();
                        if (s > ecoBar->interpolate(valueDate)*(1.0-FP_MIN))
                        {
                            UpperBarBreached = true;
                            UpperBarBreachDate = valueDate;
                        }
                        break;
                    }
                }
            }
            else
            {
                if (IntraDayMonitor
                    || valueDate.getTime() == UpperBarrier->getDates()[0].getTime()){
                    ecoBar = hasEcoUpBar ? upperEcoBarrier.get() : UpperBarrier.get();
                    if (s > ecoBar->interpolate(valueDate)*(1.0-FP_MIN))
                    {
                        UpperBarBreached = true;
                        UpperBarBreachDate = valueDate;
                    }
                }
            }
        }
        if ((LowerBarType == "KO") && (LowerBarBreached == false)
            && LowerBarrier->getDates()[0] <= valueDate
            && LowerBarrier->lastDate() >= valueDate) // barrier started and not finished.
        {
            if(LowerBarrier->getInterp() == "N")
            {
                for(i = 0; i < LowerBarrier->length(); i++)
                {
                    if (LowerBarrier->getDates()[i].equals(valueDate,isSameTime))
                    {
                        ecoBar =hasEcoLoBar ? lowerEcoBarrier.get() : LowerBarrier.get();
                        if (s < ecoBar->interpolate(valueDate)*(1.0+FP_MIN))
                        {
                            LowerBarBreached = true;
                            LowerBarBreachDate = valueDate;
                        }
                        break;
                    }
                }
            }
            else
            {
                if (IntraDayMonitor
                    || valueDate.getTime() == LowerBarrier->getDates()[0].getTime()){
                    ecoBar =hasEcoLoBar ? lowerEcoBarrier.get() : LowerBarrier.get();
                    if (s < ecoBar->interpolate(valueDate)*(1.0+FP_MIN))
                    {
                        LowerBarBreached = true;
                        LowerBarBreachDate = valueDate;
                    }
                }
            }
        }
    }

    if (isExercised && ((LowerBarBreached  && LowerBarType == "KO") ||
                        (UpperBarBreached  && UpperBarType == "KO")))
        throw ModelException(method,
                             "Can not breach KO barrier and exercise.");


    DateTime matDate = exerciseSchedule->lastDate();

    bool lowerOut = LowerBarBreached  && LowerBarType == "KO";
    bool upperOut = UpperBarBreached  && UpperBarType == "KO";

    bool isOut = ((lowerOut || upperOut) && BarrierDependence != "TWO_TOUCH")
        || (lowerOut && upperOut);

    isOut = isCalled ? true : isOut;        //Adding Callable Info

    int k;
    DateTimeArray payDates = ParticipationSchedule->getDates();
    DoubleArray part= ParticipationSchedule->getValues();
    DateTimeArray exeDates = exerciseSchedule->getDates();
    DoubleArray strikes = exerciseSchedule->getValues();

    //bool isDead = valueDate >= matDate || (isExercised && canExerciseEarly) || isOut;
    bool isDead = valueDate >= matDate || isOut;
    if (!isDead)
        return false; // not dead yet

    DateTime settlementDate = instSettle->settles(matDate, asset.get());

    // check already settled before final settlement of Strucutore
    if (isOut)
    {//find out next payment dates
        if (isCalled && (k = Neighbour(callDate, payDates, 0, payDates.size()-1, 1)) >= 0){
            settlementDate = payDates[k];
        }
        else if((k = Neighbour(lowerOut ? LowerBarBreachDate : UpperBarBreachDate,
                               payDates, 0, payDates.size()-1, 1)) >= 0){
            settlementDate = payDates[k];
        }
        //if cannot find, use settlemenDate of whole structure.
    }

    if (valueDate >= settlementDate)
    {// settled already
        results->storePrice(0.0, discount->getCcy());
        addOutputRequests(control, results, 0.0, 0.0);
        return true;
    }

    double  value = 0.0;
    // sort out ko case first
    if (isOut) // two touch case has problem choosing rebate - use lower for now
    {
        if (0 <= k && k < payDates.size() ){
            // ignore the timing, as Pyramid EOD process run as of SoD.
            if (exeDates[k].daysDiff(valueDate)>0)
                return false;       //not yet detemined final equity leg.
        }
        if (settlementDate >= valueDate) // should really be > but then can't record cash flow if no delay
        {
            if (lowerOut && !!LowerRebate)
                value = LowerRebate->interpolate(LowerBarBreachDate);
            else if (upperOut && !!UpperRebate)
                value = UpperRebate->interpolate(UpperBarBreachDate);
            else if (isCalled)
                value = CallSchedule->interpolate(callDate) * notional;
        }
        else{// already terminated.  Just return 0.0.
            results->storePrice(0.0, discount->getCcy());
            return true;
        }
    }
    else if (settlementDate < valueDate){
        results->storePrice(0.0, discount->getCcy());
        return true;
    }

    double refLevel = SpotFixing<0? initialSpot : SpotFixing;
    if ((k = Neighbour(valueDate, exeDates, 0, exeDates.size()-1, -1)) >= 0){
        if (exeDates[k] <= valueDate && valueDate < payDates[k]){
            double spot = (exeDates[k] == valueDate)? asset->getSpot() : samples[k].amount;
            double eqGain = part[k]* GetIntrinsic(spot,strikes[k],isCall, PayoffMode != "FORWARD");

            if (CapStrikeArray.size() > 0)
                eqGain -= part[k]* GetIntrinsic(spot,CapStrikeArray[k],isCall, PayoffMode != "FORWARD");
            if (AddStrikeArray1.size() > 0)
                eqGain += PartStrike1[k] * GetIntrinsic(spot,AddStrikeArray1[k], PayType1 == "C", PayType1 != "F");
            if (AddStrikeArray2.size() > 0)
                eqGain += PartStrike2[k] * GetIntrinsic(spot,AddStrikeArray2[k], PayType2 == "C", PayType2 != "F");

            eqGain *= notional/refLevel;

            //put the value to settlement date due to discount later.
            eqGain *= discount->pv(settlementDate,payDates[k]);

            value += eqGain;    //pay Equity Leg.
        }
    }

    // pv from settlement to today
    value *= discount->pv(valueDate, settlementDate);

    if (hasAvgOut && !isOut){   // when isOut, AvgOut should be KO.
        CEqGainKONote::AverageOutSP aOP = CEqGainKONote::AverageOutSP(new CEqGainKONote::AverageOut(avgOutPerf));
        DateTime avgSettle = instSettle->settles(aOP->avgOut->getLastDate(), asset.get());
        if (valueDate < avgSettle) {
            value += avgValue(aOP, refLevel, refLevel, notional, valueDate);
        }
    }

    // add Swap Leg (Pay FixLeg / Recieve Libor Leg)
    KOStubRuleSP koRule = KOStubRuleSP(new KOStubRule("B",false,"AsSchedule"));    //koStubRule,isAccrueUpToSettle,payTimingRule
    KOFixedLegSP koFix = KOFixedLegSP(new KOFixedLeg(fixedLeg,koRule,instSettle));
    KOLiborLegSP koFlt = KOLiborLegSP(new KOLiborLeg(liborLeg,koRule,instSettle,valueDate,discount));

    value -= koFlt->getKOPV(valueDate);
    value -= koFix->getKOPV(valueDate,discount.get());

    if (SwitchPlainSwap){
        value += liborLeg->getPV(valueDate, discount.get());
        value += fixedLeg->getPV(valueDate, discount.get());
    }

    // store results
    results->storePrice(value, discount->getCcy());
    addOutputRequests(control, results, value, 0.0);

    return true;
}

// **** now a copy of Vanilla, to do : add more validations
void CEqGainKONote::Validate()
{
    static const string method = "CEqGainKONote::Validate";
    int i;
    // just check the things that aren't/cannot be checked in
    // validatePop2Object

    CDblBarrier::Validate();

    if (RebateAtMat){
        throw ModelException(method, "Rebate At Mat should be false!!");
    }

    const DateTime& matDate= exerciseSchedule->lastDate();
    int numItems = liborLeg->getSize();
    if (liborLeg->AccrualDates[numItems-1] > matDate)
        throw ModelException(method, "LiborAccrualDates["+Format::toString(numItems)+"] is later than maturity");
	numItems = fixedLeg->getSize();
    if (fixedLeg->AccrueStartDates[numItems-1] > matDate)
        throw ModelException(method, "FixedAccrueStartDates["+Format::toString(numItems)+"] is later than maturity");
    if (numItems > 0){
        if (fixedLeg->isDCC())
            throw ModelException(method, "isDcc = true in FixedLeg is not allowed.");
    }

    //following validation should be relaxed if you fixed the bug in InitProd.
    //In current routine, if the accrue start date is later than previous coupon payment date,
    //It takes always 0 value on Accrue start date.  It's bug.
/*    for (i=1 ; i<numItems; i++){
        if (fixedLeg->AccrueStartDates[i]> fixedLeg->PaymentDatesArray[i-1])
            throw ModelException(method, "FixedAccrueStartDates["+Format::toString(i)+"] should be earlier than previous PaymentDatesArray");
    }
*/

    int numK1 = exerciseSchedule->length();
    int numK2 = 0;
    if (CapStrikeArray.size()>0)
    {// check the size of cap strikes, it it's given.
        numK2 = CapStrikeArray.size();
	    if (numK2 != numK1 && numK2 > 0 )
            throw ModelException(method, "CapSchedule ["+Format::toString(numK2)+"] are given, but it should be same size as excersizeSchedule ["+Format::toString(numK1)+"]. ");
    }

    int numPart = ParticipationSchedule->length();
    if (numK1 != numPart)
        throw ModelException(method, "Number of Participation Schedule ["+Format::toString(numPart)+"] is not same as number of exerciseSchedule ["+Format::toString(numK1)+"]. ");

    if (AddStrikeArray1.size()>0)
    {// check the size of cap strikes, it it's given.
        numK2 = AddStrikeArray1.size();
        numPart = PartStrike1.size();
	    if (numK2 != numK1)
            throw ModelException(method, "AddStrikeArray1 ["+Format::toString(numK2)+"] are given, but it should be same size as excersizeSchedule ["+Format::toString(numK1)+"]. ");
        if (numK2 != numPart)
            throw ModelException(method, "PartStrike1 ["+Format::toString(numPart)+"] is not same size as  AddStrikeArray1["+Format::toString(numK2)+"]. ");
    }else if (PartStrike1.size()>0){
        throw ModelException(method, "PartStrike1 should be empty when AddStrikeArray1 is empty.");
    }

    if (AddStrikeArray2.size()>0)
    {// check the size of cap strikes, it it's given.
        numK2 = AddStrikeArray2.size();
        numPart = PartStrike2.size();
	    if (numK2 != numK1)
            throw ModelException(method, "AddStrikeArray2 ["+Format::toString(numK2)+"] are given, but it should be same size as excersizeSchedule ["+Format::toString(numK1)+"]. ");
        if (numK2 != numPart)
            throw ModelException(method, "PartStrike2 ["+Format::toString(numPart)+"] is not same size as  AddStrikeArray2["+Format::toString(numK2)+"]. ");
    }else if (PartStrike2.size()>0){
        throw ModelException(method, "PartStrike2 should be empty when AddStrikeArray2 is empty.");
    }

    ///// arrange the sample event ////
    // disallow to have both of two type sampled value
//    if (samples.size()>0 && !!pastValues)
        //throw ModelException(method, "There are data in both of pastValues and samples.  Use samples, only.");

    // copy pastValues to samples
    if (samples.size()==0 && !!pastValues){
        DateTimeArray sampleDates = pastValues->getDates();
        DoubleArray sampleValues = pastValues->getValues();
        //samples = CashFlowArray(0);
        samples.resize(sampleDates.size());
        for (i = 0; i < sampleDates.size(); i++){
            samples[i] = CashFlow(sampleDates[i],sampleValues[i]);
            //samples[i].date = sampleDates[i];
            //samples[i].amount = sampleValues[i];
            //samples->push_back(CashFlow(sampleDates[i],sampleValues[i]));
        }
    }

    if (samples.size()>0){
        DateTimeArray exeDates = exerciseSchedule->getDates();
        for (i=0;i<samples.size();i++){
            // this validation should be before for loop, but not compatible existing cases...
            // as some inst doesn't have sample if it's in future.
            if (i >= exeDates.size()){
                throw ModelException(method, "Samples's size ["+Format::toString(samples.size())+
                                             "] should be equal to exerciseSchedule's size ["
                                             +Format::toString(exeDates.size())+"] .");
            }
            if (exeDates[i] < valueDate && exeDates[i] != samples[i].date ){
                throw ModelException(method, "the exersiceSchedule ("+
                                    (exeDates[i].toString()) +
                                     ") is not same as the samples dates (" +
                                     (samples[i].date.toString()) + ")");
            }
        }
    }

    if (PayoffMode == "BINARY")
        throw ModelException(method, "BINARY payoff are not supported.  Please use call-spread by setting cap level.");

    if (valueDate > exerciseSchedule->firstDate() && samples.size() == 0)
        throw ModelException(method, "valueDate is later than first exercise Date.  You need samples.");

    if (isCalled && UpperBarBreached)
        throw ModelException(method, "isCalled and UpperBarBreached are true.  It's not permitted.");
    if (isCalled && LowerBarBreached)
        throw ModelException(method, "isCalled and LowerBarBreached are true.  It's not permitted.");

    if (UpperBarType == "KI" || LowerBarType == "KI")
        throw ModelException(method, "Knock-In barrier is not supported.");

    if (!!CallSchedule){
        if (CallSchedule->getInterp() != "N")
            throw ModelException(method, "Interp method of CallSchedule is only available 'N' only.");
        if (CallSchedule->getValues().size() <=0 && isCallable)
            throw ModelException(method, "isCallable is true, but no CallSchedule found");
    }

    if (!!UpperBarrier){
        if (UpperBarrier->getInterp() != "N")
            throw ModelException(method, "Interp method of UpperBarrier is only available 'N' only.");
    }

    if (!!LowerBarrier){
        if (LowerBarrier->getInterp() != "N")
            throw ModelException(method, "Interp method of LowerBarrier is only available 'N' only.");
    }

    if (hasAvgOut){
        if(!avgOutPerf.get()){
            throw ModelException(method, "hasAvgOut = true, but no Average Out data are given");
        }
        DateTime firstAvgDate = avgOutPerf->avgSamples[0].date;
        if (!!UpperBarrier){
            if (UpperBarrier->length()>0 && UpperBarType != "NA"){
                if (firstAvgDate.daysDiff(UpperBarrier->lastDate()) < 2){
                    throw ModelException(method,
                                        "The first date of Average sample [" + firstAvgDate.toString() +
                                        "] must be at least 1 day after the last Upper Barrier Date.");
                }
            }
        }
        if (!!LowerBarrier){
            if (LowerBarrier->length()>0 && LowerBarType != "NA"){
                if (firstAvgDate.daysDiff(LowerBarrier->lastDate()) < 2){
                    throw ModelException(method,
                                        "The first date of Average sample [" + firstAvgDate.toString() +
                                        "] must be at least 1 day after the last Lower Barrier Date.");
                }
            }
        }
        if (!!CallSchedule){
            if(CallSchedule->length()>0){
                if(firstAvgDate.daysDiff(CallSchedule->lastDate()) < 2){
                    throw ModelException(method,
                                        "The first date of Average sample [" + firstAvgDate.toString() +
                                        "] must be at least 1 day after the last Call Date.");
                }
            }
        }
    }

}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CEqGainKONote::sensShift(Theta* shift)
{
    DateTime newDate = shift->rollDate(valueDate);
//    DateTime newDate = valueDate.rollDate(1);
    DateTimeArray exeDates = exerciseSchedule->getDates();

    if(samples.size()>0){
        int k = Neighbour(newDate, CashFlow::dates(samples),0,samples.size()-1,-1);
        if (k >=0){
            // When today is monitoring date, Theta Tweak requires the past sampling
            if (valueDate < samples[k].date && samples[k].date <= newDate)
                samples[k].amount = asset->getSpot();
        }
    }
    bool store_fwdStarting = fwdStarting;
    CDblBarrier::sensShift(shift);      //call parent class's sensShift
    // When theta tweak makes started option, need to scale the additional strike level, too.
    if (store_fwdStarting == true && fwdStarting == false)
        scalingStrike = initialSpot;
    else
        scalingStrike = 1.0;

    return true;
};

double CEqGainKONote::avgValue(const AverageOutSP aop, double const spot, double SpotRef, const DateTime& when) const {
    double s = Maths::max  (spot,0.00001);
    double strikeScale = SpotRef/s;
    double notl = s/SpotRef * notional;
    double refLevel = s;
    return avgValue(aop, strikeScale, refLevel, notl, when);
}
double CEqGainKONote::avgValue(const AverageOutSP aop,
                               double const strikeScale,
                               double const refLevel,
                               double const notl,
                               const DateTime& when) const {
    static const string method = "CEqGainKONote::avgValue";
    try {
        double value = 0.0;

        CClosedFormLN cfln("VolPreferred");

        // use the pastAvgOut if there is.  To DO...
        // SampleListSP sample = !aop->pastAvgOut.get() ? aop->avgOut: aop->pastAvgOut;

        for (int i = 0; i < aop->strikes.size(); i++) {
            InstrumentSP avg(Average::makeAvgSpot(aop->isCalls[i],
                                                  aop->avgOut->getLastDate(),
                                                  aop->strikes[i]*strikeScale,
                                                  aop->avgOut.get(),
                                                  instSettle.get(),
                                                  premiumSettle.get(),
                                                  asset.get(),
                                                  ccyTreatment,
                                                  discount.get(),
                                                  valueDate,
                                                  when > valueDate,
                                                  when,
                                                  false, // always notional
                                                  notl,
                                                  refLevel));

            value += cfln.calcPrice(avg.get(), 0) *
                aop->longShort[i];
        }
        value *= aop->participation;
        value /= discount->pv(when);        // return value at when, not valueDate.
        return value;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//------------------------------------------//    
// for event handling
//------------------------------------------//    
// KnownCashflows::IEventHandler interface
void CEqGainKONote::getEvents(const KnownCashflows* flows,
                IModel* model, 
                const DateTime& eventDate,
                EventResults* events) const {
    static const string method = "CallableEquityKOSwap::getEvents";
//    DateTime hitDate;
//    TDeadType deadType;
//    bool isDead = isAlreadyHit(&hitDate, deadType);
//    CashFlowArraySP cfl = getKnownCashFlow(hitDate, deadType);
    CashFlowArraySP cfl = getKnownCashFlow();
    if (!cfl->empty()) {
        events->addEvent(new KnownCashflows(eventDate, cfl, 
                                            discount->getCcy()));
    }
}


// BarrierBreach::IEventHandler interface
void CEqGainKONote::getEvents(const BarrierBreach* breach,
                IModel* model, 
                const DateTime& eventDate,
                EventResults* events) const {
    static const string method = "CEqGainKONote::getEvents";

    bool isConUp = true;
    bool isConDn = true;
    if (!IntraDayMonitor){
        if (MonitoringDependence == "BOTH")
            isConUp = false;
        else if (MonitoringDependence == "UPPER")
            isConUp = false;
        else if (MonitoringDependence == "LOWER")
            isConDn = false;
    }

    if (UpperBarType != "NA"){
        string BarType = (UpperBarType == "KO") ? BarrierBreach::KNOCK_OUT : BarrierBreach::KNOCK_IN;
        Barrier1F bar1F = Barrier1F(UpperBarrier,true,IntraDayMonitor,UpperBarBreachDate, 
                                    UpperBarBreached, asset,initialSpot);
        bar1F.makeBarrierEvent(events, eventDate, "EGK UpperBarrier",BarType);
    }
    if (LowerBarType != "NA"){
        string BarType = (LowerBarType == "KO") ? BarrierBreach::KNOCK_OUT : BarrierBreach::KNOCK_IN;
        Barrier1F bar1F = Barrier1F(LowerBarrier,false,IntraDayMonitor,LowerBarBreachDate, 
                                    LowerBarBreached, asset,initialSpot);
        bar1F.makeBarrierEvent(events, eventDate, "EGK LowerBarrier",BarType);
    }
}



/********************************************************************************************************************************/
/////////////////////////////////////////////////////////
//           tree1f/fd product
/////////////////////////////////////////////////////////
/** EqGainKONote product payoff for a tree */
class CEqGainKONoteFDProd:  public DblBarrierFDProd
{
public:

    CEqGainKONoteFDProd(const CEqGainKONote* EGK, FDModel* model):DblBarrierFDProd(EGK, model),
                                                   inst(EGK), averageOut()
    {
        numPrices = 5;

        if (inst->hasAvgOut)
        {
            // now if there's just 1 avg out date (which must be maturity)
            // we don't need to worry about pricing fwd starting avg options
            // later on as we just do the intrinsic
            DateTime firstAvgOut = inst->avgOutPerf->avgSamples[0].date;
            if (inst->avgOutPerf->avgSamples.size() == 1)
            {
                avgStart = firstAvgOut;
            }
            else
            {
                // get the end of the day before avg out starts
                avgStart = DateTime(firstAvgOut.getDate()-1, DateTime::END_OF_DAY_TIME);
            }
        }
    }

	/** initialise product specific data */
    virtual void init(CControl* control) const;

    virtual void initProd();

     /** product payoff method at maturity */
	void prod_BWD_T(const TreeSlice & spot,
								   int step,
								   int bot,
                                   int top,
							       int pStart,
							       int pEnd,
                                   const vector< TreeSliceSP > & price);

    /** product payoff method at steps earlier than maturity */
	void prod_BWD(const TreeSlice & spot,
						        int step,
						        int bot,
								int top,
								int pStart,
								int pEnd,
								const vector< TreeSliceSP > & price);

    /** extra output requests */
    virtual bool Positive() { return false; }

    /** make price refinement - combine Libor leg with step down coupons */
    double scalePremium(double P0,
						double P1,
			 			YieldCurveConstSP disc);

   /** extra output requests */
    virtual void recordOutput(Control* control,
							  YieldCurveConstSP disc,
							  Results* results);

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
	{
        return inst->ccyTreatment;
    }

	/** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    virtual void update(int& step,
                        FDProduct::UpdateType type);

protected:

    double              SettleOptionPrice;      //Store the settlement option price
    double              KOSwapPrice;            //Store the conditinal swap price

private:
    const CEqGainKONote* inst;

    vector<double>  stepStrike2;
    vector<double>  StepPart;
    BoolArray       isExeTime;
    BoolArray       isCallStep;

    vector<double>  stepCall_K;

    vector<double>  stepAddStrike1;
    vector<double>  StepAddPart1;
    vector<double>  stepAddStrike2;
    vector<double>  StepAddPart2;
    bool            isCall1;
    bool            isCall2;
    bool            isMax1;
    bool            isMax2;
    bool            isCap;

    vector<double>  stepPVFactToPayDate;        // PV factor for from payment date to a step
    vector<double>  knownEquity;                // Equity Payment already known, on payment date.
    vector<double>  koFixedValue;        // The value after KO at each tree step.
    vector<double>  koFloatValue;        // The value after KO at each tree step.

    string          koStubRule;         //KO StubRule.

    inline void calcBarrier(const double* s,
                                  int step,
                                  int j,
                                  int currIdx,
                                  double koValue,
                                  const vector< double * > & p);

    inline void calcCall(const double* s,
                               int step,
                               int j,
                               int currIdx,
                               double koValue,
                               const vector< double * > & p);

    inline double eqGainIntrinsic(double part,
                                  double partCap,
                                  double partAdd1,
                                  double partAdd2,
                                  double s,
                                  double k1,
                                  double k2,
                                  double k3,
                                  double k4);

    CEqGainKONote::AverageOutSP      averageOut;         // a copy of avgOutPerf for internal use.
    double          SpotRef;
    int             stepNumofFirstAvg;  // number of step, which is the first average out.
    DateTime        avgStart;           // average Start Date.

};

/** create a fd payoff product **/
FDProductSP CEqGainKONote::createProduct(FDModel* model) const
{
    return FDProductSP( new CEqGainKONoteFDProd(this, model) );
}

/** initialise tree1f
    Based on CDblBarrier InitTreee
    Need to add critdate following
        - exerciseSchedule
        - Barrier Monitoring Schedule
        - accrue start date
        - Call Schedule
*/
void CEqGainKONoteFDProd::init(CControl* control) const
{
    static const string method = "CEqGainKONoteFDProd::init()";
    try
    {
        int i;
        if(tree1f)
        {
            tree1f->SetDivAmountTreatment(false);

            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode = numIns;
        }

        DateTimeArray critDates = inst->exerciseSchedule->getDates();

        // remove exercise date from crit date
        critDates.erase(critDates.end()-1);

        // add barrier dates to critical dates)
        DateTime mat = inst->exerciseSchedule->lastDate();
        if (inst->UpperBarType == "KO" || inst->UpperBarType == "KI")
            CDblBarrier::addCritBarDates(inst->UpperBarrier, inst->getValueDate(), mat, tree1f, critDates);
        if (inst->LowerBarType == "KO" || inst->LowerBarType == "KI")
            CDblBarrier::addCritBarDates(inst->LowerBarrier, inst->getValueDate(), mat, tree1f,critDates);

        // Add Call Dates into critDates
        if(inst->isCallable)
        {
            DateTimeArray callDates = inst->CallSchedule->getDates();
            for (i = 0; i<callDates.size(); i++)
                critDates.push_back(callDates[i]);
        }
        // Add Fixed Accrue Start Dates into critDates  + Final Accrue End Dates
        for (i = 0; i<inst->fixedLeg->getSize(); i++)
            critDates.push_back(inst->fixedLeg->AccrueStartDates[i]);
        critDates.push_back(inst->fixedLeg->AccrueEndDates[inst->fixedLeg->getSize()-1]);

        // Add libor Accrue Dates into critDates
        for (i = 0; i<=inst->liborLeg->getSize(); i++)
            critDates.push_back(inst->liborLeg->AccrualDates[i]);

        // Add first Average Out Date into critDates
        if (inst->hasAvgOut)
            critDates.push_back(avgStart);

        DateTimeArray segDates(1);

        // this needs change if fwd start tree has to start today !!!
        if (inst->fwdStarting && inst->startDate>inst->valueDate)
        {
            segDates[0] = getStartDate();  // to do: this needs to be checked !!!
            tree1f->controlSameGridFwdStart(inst->ccyTreatment);
        }
        else
            segDates[0] = inst->valueDate;


        // set up segDensity.  Currently, looks only Up-Out Barrier.
        // add segSensity from 7day before to 7 day after of next barrier date.
        // when 
        int addSeg = 0;
        if (inst->segDensity > 0){
            if(tree1f)
                tree1f->equalTime = true;
            else
                ModelException(method, "non-zero segDensity is allowed only for Tree1f engine.");
            if (inst->UpperBarrier->getInterp() == "N" && inst->UpperBarType == "KO"){
                DateTimeArray upBarDates = inst->UpperBarrier->getDates();
                for (i=0; i<upBarDates.size(); i++) {
                    DateTime preBar = upBarDates[i].rollDate(-7);
                    DateTime aftBar = upBarDates[i].rollDate(7);
                    if (aftBar > mat){
                        break;
                    }else{
                        if (inst->valueDate < preBar) {
                            segDates.push_back(preBar);
                            segDates.push_back(aftBar);
                            break;
                        }else if (inst->valueDate < aftBar){
                            segDates.push_back(aftBar);
                            break;
                        }
                    }
                }
            }
            else 
                ModelException(method, "non-zero segDensity is not allowed except for "
                                       "Up-Out Bermudan case (i.e. Type = KO & interp=N).");
        }

        IntArray density(segDates.size(),1);
        if (inst->segDensity>0){
            switch (segDates.size()){
                case 3:
                    density[1] = inst->segDensity;
                    break;
                case 2:
                    density[0] = inst->segDensity;
                    break;
                default:
                    ModelException(method, "Internal Error.  segDates.size [" 
                                            +Format::toString(inst->segDensity)+
                                           "] should be 2 or 3 for segDenisty>0.");
            }
        }

        // make sure the last date is spot on
        segDates.push_back(mat);
        // check number of segDates & Density
        if (inst->segDensity>0){
            ASSERT(segDates.size() == density.size()+1);
        }else{
            ASSERT(segDates.size() == 2);
            ASSERT(density.size() == 1);
        }
        
        // add critical dates
        model->addCritDates( critDates );

        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

void CEqGainKONoteFDProd::initProd()
{
    static const string method = "CEqGainKONoteFDProd::initProd";
    try{

        DblBarrierFDProd::initProd();

        int callIdx, exeIdx;
        int i =0;
        int numStep = model->getLastStep();
        double fwdAtStart = inst->scalingStrike, scaleFactor = 1.0;
        SpotRef = inst->SpotFixing;

        if (inst->SpotFixing < 0.0)
            SpotRef = inst->initialSpot;

        // get spot at start if needed
        if (inst->fwdStarting && model->getDate(0)>inst->valueDate)
        {
            fwdAtStart = inst->asset->fwdValue(inst->startDate);
            if (!inst->RebateNotScaled)
                scaleFactor = fwdAtStart;
            SpotRef = fwdAtStart;
        }

        // scale adjust of gammaThreshold
        if (tree1f->gammaNodesInterval>1){
            // on tree, price is always for one notional in this product
            tree1f->SetGammaThresholdScaled(SpotRef, 1.0);
        }

        // Tree doesn't call Init when tweaking....
        if (inst->hasAvgOut)
        {
            // now if there's just 1 avg out date (which must be maturity)
            // we don't need to worry about pricing fwd starting avg options
            // later on as we just do the intrinsic
            DateTime firstAvgOut = inst->avgOutPerf->avgSamples[0].date;
            if (inst->avgOutPerf->avgSamples.size() == 1)
            {
                avgStart = firstAvgOut;
            }
            else
            {
                // get the end of the day before avg out starts
                avgStart = DateTime(firstAvgOut.getDate()-1, DateTime::END_OF_DAY_TIME);
            }
        }

        // get participation
        StepPart.resize(numStep+1);
        stepStrike2.resize(numStep+1);
        isExeTime.resize(numStep+1);

        stepCall_K.resize(numStep+1);
        isCallStep.resize(numStep+1);

        StepAddPart1.resize(numStep+1);
        stepAddStrike1.resize(numStep+1);
        StepAddPart2.resize(numStep+1);
        stepAddStrike2.resize(numStep+1);
        stepPVFactToPayDate.resize(numStep+1);
        knownEquity.resize(numStep+1);

        koFixedValue.resize(numStep+1);
        koFloatValue.resize(numStep+1);
        koStubRule  = "B";

        stepNumofFirstAvg = numStep+100;    // If it's not found, never comes in to the tree steps.

        double part=0.0;
        DateTimeArray exeDates = inst->exerciseSchedule->getDates();
        DateTimeArray payDates = inst->ParticipationSchedule->getDates();
        DateTimeArray callDates;
        if (inst->isCallable)
            callDates = inst->CallSchedule->getDates();


        isCap = (inst->CapStrikeArray.size() == inst->exerciseSchedule->length());
        DoubleArray tmpArray;
        if (isCap)
            tmpArray = inst->CapStrikeArray;
        else{
            tmpArray = DoubleArray(0);    // dummy.
            for (i=0;i<exeDates.size();i++)
                tmpArray.push_back(-FP_MIN * 2.0);
        }
        Schedule secondExeLvl(inst->exerciseSchedule->getDates(),
                              tmpArray,
                              inst->exerciseSchedule->getInterp());


        bool isAdd1 = (inst->AddStrikeArray1.size() == inst->exerciseSchedule->length());
        DoubleArray tmpArray1;
        if (isAdd1)
            tmpArray1 = inst->AddStrikeArray1;
        else
            tmpArray1 = tmpArray;    // dummy. not used
        Schedule addExeLvl1(inst->exerciseSchedule->getDates(),
                            tmpArray1,
                            inst->exerciseSchedule->getInterp());
        bool isAdd2 = (inst->AddStrikeArray2.size() == inst->exerciseSchedule->length());
        DoubleArray tmpArray2;
        if (isAdd2)
            tmpArray2 = inst->AddStrikeArray2;
        else
            tmpArray2 = tmpArray;    // dummy. not used
        Schedule addExeLvl2(inst->exerciseSchedule->getDates(),
                            tmpArray2,
                            inst->exerciseSchedule->getInterp());

        //set option type for additional strikes
        if(inst->PayType1 == "N")
            isCall1 = inst->isCall;
        else if (inst->PayType1 == "C")
            isCall1 = true;
        else if (inst->PayType1 == "P")
            isCall1 = false;
        else if (inst->PayType1 == "F")
            isCall1 = true;

        if(inst->PayType2 == "N")
            isCall2 = inst->isCall;
        else if (inst->PayType2 == "C")
            isCall2 = true;
        else if (inst->PayType2 == "P")
            isCall2 = false;
        else if (inst->PayType2 == "F")
            isCall2 = true;

        isMax1 = isMax2 = true;
        if(inst->PayType1 == "F" ||
           (inst->PayType1 == "N" && PayMode == FORWARD))
            isMax1 = false;
        if(inst->PayType2 == "F" ||
           (inst->PayType2 == "N" && PayMode == FORWARD))
            isMax2 = false;

        DoubleArray pastLvl = DoubleArray(0);
        if (inst->samples.size() == 0)
        {
            for (i=0; i<exeDates.size(); i++)
                pastLvl.push_back(-1.0); //negative means not observed yet.
        }
        else
        {
            for (i=0; i<inst->samples.size(); i++)
                pastLvl.push_back(inst->samples[i].amount);
        }

        // find the first exeDates after valuation date.
        for (exeIdx=0; exeIdx< exeDates.size(); exeIdx++)
        {
            if (exeDates[exeIdx] >= model->getDate(0))
                break;
        }

        // find the first callDates after valuation date.
        if (inst->isCallable)
        {
            for (callIdx=0; callIdx < callDates.size(); callIdx++)
            {
                if (callDates[callIdx] >= model->getDate(0))
                    break;
            }
            if (callIdx >= callDates.size())
                callIdx --; //If all calldates are past, need to avoid access memory error.
        }
        // set up average out
        if (inst->hasAvgOut)
        {
            averageOut = CEqGainKONote::AverageOutSP(new CEqGainKONote::AverageOut(inst->avgOutPerf));
        }

        KOStubRuleSP koRule = KOStubRuleSP(new KOStubRule("B",false,"AsSchedule"));    //koStubRule,isAccrueUpToSettle,payTimingRule
        KOFixedLegSP koFix = KOFixedLegSP(new KOFixedLeg(inst->fixedLeg,koRule,inst->instSettle));
        KOLiborLegSP koFlt = KOLiborLegSP(new KOLiborLeg(inst->liborLeg,koRule,inst->instSettle,inst->valueDate,inst->discount));
        double pv_scaleFactor;
        for (i=0; i<=numStep; i++)
        {
            DateTime stepDate = model->getDate(i);

            stepPVFactToPayDate[i] = inst->discount->pv(stepDate, payDates[exeIdx]);
            pv_scaleFactor = stepPVFactToPayDate[i]  * inst->notional/SpotRef;

            // Calculate Known Coupon
            // All coupon is known to be paid on the next day to accrue starting date.
            // Thus, make the flag on the date of accrue starting date.
            // On this date, accrue started coupon is added to price stream.
            // if touch the barrier or called, rebete or callable strike compensate those started coupons.
            // Also, the maturities for each leg is marked as "isAccrue" flag.

            knownEquity[i] = 0.0;		//initilizing.

            koFixedValue[i] = 0.0;  //initilizing.
            koFloatValue[i] = 0.0;  //initilizing.

            // libor leg
            //koFloatValue[i] = inst->liborLeg->getKOPV(inst->valueDate,stepDate,inst->discount.get(),koStubRule);
            //koFixedValue[i] = inst->fixedLeg->getKOPV(stepDate,inst->discount.get(),koStubRule);
            koFloatValue[i] = koFlt->getKOPV(stepDate);
            koFixedValue[i] = koFix->getKOPV(stepDate,inst->discount.get());


            if (exeDates[exeIdx] == stepDate)
            {
                part = inst->ParticipationSchedule->interpolate(payDates[exeIdx]);
                StepPart[i] = part * pv_scaleFactor;
                if (isCap)
                    stepStrike2[i] = fwdAtStart * secondExeLvl.interpolate(stepDate);
                else
                    stepStrike2[i] = -1.0;

                if (isAdd1)
                {
                    stepAddStrike1[i] = fwdAtStart * addExeLvl1.interpolate(stepDate);
                    StepAddPart1[i] = inst->PartStrike1[exeIdx] * pv_scaleFactor;
                }
                else
                    stepAddStrike1[i] = -1.0;

                if (isAdd2)
                {
                    stepAddStrike2[i] = fwdAtStart * addExeLvl2.interpolate(stepDate);
                    StepAddPart2[i] = inst->PartStrike2[exeIdx] * pv_scaleFactor;
                }
                else
                    stepAddStrike2[i] = -1.0;

                isExeTime[i] = true;
                exeIdx++;
            }
            else    //no schedule
            {
                StepPart[i] = 0.0;
                stepStrike2[i] = -1.0;  // not used
                isExeTime[i] = false;

                StepAddPart1[i] = StepAddPart2[i] = 0.0;
                stepAddStrike1[i] = stepAddStrike2[i] = -1.0;

            }
            if (inst->isCallable)
            {
                if( callDates[callIdx] == stepDate)
                {
                    // Call Strike Level is adjusted by DiscFactor to payment date.
                    stepCall_K[i] = pv_scaleFactor * SpotRef * inst->CallSchedule->interpolate(stepDate);
                    isCallStep[i] = true;
                    callIdx++;
                    if (callIdx >= callDates.size())
                        callIdx --;     //avoid Uninitilized Memory Read.
                }
                else
                {
                    stepCall_K[i] = 99999999999999999.9;    //dummy value and not used.
                    isCallStep[i] = false;
                }
            }
            if (inst->hasAvgOut)
            {
                if (stepDate == avgStart)
                    stepNumofFirstAvg = i;
            }

			// equity payment
			// Only if exeDate is past and not yet pay case.
			// Currently, KO or Call never cancel those payment.  If you change this rule, you should look here, again.
			if (exeIdx > 0)
            {
				if (exeDates[exeIdx-1] < inst->valueDate  && inst->valueDate < payDates[exeIdx-1])
                {     // if PayDateTime = ValTime, it drops known value.
                    bool isKnown = false;
                    if (i==numStep)
                    {
                        if (stepDate <= payDates[exeIdx-1])
                            isKnown = true;
                    }
                    else
                    {
                        if (stepDate <= payDates[exeIdx-1] && payDates[exeIdx-1] < model->getDate(i+1))	// in future, paymentDate should be crit dates.
				            isKnown = true;
                    }
                    if (isKnown)
                    {
                        double partMain = inst->ParticipationSchedule->interpolate(payDates[exeIdx-1]);
                        double partCap  = isCap ? partMain : 0.0;
                        double partAdd1 = (isAdd1 ? inst->PartStrike1[exeIdx-1]:0);
                        double partAdd2 = (isAdd2 ? inst->PartStrike2[exeIdx-1]:0);
                        knownEquity[i] = eqGainIntrinsic(partMain, partCap, partAdd1, partAdd2, pastLvl[exeIdx-1],
                                            inst->exerciseSchedule->interpolate(exeDates[exeIdx-1]),
                                            secondExeLvl.interpolate(exeDates[exeIdx-1]),
                                            addExeLvl1.interpolate(exeDates[exeIdx-1]),
                                            addExeLvl1.interpolate(exeDates[exeIdx-1]));
                        knownEquity[i] *= inst->discount->pv(stepDate, payDates[exeIdx-1]) * inst->notional/SpotRef;
                    }
                }
			}

        }
        // When first average out is past, add the expected value of average to the price.
        if (inst->hasAvgOut && stepNumofFirstAvg > numStep)
        {
            if (avgStart < inst->valueDate)
                knownEquity[numStep] += inst->avgValue(averageOut, SpotRef, SpotRef, inst->notional, avgStart);
            else
                throw ModelException(method, "The first date of Average sample cannot found on the Tree Time Steps.");
        }

    }
    catch(exception& e)
    {// continue if indicative vol fail - non BS vol
        throw ModelException(e, method,"Failed at InitProd.");
    }
}

/*------------------------------------------------------------------------------
*   Name         :  CEqKONote::TreePayoff
*
*   Description  :	calc payoff at exercise and provide tree end smoothing.
*					all CalPayoff implementation must provide for T-1 and any step before that.
*					tree assumes price available from T-1.
*                   <Fixed/Libor Coupon Part>
*					Price[][1][] : equity option price
*					Price[][0][] : coupon (digital) option price
*                   <Equity Link Coupon Part>
*					Price[][2][] : euroPrice of coupon (digital) option
*					Price[][3][] : euroPrice of equity option price (with KO)
*                   Price[][4][] : euroPrice of each equity option price (w/o KO).  Need to be refresh at each exercise date.
*
*   Returns      :  none
*------------------------------------------------------------------------------*/
void CEqGainKONoteFDProd::prod_BWD_T(const TreeSlice & spot,
								     int step,
								     int bot,
                                     int top,
							         int pStart,
							         int pEnd,
                                     const vector< TreeSliceSP > & price)
{
    static const string method = "CEqGainKONoteFDProd::prod_BWD_T";

    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    int j;
    int currIdx = tree1f->GetSliceIdx();
    double eqGain = 0.0;

    double koValue = 0.0;
    DateTime stepDate = model->getDate(step);

    // calc terminal price at T first
    int numOfSteps = model->getLastStep();
    double strikeT = stepStrike[numOfSteps]; //  Suppose this is Strike.Exercise.Interpolate(TradeTime[NumOfStep]);
	for (j = bot; j<=top; j++)
	{// floor and ceiling calculated


        eqGain = eqGainIntrinsic(StepPart[numOfSteps],
                                 isCap?StepPart[numOfSteps]:0.0,
                                 StepAddPart1[numOfSteps],
                                 StepAddPart2[numOfSteps],
                                 s[j],
                                 strikeT,
                                 stepStrike2[numOfSteps],
                                 stepAddStrike1[numOfSteps],
                                 stepAddStrike2[numOfSteps]);

        (p[4])[j] = eqGain;

        if (inst->hasAvgOut && stepNumofFirstAvg == step)
        {
            (p[4])[j] += inst->avgValue(averageOut, s[j], SpotRef, stepDate);
        }

        (p[0])[j] = (p[2])[j] = 0.0;        
        (p[1])[j] = (p[3])[j] = knownEquity[step];

        koValue = -koFixedValue[step] - koFloatValue[step];

        // Calculate Barrier Event
        calcBarrier(s, numOfSteps, j, currIdx, koValue , p);

        if (isCallStep[step])
        {//Check Callable
            calcCall(s, numOfSteps, j, currIdx, koValue , p);
        }
    }
}

/** product payoff method at steps earlier than maturity */
void CEqGainKONoteFDProd::prod_BWD(const TreeSlice & spot,
						        int step,
						        int bot,
								int top,
								int pStart,
								int pEnd,
								const vector< TreeSliceSP > & price)
{
    static const string method = "CEqGainKONoteFDProd::prod_BWD";

    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    int currIdx = tree1f->GetSliceIdx();
    double eqGain = 0.0;
    double firstEqGain = 0.0;
    double koValue = 0.0;
    DateTime stepDate = model->getDate(step);

    for (int j = bot; j<=top; j++)
	{
        eqGain = eqGainIntrinsic(StepPart[step],
                                 isCap?StepPart[step]:0.0,
                                 StepAddPart1[step],
                                 StepAddPart2[step],
                                 s[j],
                                 stepStrike[step],
                                 stepStrike2[step],
                                 stepAddStrike1[step],
                                 stepAddStrike2[step]);

        // add Average option
        if (inst->hasAvgOut && stepNumofFirstAvg == step)
        {
            (p[4])[j] += inst->avgValue(averageOut, s[j], SpotRef, stepDate);
        }

        // Added already known payment.  If Callable/KO will affect knownEquity, it should be re-considered.
        (p[3])[j] += knownEquity[step];

        // Prepare Plain Value. (KO is achived later).
        if (isExeTime[step]) // || step == 0)
        {
            (p[3])[j] += (p[4])[j];
            (p[4])[j] = eqGain;
        }

        // if KO or called, following amount is going to be pay ("-" means cancelled).
        // NB : Plain Swap value would be added if SwithPlainSwap = true.
        koValue = -koFixedValue[step] - koFloatValue[step];

        // Calculate Barrier Event
        calcBarrier(s, step, j, currIdx, koValue, p);

        if (isCallStep[step])
        {//Check Callable
            calcCall(s, step, j, currIdx, koValue, p);
        }

        //if(isExeTime[step] && step ==0){
        if(step ==0){
            // To Do :  Look at the past and add only when the KO is not breached in past.
            // Currently, rely on the Breached or not.
            // This should be done after calcCall, because this value should be paid in anycase.
            (p[3])[j] += (p[4])[j];  //eqGain;
            (p[1])[j] += (p[4])[j];  //eqGain;
        }

    }
}

inline void CEqGainKONoteFDProd::calcBarrier(const double* s,
                                                   int step,
                                                   int j,
                                                   int currIdx,
                                                   double koValue,
                                                   const vector< double * > & p)
{
    if (s[j]<LowerBar * (1.0+FP_MIN))
    {
        if (LType == KO)
        {
            (p[0])[j] = (p[2])[j] = 0.0;
            (p[1])[j] = (p[3])[j] = stepPVFactToPayDate[step] * LBarRebate[step] + koValue  ;
        }
        else if (LType == KI || UType !=KI)
        {//Knock-In  (or No KI Barrier)
            (p[0])[j] = (p[2])[j];
            (p[1])[j] = (p[3])[j];
        }
        //else  Just Roll (No Add Cpn/EqGain to (p[0]) and [1])
    }
    else if (s[j] > UpperBar * (1.0-FP_MIN))
    {
        if (UType == KO)
        {
            (p[0])[j] = (p[2])[j] = 0.0;
            (p[1])[j] = (p[3])[j] = stepPVFactToPayDate[step] * UBarRebate[step] + koValue ;
        }
        else if (UType == KI || LType != KI)
        {//Knock-In (Or There is no KI Barrier)
            (p[0])[j] = (p[2])[j];
            (p[1])[j] = (p[3])[j];
        }
        //else  Just Roll (No Add Cpn/EqGain to (p[0]) and [1])
    }
    else //between Low and Up barrier.
    {
        if (UType == KI || LType == KI)
        {//If KI Barrier, eqGain Or Coupon are not added.
    		// Just Roll
        }
        else
        {
    		(p[0])[j] = (p[2])[j];
            (p[1])[j] = (p[3])[j];
        }
    }
}



inline void CEqGainKONoteFDProd::calcCall(const double* s,
                                                int step,
                                                int j,
                                                int currIdx,
                                                double koValue,
                                                const vector< double * > & p)
{
    if ((p[0])[j] + (p[1])[j] > stepCall_K[step] + koValue)
    {//Call : Plain Swap Value is excluded.
     //koValue is almost "Rebate", i.e. the value lost/get when it's called.
        (p[0])[j] = (p[2])[j] = stepCall_K[step] + koValue ;
        (p[1])[j] = (p[3])[j] = 0.0;
    }
    else
    {
        //Non Call
    }
}

inline double CEqGainKONoteFDProd::eqGainIntrinsic(double part,
                                                   double partCap,
                                                   double partAdd1,
                                                   double partAdd2,
                                                   double s,
                                                   double k1,
                                                   double k2,
                                                   double k3,
                                                   double k4)
{
        double value = part*GetIntrinsic(s,k1,inst->isCall, PayMode!=FORWARD);
        if (k2 > -FP_MIN)
            value -= partCap*GetIntrinsic(s,k2,inst->isCall, PayMode!=FORWARD);
        if (k3 > -FP_MIN)
            value += partAdd1*GetIntrinsic(s,k3,isCall1, isMax1);
        if (k4 > -FP_MIN)
            value += partAdd2*GetIntrinsic(s,k4,isCall2, isMax2);
        return value;
}

void CEqGainKONoteFDProd::recordOutput(Control* control,
                                       YieldCurveConstSP disc,
                                       Results* results)
{
    // get prices at t=0
    // save price
    double price = scalePremium(model->getPrice0( *slices[0] ), model->getPrice0( *slices[1] ), disc);
    results->storePrice(price, disc->getCcy());

    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate= inst->exerciseSchedule->lastDate();
        double         indVol;
        // calculate indicative vol
        try{
            if ( matDate.isGreater(inst->valueDate) )
            {

                DateTime imntStartDate = inst->fwdStarting?
                                 inst->startDate: inst->valueDate;

                // get vol request
                CVolRequestConstSP lnVolRequest = GetLNRequest();

                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!vol){
                    throw ModelException("CEqGainKONoteFDProd::recordOutput", 
                                         "No Black Scholes Vol");
                }

                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            else
            {
                indVol = 0.0;
            }
        }
        catch(exception&)
        {// continue if indicative vol fail - non BS vol
           indVol = 0.0;
        }

        // Store Settlement Option Price
        results->storeScalarGreek(SettleOptionPrice, Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName("SettleOption")));
        // Store Conditional Leg
        results->storeScalarGreek(KOSwapPrice, Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName("KOSwap")));

        inst->addOutputRequests(control,
                                results,
                                price,
                                indVol);

        // for the client valuation purpose.  The Libor funding cost (spreads) is tweaked
        OutputRequest* request = control->requestsOutput(OutputRequest::LIBOR_FUNDING_TWEAK);
        if (request && !!inst->liborLeg.get()) {
            CControlSP ctrl(Control::makeFromFlags("", 0.0));   // to avoid iterative calculation?
            IModelSP tree(copy(model));                         // better to use a copy of it?
            inst->liborLeg->tweakSpreads(true);
            InstrumentSP imnt(copy(inst));
            double price_tweak = tree->calcPrice(imnt.get(), ctrl.get());
            results->storeRequestResult(request, price_tweak - price); 
            inst->liborLeg->tweakSpreads(false);
        }
    }
}

void CEqGainKONoteFDProd::update(int& step, FDProduct::UpdateType type)
{
	// we assume just need one und level for spot here
    const TreeSlice & s = payoffIndex->getValue( step );
    int bot, top;
    s.getCalcRange( bot, top );

    const vector< TreeSliceSP > & price = slices;
    int pStart = 0, pEnd = price.size() - 1;

    if (type == FDProduct::BWD_T)
	{
		prod_BWD_T(s,
				   step,
				   bot,
				   top,
                   pStart,
				   pEnd,
				   price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
		{
		    prod_BWD_T(   *insNodes,
					      step,
					      0,
					      tree1f->NumOfInsertNode-1,
                          pStart,
					      pEnd,
					      *insPrices);
        }

    }
    else if(type == FDProduct::BWD)
	{
		prod_BWD( s,
				  step,
				  bot,
				  top,
                  pStart,
				  pEnd,
				  price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
		{
		    prod_BWD(     *insNodes,
					      step,
					      0,
					      tree1f->NumOfInsertNode-1,
                          pStart,
					      pEnd,
					      *insPrices);
        }
    }
}

double CEqGainKONoteFDProd::scalePremium(double P0,
										 double P1,
										 YieldCurveConstSP disc)
{
	static const string method = "CEqGainKONoteFDProd::ScalePremium";

    double fwdStartDF = 1.0;
        if (inst->fwdStarting)
            fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));

    SettleOptionPrice = fwdStartDF * P1;     //Store the settlement option price
    KOSwapPrice = fwdStartDF * P0;           //Store the conditinal swap price

	// combine results if needed
	if (!inst->RemoveSettlementOption) // combine with option
        P0 += P1;

	// substract PlainSwap
    if (inst->SwitchPlainSwap)
    {
        P0 += inst->liborLeg->getPV(inst->valueDate, inst->discount.get())/fwdStartDF;
        P0 += inst->fixedLeg->getPV(inst->valueDate, inst->discount.get())/fwdStartDF;
    }

	return fwdStartDF * P0;
}

DRLIB_END_NAMESPACE
