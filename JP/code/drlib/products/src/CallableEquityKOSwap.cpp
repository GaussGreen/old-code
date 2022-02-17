//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CallableEquityKOSwap.cpp
//
//   Description   Callable Equity KO Swap 
//
//
//   $Log: CallableEquityKOSwap.cpp,v $
//   Revision 1.3  2005/04/26 09:37:09  Kkitazaw
//   Fixing many bugs.
//
//   Revision 1.2  2005/04/22 09:43:53  Kkitazaw
//   Fixed around FixedLeg, so as to use DCC.
//   Multiplying settlePV for koELC to use instrument settlement.
//
//   Revision 1.1  2005/04/20 08:35:54  Kkitazaw
//   General Equity KO Swap with Equity Linked Coupons.
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/AssetUtil.hpp"

#include "edginc/SwapLegIntFace.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"

#include "edginc/Average.hpp"
#include "edginc/CallableEquityKOSwap.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/EqLinkCashFlow.hpp"
#include "edginc/LatticeProdEDR.hpp"

#include "edginc/Tree1fLV.hpp"
#include "edginc/FD1DLV.hpp"

DRLIB_BEGIN_NAMESPACE

// given exercise schedule, isAmer flag and underlying, set exercise flag on each time step.
// Comparing to SetStepExercise, the last step is not always.  The end of sched could be earlier than
// last time step.
void SetStepEvents(vector<bool>& stepExercise, ScheduleConstSP sched, CAssetConstSP  Underlier, DateTimeArray StepDates) 
{
//    static const string method = "Model1F::SetStepEvents";
//    try {
        //vector<bool> stepExercise;
        int NumOfStep = StepDates.size()-1;
        DateTimeArray   schedDates(sched->getDates());
        DateTime        exerStart = schedDates[0];

        int iStep;

        stepExercise.resize(NumOfStep+1);

        // all init to false first
        for (iStep =0; iStep <= NumOfStep; iStep++)
            stepExercise[iStep] = false;

        // find first time point inside the sched
        iStep = 0;
        if ( !(sched->length() == 1) ) {
            while (iStep < NumOfStep && StepDates[iStep].getDate() < exerStart.getDate()) {
                iStep++;
            }
        }

        if (iStep >= NumOfStep)
            return; // There is date in sched.

        // treat bermudan first
        if (sched->getInterp() == Schedule::INTERP_NONE)
        {
            // init Bermudan index
            int iBermudan = 0;
            // locate first exercisable date
            int nBermudan = schedDates.size();
            while (iBermudan < nBermudan-1)
            {
                if (schedDates[iBermudan] < StepDates[0])
                    iBermudan++;
                else
                    break;
            }
            // loop through steps
            for (; iStep <= NumOfStep; iStep++)
            {
                if (iBermudan >= nBermudan)
                    break;
                if (StepDates[iStep].equals(schedDates[iBermudan], false))
                    // && Underlier->getHoliday()->isBusinessDay(StepDates[iStep]) this can rely on external validation ?!
                {
                    stepExercise[iStep] = true;

                    // exaust same day 
                    while (iStep + 1 <NumOfStep && StepDates[iStep + 1].equals(schedDates[iBermudan], false))
                    {
                        iStep++;
                        stepExercise[iStep] = true;
                    }
                    iBermudan++;
                }
            }
        }
        else
        {// now American style
            // all points starting with iStep are potentially exercisable (ie after first sched date)
            // the simple rule is: a timepoint is exercisable if there is a chance to exercise any time after the previous point
            if (AssetUtil::getHoliday(Underlier.get())->isBusinessDay(StepDates[iStep]))
            {
                stepExercise[iStep] = true;
            }
            else
            {
                stepExercise[iStep] = false;
            }

            // loop through steps
            DateTime lastSched = sched->lastDate();
            for (iStep = iStep+1; iStep <= NumOfStep; iStep++)
            {
                if (AssetUtil::hasTradingTime(Underlier,
                                              StepDates[iStep-1],
                                              StepDates[iStep]) 
                  && StepDates[iStep] <= lastSched)
                {
                    stepExercise[iStep] = true;
                }
            }
        }
//    return stepExercise;
//    }
//    catch (exception& e) {
//        throw ModelException(e, method);
//    }
}



// Applying smoothed curve on price array.
// when peakDelta is too small, it returns false with least delta.  
bool SmoothPrice(const double* s, const double peakDelta, 
                 const bool isCheaper,  // -1 : down, 1 : upside
                 const bool isUp,
                 const int iNode,       // first node inside upTrigger or last node inside lower Trigger
                 const int bot, 
                 const int top,
                 const double barLevel,
                 const double hitPrice,
                 const double preHitPrice,
                 double*price,
                 bool &needOWInsPrice,
                 double &lowLevel,
                 double &upLevel,
                 double &lowPrice,
                 double &upPrice, 
                 double &orgDelta)      // (O) return original delta around bar.  When smoothing fails, return least delta.
{
    int i;
    int iEnd = iNode;
    int smoothSide;
    double leastDelta = 1e99;

    // first, define the delction,
    int iDnNode = isUp ? iNode-1 : iNode;
    int iHiNode = isUp ? iNode : iNode +1;

    double startPrice;
    int iShift = 0;

    lowPrice = isUp ? preHitPrice : hitPrice;
    upPrice = isUp ? hitPrice : preHitPrice;
                                                //         ______
    if (lowPrice < upPrice){                    //    ____|         (Original)
        if (isCheaper){                         //          _____    
            smoothSide = 1;                     // => ____|/
            startPrice = lowPrice;
            needOWInsPrice = isUp;
        }
        else {                                  //         ______    
            smoothSide = -1;                    // => ___/|
            startPrice = upPrice;
            needOWInsPrice = !isUp;
        }
    }                                           //    ____
    else {                                      //        |______   (Original)
        if (isCheaper) {                        //    ___
            smoothSide = -1;                    // =>    \|______
            startPrice = upPrice;
            needOWInsPrice = !isUp;
        }
        else {                                  //    ____
            smoothSide = 1;                     // =>     |\_____
            startPrice = lowPrice;
            needOWInsPrice = isUp;
        }
    }


    // for the case no smoothing, store the original status, first.
    lowLevel = smoothSide == 1 ?       barLevel : s[iDnNode];
    upLevel  = smoothSide == 1 ?     s[iHiNode] : barLevel;
    lowPrice = smoothSide == 1 ?     startPrice : price[iDnNode];
    upPrice  = smoothSide == 1 ? price[iHiNode] : startPrice;

    double preDelta = 0.0;
    double afterDelta = 0.0;
    
    // up and smooth up side or down and smooth down side should smooth iNode.
    if (isUp && smoothSide == 1)        
        iShift = -1;
    else if (!isUp && smoothSide == -1)
        iShift = +1;

    // find out the smoothing range, so as to reduce the delta to peakDelta
    // the maximum delta from smoothing function  = 15/16 * dP/dx = 15/16 * dP / dS * 2.0
    // 2.0 comes from x takes from -1 to 1.  (dX/dS = 2.0 / (S(h)-S(l)));
    // ignoreing 15/16, here.
    for (i= iNode + iShift + smoothSide; bot<=i && i<=top ; i += smoothSide){
        afterDelta = (startPrice - price[i]) / (barLevel - s[i]) * (barLevel+s[i]); 
        // keep the least delta for just output purpose.
        if (i == iNode + iShift + smoothSide){
            leastDelta = afterDelta;
            orgDelta = afterDelta;
        } else if (fabs(leastDelta) > fabs(afterDelta))
            leastDelta = afterDelta;
        //
        if ( fabs(afterDelta) < peakDelta ){
            iEnd = i;
            break;
        }
        preDelta = afterDelta;
    }
    if (i < bot || i > top){
        orgDelta = leastDelta;  // return least Delta
        return false; 
    }
    else if (iEnd == iNode+iShift+smoothSide){   

        //no need smoothing as delta is not huge.
        return true;
    }
    
    int iBefore = iEnd - smoothSide;
    double targetDelta = preDelta < 0? -peakDelta : peakDelta;  //peakDelta is absolute, so add sign.
    // here, search for spot level which gives peak Delta by linear interpolation.
    // Already I tried spline, but it didn't work, because delta curve is not smoothed at all.
    double limitLevel = s[iBefore] + (targetDelta - preDelta) * (s[iEnd]-s[iBefore]) / (afterDelta-preDelta);
    double limitPrice = price[iBefore] + (limitLevel - s[iBefore]) * (price[iEnd]-price[iBefore]) / (s[iEnd] - s[iBefore]);
    if (smoothSide ==1){
        lowPrice = startPrice; 
        lowLevel = barLevel;
        upLevel = limitLevel;
        upPrice = limitPrice;
        //upPrice  = price[iEnd];
        //upLevel  = s[iEnd];
        //// Shit Price on Smoothed function ////
        for (i=iNode ; i < iEnd ; i += smoothSide){        
            if (isCheaper)
                price[i] = Maths::min(SmoothValue(upPrice, upLevel, lowPrice, lowLevel, s[i]), price[i]);
            else
                price[i] = Maths::max(SmoothValue(upPrice, upLevel, lowPrice, lowLevel, s[i]), price[i]);
        }    
    }
    else {
        upPrice  = startPrice; 
        upLevel  = barLevel;
        lowPrice = limitPrice;
        lowLevel = limitLevel;
        //lowPrice = price[iEnd];
        //lowLevel = s[iEnd];
        for (i=iNode ; i > iEnd ; i += smoothSide){        
            if (isCheaper)
                price[i] = Maths::min(SmoothValue(upPrice, upLevel, lowPrice, lowLevel, s[i]), price[i]);
            else
                price[i] = Maths::max(SmoothValue(upPrice, upLevel, lowPrice, lowLevel, s[i]), price[i]);
        }    
    }
    return true;
}

double SplinePrice(const double* s, const double* price, const int bot, const int top, 
                   const int jNode, const double interpLvl) {
    int i, tmpBot, tmpTop, tmpSize;
    double y1 = 2e30; // default to natural spline
    double yn = 2e30;
    double result;
    // make a copy of part of stock price array
    //tmpBot = jNode - iWidth < -bot ? -bot : jNode - iWidth;
    //tmpTop = jNode + iWidth >  top ?  top : jNode + iWidth;
    tmpBot = -bot;
    tmpTop =  top;
    tmpSize = tmpTop - tmpBot + 1;
    vector<double> s_copy(tmpSize);
    vector<double> p_copy(tmpSize);
    for (i=0; i<tmpSize; i++) {// copy all nodes
        s_copy[i] = s[tmpBot+i];
        p_copy[i] = price[tmpBot+i];
    }
    vector<double> y2(tmpSize);
    // call cubic spline routines
    spline(&*s_copy.begin()-1, &*p_copy.begin()-1, tmpSize, y1, yn,&*y2.begin()-1);
    splint(&*s_copy.begin()-1, &*p_copy.begin()-1, &*y2.begin()-1, tmpSize, interpLvl, &result);

    return result;
}

void CallableEquityKOSwap::Validate(){
    static const string method = "CallableEquityKOSwap::Validate";

    // fill in any samples
    if (isAvgIn){
        // make SampleList for Product Class.
        DateTimeArray avgDates = CashFlow::dates(avgIn);
        int nSize = avgDates.size();
        DoubleArray pastLevels(nSize,0.0);
        for (int iAvg = 0; iAvg < nSize; iAvg++){
            pastLevels[iAvg] = avgIn[iAvg].amount;
        }
        DoubleArray weights(nSize, 1.0/nSize);
        avgInSL = SampleListSP(new SampleList(avgDates, 
                                             pastLevels,
                                             weights));
    }

    if (hasEqLeg) {
        if (koELC->getObservDates().size()<=0)
            throw ModelException(method, "koELC's observation date is empty.");
    }

    DateTime lastEqStep = DateTime(0,0);
    if (hasUpTrigger) {
        if (!upTrigger.barrier.get())
            throw ModelException(method, "Up Trigger barrier schedule is empty.");
        if (!upTrigger.ecoBarrier.get())
            throw ModelException(method, "Up Trigger barrier schedule is empty.");
        if (upTrigger.getSmoothType() != NO_SMOOTHING && upTrigger.peakDelta < 0)
            throw ModelException(method, "Up Trigger's peakDelta should be positive when it smooth");
        DateTime last = upTrigger.barrier->getDates().back();
        lastEqStep = lastEqStep < last ? last : lastEqStep;
    }

    if (hasLowTrigger) {
        if (!lowTrigger.barrier.get())
            throw ModelException(method, "Low Trigger barrier schedule is empty.");
        if (!lowTrigger.ecoBarrier.get())
            throw ModelException(method, "Low Trigger barrier schedule is empty.");
        if (lowTrigger.getSmoothType() != NO_SMOOTHING && lowTrigger.peakDelta < 0)
            throw ModelException(method, "Low Trigger's peakDelta should be positive when it smooth");
        DateTime last = lowTrigger.barrier->getDates().back();
        lastEqStep = lastEqStep < last ? last : lastEqStep;
    }

    // validation of matDate
    if (hasUpTrigger || hasLowTrigger){
        DateTime endTreeDate = getLastSimDate();
        if (endTreeDate < lastEqStep){
            throw ModelException(method, "maturity of the finalPerf [" + endTreeDate.toString() 
                                        + "] should not be earlier than \n [" 
                                        + lastEqStep.toString() 
                                        + "], which is the last of any earlly termination dates.");
        }
    }

    // set up fwdStarting
    // when isAvgIn, it's fwdStarting, and fwdStarting date is last Avg-In date,
    // as model build tree from last avg-in date, by replacing the fwd level by avg-in fwd.
    if (isAvgIn){
        DateTimeArray avgDates = CashFlow::dates(avgIn);
        int nSize = avgDates.size();

        // fill in any samples
        // make SampleList for Product Class.
        DoubleArray pastLevels(nSize,0.0);
        for (int iAvg = 0; iAvg < nSize; iAvg++){
            pastLevels[iAvg] = avgIn[iAvg].amount;
        }
        DoubleArray weights(nSize, 1.0/nSize);
        avgInSL = SampleListSP(new SampleList(avgDates, 
                                             pastLevels,
                                             weights));
        if (nSize > 1) {
            throw ModelException(method, "AvgIn with multi sample dates is not allowed, yet.  Number of sampling should be one (i.e. FwdStarting).");
            AvgIn avgInGen(avgInSL,asset.get());
            double fwd, bsVar;
            avgInGen.getFwdAndBSVar(valueDate,nSize-1,fwd,bsVar);
            initialSpot = fwd * exp(-0.5*bsVar);
        }


        startDate = getStartDate();
        if (hasEqLeg){
            DateTimeArray obs = koELC->getObservDates();
            if (startDate > obs[0]){
                throw ModelException(method, "Last Average-In dates[" +
                                            startDate.toString() + 
                                            "] should be earlier than 1st koELC observatoin date[" + 
                                            obs[0].toString() + "].");
            }
        }
        if (hasUpTrigger) {
            DateTime fstBarDt = upTrigger.barrier->firstDate();
            if (startDate > fstBarDt){
                throw ModelException(method, "Last Average-In dates[" +
                                            startDate.toString() + 
                                            "] should be earlier than 1st upper barrier date[" + 
                                            fstBarDt.toString() + "].");
            }
        }
        if (hasLowTrigger) {
            DateTime fstBarDt = lowTrigger.barrier->firstDate();
            if (startDate > fstBarDt){
                throw ModelException(method, "Last Average-In dates[" +
                                            startDate.toString() + 
                                            "] should be earlier than 1st lower barrier date[" + 
                                            fstBarDt.toString() + "].");
            }
        }

        fwdStarting  = startDate > valueDate;
        if (!fwdStarting){
            if (!useIniSpotAsRefLvl){
                double sum = 0.0;
                for (int i=0; i<avgDates.size();i++){
                    sum += avgIn[i].amount;
                }
                if (sum < 1.0e-15){
                    throw ModelException(method, "sampled of average-in should not be zero!!");
                }
                this->initialSpot = sum / avgDates.size();
            }
            else{
                if (this->initialSpot < 1.0e-10){
                    throw ModelException(method, "initialSpot should be > 1.0e-10.  Or turn-off useIniSpotAsRefLvl flag.");
                }
            }
        } else {
            // validate so as to avoid too much "average-in" product....
            DateTime firstAvgIn = avgDates[0];
            DateTime lastAvgIn = avgDates[avgDates.size()-1];
            DateTime firstAvgOut = finalPerfMaker->getMatStartDate();
            int daysOfTree = firstAvgOut.daysDiff(lastAvgIn);
            int daysAvgIn  = lastAvgIn.daysDiff(firstAvgIn);
            if (daysOfTree / (daysAvgIn+1) < 20 ){
                throw ModelException(method, "AvgIn period ("+ Format::toString(daysAvgIn) + 
                                             " days) is too long comparing to the option life days (" +
                                             Format::toString(daysOfTree) + 
                                             " days). \n  (AvgIn period) / (option life days) should be greater than 20. ");
            }
        }

    } else {
        fwdStarting = false;
        this->initialSpot = initialSpot;
    }
    // Call Generic1Factor validate()
    validate();
}

//------------------------------------------//
//  Class for Callable                 //
//------------------------------------------//
CallableEquityKOSwap::CallSchedule::CallSchedule():CObject(TYPE),isPuttable(false) {}; 

void CallableEquityKOSwap::CallSchedule::validatePop2Object() {
    static const string method("CallableEquityKOSwap::CallSchedule::validatePop2Object");
    try {
        if (!!callLevels){
            if (callLevels->getInterp() != Schedule::INTERP_NONE)
                throw ModelException(method, "only interp = N (on dates) is allowed.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void CallableEquityKOSwap::CallSchedule::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Callable Deposit Call Schedule");
    REGISTER(CallableEquityKOSwap::CallSchedule, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCallSchedule);
    FIELD(isPuttable, "is puttable");
    FIELD_MAKE_OPTIONAL(isPuttable);
    FIELD(callLevels, "schedule");
    FIELD_MAKE_OPTIONAL(callLevels);
    //FIELD(notification, "notification dates");
    //FIELD_MAKE_OPTIONAL(notification);
    //FIELD(schedule, "schedule");
    //FIELD_MAKE_OPTIONAL(schedule);
}

IObject* CallableEquityKOSwap::CallSchedule::defaultCallSchedule(){
    return new CallableEquityKOSwap::CallSchedule();
}

//------------------------------------------//
//  Class for Barrier (Trigger)        //
//------------------------------------------//
CallableEquityKOSwap::Trigger::Trigger():CObject(TYPE) {
    peakDelta = 0.0;
}; 

CallableEquityKOSwap::Trigger::Trigger(const bool isUp,const ScheduleSP barrier,
                                        const ScheduleSP ecoBarrier, const bool intraDayMonitor,
                                        const bool hasRebate, const ScheduleSP  rebate,
                                        const string smoothingType, const double peakDelta,
                                        DateTimeArray paymentDates):CObject(TYPE),
                                        isUp(isUp),barrier(barrier),ecoBarrier(ecoBarrier),
                                        intraDayMonitor(intraDayMonitor),hasRebate(hasRebate),rebate(rebate),
                                        smoothingType(smoothingType), peakDelta(peakDelta), paymentDates(paymentDates){
    static const string method("CallableEquityKOSwap::Trigger::Trigger()");
    try {
        validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void CallableEquityKOSwap::Trigger::validatePop2Object() {
    static const string method("CallableEquityKOSwap::Trigger::validatePop2Object");
    try {
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// convert input schedule to the array along to the time line
void CallableEquityKOSwap::Trigger::setStepRebatePVed(const DateTimeArray timeLine,       // (I) Date Array to be based on.
                                                  const vector<bool>& isKO,           // (I) is the step KO.
                                                  const YieldCurveConstSP discount,
                                                  const InstrumentSettlementSP  instSettle,
                                                  const CAssetWrapper asset,
                                                  vector<double>& reb)          // (O) the rebete value
                                                  const {
    static const string method("CallableEquityKOSwap::Trigger::setStepRebate");
    try {
        int iStep;
        int numSteps = isKO.size();
        if (numSteps != timeLine.size())
            throw ModelException(method, "isKO array and timeLine should be same size");

        reb.clear();
        reb.resize(numSteps,0.0);
        int jCpn = 0;

        bool overwritePayDates = paymentDates.size() > 0;
        if (hasRebate){
            for (iStep = 0; iStep < numSteps; iStep++) {
                DateTime stepDate = timeLine[iStep];
                if (isKO[iStep]){
                    reb[iStep] = rebate->interpolate(stepDate);
                    if (overwritePayDates){
                        // find the earliest (or equal) date in paymentDates
                        jCpn = Neighbour(stepDate,paymentDates,jCpn,paymentDates.size()-1,1);
                        if (jCpn >= 0)
                            reb[iStep] *= discount->pv(stepDate, paymentDates[jCpn]);
                        else
                            throw ModelException(method, "not found the rebate's payment date for tree step["
                                                        + stepDate.toString() + "]. ");
                    }else
                        reb[iStep] *= instSettle->pvAdjust(stepDate,discount.get(), asset.get());
                }
            }             
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// convert barrier schedule to absolute level along with the time line
void CallableEquityKOSwap::Trigger::setStepBarrier(const DateTimeArray timeLine,       // (I) Date Array to be based on.
                                                  const bool          isUP,           // (I) is Up Barrier or not.  
                                                  const double        scale,          // (I) scale factor (refLevel)
                                                  const vector<bool>& isKO,           // (I) is the step KO.
                                                  vector<double>& bar)          // (O) the tribber level of each step.                                                  
                                                  const {
    static const string method("CallableEquityKOSwap::Trigger::setStepBarrier");
    try {
        int iStep;
        double EXTREME = isUp ? 9999999.99 : 0.0;
        int numSteps = isKO.size();
        if (numSteps != timeLine.size())
            throw ModelException(method, "isKO array and timeLine should be same size");

        bar.resize(numSteps,EXTREME);

        for (iStep = 0; iStep < numSteps; iStep++) {
            DateTime stepDate = timeLine[iStep];
            if (isKO[iStep]){
                bar[iStep] = barrier->interpolate(stepDate)*scale;
            }
        }                      
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// check is already touched before.
bool CallableEquityKOSwap::Trigger::isTriggered(const CashFlowArray samples,
                                                const DateTime valDate,
                                                double spot,
                                                double refLevel,
                                                DateTime* hitDate)  const{
    static const string method("CallableEquityKOSwap::Trigger::isTriggered");
    try {
        bool isHit = false;            
        DateTimeArray smplDates = CashFlow::dates(samples);
        DateTimeArray barDates = barrier->getDates();
        if (smplDates.size()>=barDates.size()){
            vector<int> smplIdx = DateTime::getIndexes(smplDates, barDates);
            if (barrier->getInterp() != Schedule::INTERP_NONE){
                // following is not sufficient.  Need to consider how monitor the continou bar....
                //double barLvl = barrier->interpolate(valDate);
                //hitDate = DateTime(valDate);
                //isHit = isUp ? perf > barLvl : perf < barLvl;
            }
            else {
                DoubleArray barLevel = (ecoBarrier.get() && !ecoBarrier->getDates().empty()) ?
                    ecoBarrier->getValues() : barrier->getValues();
                for (int i=0; i<barDates.size(); i++){
                    if (valDate >= smplDates[smplIdx[i]]){
                        if (isUp && samples[smplIdx[i]].amount/refLevel > barLevel[i] ||
                            !isUp && samples[smplIdx[i]].amount/refLevel < barLevel[i]){
                            isHit = true;
                            *hitDate = barDates[i];
                            break;
                        }
                    }                
                }
            }
        }
        return isHit;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return trigger type by enum format.
// formally, I designed that trigger object has enum member, but
// the enum memeber could not be copied (or initiazlied) at greeks calculation,
// especially for LEGAL_TARM_FV.  So, instead of having member enum variable,
// I made the function to ask the type.
TSmoothingType CallableEquityKOSwap::Trigger::getSmoothType() const
{
    TSmoothingType smoothType;
    if (smoothingType == "DECREASE")
        smoothType = DECREASE;
    else if (smoothingType == "INCREASE")
        smoothType = INCREASE;
    else
        smoothType = NO_SMOOTHING;
    return smoothType;
}

// overwrite the barrier by eco barrier for LegalTerms sens
bool CallableEquityKOSwap::Trigger::useEcoBar(){
    if (ecoBarrier.get() && barrier.get()){
        barrier = ecoBarrier;
        smoothingType = "NO_SMOOTHING";
        return true;
    }
    else
        return false;
}

// for event handling
// this should be called after deteced breach Date using Samples.
void CallableEquityKOSwap::Trigger::makeBarEvnt(DateTime eventDate,
                                                DateTime breachDate,
                                                EventResults* events,
                                                CAssetWrapper asset,
                                                double initialSpot) const{
    bool isBreached = true; // already detected.
    Barrier1F bar1F = Barrier1F(barrier,isUp,intraDayMonitor,breachDate, isBreached, 
                                asset,initialSpot);
    bar1F.makeBarrierEvent(events, eventDate, "CKS Trigger",BarrierBreach::KNOCK_OUT);
}

/* get Rebate
bool CallableEquityKOSwap::Trigger::getRebate(const DateTime hitDate)const{
    static const string method("CallableEquityKOSwap::Trigger::getRebate");
    try {
        if (barrier->getInterp() != Schedule::INTERP_NONE){
            // following is not sufficient.  Need to consider how monitor the continou bar....
            //double barLvl = barrier->interpolate(valDate);
            //hitDate = DateTime(valDate);
            //isHit = isUp ? perf > barLvl : perf < barLvl;
        }
        else {
            DoubleArray barLevel = barrier->getValues();
            for (int i=0; i<barDates.size(); i++){
                if (valDate >= smplDates[smplIdx[i]]){
                    if (isUp && samples[smplIdx[i]].amount/refLevel > barLevel[i] ||
                        !isUp && samples[smplIdx[i]].amount/refLevel < barLevel[i]){
                        isHit = true;
                        hitDate = &barDates[i];
                        break;
                    }
                }                
            }
        }
        return isHit;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
*/

/** Invoked when Class is 'loaded' */
void CallableEquityKOSwap::Trigger::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Callable Deposit Call Schedule");
    REGISTER(CallableEquityKOSwap::Trigger, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultTrigger);
    FIELD(isUp, "is UpOut");
    FIELD_MAKE_OPTIONAL(isUp);
    FIELD(intraDayMonitor, "true for contionus monitoring.");
    FIELD(barrier, "Trigger Schedule");
    FIELD_MAKE_OPTIONAL(barrier);
    FIELD(ecoBarrier, "Trigger Schedule");
    FIELD_MAKE_OPTIONAL(ecoBarrier);
    FIELD(hasRebate, "has Rebate");
    FIELD(rebate,  "Rebate for this Trigger");
    FIELD_MAKE_OPTIONAL(rebate);
    FIELD(smoothingType, "smoothing type of trigger, NONE, UPPER or LOWER");
    FIELD_MAKE_OPTIONAL(smoothingType);
    FIELD(peakDelta, "smoothing size.  model will smooth so as to reduce the peak delta by this number");
    FIELD_MAKE_OPTIONAL(peakDelta);
    FIELD(paymentDates, "payment date overwrite for rebate.  empty => use the instrument's settle");
    FIELD_MAKE_OPTIONAL(paymentDates);
}
    
IObject* CallableEquityKOSwap::Trigger::defaultTrigger(){
    return new CallableEquityKOSwap::Trigger();
}


DateTime CallableEquityKOSwap::getStartDate() const
{//  return startDate, by looking at isAvgIn or not.
    static const string method("CallableEquityKOSwap::getStartDate");
    try {
        DateTime stDate;
        if (isAvgIn){
            DateTimeArray avgDates = CashFlow::dates(avgIn);
            int nSize = avgDates.size();
            //if (nSize > 1) 
                //throw ModelException(method, "AvgIn is not allowed, yet.  Number of sampling should be one to be FwdStarting.");

            stDate = avgDates[nSize-1];
        }
        else 
            stDate = valueDate;

        return stDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return true if it's average-in and need to multiple layer of the price array.
bool CallableEquityKOSwap::isMultiState() const
{
    if (valueDate < getStartDate()){
        if (isAvgIn){
            DateTimeArray avgDates = CashFlow::dates(avgIn);
            int nSize = avgDates.size();
            if (nSize>1)
                return true;
        }
    }
    return false;
}

CallableEquityKOSwap::CallableEquityKOSwap():Generic1Factor(TYPE), 
                        hasEqLeg(false),
                        hasFixedLeg(false), hasFloater(false), 
                        hasUpTrigger(false),hasLowTrigger(false), 
                        hasCallSchedule(false), isAvgIn(false)
#ifdef  TREE_THETA_CAP
                        ,
                        samples(CashFlowArray(0)), numOfAvgInSV(1),
                        isPositiveThetaSmooth(false),thetaSmoothThrehold(0.0)
#endif
                        {

};

CallableEquityKOSwap::CallableEquityKOSwap(
                        //Generic1Fcator
                        const DateTime                valueDate,
                        const bool                    oneContract,
                        const double                  notional,
                        const double                  initialSpot,
                        const string                  ccyTreatment,
                        const InstrumentSettlementSP  instSettle,
                        const InstrumentSettlementSP  premiumSettle,
                        const CAssetWrapper           asset,
                        const YieldCurveWrapper       discount,
                        // for CKS
                        const bool hasEqLeg,
                        EqCpnKOMakerSP koELC,
                        IFinalPerfMakerSP  finalPerfMaker,         
                        const double       redemption,    
                        const bool hasFixedLeg,
                        const bool hasFloater,
                        KOFixedLegMaker  koFixedLeg,
                        KOLiborLegMaker  koFloater,
                        const bool hasUpTrigger,
                        const bool hasLowTrigger,
                        Trigger     upTrigger,
                        Trigger     lowTrigger,
                        const bool hasCallSchedule,
                        CallSchedule callSchedule,    
                        const bool          isAvgIn,
                        CashFlowArray avgIn,
                        const bool useIniSpotAsRefLvl,
                        CashFlowArray samples,
                        int numOfAvgInSV,
                        InstrumentSettlementSP instSettleForFP
#ifdef  TREE_THETA_CAP
                        ,
                        const double thetaSmoothThrehold,
                        const bool isPositiveThetaSmooth
#endif
                        ):Generic1Factor(TYPE), 
                            hasEqLeg(hasEqLeg),koELC(koELC),finalPerfMaker(finalPerfMaker),
                            redemption(redemption),hasFixedLeg(hasFixedLeg),hasFloater(hasFloater), 
                            koFixedLeg(koFixedLeg),koFloater(koFloater),hasUpTrigger(hasUpTrigger),
                            hasLowTrigger(hasLowTrigger),upTrigger(upTrigger),lowTrigger(lowTrigger),
                            hasCallSchedule(hasCallSchedule), callSchedule(callSchedule),
                            isAvgIn(isAvgIn), avgIn(avgIn), useIniSpotAsRefLvl(useIniSpotAsRefLvl),
                            samples(samples), numOfAvgInSV(numOfAvgInSV), instSettleForFP(instSettleForFP)
#ifdef  TREE_THETA_CAP
                            ,
                            thetaSmoothThrehold(thetaSmoothThrehold),isPositiveThetaSmooth(isPositiveThetaSmooth)
#endif
                            {
    this->valueDate=valueDate; 
    this->oneContract=oneContract; 
    this->notional=notional; 
    this->ccyTreatment=ccyTreatment; 
    this->instSettle=instSettle; 
    this->premiumSettle=premiumSettle; 
    this->asset=asset; 
    this->discount=discount; 
    this->initialSpot = initialSpot;

    // fill in any samples
    if (isAvgIn){
        // make SampleList for Product Class.
        DateTimeArray avgDates = CashFlow::dates(avgIn);
        int nSize = avgDates.size();
        DoubleArray pastLevels(nSize,0.0);
        for (int iAvg = 0; iAvg < nSize; iAvg++){
            pastLevels[iAvg] = avgIn[iAvg].amount;
        }
        DoubleArray weights(nSize, 1.0/nSize);
        avgInSL = SampleListSP(new SampleList(avgDates, 
                                             pastLevels,
                                             weights));
    }
};


/** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
bool CallableEquityKOSwap::avoidVegaMatrix(const Model *model){
    if (CTree1fLV::TYPE->isInstance(model)) {
        return true; // do pointwise instead
    }
    return false;
}

// make a LN vol request - not in real use as prefer LV tree
CVolRequestLN* CallableEquityKOSwap::makeVolRequest() const {
    static const string method("CallableEquityKOSwap::makeVolRequest");
    try {
        DateTime imntStartDate = fwdStarting ? startDate : valueDate;
        //DateTime maturity = avgOutPerf->getFirstDate();
        DateTimeArray dts = finalPerfMaker->getMatDates();
        DateTime maturity = dts[0];
        double strike = 1.0;        // always ATM
        if (fwdStarting == false)
            strike *= initialSpot;
        return new LinearStrikeVolRequest(strike,
                                          imntStartDate,
                                          maturity,
                                          fwdStarting);                                          
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns all strikes the CallableEquityKOSwap is sensitive to  */
DoubleArraySP CallableEquityKOSwap::getSensitiveStrikes(OutputNameConstSP outputName,
                                                        const Model*      model){
    static const string method("CallableEquityKOSwap::getSensitiveStrikes");
    try {
        DoubleArraySP sensStrikes(new DoubleArray(0));
        if (avoidVegaMatrix(model)) {
            throw ModelException(method, 
                                 "VEGA_MATRIX is not valid for this "
                                 "instrument");
        }
        CVolRequestConstSP volRequest(makeVolRequest());

        SensitiveStrikeDescriptor sensStrikeDesc;
        sensStrikeDesc.forwardOnly = false;
        asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                   sensStrikeDesc, sensStrikes);

        return sensStrikes;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** when to stop tweaking */
DateTime CallableEquityKOSwap::endDate(const Sensitivity* sensControl) const{
    //DateTime maturity = avgOutPerf->getLastDate();    
    DateTimeArray dts = finalPerfMaker->getMatDates();    
    DateTime maturity = dts[dts.size()-1];    
    DateTime instEnd  = instSettle->settles(maturity, asset.get());
    DateTime assetEnd = asset->settleDate(maturity);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;    
}

// return final simulation date. 
DateTime CallableEquityKOSwap::getLastSimDate() const{
    static const string method = "CallableEquityKOSwap::getLastSimDate";
    try {
        //DateTime avgOutStart = getAvgOutStart();
        DateTime avgOutStart = finalPerfMaker->getMatStartDate();
        DateTime maturity = avgOutStart;
        DateTimeArray avgDates = finalPerfMaker->getMatDates();

        if (hasEqLeg){
            DateTimeArray obs = koELC->getObservDates();
            if(avgOutStart < obs[obs.size()-1])
                maturity = obs[obs.size()-1];            
        }

        return maturity;
    } catch (exception& e) {
        throw ModelException(e, method);
    }    
}

bool CallableEquityKOSwap::sensShift(Theta* shift) {
    static const string method = "CallableEquityKOSwap::sensShift";
    try  {
        DateTime newDate = shift->rollDate(valueDate);

        finalPerfMaker->sensShift(shift, asset.get(), valueDate);

        if (isAvgIn){
            if (!avgInSL.get()){
                // make SampleList for Product Class.
                DateTimeArray avgDates = CashFlow::dates(avgIn);
                int nSize = avgDates.size();
                DoubleArray pastLevels(nSize,0.0);
                for (int iAvg = 0; iAvg < nSize; iAvg++){
                    pastLevels[iAvg] = avgIn[iAvg].amount;
                }
                DoubleArray weights(nSize, 1.0/nSize);
                avgInSL = SampleListSP(new SampleList(avgDates, 
                                                     pastLevels,
                                                     weights));
            }
            else
                avgInSL->roll(shift->getUtil(valueDate), 0, asset.get());
        }

        // and fixings
        if (hasFloater) {
            koFloater.setFixingforThetaShift(valueDate,
                                             discount.get(),
                                             newDate);
        }

        // looking for any samples which is between valueDate and newDates
        for (int i=0; i<samples.size();i++){
            if (valueDate <= samples[i].date && samples[i].date <= newDate){
                if (samples[i].amount < 1.0E-15)
                    samples[i].amount = asset->getThetaSpotOnDate(shift, samples[i].date);
            }
        }
         
        startDate = getStartDate();
        fwdStarting  = startDate > valueDate;
        if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
            startDate.isGreaterOrEqual(valueDate))
        {
            fwdStarting = false;
            initialSpot = asset->getThetaSpotOnDate(shift, startDate);
        }

        // roll the parent (updates value date etc)
        Generic1Factor::sensShift(shift);
        
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
    return true; // our components have theta type sensitivity
}
  
bool CallableEquityKOSwap::sensShift(LegalTerms* shift) {
    // 1. upperEcoBarrier -> UpperBarrier
    // 2. lowerEcoBarrier -> LowerBarrier
    // 3. Anything else?
    bool isSucceed = true;
    if (hasUpTrigger)
        isSucceed = upTrigger.useEcoBar();
    if (hasLowTrigger)
        isSucceed = lowTrigger.useEcoBar();
    if (finalPerfMaker->hasBarrier())
        isSucceed = finalPerfMaker->useEcoBar();
    return isSucceed; 
}

bool CallableEquityKOSwap::priceDeadInstrument(CControl* control, CResults* results) const{
    static const string method = "CallableEquityKOSwap::priceDeadInstrument";
    try  {
        if (fwdStarting)
            return false;

        int i;
        bool isDead = false;
        bool isEnd = false;
        double value = 0.0;
        DateTime hitDate = valueDate;

        TDeadType deadType;
        // check the past samples    
        isDead = isAlreadyHit(&hitDate, deadType);

        // check past maturity or not.
        // Finding out the end of tree.  
        // if there's just 1 avg out date (which must be maturity).
        // Otherwise, set the end of the day before avg out starts
        isEnd = valueDate >= getLastSimDate();

        if (isEnd||isDead){
            // get Results
            CashFlowArraySP koCFL = getKnownCashFlow(hitDate, deadType);
            for (i=0;i<koCFL->size();i++){
                value += (*koCFL)[i].date > valueDate ? 
                         (*koCFL)[i].amount * discount.get()->pv(valueDate, (*koCFL)[i].date) : 0.0;
            }

            value *= notional;
            results->storePrice(value, discount->getCcy());
            if (control && control->isPricing()) {
                recordOutputRequests(control, results);
            }

        }
        return isDead || isEnd;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }            
}

// check already hit the trigger or not.  Return hit Date and also which barrier is hitted.
bool CallableEquityKOSwap::isAlreadyHit(DateTime* hitDate, TDeadType& deadType) const{
    static const string method = "CallableEquityKOSwap::isAlreadyHit";
    try{
        TDeadType test = Alive;
        deadType = test;
        deadType = Alive;
        deadType = Alive;
        bool isDead = false;
        // check the past samples    
        if (samples.size()>0 && valueDate >= samples[0].date){ 
            bool isUpHit = false;
            bool isDnHit = false;
            DateTime upHitDate = DateTime(0,0);
            DateTime dnHitDate = DateTime(0,0);
            if (hasUpTrigger){
                isUpHit = upTrigger.isTriggered(samples, valueDate, asset->getSpot(), initialSpot, &upHitDate);
            }
            if (hasLowTrigger){
                isDnHit = lowTrigger.isTriggered(samples, valueDate, asset->getSpot(), initialSpot, &dnHitDate);
            }
            isDead = isUpHit || isDnHit;
            if(isUpHit){
                // take the earlier hit date.
                if (isUpHit && isDnHit && upHitDate > dnHitDate) {
                    *hitDate = dnHitDate;   
                    deadType = LowHit;
                } else {
                    *hitDate = upHitDate;
                    deadType = UpHit;
                }
            } else if(isDnHit) {
                *hitDate = dnHitDate;
                deadType = LowHit;
            }
            if (hasCallSchedule){
                // to do ?
            }
        }
        return isDead;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }            
}

// Return the critical Date //
DateTimeArray CallableEquityKOSwap::getCritDates() const{
    static const string method = "CallableEquityKOSwap::getCritDates";
    try{
        int i;
        DateTimeArray critDates;
        // 1. EQL
        if (hasEqLeg){
            critDates = koELC->getObservDates();
            // avg out
            critDates.push_back(finalPerfMaker->getMatStartDate());
        }
        // 2. Fixed/Floater
        if (hasFixedLeg){
            KOFixedLegSP koFix(koFixedLeg.makeKOFixedLeg(instSettle));
            DateTimeArray critFix = koFix->getCritDates();
            critDates = DateTime::merge(critDates,critFix);
        }
        if (hasFloater){
            KOLiborLegSP koFlt(koFloater.makeKOLiborLeg(instSettle, valueDate));
            DateTimeArray critFlt = koFlt->getCritDates();
            critDates = DateTime::merge(critDates,critFlt);
        }
        // 3. Trigger Dates
        if (hasUpTrigger){
            const DateTimeArray& barDates = upTrigger.barrier->getDates(); 
            for (i = 0; i<barDates.size(); i++){
                critDates.push_back(barDates[i]);
                critDates.push_back(barDates[i].rollDate(-2));
            }
        }
        if (hasLowTrigger){
            const DateTimeArray& barDates = lowTrigger.barrier->getDates(); 
            for (i = 0; i<barDates.size(); i++){
                critDates.push_back(barDates[i]);
                critDates.push_back(barDates[i].rollDate(-2));
            }
        }
        if (hasCallSchedule){
            const DateTimeArray& callDates = callSchedule.callLevels->getDates(); 
            for (i = 0; i<callDates.size(); i++){
                critDates.push_back(callDates[i]);
                critDates.push_back(callDates[i].rollDate(-2));
                //critDates.push_back(callSchedule.notification[i]);
                //critDates.push_back(callSchedule.notification[i].rollDate(-2));
            }
        }
        // 4. final Perf
        DateTimeArray critFinal;
        if (finalPerfMaker->hasCritDates(critFinal)){
            critDates = DateTime::merge(critDates,critFinal);
        }
        return critDates;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }            
}

CashFlowArraySP CallableEquityKOSwap::getKnownCashFlow(DateTime hitDate, TDeadType deadType) const{
    static const string method = "CallableEquityKOSwap::getKnownCashFlow";
    try{
        bool isDead = deadType != Alive;
        CashFlowArray eqCFL = CashFlowArray(0);
        CashFlowArraySP cflsp = CashFlowArraySP(new CashFlowArray(0));

        // make up known CashFlow
        if (hasFloater){
            KOLiborLegSP koFlt(koFloater.makeKOLiborLeg(instSettle, valueDate));
            cflsp = koFlt->makeKnownCashFlows(isDead? hitDate : valueDate, isDead);
        }

        if (hasEqLeg){
            DateTime vlDt = valueDate;
            if (isDead){
                bool found = koELC->getPayDate(hitDate, instSettle, &vlDt);
                if (!found || vlDt > valueDate)
                    vlDt = valueDate;   
            }
            eqCFL = koELC->getKnownCashFlow(samples, vlDt, initialSpot, instSettle, asset);
            cflsp = CashFlow::merge(CashFlowArraySP(new CashFlowArray(eqCFL)),cflsp);
        }
        
        if (isDead){
            // Already Knocked Out Case
            CashFlow fixCPN;
            if (hasFixedLeg){
                KOFixedLegSP koFix = KOFixedLegSP(koFixedLeg.makeKOFixedLeg(instSettle));
                fixCPN = koFix->getKOCashFlow(hitDate);                                
                cflsp = CashFlow::merge(CashFlowArraySP(new CashFlowArray(1, fixCPN)),cflsp);
            }
            bool hasReb = true;
            CashFlow reb;
            switch (deadType) {
            case UpHit :
                reb = CashFlow(instSettle->settles(hitDate, asset.get()),upTrigger.rebate->interpolate(hitDate));
                break;
            case LowHit :
                reb = CashFlow(instSettle->settles(hitDate, asset.get()),lowTrigger.rebate->interpolate(hitDate));
                break;
            case Called:
                // need to do
                break;
            default:
                hasReb = false;
            }
            if (hasReb){
                // need to make it array to use merge
                cflsp = CashFlow::merge(cflsp, CashFlowArraySP(new CashFlowArray(1, reb)));
            }
        } 
        else {
            // Not terminated yet....

            // check past maturity or not.
            bool isEnd = valueDate >= getLastSimDate();

            if (isEnd) {
                DateTimeArray dt(1);
                dt[0] = finalPerfMaker->getMatStartDate();
                DateTimeArray dts = finalPerfMaker->getMatDates();
                InstrumentSettlementSP myInstSettle;
                if (!instSettleForFP.get())
                    myInstSettle = instSettle;
                else
                    myInstSettle = instSettleForFP;
                DateTime payDateFP = myInstSettle->settles(dts[dts.size()-1], asset.get());
                IFinalPerfSP fp = finalPerfMaker->getFinalPerf(asset.get(), valueDate, initialSpot, discount.getSP(), 
                                                               myInstSettle, premiumSettle, ccyTreatment, dts);    //dts is dummy
                CashFlowArraySP fnlCfl = CashFlowArraySP(new CashFlowArray(0));
                // finalVal should be intrinsic, not pv adjusted.
                double spot = (valueDate == samples[samples.size()-1].date) ?
                                    asset->getSpot() : samples[samples.size()-1].amount;
                double finalVal = fp->getIntrinsic(spot,initialSpot); 
                finalVal += redemption;                
                fnlCfl->push_back(CashFlow(payDateFP, finalVal));
                cflsp = CashFlow::merge(cflsp, fnlCfl);
            }
            else {
                // not terminated, but should be some known cash flow from EQ and Swap.
            }
            // even the fixed Leg are not Known to be settled, it return the all future cashflow
            // and past cash flows, by request from MO.  It's same in EGK or EDS.
            // fixed Leg - maybe would like to implement in SwapLegIntFace.cpp
            CashFlowArraySP fixCFLsp = CashFlowArraySP(new CashFlowArray(0));
            if (hasFixedLeg){
                KOFixedLegSP koFix = KOFixedLegSP(koFixedLeg.makeKOFixedLeg(instSettle));
                DateTime zero(0,0);
                fixCFLsp = koFix->getKnownCashFlows(zero);
                cflsp = CashFlow::merge(cflsp, fixCFLsp);
            }

        }

        return cflsp;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }            
}
/** extra output requests */
void CallableEquityKOSwap::recordOutputRequests(Control* control, 
                          Results* results) const {
    try {
        int i;

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       getLastSimDate(),
                                       valueDate,
                                       asset.get());


        // KNOWN_CASHFLOWS
        if (control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS)) {
            DateTime hitDate;
            TDeadType deadType;
            bool isDead = isAlreadyHit(&hitDate, deadType);
            CashFlowArraySP knownCFL = getKnownCashFlow(hitDate, deadType);

            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    knownCFL.get());   

        }

        // PAYMENT_DATES
        if (control->requestsOutput(OutputRequest::PAYMENT_DATES)) {
            // add all eq coupon paydates.
            DateTimeArray dates = koELC->getObservDates();
            DateTimeArray payDates = koELC->getPaymentDates();
            if (payDates.size()>0){ 
                for (i=0; i<dates.size(); i++){
                    dates[i] = payDates[i];
                }
            }else{// if no paymentDates are given, use instSettle.
                for (i=0; i<dates.size(); i++){
                    dates[i] = instSettle->settles(dates[i], asset.get());
                }
            }
            if (hasFloater){
                DateTimeArray paydts = koFloater.getPayDates();
                for (i=0; i<paydts.size(); i++)
                    dates.push_back(paydts[i]);
            }
            if (hasFixedLeg){
                DateTimeArray paydts = koFixedLeg.getPayDates();
                for (i=0; i<paydts.size(); i++)
                    dates.push_back(paydts[i]);
            }
            // final redemption
            DateTimeArray dts = finalPerfMaker->getMatDates();
            DateTime payDateFP;
            if (!instSettleForFP.get())
                payDateFP = instSettle->settles(dts[dts.size()-1], asset.get());
            else
                payDateFP = instSettleForFP->settles(dts[dts.size()-1], asset.get());
            dates.push_back(payDateFP); 
            OutputRequestUtil::recordPaymentDates(control,results,&dates); 
        }

        // BARRIER_LEVEL
        if (control->requestsOutput(OutputRequest::BARRIER_LEVEL)) {
            if (hasUpTrigger || hasLowTrigger) {
            
                // report barrier levels over a date range
                DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
                BarrierLevelArraySP levels(new BarrierLevelArray(0));
                // use economic barrier (if it exists)
                
                Schedule* s;
                CashFlowArraySP subset;
                if (hasLowTrigger){
                    s = lowTrigger.ecoBarrier.get(); // always use ecoBarrier.  -> validated this already.
                    subset = CashFlowArraySP(s->subset(valueDate, upperDate));
                    if (!subset->empty()) {
                        for (int i = 0; i < subset->size(); i++) {
                            BarrierLevel bl(lowTrigger.isUp,(*subset)[i].date,
                                            (*subset)[i].amount*initialSpot, lowTrigger.intraDayMonitor);
                            levels->push_back(bl);
                        }
                    }
                }
                if (hasUpTrigger){
                    s = upTrigger.ecoBarrier.get(); // always use ecoBarrier.  -> validated this already.
                    subset = CashFlowArraySP(s->subset(valueDate, upperDate));
                    if (!subset->empty()) {
                        for (int i = 0; i < subset->size(); i++) {
                            BarrierLevel bl(upTrigger.isUp,(*subset)[i].date,
                                            (*subset)[i].amount*initialSpot, upTrigger.intraDayMonitor);
                            levels->push_back(bl);
                        }
                    }
                }

                // add record from finalPerfMaker.  Currently only for BARRIER_LEVEL of KnockInPerf...
                bool isCont;                
                ScheduleSP ki_bar = finalPerfMaker->getBarrier(true, isCont);  // use economic barrier (if it exists)
                // report barrier levels over a date range
                if (ki_bar.get()){
                    DateTime barDt = BarrierLevel::barrierWindow(valueDate);
                    subset = CashFlowArraySP(ki_bar->subset(valueDate, barDt));
                    for (int i = 0; i < subset->size(); i++) {
                        BarrierLevel bl(false,(*subset)[i].date,(*subset)[i].amount,isCont);
                        levels->push_back(bl);
                    }
                }
                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       levels.get());
            }
        } 
    }
    catch (exception&) {
        // don't die if any of these fail
    }
}  

/** Get the asset and discount market data */
void CallableEquityKOSwap::GetMarket(const IModel*          model, 
                                     const CMarketDataSP    market) {
    // parent
    Generic1Factor::GetMarket(model, market);
    if(hasFloater){
        koFloater.getMarket(model, market.get(), discount);
    }
}


//------------------------------------------//    
// for event handling
//------------------------------------------//    
// BarrierBreach::IEventHandler interface
void CallableEquityKOSwap::getEvents(const BarrierBreach* breach,
                IModel* model, 
                const DateTime& eventDate,
                EventResults* events) const {
    static const string method = "CallableEquityKOSwap::getEvents";

    DateTime hitDate;
    TDeadType deadType;
    bool isDead = isAlreadyHit(&hitDate, deadType);
    
    double spot = asset->getSpot();        
    double barLvl = 0.0;        
    if (isDead){
        if (deadType == UpHit)
            upTrigger.makeBarEvnt(eventDate, hitDate, events, asset, initialSpot);
        else if (deadType == LowHit)
            lowTrigger.makeBarEvnt(eventDate, hitDate, events, asset, initialSpot);
    }
    else{
        //check the IFinalPerf (KnockIn Perf)
        bool isCont;
        // use barrier, which would be tweaked to be eco barrier.
        ScheduleSP ki_bar = finalPerfMaker->getBarrier(false, isCont);  
        DateTime endDate = eventDate;
        if (finalPerfMaker->isEnd(eventDate, spot/initialSpot, &endDate, barLvl)){
            // return, it's always isUp=false (Down In)
            Barrier1F bar1F = Barrier1F(ki_bar,false,isCont,endDate, true, 
                                        asset,initialSpot);
            bar1F.makeBarrierEvent(events, eventDate, 
                            "CKS KnockInPerf",BarrierBreach::KNOCK_IN);
        }
    }
}

// KnownCashflows::IEventHandler interface
void CallableEquityKOSwap::getEvents(const KnownCashflows* flows,
                IModel* model, 
                const DateTime& eventDate,
                EventResults* events) const {
    static const string method = "CallableEquityKOSwap::getEvents";
    DateTime hitDate;
    TDeadType deadType;
    bool isDead = isAlreadyHit(&hitDate, deadType);
    CashFlowArraySP cfl = getKnownCashFlow(hitDate, deadType);
    if (!cfl->empty()) {
        events->addEvent(new KnownCashflows(eventDate, cfl, 
                                            discount->getCcy()));
    }
}

// -------------------------------------------------------------------------- //
//                              product class
// -------------------------------------------------------------------------- //
CallableEquityKOSwapFDProd::CallableEquityKOSwapFDProd(const CallableEquityKOSwap* cks, FDModel* m) :
    LatticeProdEDRIns(m,0,0),
    cks(cks), eqPerfs(cks->hasEqLeg ? cks->koELC->getObservDates().size() : 0,0.0),
    numState(0)
{
    // initialization
    stepIsUpKO.clear();     
    stepIsDnKO.clear();     
    stepUpRebate.clear();
    stepDnRebate.clear();
    stepCallLevel.clear();
    stepIsCallable.clear();
    stepFixedFee.clear();
    stepFloatFee.clear();
    stepKOCpns.clear();
    stepEqCpnFactor.clear();
    State.clear();
    stepPerfIndex.clear();


//    if( ! tree1f )
//    {
//        throw ModelException( "CallableEquityKOSwapFDProd::CallableEquityKOSwapFDProd", "Instrument of type "+
//                             cks->getClass()->getName() +
//                             " can be priced by CTree1f only" );
//    }

    // first: set discount curve
    if( tree1f )
        tree1f->setDiscountCurve( cks->discount.getSP() );

    // second: create spot payoff
    payoffIndex = model->createProduct( IProdCreatorSP( new
        IndexSpecEQ( cks->asset.getName(), cks->asset, cks->ccyTreatment ) ) );

    //avgOutStart = cks->getAvgOutStart();
    avgOutStart = cks->finalPerfMaker->getMatStartDate();
    matDate = cks->getLastSimDate();

    baseNP = 2 + cks->finalPerfMaker->addNumPrices();        
    numPrices = cks->isAvgIn ? cks->numOfAvgInSV*baseNP : baseNP;  // # prices at root
    numPrices += cks->finalPerfMaker->addNumPrices();

    if (tree1f) 
    {
        if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) 
        {
            tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
        }

        if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION)
        {
            if (cks->hasUpTrigger)
                numIns ++;
            if (cks->hasLowTrigger)
                numIns ++;
            numIns += cks->finalPerfMaker->addInsNode();
        }
        tree1f->NumOfPrice = numPrices;	
        tree1f->NumOfInsertNode = numIns;         
    }else if (fd1dRet){
        if (cks->hasUpTrigger)
            numIns ++;
        if (cks->hasLowTrigger)
            numIns ++;
        numIns += cks->finalPerfMaker->addInsNode();
    }
    
    if (cks->hasEqLeg)
        koELC = cks->koELC;

//    zeroProd = DYNAMIC_POINTER_CAST<ZeroBondProd>(model->createProduct(IProdCreatorSP(new ZeroBond(
//        sched.resetEff, sched.pay, discYCName))));

        
}


// Destruct
CallableEquityKOSwapFDProd::~CallableEquityKOSwapFDProd()
{}

/** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
    isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
void CallableEquityKOSwapFDProd::update(int& step, 
                                        FDProduct::UpdateType type)
{
    const TreeSlice & s = payoffIndex->getValue( step );
    const vector< TreeSliceSP > & price = slices;

    if(type == BWD_NODE_INSERTION)
    {           
        ( this->*calcPayoff )(step, *insNodes, *insPrices);

    }else{
        ( this->*calcPayoff )(step, s, price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)    {
            ( this->*calcPayoff )(step, *insNodes, *insPrices);
        }

        if(type == FDProduct::BWD){
            postCalc(step, s, price);                                               
        }
    }
}

// Set Boolean Array along time steps, which indicate the time step is 
// corresponding to the event date of schedule.
// Jus using SetStepExercise is not good 
// because maturity date automatically true or last event date will be missed when it's not = lastNode. 
void CallableEquityKOSwapFDProd::SetStepEvent(vector<bool>& stepIsEvent, 
                                              ScheduleConstSP sched) 
{
    int lastStep = stepIsEvent.size()-1;
    DateTime lastSched = sched->lastDate();
    SetStepEvents(stepIsEvent,
                  sched, 
                  cks->asset.getSP(),
                  model->getDates());
}

void CallableEquityKOSwapFDProd::SetPerfIndex(const DateTimeArray timeLine, 
                                              const DateTimeArray observDates)
{
    static const string method("CallableEquityKOSwapProd::EqCpnKOLeg::SetPerfIndex");
    try 
    {        
        int numSize = timeLine.size();
        int numObs = observDates.size();
        stepPerfIndex.resize(timeLine.size(), -1);
        // looking for the next observdates in future.
        int j = Neighbour(timeLine[0], observDates, 0, numObs-1,1);
        for (int i=0; i<numSize && 0<=j && j<numObs ; i++)
        {
            if (observDates[j] == timeLine[i])
            {
                stepPerfIndex[i] = j;
                j++;
            }
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** initialise tree1f - allow product customisation */
void CallableEquityKOSwapFDProd::init(CControl* control) const
{
    static const string method = "CallableEquityKOSwapFDProd::init()";
    try
    {
        if (fd1dRet){
            InitFD1D(control);
        }else{
            if (tree1f)
            {
                if ( cks->fwdStarting ) 
                    tree1f->controlSameGridFwdStart(cks->ccyTreatment);			

#ifdef  TREE_THETA_CAP
            // set up theta smoothing
            tree1f->activateThetaSmooth(cks->isPositiveThetaSmooth, cks->thetaSmoothThrehold);
#endif

                tree1f->SetDivAmountTreatment(false);

                tree1f->DEBUG_UseCtrlVar = false;
		    }

            // compile list of critical dates
            DateTimeArray critDates = cks->getCritDates();
            if (cks->isMultiState())
                critDates.push_back(cks->avgInSL->getLastDate());

            // add critical dates
            model->addCritDates( critDates );

            DateTimeArray segDates(2);
            segDates.resize(2);
            
            segDates[0] = cks->fwdStarting ? cks->startDate : cks->valueDate; 
            segDates[1] = matDate;

            IntArray density( 1, 1 );

            // prepare timeline set up
            model->initSegments( segDates, density );
        }
    }
    catch (exception& e) 
	{
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void CallableEquityKOSwapFDProd::initProd() 
{
    static const string method = "CallableEquityKOSwapFDProd::initProd()";
    try 
    { 
        initSlices( numPrices );
        initInsertNode();

        int iStep;
        int k;  // number of Average-In State.
        int numSteps = model->getLastStep();
        iStepAvgInEnd = -1;
        InstrumentSettlementSP myInstSettle;
        if (!cks->instSettleForFP.get())
            myInstSettle = cks->instSettle;
        else
            myInstSettle = cks->instSettleForFP;
        State.clear();  // initialize State
        if (cks->isMultiState())
        {
            numState = cks->numOfAvgInSV;
            State.resize(numState);

            AvgIn avgInGen(cks->avgInSL,cks->asset.get());
            DateTimeArray avgDates = CashFlow::dates(cks->avgIn);
            int nSize = avgDates.size()-1;

            pdf_AvgIn.resize(numState);

            double interval  = 0.05;    // need to review about interval.
            double spot = cks->asset->getSpot();
            double sEnd;
            double varEnd;
            avgInGen.getFwdAndBSVar(cks->valueDate, nSize, sEnd, varEnd);
            double sCenter = sEnd*exp(-0.5*varEnd);  //using peak level, not fwd.
            double mid = double(numState-1)/2.0;
            double sum_pdf=0.0;
            for (k=0;k<numState;k++)
            {                
                double x = interval * (double)(k-mid);
                sEnd = sCenter * exp(x*varEnd);
                State[k].refLevel = avgInGen.getMeanCondAvgIn(cks->valueDate, spot,sEnd,varEnd);
                //pdf_AvgIn[k] = exp(-(x-mean)*(x-mean)/varEnd/2.0)/(2.0*Maths::PI*sqrt(varEnd));
                sum_pdf += pdf_AvgIn[k];

                // set up finalPerf
                State[k].finalPerf = cks->finalPerfMaker->getFinalPerf(cks->asset.get(), cks->valueDate, 
                                                                        State[k].refLevel, cks->discount.getSP(), 
                                                                        myInstSettle, cks->premiumSettle, cks->ccyTreatment, 
                                                                        model->getDates());
            }
            for (k=0;k<numState;k++)
                pdf_AvgIn[k] = 1.0;
            
            iStepAvgInEnd = Neighbour(cks->avgInSL->getLastDate(),model->getDates(),0,numSteps,0);
        }
        else
        {
            numState = 1;
            State.resize(numState);
            State[0].refLevel = cks->fwdStarting ? cks->asset->fwdValue(cks->startDate) : cks->initialSpot;            
            State[0].finalPerf = cks->finalPerfMaker->getFinalPerf(cks->asset.get(), cks->valueDate, 
                                                                    State[0].refLevel, cks->discount.getSP(), 
                                                                    myInstSettle, cks->premiumSettle, cks->ccyTreatment, 
                                                                    model->getDates());
        }

        if (tree1f){
            // scale adjust of gammaThreshold
            if (tree1f->gammaNodesInterval>1){
                // on tree, price is always for one notional in this product
                // need to consider for nuState > 1 case.....
                tree1f->SetGammaThresholdScaled(State[0].refLevel, 1.0);
            }

#ifdef  TREE_THETA_CAP
            // scale adjust for theta smoothing
            if (!Maths::isZero(cks->thetaSmoothThrehold)){
                tree1f->ScaleThetaSmoothThreshold(1.0, 1.0);    // node price are as of bond value.
            }
#endif

        }
        // tree doesnt' call Init, so need to redefine.
        //avgOutStart = cks->getAvgOutStart();
        avgOutStart = cks->finalPerfMaker->getMatStartDate();

        notPaidValues = 0.0;

        for (k=0; k<numState; k++)
        {        
            double UP_EXTREME = State[k].refLevel/FP_MIN;
            State[k].stepUpBarrier.resize(numSteps+1, UP_EXTREME);
            State[k].stepDnBarrier.resize(numSteps+1, -1.0/UP_EXTREME);
        }
        stepUpRebate.resize(numSteps+1, 0.0);
        stepDnRebate.resize(numSteps+1, 0.0);
        stepPerfIndex.resize(numSteps+1, -1);            

        
        // 1.  Set Up Barrier Information     
        int upBarTime, dnBarTime;
        stepIsUpKO.resize(numSteps+1,false);
        stepIsDnKO.resize(numSteps+1,false);                        
        if (cks->hasUpTrigger)
        {
            // check same time schedule
            upBarTime = cks->upTrigger.barrier->lastDate().getTime();
            if (!cks->upTrigger.intraDayMonitor){
                if (! cks->upTrigger.barrier->timesAreAll(upBarTime)){
                    throw ModelException(
                        "All monitoring dates on barrier schedule must have identical "
                        "time of day for intraDayMonitor = false case.");
                }
            }
        
            SetStepEvent(stepIsUpKO, cks->upTrigger.barrier);
            cks->upTrigger.setStepRebatePVed(model->getDates(), 
                                         stepIsUpKO,
                                         cks->discount.getSP(),
                                         cks->instSettle,
                                         cks->asset,
                                         stepUpRebate);
            for (k=0; k<numState; k++)
            {        
                cks->upTrigger.setStepBarrier(model->getDates(), true, 
                                              State[k].refLevel, stepIsUpKO, State[k].stepUpBarrier);
            }
            if (cks->upTrigger.getSmoothType() != NO_SMOOTHING)
            {
                if (!tree1f)
                    throw ModelException(method, "smoothingType is available only for Tree1f.");
                if (tree1f->GetSmoothMethod() != CTree1f::NODE_INSERTION)
                {
                    throw ModelException(method, "Tree smoothing Method should be DEFAULT (NODE_INSERTION), when smoothingType is not NO_SMOOTHING");
                }
                db_UpBarSLow = DoubleArray(0);
                db_UpBarSHi = DoubleArray(0);
                db_UpBarPLow = DoubleArray(0);
                db_UpBarPHi = DoubleArray(0);
                db_UpBarDelta = CashFlowArraySP(new CashFlowArray(0));
            }
        }
        if (cks->hasLowTrigger)
        {
            // check same time schedule
            dnBarTime = cks->lowTrigger.barrier->lastDate().getTime();
            if (!cks->lowTrigger.intraDayMonitor){
                if (! cks->lowTrigger.barrier->timesAreAll(dnBarTime)){
                    throw ModelException(
                        "All monitoring dates on barrier schedule must have identical "
                        "time of day for intraDayMonitor = false case.");
                }
            }

            SetStepEvent(stepIsDnKO,cks->lowTrigger.barrier);
            cks->lowTrigger.setStepRebatePVed(model->getDates(), 
                                         stepIsDnKO, 
                                         cks->discount.getSP(),
                                         cks->instSettle,
                                         cks->asset,
                                         stepDnRebate);
            for (k=0; k<numState; k++)
            {        
                cks->lowTrigger.setStepBarrier(model->getDates(), false, 
                                              State[k].refLevel, stepIsDnKO, State[k].stepDnBarrier);
            }
            if (cks->lowTrigger.getSmoothType() != NO_SMOOTHING)
            {
                if (!tree1f)
                    throw ModelException(method, "smoothingType is available only for Tree1f.");
                if (tree1f->GetSmoothMethod() != CTree1f::NODE_INSERTION)
                {
                    throw ModelException(method, "Tree smoothing Method should be DEFAULT (NODE_INSERTION), when smoothingType is not NO_SMOOTHING");
                }
                db_DnBarSLow = DoubleArray(0);
                db_DnBarSHi = DoubleArray(0);
                db_DnBarPLow = DoubleArray(0);
                db_DnBarPHi = DoubleArray(0);
                db_DnBarDelta = CashFlowArraySP(new CashFlowArray(0));
            }
        }
        upSmoothType = cks->upTrigger.getSmoothType();
        lowSmoothType = cks->lowTrigger.getSmoothType();

        
        // 2.  Callable Information    
        // do some validation, here.  
        // I think I cannot move those to class's validatePop2Object
        // because the class itself could be empty if it's not used.
        if (cks->hasCallSchedule) 
        {
            stepIsCallable.resize(numSteps+1,false);
            stepCallLevel.resize(numSteps+1);
            SetStepEvent(stepIsCallable,cks->callSchedule.callLevels);
            for (iStep = 0; iStep < numSteps; iStep++) 
            {
                DateTime stepDate = model->getDate(iStep);
                if (stepIsCallable[iStep])
                    stepCallLevel[iStep] = cks->callSchedule.callLevels->interpolate(stepDate);                    
            }
            DateTime lastCall = cks->callSchedule.callLevels->lastDate();
            if (lastCall > matDate)
            {
                throw ModelException(method, "The last Call Date ("+lastCall.toString() + 
                                              ") should be equal or earlier than maturity ("+
                                              matDate.toString() + ").");
            }
        }

        // 3.  Set the payment date flag.
        stepFixedFee.clear();
        stepFloatFee.clear();
        stepFixedFee.resize(numSteps+1,0.0);
        stepFloatFee.resize(numSteps+1,0.0);
        int i,j;
        if (cks->hasFixedLeg)
        {
            DateTimeArray paydates = cks->koFixedLeg.getPayDates();
            koFix = KOFixedLegSP(cks->koFixedLeg.makeKOFixedLeg(cks->instSettle));
            CashFlowArrayConstSP fixCFL = koFix->getCashFlows();
            // search the next coming paydates
            j = Neighbour(model->getDate(0), paydates,0,paydates.size()-1, 1);
            for (i=0; i<=numSteps && j<paydates.size() ; i++){
                if (model->getDate(i) == paydates[j])
                {
                    stepFixedFee[i] = (*fixCFL)[j].amount;                    
                    j++;
                }
            }
        }
        if (cks->hasFloater)
        {
            koFlt = KOLiborLegSP(cks->koFloater.makeKOLiborLeg(cks->instSettle,cks->valueDate));
            // prepare the deteministic cash flow array at here
            CashFlowArrayConstSP cfl = koFlt->getCashFlows(cks->valueDate, cks->discount.get());
            DateTimeArray paydates = CashFlow::dates(*cfl);
            // search the next coming paydates
            // Skip at last step.  At last, remainedCpn will handle the remained coupon CalcPayoff. 
            j = Neighbour(model->getDate(0), paydates,0,paydates.size()-1, 1);
            for (i=0; i<numSteps && j<paydates.size(); i++){
                if (model->getDate(i) == (*cfl)[j].date){
                    // when T(0) is payment date, it should drop the value.
                    stepFloatFee[i] = (i==0) ? 0.0 : (*cfl)[j].amount;
                    j++;
                }
            }
        }


        // 4. Set up equity linked coupon.
        stepEqCpnFactor.resize(numSteps+1,0.0);
        iStepAvgOutStart = numSteps;   //if no EqLeg, mat is avgOutStart.
        if (cks->hasEqLeg)
        {
            DateTimeArray obs = cks->koELC->getObservDates();
            eqCpnKOLeg = EqCpnKOLeg(cks->koELC,eqPerfs,cks->valueDate);
            SetPerfIndex(model->getDates(), koELC->getObservDates());
            eqCpnKOLeg.setMonStep(model->getDates());
            // currently looking at only Upper KO.  (To Do) for down KO!!
            eqCpnKOLeg.setKOFactor(model->getDates(),stepIsUpKO,stepEqCpnFactor);
            iStepAvgOutStart = Neighbour(avgOutStart,model->getDates(),0,numSteps,0);                
            eqCpnKOLeg.setDiscFacts(model->getDates(), cks->discount.getSP(), cks->instSettle, cks->asset);
            // collectiong already known payoff, but not paid.
            if (cks->samples.size()>0)
            {
                CashFlowArray cfl = koELC->getKnownCashFlow(cks->samples, cks->valueDate, State[0].refLevel, cks->instSettle, cks->asset, true);
                for (i=0; i<cfl.size(); i++){
                    if (cks->valueDate<cfl[i].date)
                        notPaidValues += cfl[i].amount * cks->discount.get()->pv(cks->valueDate, cfl[i].date);
                }
            }
            else if(cks->valueDate > obs[0])
            {
                throw ModelException(method, "missing the samples for the first equity link coupon observ date[" 
                        + obs[0].toString() + "] ");
            }
        }

        // set up coupons when KO happens on the corresponding date.
        stepKOCpns.clear();
        stepKOCpns.resize(numSteps+1,0.0);  // initialized by 0 value!!
        bool isUseKOCpnEoD = false;
        double koCpnEoD = 0.0;
        DateTimeArray fltAccDts = cks->hasFloater ? koFlt->getAccrueDates() : DateTimeArray(0);
        DateTimeArray fixAccDts = cks->hasFixedLeg ? koFix->getAccrueDates() : DateTimeArray(0);
        int iFltAcc = fltAccDts.size()-1;
        int iFixAcc = fixAccDts.size()-1;
        for (i=numSteps; i>=0; i--){
            DateTime stepDate = model->getDate(i);            
            DateTime koDateFlt = stepDate;
            DateTime koDateFix = stepDate;
        
            // A trick of barrier step in Tree Treatment.
            // For the case that When trigger is "Closing" sampling (EOD), 
            // and Accrue Start is on the barrier date's SOD,
            // i.e. the coupon is already accrue started at EOD.  In such a case the behavior must be....
            // For "B" case, it won't be cancel.
            // For "S" case, a few day's accrue fee would be paid.
            // For "N" case, no need to worry...
            // Headache of tree settting is that the barrier is active for all node in the same date.
            // Thus, the SOD nodes shouldn't have barrier but it actually has.  
            // To avoid this, use the timing of barrier even for SOD or any other timing nodes.  

            // clearly, the double KO doesn't work!!
            if (iFltAcc>=0){
                if (fltAccDts[iFltAcc].equals(stepDate,false)){
                    if (stepIsUpKO[i] && !cks->upTrigger.intraDayMonitor){
                        if (stepDate.getTime() < upBarTime && fltAccDts[iFltAcc].getTime() < upBarTime){
                            koDateFlt = DateTime(stepDate.getDate(),upBarTime);
                        }
                    }
                    if (stepIsDnKO[i] && !cks->lowTrigger.intraDayMonitor){
                        if (stepDate.getTime() < dnBarTime && fltAccDts[iFltAcc].getTime() < dnBarTime){
                            koDateFlt = DateTime(stepDate.getDate(),dnBarTime);
                        }
                    }
                    if (i>0 && !fltAccDts[iFltAcc].equals(model->getDate(i-1),false))
                        iFltAcc--;
                }
            }
            if (iFixAcc>=0){
                if (fixAccDts[iFixAcc].equals(stepDate,false)){
                    if (stepIsUpKO[i] && !cks->upTrigger.intraDayMonitor){
                        if (stepDate.getTime() < upBarTime && fixAccDts[iFixAcc].getTime() < upBarTime){
                            koDateFix = DateTime(stepDate.getDate(),upBarTime);
                        }
                    }
                    if (stepIsDnKO[i] && !cks->lowTrigger.intraDayMonitor){
                        if (stepDate.getTime() < dnBarTime && fixAccDts[iFixAcc].getTime() < dnBarTime){
                            koDateFix = DateTime(stepDate.getDate(),dnBarTime);
                        }
                    }
                    if (i>0 && !fixAccDts[iFixAcc].equals(model->getDate(i-1),false))
                        iFixAcc--;
                }
            }

            stepKOCpns[i]  = cks->hasFloater  ? koFlt->getKOValue(koDateFlt, cks->discount.get()) : 0.0;
            stepKOCpns[i] += cks->hasFixedLeg ? koFix->getKOValue(koDateFix, cks->discount.get()) : 0.0; 
        }


        
        //extra init prod data
        if (fd1dRet)
        {
            //to change
            InitProdFD1D();
        }

        // if tree1f engine then use original way (for performance) otherwise slice operators
        bool doNotUseSliceOper = (fd1dRet || tree1f);
        if( doNotUseSliceOper )
        {
            calcPayoff = &CallableEquityKOSwapFDProd::CalcPayoff;
        }
        else
        {
            if (upSmoothType != NO_SMOOTHING || lowSmoothType != NO_SMOOTHING)
                throw ModelException(method, "smoothing is not available with _oper");
            calcPayoff = &CallableEquityKOSwapFDProd::CalcPayoff_oper;
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

//output results 
void CallableEquityKOSwapFDProd::recordOutput(Control* control, 
                                              YieldCurveConstSP disc, Results* results)
{
    // get prices at t=0
    double price0 = model->getPrice0( *slices[0] );

    // scaling and adding know values
    double fwdStartDF = 1.0;
    if (cks->fwdStarting && model->getDate(0).getDate() > cks->valueDate.getDate())
        fwdStartDF = disc->pv(cks->valueDate, model->getDate(0));
 
    double price = cks->notional*fwdStartDF*(price0+notPaidValues);

    // save price
    results->storePrice(price, disc->getCcy());
    cks->recordOutputRequests(control, results); 

    // for the client valuation purpose.  The Libor funding cost (spreads) is tweaked
    OutputRequest* request = control->requestsOutput(OutputRequest::LIBOR_FUNDING_TWEAK);
    if (request && cks->hasFloater) {
        CControlSP ctrl(Control::makeFromFlags("", 0.0));   // to avoid iterative calculation?
        IModelSP tree(copy(model));                         // better to use a copy of it?
        koFlt->tweakSpreads(true);
        InstrumentSP imnt(copy(cks));
        double price_tweak = tree->calcPrice(imnt.get(), ctrl.get());
        results->storeRequestResult(request, price_tweak - price); 
        koFlt->tweakSpreads(false);
    }

    // ---debug output---
    if (control->isPricing()){
        int i;
        int numSmth = db_DnBarSLow.size();
        CashFlowArraySP tmp(new CashFlowArray(0));            
        if (numSmth > 0){
            CDoubleMatrixSP smthMatrix(new DoubleMatrix(numSmth, 4));
            for (i=0; i<numSmth; i++){
                (*smthMatrix)[numSmth-1-i][0] = db_DnBarSLow[i];
                (*smthMatrix)[numSmth-1-i][1] = db_DnBarSHi[i];
                (*smthMatrix)[numSmth-1-i][2] = db_DnBarPLow[i];
                (*smthMatrix)[numSmth-1-i][3] = db_DnBarPHi[i];
                tmp->push_back((*db_DnBarDelta)[numSmth-1-i]);
            }
            results->storeGreek(smthMatrix, Results::DEBUG_PACKET, OutputNameSP(new OutputName("LowSmoothLevels")));
            results->storeGreek(db_DnBarDelta, Results::DEBUG_PACKET, OutputNameSP(new OutputName("LowBarDelta")));
        }
        numSmth = db_UpBarSLow.size();
        if (numSmth > 0){
            CDoubleMatrixSP smthMatrix(new DoubleMatrix(numSmth, 4));
            for (i=0; i<numSmth; i++){
                (*smthMatrix)[numSmth-1-i][0] = db_UpBarSLow[i];
                (*smthMatrix)[numSmth-1-i][1] = db_UpBarSHi[i];
                (*smthMatrix)[numSmth-1-i][2] = db_UpBarPLow[i];
                (*smthMatrix)[numSmth-1-i][3] = db_UpBarPHi[i];
                tmp->push_back((*db_UpBarDelta)[numSmth-1-i]);
            }
            results->storeGreek(smthMatrix, Results::DEBUG_PACKET, OutputNameSP(new OutputName("UpSmoothLevels")));
            results->storeGreek(tmp, Results::DEBUG_PACKET, OutputNameSP(new OutputName("UpBarDelta")));
        }
    }
    // ------------------
}

/** called before PayoffAtMat, PayoffBeforeMat and tree roll() */
// calculate barrier and place barrier at inserted node if needed
void CallableEquityKOSwapFDProd::preCalc(int step) 
{
    static const string method("CallableEquityKOSwapFDProd::preCalc");

    try 
    {
        int k;            
        bool doubleBar = cks->hasLowTrigger && cks->hasUpTrigger;            
        if (tree1f){
            int idx = tree1f->getSliceIndex(step);

            for (k=0; k < numState; k++)
            {
                State[k].adjUpBarrier = State[k].stepUpBarrier[step];
                State[k].adjDnBarrier = State[k].stepDnBarrier[step];

                int insNodeIdx = doubleBar ? 2*k : k;
                State[k].finalPerf->upDatePerf(tree1f, step);
                double insLvl;
                if (State[k].finalPerf->insNodeLevel(step, insLvl))
                    tree1f->SetInsertNode(idx, doubleBar ? insNodeIdx+2 : insNodeIdx+1, insLvl, 0); // 0 mean KI

                // taking care of discrete monitoring, treat as daily monitor
                if (cks->hasUpTrigger && !cks->upTrigger.intraDayMonitor)
                {
                    vector<double> vol;            
                    tree1f->GetStepVol(step, vol, &State[k].adjUpBarrier, 0, 0);                                
                    Barrier::BarrierAdjustment(vol[0], true, State[k].adjUpBarrier);
                }
                if (cks->hasLowTrigger && !cks->lowTrigger.intraDayMonitor)
                {
                    vector<double> vol;            
                    tree1f->GetStepVol(step, vol, &State[k].adjDnBarrier, 0, 0);                                
                    Barrier::BarrierAdjustment(vol[0], false, State[k].adjDnBarrier);
                }
                if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
                { 
                    // apparently last param = 1 means KO
                    if (cks->hasLowTrigger)
                        tree1f->SetInsertNode(idx, insNodeIdx, State[k].adjDnBarrier, 1); 
                    if (cks->hasUpTrigger)
                        tree1f->SetInsertNode(idx, doubleBar ? insNodeIdx+1 : insNodeIdx, State[k].adjUpBarrier, 1); 
                }

                // need to this this, here.
                State[k].foundUpLvl = false;
                State[k].foundDnLvl = false;
            }
        }else if(fd1dRet)
        {
            int idx = fd1dRet->getSliceIndex(step);

            for (k=0; k < numState; k++)
            {
                State[k].adjUpBarrier = State[k].stepUpBarrier[step];
                State[k].adjDnBarrier = State[k].stepDnBarrier[step];

                if (cks->hasUpTrigger && !cks->upTrigger.intraDayMonitor)
                {
                    vector<double> vol;
                    fd1dRet->GetStepVol(step, vol, &State[k].adjUpBarrier, 0, 0); // get vol at barrier
                    Barrier::BarrierAdjustment(vol[0], true, State[k].adjUpBarrier);
                }
                if (cks->hasLowTrigger && !cks->lowTrigger.intraDayMonitor)
                {
                    vector<double> vol;
                    fd1dRet->GetStepVol(step, vol, &State[k].adjDnBarrier, 0, 0); // get vol at barrier
                    Barrier::BarrierAdjustment(vol[0], false, State[k].adjDnBarrier);
                }

                // currently, only KO available.  Need to work with KI.
                if (cks->hasUpTrigger) 
                {
                    fd1dRet->getUpBarrier( 0, State[k].adjUpBarrier);
                } 

                if (cks->hasLowTrigger) 
                {
                    fd1dRet->getDownBarrier( 0, State[k].adjDnBarrier);
                }
            }
        }else{
            const FD1DLV * modelLV = dynamic_cast< const FD1DLV * >( model );
            const TreeSlice & s = payoffIndex->getValue( step );
            // adjust barrier if needed
            vector< double > vol;
            for (k=0; k < numState; k++)
            {
                int insNodeIdx = doubleBar ? 2*k : k;
                State[k].adjUpBarrier = State[k].stepUpBarrier[step];
                State[k].adjDnBarrier = State[k].stepDnBarrier[step];

                if (modelLV){
                    if (!cks->upTrigger.intraDayMonitor ){
                        modelLV->GetStepVol( step, vol, &State[k].adjUpBarrier, 0, 0 ); // get vol at barrier
                        Barrier::BarrierAdjustment( vol[0], true, State[k].adjUpBarrier );
                    }
                    if (!cks->lowTrigger.intraDayMonitor ){
                        modelLV->GetStepVol( step, vol, &State[k].adjDnBarrier, 0, 0 ); // get vol at barrier
                        Barrier::BarrierAdjustment( vol[0], false, State[k].adjDnBarrier );
                    }
                }
                
                if (cks->hasUpTrigger) {
                    double koValue = stepKOCpns[step] + stepUpRebate[step];
                                    // TO DO get eqCpn Value!
                                    // + stepEqCpnFactor[step] * ( *slices[insNodeIdx+1] );  
                    model->addCriticalLevel(
                        step, s, State[k].adjUpBarrier, *slices[ insNodeIdx ], 
                        FDModel::LEVEL_BARRIER_UP, koValue );
                }
                if (cks->hasLowTrigger){
                    double koValue = stepKOCpns[step] + stepDnRebate[step];
                                    // TO DO get eqCpn Value!
                                    //+ stepEqCpnFactor[step] * ( *slices[insNodeIdx + 1 + doubleBar?1:0 ] ); 
                    model->addCriticalLevel(
                        step, s, State[k].adjDnBarrier, *slices[ insNodeIdx + doubleBar?1:0 ], 
                        FDModel::LEVEL_BARRIER_DOWN, koValue );
                }
            }
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

///** initialise tree1f - allow product customisation */
//copy it from inittree
//to rewrite after new interface
void CallableEquityKOSwapFDProd::InitFD1D(Control*    control) const {

    //set up segment info.
    //set bar dates as seg dates
      IntArray isAddedSeg ;
    DateTimeArray segDates;
    DateTime t0;
    TimeMetricConstSP metric = model->getTimeMetric();
    const DateTime& matDate= cks->getLastSimDate();

    if (cks->fwdStarting && cks->startDate>cks->valueDate)
    {
        t0 = cks->startDate;  // to do: this needs to be checked !!!
    }
    else{
        t0 = cks->valueDate;
    }

    ScheduleSP tmpUpBar, tmpDnBar;
    string upBarType, dnBarType;
    if (cks->hasUpTrigger) {
        tmpUpBar = cks->upTrigger.barrier;
        upBarType = "KO";
    }
    else{
        tmpUpBar = ScheduleSP(   );
        upBarType = "NA";
    }
    if (cks->hasLowTrigger) {
        tmpDnBar = cks->lowTrigger.barrier;
        dnBarType = "KO";
    }
    else{
        tmpDnBar = ScheduleSP(   );
        dnBarType = "NA";
    }


    FDUtils::SetSegDates(t0,
            matDate,
            metric,
            tmpUpBar, 
            tmpDnBar, 
            upBarType,
            dnBarType,
            segDates, &isAddedSeg);

    //end with seg


    // all exercise dates are copied to critical dates
    
    DateTimeArray critDates = cks->getCritDates();
    // remove exercise date from crit date
    // critDates.erase(critDates.end()-1);
    // add div event dates if needed
    EventAssetMove divEvent;

    //treeif, FD1F and FD1D model names are diff. and set up take diff param.
    //review, chg fd's name based on treeif's structure,aaaaaa,
    //FD1D set up

    //using var grid
    DoubleArray outCritSpacePts;

    outCritSpacePts.resize(1);
    //get roughly the barriers level (average) from products 
//    if (inNeedVariableGrid == true){
        FDUtils::setCritSpacePtsAll(
            t0,
            matDate,
            tmpUpBar, 
            tmpDnBar, 
            upBarType,
            dnBarType,
            outCritSpacePts);
//    }

    fd1dRet->setNumOfInsertNode( numIns );
    fd1dRet->setNbOfProd(1);
    fd1dRet->setMaxNumOfValue( numPrices );

    // add critical dates
    model->addCritDates( critDates );

    IntArray density(1, 1);
      
    // prepare model set up
    //dblBarrier need extra data, call fd1dRet->addInitData to fill iniDataExtra
    //otherwise, just call model->addInitData to fill initData.
    fd1dRet->initSegments(segDates, density, outCritSpacePts, &isAddedSeg);
}

void CallableEquityKOSwapFDProd::InitProdFD1D(){
    //InitProd();

    /* not yet KI
    if (UType == KI || LType == KI){

        if (LType != NA){ //Ko, KI
            fd1dRet->barrier[0]->needSpecialFDDownBar = true; //KIKO
            if (LType == KO){
                fd1dRet->barrier[1]->needSpecialFDDownBar = true; //KIKO
            }
        }

        if (UType != NA){//KO, KI
            fd1dRet->barrier[0]->needSpecialFDUpBar = true; //KIKO
            if(UType == KO){
                fd1dRet->barrier[1]->needSpecialFDUpBar = true; //KIKO
            }
        }
    }else{        */
        if (cks->hasUpTrigger){
            fd1dRet->barrier[0]->needSpecialFDUpBar = true; 
        }

        if(cks->hasLowTrigger){
            fd1dRet->barrier[0]->needSpecialFDDownBar = true; 
        }
    //}
}


// 1.  Add Equity Cpn and prepare the protected amount.
// 2.  Overwrite the prices by the KO values
// 3.  Overwrite the prices by eary exerercise.
// 4.  Adding fixed/floater coupon.  
// NB :  All cash flow are designed "InArears" = Y.  The level fixed at beggining, and 
//       paid later.  Equity Coupon is oppsite.  
//       Also, Accrue Date are "Include" start and "Exclude" end.
//       While Equity's Accrue Date are "Include" end (observ date) and "Exclude" start (previous observ date).

void CallableEquityKOSwapFDProd::CalcPayoff(int step, 
                                            const TreeSlice & spot, 
                                            const vector< TreeSliceSP > & price)
{
    //int step, int idx, double upKOCM, double dnKOCM, double koCpns, 
    //DateTime stepDate, double s, double& value)
    int bot, top;
    spot.getCalcRange( bot, top );    
    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );
    int pEnd = price.size()-1;

    bool isAtMat = (step == model->getLastStep());
        
    DateTime stepDate = model->getDate(step);

    int numFP = cks->finalPerfMaker->addNumPrices();
    
    double remainedCpn = 0.0;
    if (isAtMat)
    {
        // notional remaining at stake (if breached dates)
        double settlePV;
        if (!cks->instSettleForFP.get()){
            settlePV = cks->instSettle->pvAdjust(stepDate,
                                            cks->discount.get(), 
                                            cks->asset.get());
        }else{
            settlePV = cks->instSettleForFP->pvAdjust(stepDate,
                                                      cks->discount.get(), 
                                                      cks->asset.get());
        }    
        double redemption = cks->redemption * settlePV;

        double fixedRemain = cks->hasFixedLeg ? koFix->getPV(stepDate, cks->discount.get()) : 0.0;
        double floatRemain = 0.0;
        if (cks->hasFloater) 
        {
            floatRemain = koFlt->getPV(stepDate, cks->discount.get());
        }                                                  
        remainedCpn = redemption + fixedRemain + floatRemain;
    }


    int idx = stepPerfIndex[step];  // idx<0 means it's not payoff step.

    double koCpns = stepKOCpns[step];

    int k;
    int jNode;
    for (k=0; k<numState; k++){            
        if (isAtMat){
            for (jNode = bot; jNode<=top; jNode++)// initialize
                    (p[pEnd])[jNode] = 0.0;
        }
        // barriers for continuous/daily monitoring and rebates
        double upKOCM = State[k].adjUpBarrier;
        double dnKOCM = State[k].adjDnBarrier;

        int iPrice = baseNP*k; 
        
        // Store eq Link Coupon price into price[k+1] when there is payoff (idx>=0)
        if (!eqCpnKOLeg.isNull && idx >= 0){
            for (jNode=bot; jNode<=top; jNode++){                
                eqPerfs[idx] = s[jNode]/State[k].refLevel;
                (p[iPrice+1])[jNode] = eqCpnKOLeg.getValue(step);
            }
        }
        // pay-off
        if (!State[k].foundDnLvl) 
            State[k].iLstDnNode.resize(2, bot-1);
        else
            State[k].iLstDnNode.resize(0);
        if (!State[k].foundUpLvl) 
            State[k].iFstUpNode.resize(2,  top+1);
        else
            State[k].iFstUpNode.resize(0);

        vector<double> noKOPrice;
        vector<double> s_noKO;
        double koEqnCpn;
        for (jNode=bot; jNode<=top; jNode++) {
            //add the avgOut premium.
            if (isAtMat)
                (p[iPrice])[jNode] = remainedCpn;                
            if (step == iStepAvgOutStart){
                // calculate the finalPerf's value on tree.
                (p[iPrice])[jNode] += State[k].finalPerf->getValue(0, step, s[jNode], State[k].refLevel);
                for (int iFP = 1; iFP <= numFP; iFP++){
                    (p[iPrice+1+iFP])[jNode] = State[k].finalPerf->getValue(iFP, step, s[jNode], State[k].refLevel);
                }
            } 
            // add finalPerf's results on instrument prices.
            State[k].finalPerf->upDateValues(step, s[jNode], (p[iPrice+1+1])[jNode], (p[iPrice+1+numFP])[jNode]);
                
        
            // 1. Add Equity Payment on observ date.  
            if (!eqCpnKOLeg.isNull && idx >= 0)
                (p[iPrice])[jNode] += (p[iPrice+1])[jNode];

            // The equity payments, which is not cancelled, are stored to koEqnCpn 
            // so as to avoid early termination.
            koEqnCpn = stepEqCpnFactor[step] * (p[iPrice+1])[jNode];

            // 2. KO values
            if (iPrice ==0 && (upSmoothType != NO_SMOOTHING 
                            ||lowSmoothType != NO_SMOOTHING) ){//store no KO prices
                s_noKO.push_back(s[jNode]);
                noKOPrice.push_back((p[iPrice])[jNode] + stepFixedFee[step] + stepFloatFee[step]);
            }
            if (s[jNode] > upKOCM*(1.0-FP_MIN) && cks->hasUpTrigger){
                if (!State[k].foundUpLvl){   // store the price before HitBar
                    State[k].iFstUpNode[iPrice] = jNode;
                    State[k].foundUpLvl = true;
                }
                (p[iPrice])[jNode] =  koCpns + stepUpRebate[step] + koEqnCpn;
                for (int iFP = 1; iFP <= numFP; iFP++)
                    (p[iPrice+1+iFP])[jNode] = 0.0; // KO the final Payoff
            }
            if (s[jNode] < dnKOCM*(1.0+FP_MIN) && cks->hasLowTrigger) {
                if(!State[k].foundDnLvl && jNode < top && s[jNode+1] > dnKOCM*(1-FP_MIN)){// store the price before HitBar
                    State[k].iLstDnNode[iPrice] = jNode;
                    State[k].foundDnLvl = true;
                }
                (p[iPrice])[jNode] =  koCpns + stepDnRebate[step] + koEqnCpn;
                for (int iFP = 1; iFP <= numFP; iFP++)
                    (p[iPrice+1+iFP])[jNode] = 0.0; // KO the final Payoff
            } 

            // 3. early exercise
            if (cks->hasCallSchedule){
                if (stepIsCallable[step]){
                    double intrinsic = koCpns + koEqnCpn + stepCallLevel[step];
                    if ((!cks->callSchedule.isPuttable && intrinsic < (p[iPrice])[jNode]) ||
                        ( cks->callSchedule.isPuttable && intrinsic > (p[iPrice])[jNode]))
                        (p[iPrice])[jNode] = intrinsic;                        
                }
            }

            // 4. add all cash flow, which are not affected by early termination.
            (p[iPrice])[jNode] += stepFixedFee[step] + stepFloatFee[step];

            // 5. add FinalPerf
            if (step==0 && numFP>0)
                (p[iPrice])[jNode] +=  (p[iPrice+1+1])[jNode];

        }

        if (State[k].iFstUpNode.size()>0 && State[k].foundUpLvl && upSmoothType != NO_SMOOTHING){//store no KO prices
            int tmpSize = s_noKO.size();
            vector<double> y2(tmpSize);
            spline(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, tmpSize, 2e30, 2e30,&*y2.begin()-1);
            splint(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, &*y2.begin()-1, tmpSize, upKOCM, &State[k].prePriceAtUpBar);
                //prePriceAtUpBar = SplinePrice(&s_noKO[0], &noKOPrice[0], 0, s_noKO.size(), iFstUpNode[0]-bot, upKOCM);
        }
        if (State[k].iLstDnNode.size()>0 && State[k].foundDnLvl && lowSmoothType != NO_SMOOTHING){//store no KO prices
            int tmpSize = s_noKO.size();
            vector<double> y2(tmpSize);
            spline(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, tmpSize, 2e30, 2e30,&*y2.begin()-1);
            splint(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, &*y2.begin()-1, tmpSize, dnKOCM, &State[k].prePriceAtDnBar);
    //            prePriceAtDnBar = SplinePrice(&s_noKO[0], &noKOPrice[0], 0, s_noKO.size(), iLstDnNode[0]-bot, dnKOCM);
        }

        // to avoid remakr the insert level at insert node calc, make it true regardless it's found or not.
        State[k].foundUpLvl = State[k].foundDnLvl = true;
    }
}

// using Slice Operator Version (for general model, including hyb3). 
void CallableEquityKOSwapFDProd::CalcPayoff_oper(int step,const TreeSlice & spot,const vector< TreeSliceSP > & price)
{
    bool isAtMat = (step == model->getLastStep());
    int pEnd = price.size()-1;
        
    DateTime stepDate = model->getDate(step);

    int numFP = cks->finalPerfMaker->addNumPrices();
    
    double remainedCpn = 0.0;
    if (isAtMat)
    {
        double settlePV = 0.0;
        if (!cks->instSettleForFP.get()){
            settlePV = cks->instSettle->pvAdjust(stepDate,
                                            cks->discount.get(), 
                                            cks->asset.get());
        }else{
            settlePV = cks->instSettleForFP->pvAdjust(stepDate,
                                                      cks->discount.get(), 
                                                      cks->asset.get());
        }    
        double redemption = cks->redemption * settlePV;

        double fixedRemain = cks->hasFixedLeg ? koFix->getPV(stepDate, cks->discount.get()) : 0.0;
        double floatRemain = 0.0;
        if (cks->hasFloater) 
        {
            floatRemain = koFlt->getPV(stepDate, cks->discount.get());
        }                                                  
        remainedCpn = redemption + fixedRemain + floatRemain;
    }

    // notional remaining at stake (if breached dates)
    // TO DO : rate dependent slices....
//    if (!zeroProd){
//        const TreeSlice &zeroSlice  = zeroProd->getValue(step, payDate);
//
//    }else{
//
//    }
//        

    int idx = stepPerfIndex[step];  // idx<0 means it's not payoff step.

    double koCpns = stepKOCpns[step];

    int k;
    for (k=0; k<numState; k++){            
        if (isAtMat){
            // initialize
            *price[pEnd] = 0.0;
        }

        // TO DO : GetSliceIdx() !!
        // barriers for continuous/daily monitoring and rebates
        double upKOCM = State[k].adjUpBarrier;
        double dnKOCM = State[k].adjDnBarrier;

        int iPrice = baseNP*k; 
        
        // Store eq Link Coupon price into price[k+1] when there is payoff (idx>=0)
        if (!eqCpnKOLeg.isNull && idx >= 0){
            *price[iPrice+1] = EqCpnKOLeg::getValue_oper(
                                        spot, 
                                        step, 
                                        State[k].refLevel, 
                                        eqPerfs[idx], 
                                        eqCpnKOLeg);
        }
        // pay-off
//        if (!State[k].foundDnLvl) 
//            State[k].iLstDnNode.resize(2, bot-1);
//        else
//            State[k].iLstDnNode.resize(0);
//        if (!State[k].foundUpLvl) 
//            State[k].iFstUpNode.resize(2,  top+1);
//        else
//            State[k].iFstUpNode.resize(0);

        vector<double> noKOPrice;
        vector<double> s_noKO;
//        double koEqnCpn;
//        for (jNode=bot; jNode<=top; jNode++) {
//            //add the avgOut premium.
        if (isAtMat)
            *price[iPrice] = remainedCpn;                
        if (step == iStepAvgOutStart){
            // calculate the finalPerf's value on tree.
            *price[iPrice] += IFinalPerf::getValue_oper(spot, //*price[0],
                                                            State[k].finalPerf,
                                                            0,            //int iPrice, 
                                                            step,         //int simDateIdx, 
                                                            State[k].refLevel //double k,
                                                            );
        } 
        // add finalPerf's results on instrument prices.
        if (numFP>0){
            *price[iPrice+1+numFP] = IFinalPerf::upDateValues_oper(
                                                State[k].finalPerf, 
                                                step, 
                                                spot, 
                                                *price[iPrice+1+1], 
                                                *price[iPrice+1+numFP]);
        }                
    
        // 1. Add Equity Payment on observ date.  
        if (!eqCpnKOLeg.isNull && idx >= 0)
            *price[iPrice] += *price[iPrice+1];

        // The equity payments, which is not cancelled, are stored to koEqnCpn 
        // so as to avoid early termination.
        //koEqnCpn = stepEqCpnFactor[step] * ( *price[iPrice+1] );

        // 2. KO values
        #define up_payoff( x )                                  \
            IF( spot > upKOCM*(1.0-FP_MIN))                     \
                koCpns + stepUpRebate[step]                     \
                + stepEqCpnFactor[step] * ( *price[iPrice+1] )  \
            ELSE                                                \
                x                                               \
            ENDIF                                               \


        #define dn_payoff( x )                                  \
            IF( spot < dnKOCM*(1.0+FP_MIN))                     \
                koCpns + stepDnRebate[step]                     \
                + stepEqCpnFactor[step] * ( *price[iPrice+1] )  \
            ELSE                                                \
                x                                               \
            ENDIF                                               \

        
        #define up_fp(x)                                            \
            IF( spot > upKOCM*(1.0-FP_MIN) && cks->hasUpTrigger)    \
                0.0                                                 \
            ELSE                                                    \
                x                                                   \
            ENDIF                                                   \

        #define dn_fp(x)                                            \
            IF( spot < dnKOCM*(1.0+FP_MIN) && cks->hasLowTrigger)   \
                0.0                                                 \
            ELSE                                                    \
                x                                                   \
            ENDIF                                                   \

        if (cks->hasUpTrigger)
            *price[iPrice] = up_payoff( *price[iPrice] );
        if (cks->hasLowTrigger)
            *price[iPrice] = dn_payoff( *price[iPrice] );

        for (int iFP = 1; iFP <= numFP; iFP++){
            if (cks->hasUpTrigger)
                *price[iPrice+1+iFP] = up_fp( spot ); // KO the final Payoff
            if (cks->hasLowTrigger)
                *price[iPrice+1+iFP] = up_fp( spot ); // KO the final Payoff
        }

        #undef up_payoff
        #undef dn_payoff
        #undef up_fp
        #undef dn_fp

        // 3. early exercise
        if (cks->hasCallSchedule){
            double callOrPut = 0.0;     // cannot declare in macro
            double buff = 0.0;
            //double intrinsic = koCpns + koEqnCpn + stepCallLevel[step];

            #define callable(x , y)                                            \
                IF( callOrPut * ( x - ( stepEqCpnFactor[step] * y +  buff) ) > 0 )  \
                    stepEqCpnFactor[step] * y + buff                    \
                ELSE                                                    \
                    x                                                   \
                ENDIF                                                   \

            if (stepIsCallable[step]){
                buff = koCpns + stepCallLevel[step];
                callOrPut = cks->callSchedule.isPuttable ? 1.0 : -1.0;
                *price[iPrice] = callable(*price[iPrice] , *price[iPrice+1]);
            }
            #undef callable
        }

        // 4. add all cash flow, which are not affected by early termination.
        *price[iPrice] += stepFixedFee[step] + stepFloatFee[step];

        // 5. add FinalPerf
        if (step==0 && numFP>0)
            *price[iPrice] += *price[iPrice+1+1];

        if (State[k].iFstUpNode.size()>0 && State[k].foundUpLvl && upSmoothType != NO_SMOOTHING){//store no KO prices
//            int tmpSize = s_noKO.size();
//            vector<double> y2(tmpSize);
//            spline(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, tmpSize, 2e30, 2e30,&*y2.begin()-1);
//            splint(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, &*y2.begin()-1, tmpSize, upKOCM, &State[k].prePriceAtUpBar);
//                //prePriceAtUpBar = SplinePrice(&s_noKO[0], &noKOPrice[0], 0, s_noKO.size(), iFstUpNode[0]-bot, upKOCM);
        }
        if (State[k].iLstDnNode.size()>0 && State[k].foundDnLvl && lowSmoothType != NO_SMOOTHING){//store no KO prices
//            int tmpSize = s_noKO.size();
//            vector<double> y2(tmpSize);
//            spline(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, tmpSize, 2e30, 2e30,&*y2.begin()-1);
//            splint(&*s_noKO.begin()-1, &*noKOPrice.begin()-1, &*y2.begin()-1, tmpSize, dnKOCM, &State[k].prePriceAtDnBar);
//    //            prePriceAtDnBar = SplinePrice(&s_noKO[0], &noKOPrice[0], 0, s_noKO.size(), iLstDnNode[0]-bot, dnKOCM);
        }

        // to avoid remakr the insert level at insert node calc, make it true regardless it's found or not.
        State[k].foundUpLvl = State[k].foundDnLvl = true;
    }
}
void CallableEquityKOSwapFDProd::postCalc(int step, const TreeSlice & spot, const vector< TreeSliceSP > & price) 
{
    static const string method("CallableEquityKOSwapFDProd::postCalc");
    try 
    {
        if (tree1f){
            int bot, top;
            spot.getCalcRange( bot, top );    

            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            int idx = tree1f->GetSliceIdx();
            bool isCheaper;
            double insPrice;                
            int iPrice;
            bool needOWInsPrice = false;
            bool isSuccess;
            int k;
            // smoothing price curve.  Skip pEnd as it's EqGain performance (plain vanilla).
            if (stepIsDnKO[step])
            {
                DateTime stepDate = model->getDate(step);
                double sLow,sHigh,pLow,pHigh,orgDelta;
                isCheaper = cks->lowTrigger.getSmoothType() == DECREASE;                
                for (k = 0; k<numState; k++)
                {
                    if (cks->lowTrigger.getSmoothType() != NO_SMOOTHING && State[k].adjDnBarrier > s[-bot])
                    {
                        iPrice = 2*k;
                        tree1f->getInsertNodePrice(idx, iPrice, 0, &insPrice);
                        isSuccess = SmoothPrice(s, cks->lowTrigger.peakDelta, isCheaper, false,  
                                    State[k].iLstDnNode[iPrice], bot, top, State[k].adjDnBarrier, insPrice, 
                                    State[k].prePriceAtDnBar, (p[iPrice]), needOWInsPrice,sLow,sHigh,pLow,pHigh,orgDelta);
                        if (!isSuccess)
                        {
                            throw ModelException(method, "The peakDelta["+ Format::toString(cks->lowTrigger.peakDelta) + 
                                "] is too small.  \n The least delta at this step["+ stepDate.toString() +
                                "] = " + Format::toString(orgDelta) + 
                                ".  peakDelta should be enough bigger. ");
                        }
                        // need to overwrite insert node price...
                        if (needOWInsPrice)
                        {
                            tree1f->SetInsertNodeAndPrice(idx, 
                                                        0, 
                                                        State[k].adjDnBarrier, 
                                                        0, 
                                                        iPrice, 
                                                        iPrice, 
                                                        State[k].prePriceAtDnBar);
                        }
                        if (iPrice == 0){//keep debug record
                            db_DnBarSLow.push_back(sLow/State[0].refLevel);
                            db_DnBarSHi.push_back(sHigh/State[0].refLevel);
                            db_DnBarPLow.push_back(pLow);
                            db_DnBarPHi.push_back(pHigh);
                            db_DnBarDelta->push_back(CashFlow(stepDate, orgDelta));
                        }
                    }                    
                }
            }
            if (stepIsUpKO[step])
            {
                DateTime stepDate = model->getDate(step);
                isCheaper = cks->upTrigger.getSmoothType() == DECREASE;
                double sLow,sHigh,pLow,pHigh,orgDelta;
                for (k = 0; k<numState; k++)
                {
                    if (cks->upTrigger.getSmoothType() != NO_SMOOTHING && State[k].adjUpBarrier < s[top])
                    {
                        iPrice = 2*k;
                        tree1f->getInsertNodePrice(idx, iPrice, cks->hasLowTrigger ? 1 : 0, &insPrice);
                        isSuccess = SmoothPrice(s, cks->upTrigger.peakDelta, isCheaper,  true, 
                                    State[k].iFstUpNode[iPrice], bot, top, State[k].adjUpBarrier, insPrice, 
                                    State[k].prePriceAtUpBar,(p[iPrice]), needOWInsPrice,sLow,sHigh,pLow,pHigh,orgDelta);
                        if (!isSuccess)
                        {
                            throw ModelException(method, "The peakDelta["+ Format::toString(cks->upTrigger.peakDelta) + 
                                "] is too small.  \n The least delta at this step["+ stepDate.toString() +
                                "] = " + Format::toString(orgDelta) + 
                                ".  peakDelta should be enough bigger. ");
                        }
                        // need to overwrite insert node price...
                        if (needOWInsPrice){
                            tree1f->SetInsertNodeAndPrice(idx, cks->hasLowTrigger ? 1 : 0, State[k].adjUpBarrier, 
                                                        0, iPrice, iPrice, State[k].prePriceAtUpBar);
                        }                                                
                        if (iPrice == 0){//keep debug record
                            db_UpBarSLow.push_back(sLow/State[0].refLevel);
                            db_UpBarSHi.push_back(sHigh/State[0].refLevel);
                            db_UpBarPLow.push_back(pLow);
                            db_UpBarPHi.push_back(pHigh);
                            db_UpBarDelta->push_back(CashFlow(stepDate, orgDelta));
                        }
                    }
                }
            }

            // migrate the various state results.
            if (step == iStepAvgInEnd)
            {
                for (int jNode=bot; jNode<=top; jNode++){
                    // Integrate the state variable by pdf.
                    double sum = 0.0;
                    for (k=0;k<numState;k++)
                    {
                        sum += pdf_AvgIn[k] * (p[baseNP*k])[jNode];
                    }
                    (p[0])[jNode] = sum;
                }
                // no more need state variable.  But not speed up much.....  Can remove this.
                tree1f->ResetPriceArr(step, tree1f->GetSliceIdx(), 0, 0);
            }
        }        
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

CVolRequestConstSP CallableEquityKOSwapFDProd::GetLNRequest() const
{
    // get strike and maturity date from instrument
    DateTime        matDate = cks->getLastSimDate();

    double volStrike  = cks->initialSpot;

    DateTime imntStartDate = cks->fwdStarting? 
                         cks->startDate: cks->valueDate;

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                   matDate, cks->fwdStarting));
    return volRequest;
}

// return true if already dead instrument.
bool CallableEquityKOSwapFDProd::priceDeadInst(Control * control,Results * results)
{
    return cks->priceDeadInstrument(control, results);
}

/** create a fd payoff product */
FDProductSP CallableEquityKOSwap::createProduct(FDModel* model) const
{
    return FDProductSP( new CallableEquityKOSwapFDProd(this, model) );
}

/****************************************************************************************************************/
class CallableEquityKOSwapHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("CallableEquityKOSwap");
        REGISTER(CallableEquityKOSwap, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        IMPLEMENTS(KnownCashflows::IEventHandler);
        FIELD(finalPerfMaker, "option at maturity");
        FIELD(redemption, "final redemption amount");
        FIELD(instSettleForFP, "instrument settlement rule, applied only for Final Perf.  Empty => use instSettle");
        FIELD_MAKE_OPTIONAL(instSettleForFP);

        FIELD(hasEqLeg, "has Equity Linked Coupon Leg?");
        FIELD(koELC, "Equity Linked Cash Flow");
        FIELD_MAKE_OPTIONAL(koELC);
        
        FIELD(hasFloater, "has floater?");
        FIELD(koFloater, "libor leg with KO rule");
        FIELD_MAKE_OPTIONAL(koFloater);

        FIELD(hasFixedLeg, "has fixed cash flow?");
        FIELD(koFixedLeg,  "fixed cash flow with KO rule");
        FIELD_MAKE_OPTIONAL(koFixedLeg);
                
        FIELD(hasUpTrigger, " ");
        FIELD(upTrigger, "Upper Trigger");
        FIELD_MAKE_OPTIONAL(upTrigger);
        
        FIELD(hasLowTrigger, " ");
        FIELD(lowTrigger, "Lower Trigger");
        FIELD_MAKE_OPTIONAL(lowTrigger);        
            
        FIELD(hasCallSchedule, " ");    
        FIELD(callSchedule, "Call Schedule");
        FIELD_MAKE_OPTIONAL(callSchedule);

        FIELD(isAvgIn,  "is average-in?");
        FIELD_MAKE_OPTIONAL(isAvgIn);
        FIELD(avgIn,    "average-in schedule");
        FIELD_MAKE_OPTIONAL(avgIn);
        FIELD(useIniSpotAsRefLvl, "true : ignore the sample of avg-in, and use Spot At Start as reference level.");
        FIELD_MAKE_OPTIONAL(useIniSpotAsRefLvl);

        FIELD(samples,  "past sampling results");
        FIELD_MAKE_OPTIONAL(samples);

        FIELD(avgInSL, "internal SampleList for AvgIn");
        FIELD_MAKE_TRANSIENT(avgInSL);
        FIELD(numOfAvgInSV, "number of state variable for Average-In.");
        FIELD_MAKE_OPTIONAL(numOfAvgInSV);

#ifdef  TREE_THETA_CAP
        FIELD(thetaSmoothThrehold, "The threshold for theta smoothing.");
        FIELD_MAKE_OPTIONAL(thetaSmoothThrehold);
        FIELD(isPositiveThetaSmooth, "true for smoothing so as to FV positively increase.");
        FIELD_MAKE_OPTIONAL(isPositiveThetaSmooth);
#endif
        EMPTY_SHELL_METHOD(defaultCallableEquityKOSwap);
    }

    static IObject* defaultCallableEquityKOSwap(){
        return new CallableEquityKOSwap();
    }
};

CClassConstSP const CallableEquityKOSwap::TYPE = CClass::registerClassLoadMethod(
    "CallableEquityKOSwap", typeid(CallableEquityKOSwap), CallableEquityKOSwapHelper::load);

CClassConstSP const CallableEquityKOSwap::CallSchedule::TYPE = CClass::registerClassLoadMethod(
    "CallableEquityKOSwap::CallSchedule", typeid(CallableEquityKOSwap::CallSchedule), CallableEquityKOSwap::CallSchedule::load);

CClassConstSP const CallableEquityKOSwap::Trigger::TYPE = CClass::registerClassLoadMethod(
    "CallableEquityKOSwap::Trigger", typeid(CallableEquityKOSwap::Trigger), CallableEquityKOSwap::Trigger::load);

/* for class loading */
bool CallableEquityKOSwapLoad() {
    return (CallableEquityKOSwap::TYPE != 0);
}

///////////////////////////////////////////////////////
////   Barrier1F class                   //////////////
///////////////////////////////////////////////////////
//Barrier1F::Barrier1F(){
//        barrier = ScheduleSP(0);
//        isUp = false;
//        isCont = false;
//        breachDate = DateTime(0);
//        isBreached = false;
//        asset = CAssetWrapper(0);
//        initialSpot = 0;
//}

Barrier1F::Barrier1F(ScheduleSP barrier,
                    bool  isUp,
                    bool isCont,
                    DateTime breachDate,                                  
                    bool isBreached,
                    CAssetWrapper asset,
                    double initialSpot) :
        barrier(barrier),isUp(isUp),isCont(isCont), breachDate(breachDate),
        isBreached(isBreached),asset(asset),initialSpot(initialSpot){}

void Barrier1F::makeBarrierEvent(EventResults* events,
                        const DateTime& eventDate,
                        string barrDesc,
                        string barrierType){
            bool hasEvent = false;
            DateTime eDate = breachDate;
            double level = 0; 
            double barrLvl = 0;  
            // no barrier events unless the barrier is active
            try {
                barrLvl = barrier->interpolate(eventDate) * initialSpot;
            }
            catch(exception& e){
                // cannot find the interp level.  No need to return Barrier Event
                e;  // just to avoid warning
                return;
            }
                
            // if breached flag is set only report event if it's now
            if (isBreached && breachDate == eventDate) {
                    level = asset->fwdValue(eventDate);
                    eDate = breachDate;
                    hasEvent = true;                
            }else {
                // if not doing cts monitoring we can only have event
                // if time of day is right
                if (isCont ||
                    eventDate.getTime() == barrier->firstDate().getTime()) {
                    // we are monitoring, we haven't breached yet
                    // let's see if it's breaching right now
                    level = asset->fwdValue(eventDate);
                    if (isUp ? (level > barrLvl) : (level < barrLvl)) 
                        hasEvent = true;
                }
            }

            // convert to Array
            if (hasEvent){
                StringArraySP assetNames(new StringArray(1, asset->getTrueName()));
                DoubleArraySP assetLevels(new DoubleArray(1, level));
                DoubleArraySP barrLevels(new DoubleArray(1, barrLvl));
                string monType = isCont ? BarrierBreach::CONTINUOUS :
                                          barrier->getInterp() == Schedule::INTERP_NONE ?
                                              BarrierBreach::EUROPEAN :
                                              BarrierBreach::DAILY;
                events->addEvent(new BarrierBreach(breachDate, barrDesc,
                                                   monType,
                                                   barrierType, isUp, 0,
                                                   assetNames, assetLevels, barrLevels));
            }
}


DRLIB_END_NAMESPACE
