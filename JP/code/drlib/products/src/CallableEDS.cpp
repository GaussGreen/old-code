//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CallableEDS.cpp
//
//   Description : Callable Equity Stability Swap - double out version of EDS
//                 with cal feature
//
// Client 
//   o PAYS a fee until earliest of maturity/knockout.
//   o On knockout, RECEIVES rebate, pays fee accrued. 
//   o KO is either up or down (i.e double out)
//   o Can call - pay call amount + accrued fee
// 
//   Author      : Andrew J Swain
//
//   Date        : 12 June 2003
//
//----------------------------------------------------------------------------
//
//  IMPORTANT REMARK:
//
/*  There is a slight error concerning the payment of fees:
    Fee payments should always be settled T+0, while payments concerning
    rebates and call amounts usually settle T+x.
    In this code, the payment of fees is implemented as follows:
        - If a fee payment has to be made and if on that day neither
            a rebate nor a call amount have to be paid, then the fee 
            is settled T+0.
        - If a fee payment has to be made and if on that day either 
            a rebate or a call amount (or both) have to be paid, then the 
            fee is settled in the same way as  the rebate resp the call 
            amount are settled, namely T+x.
    In reality, this is not really the case. However, implementing the payment
    of a fee in this way yields consistency with CallableEDS.cpp and 
    EquityStabilitySwap.cpp (thus the old test cases can be used as benchmarks). 
    For details, either ask Yi Gu or Eva Strasser.
*/


#include "edginc/config.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/AssetUtil.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class CallableEDS: public Generic1Factor,
                   virtual public LegalTerms::Shift,
                   virtual public FDModel::IIntoProduct, 
                   virtual public LastSensDate {
public:
    static CClassConstSP const TYPE;  

    static const string UP;
    static const string DOWN;

    // useful data holder
    class CallSchedule: public CObject {
    public:
        static CClassConstSP const TYPE;  

        // fields
        bool          isIssuerCall;         // isIssuerCall=True => issuer has right to call
        DateTimeArray notification;         // list of notification dates
        CashFlowArray schedule;             // list of call dates and corresponding call amounts
        bool          isCalled;             // is it called ?
        DateTime      callDate;             // if so, what was the NOTIFICATION date ?
        
        // transient to save repeated lookups
        int           callidx;              // if isCalled=True, then callidx is index for schedules
        
        CallSchedule():CObject(TYPE), isIssuerCall(false), isCalled(false), 
            callidx(0) {}; 

        void validatePop2Object() {
            static const string method(
                "CallableEDS::CallSchedule::validatePop2Object");
            try {
                if (notification.size() != schedule.size()) {
                    throw ModelException(method, "notification dates and call "
                                         "schedule are different lengths");
                } 

                CashFlow::ensureDatesIncreasing(schedule,
                                                "call schedule", false);
                DateTime::ensureIncreasing(notification, 
                                           "call notification", false);
                
                // can't call before notification
                int i;
                for (i = 0; i < notification.size(); i++) {
                    if (notification[i] > schedule[i].date) {
                        throw ModelException(method, 
                                             "call notification date (" + 
                                             notification[i].toString()+
                                             ") must come before call date (" + 
                                             schedule[i].date.toString() + 
                                             ")");
                    }              
                }    

                // check that notification dates and call dates are always EOD
                for (i=0; i < notification.size(); i++) {
                    if (notification[i].getTime() != DateTime::END_OF_DAY_TIME) {
                        throw ModelException(method,
                            "notification dates must always be EOD (" + 
                            notification[i].toString() + 
                            ")");
                    }
                    if (schedule[i].date.getTime() != DateTime::END_OF_DAY_TIME) {
                        throw ModelException(method,
                            "call dates must always be EOD (" + 
                            schedule[i].date.toString() +
                            ")");
                    }
                }

                if (isCalled) {
                    // see if the call date is for real
                    bool found = false;
                    for (i = 0; i < notification.size() && !found ; i++) {
                        if (callDate == notification[i]) {
                            found = true;
                            callidx = i;
                        }
                    }
                    if (!found) {
                        throw ModelException(method,
                            "structure is called, but call "
                            "date ("+callDate.toString() + 
                            ") is not in notification dates");
                    }
                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }


    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("Callable ESS Call Schedule");
            REGISTER(CallableEDS::CallSchedule, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultCallSchedule);
            FIELD(isIssuerCall, "issuer can call");
            FIELD_MAKE_OPTIONAL(isIssuerCall);
            FIELD(notification, "notification dates");
            FIELD(schedule, "schedule");
            FIELD(isCalled, "is it called?");
            FIELD(callDate, "notification date when called");
            FIELD_MAKE_OPTIONAL(callDate);
            FIELD_NO_DESC(callidx);
            FIELD_MAKE_TRANSIENT(callidx);
        }
        
        static IObject* defaultCallSchedule(){
            return new CallableEDS::CallSchedule();
        }
    };

    virtual void Validate(){
        static const string method = "CallableEDS::Validate";
        
        // Call Generic1Factor validate()
        validate();

        // make old CallableEDS and new CallableEDS compatible, so that 
        // old test cases can be used
        if (pcNotionalTranche.size()==0) {
            pcNotionalTranche.resize(1);
            pcNotionalTranche[0] = 1.0;
        }
        if (isBreachedUp) {
            breachDates.resize(1);
            breachDates[0] = breachDate;
            breachUpOrDown.resize(1);
            breachUpOrDown[0] = UP;
        }
        if (isBreachedDown) {
            breachDates.resize(1);
            breachDates[0] = breachDate;
            breachUpOrDown.resize(1);
            breachUpOrDown[0] = DOWN;
        }

        if (oneContract) {
            throw ModelException(method, 
                                 "Callable ESS with Tranches must be "
                                 "notionally booked");
        }

        // sum of tranches = 100%
        double sumNotionalTranche = 0.0;
        for (int iTranche=0; iTranche<pcNotionalTranche.size(); iTranche++) {
            sumNotionalTranche += pcNotionalTranche[iTranche];
        }
        if (!Maths::areEqualWithinTol(sumNotionalTranche, 1.0, 0.000001)) {
            throw ModelException(method,
                "the sum of percentages represented by each tranche "
                "has to be 1.0");
        }

        // intraday monitoring only possible for single tranche products
        if (intraDayMonitor){
            if (pcNotionalTranche.size()!=1) {
                throw ModelException(method,
                    "intraday monitoring only allowed for single tranche products");
            }
        }

        if (fwdStarting && startDate <= valueDate) {
            throw ModelException(method, 
                                 "instrument is fwd starting but start date ("+
                                 startDate.toString() + ") is <= today ("+
                                 valueDate.toString());
        }

        // let's not go there shall we
        if (instSettle->isMargin()) {
            throw ModelException(method, "margin settlement not supported");
        }
        if (instSettle->isPhysical()) {
            throw ModelException(method, "physical settlement not supported");
        }

        if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
            throw ModelException(method, "ccy struck not supported");
        }
          
        if (upBarrier->getInterp() == Schedule::INTERP_NONE) {
            throw ModelException(method, 
                                 "'dates only' up barrier not supported");
        }
        if (upRebate->getInterp() == Schedule::INTERP_NONE) {
            throw ModelException(method, 
                                 "'dates only' up rebate not supported");
        }
        if (downBarrier->getInterp() == Schedule::INTERP_NONE) {
            throw ModelException(method, 
                                 "'dates only' down barrier not supported");
        }
        if (downRebate->getInterp() == Schedule::INTERP_NONE) {
            throw ModelException(method, 
                                 "'dates only' down rebate not supported");
        }             

        if (downBarrier->length() == 0 || upBarrier->length() == 0) {
            throw ModelException(method, 
                                 "barriers must be of non-zero length");
        }             
        if (downRebate->length() == 0 || upRebate->length() == 0) {
            throw ModelException(method, 
                                 "rebates must be of non-zero length");
        }             


        // fee dates in strict order
        DateTime::ensureStrictlyIncreasing(feeDates, "fee payment dates", true);

        // accrue dates if necessary
        if (useAccrueDates) {
            // in order
            DateTime::ensureStrictlyIncreasing(accrueDates, "accrue dates", true);

            // one fee payment date for each accrue date
            if (accrueDates.size() != feeDates.size()) {
                throw ModelException(method,  
                                     "The number of accrue dates (" + 
                                     Format::toString(accrueDates.size()) + 
                                     ") must be equal to the number of fee dates (" + 
                                     Format::toString(feeDates.size()) + ")");
            }

            // fee payment date after accrue date
            int i;
            for(i=0; i<accrueDates.size();i++) {
                if (accrueDates[i]>feeDates[i]) {
                    throw ModelException(method,  
                                         "The " + Format::toString(i) + 
                                         " accrue date (" + accrueDates[i].toString() +
                                         ") must be equal or before the " + 
                                         Format::toString(i) + " fee date (" + 
                                         feeDates[i].toString() +")");
                }
            }
        }
            

        // VALIDATE UP & DOWN BARRIERS
        // last fee and last ko/rebate date are co-incident
        DateTime lastFee = feeDates[feeDates.size()-1];
        DateTime lastKO  = downBarrier->lastDate();
        DateTime lastReb = downRebate->lastDate();

        if (lastKO != lastReb) {
            throw ModelException(method,  
                                 "last barrier (" + lastKO.toString() + 
                                 ") and last rebate (" + lastReb.toString() + 
                                 ") must co-incide");
        }


        // first fee/ko/rebate not before initial accrue
        DateTime firstFee = feeDates[0];
        DateTime firstKO  = downBarrier->firstDate();
        DateTime firstReb = downRebate->firstDate();
       
        if (firstFee <= initialAccrueDate) {
            throw ModelException(method, "first fee (" + firstFee.toString() + 
                                 ") must be after initial accrue date (" +
                                 initialAccrueDate.toString() + ")");
        }

        // cross-check up/down barrier/rebate
        DateTime lastUpKO   = upBarrier->lastDate();
        DateTime lastUpReb  = upRebate->lastDate();
        DateTime firstUpKO  = upBarrier->firstDate();
        DateTime firstUpReb = upRebate->firstDate();

        if (lastUpKO != lastKO) {
            throw ModelException(method, 
                                 "last up barrier (" +lastUpKO .toString() + 
                                 ") and last down barrier (" +
                                 lastKO.toString() + 
                                 ") must co-incide"); 
        }
        if (lastUpReb != lastReb) {
            throw ModelException(method, 
                                 "last up rebate (" +lastUpReb .toString() + 
                                 ") and last down rebate (" +
                                 lastReb.toString() + 
                                 ") must co-incide"); 
        }            
        if (firstUpKO != firstKO) {
            throw ModelException(method, 
                                 "first up barrier (" +firstUpKO .toString() + 
                                 ") and first down barrier (" +
                                 firstKO.toString() + 
                                 ") must co-incide"); 
        }
        if (firstUpReb != firstReb) {
            throw ModelException(method, 
                                 "first up rebate (" +firstUpReb .toString() + 
                                 ") and first down rebate (" +
                                 firstReb.toString() + 
                                 ") must co-incide"); 
        }            

      
        // are there any entries in the breached schedule (either time or type)?
        bool breached = (breachDates.size()>0 || breachUpOrDown.size()>0);

        if(breached) {
            // breach dates in order
            DateTime::ensureIncreasing(breachDates, "breached dates", true);
            // # breached dates inferior or euqal to # tranches of notional
            if (breachDates.size()>pcNotionalTranche.size()) {
                throw ModelException(method, 
                    "number of breached dates must be inferior "
                    "or equal to the number of tranches of notional");
            }
            // # breached dates equals # number of breachUpOrDown
            if (breachDates.size() != breachUpOrDown.size()){
                throw ModelException(method, 
                    "number of breached dates must be consistent "
                    "with the inputs for up/down");
            }
            // content of breachUpOrDown
            bool val=true;
            for (int breach=0; breach<breachUpOrDown.size(); breach++) {
                if (breachUpOrDown[breach] != UP && breachUpOrDown[breach] != DOWN) {
                    val = false;
                }
            }
            if (!val) {
                throw ModelException(method,
                    "Up/Down must only be filled with Up or Down");
            }
        }
         

        if (fwdStarting && breached) {
            throw ModelException(method, 
                                 "instrument is fwd starting and breached");
        }

        if (breached && breachDates.back() > valueDate) {
            throw ModelException(method, "instrument is knocked out but ko "
                                 "date (" + breachDates.back().toString() + ") is "
                                 "after today (" + valueDate.toString() + ")");
        }          

        
        // if there's an economic barrier, check it lines up with risk barrier
        if (downEcoBarrier.get()) {
            DateTime firstEco = downEcoBarrier->firstDate();
            DateTime lastEco  = downEcoBarrier->lastDate();
            if (firstEco != firstKO || lastEco != lastKO) {
                throw ModelException(method, "mis-matched dates between "
                                     "down risk barrier ("+firstKO.toString()+
                                     ", " + lastKO.toString() + ") and "
                                     "down economic barrier ("+
                                     firstEco.toString()+
                                     ", " + lastEco.toString() + ")");
            }
        }
        if (upEcoBarrier.get()) {
            DateTime firstEco = upEcoBarrier->firstDate();
            DateTime lastEco  = upEcoBarrier->lastDate();
            if (firstEco != firstUpKO || lastEco != lastUpKO) {
                throw ModelException(method, "mis-matched dates between "
                                     "up risk barrier ("+firstUpKO.toString()+
                                     ", " + lastUpKO.toString() + ") and "
                                     "up economic barrier ("+
                                     firstEco.toString()+
                                     ", " + lastEco.toString() + ")");
            }
        }

        // check call schedule
        if (!callSchedule.schedule.empty()) {
            int numCalls = callSchedule.schedule.size();
            DateTime& lastCall = callSchedule.schedule[numCalls-1].date;
            if (lastCall > lastKO) {
                throw ModelException(method, 
                                     "last call date (" + lastCall.toString() + 
                                     ") is after maturity ("+lastKO.toString()+")");
            }     
            
            DateTime& firstNotif = callSchedule.notification[0];

            if (firstNotif < initialAccrueDate) {
                throw ModelException(method, 
                                     "first notification date (" +
                                     firstNotif.toString() + 
                                     ") is before initial accrue date (" +
                                     initialAccrueDate.toString() + ")"); 
            }     

            if (fwdStarting && startDate > firstNotif) {
                throw ModelException(method, 
                                     "start date (" + startDate.toString() + 
                                     ") is after first call notice date ("+
                                     firstNotif.toString()+")");                
            }    

        }   

        if (callSchedule.isCalled && callSchedule.callDate > valueDate) {
                throw ModelException(method, 
                                     "structure is called, but notification date ("+
                                     callSchedule.callDate.toString() + 
                                     ") is after today (" + 
                                     valueDate.toString() +")");
        }
        if (callSchedule.isCalled) {
            DateTime callDate = callSchedule.schedule[callSchedule.callidx].date;
            if(breachDates.size()>0 && callDate < breachDates.back()) {
                throw ModelException(method,
                    "structure is called and there are breach dates "
                    "afterwards, e.g. " + callDate.toString() + ".");
            }

        }
       
        // if there is a breach on the valueDate, make sure that you also 
        // include it in the breached schedule
        if (upBarrier->coversDateRange(valueDate, valueDate, false) && 
            (valueDate.getTime() == DateTime::END_OF_DAY_TIME) ) {
            double refLevel2 = fwdStarting ? 
                asset->fwdValue(startDate) : initialSpot;
            if (asset->getSpot() > upBarrier->lastValue()*refLevel2-FP_MIN) {
                if (breachDates.empty() || 
                    (!breachDates.empty() && breachDates.back() != valueDate)) {
                    throw ModelException(method,
                        "today, " + valueDate.toString() +", is monitoring date and"
                        " and there is a up breach, include this into breach schedule");
                }
            } else if (asset->getSpot() < downBarrier->lastValue()*refLevel2+FP_MIN) {
                if (breachDates.empty() || 
                    (!breachDates.empty() && breachDates.back() != valueDate)) {
                    throw ModelException(method,
                        "today, " + valueDate.toString() +", is monitoring date and"
                        " and there is a down breach, include this into breach schedule");
                }
            }
        }

        // and finally
        daycount = DayCountConventionSP(DayCountConventionFactory::make(dcc));
        if (breachDates.size()>0) {
            breachSettleDates.resize(breachDates.size());
            for (int i=0; i<breachDates.size();i++) {
                breachSettleDates[i] = instSettle->settles(breachDates[i], asset.get());
            }
        }
    }

private:
    friend class CallableEDSHelper;
    friend class CEDFDProd;  

    CallableEDS():Generic1Factor(TYPE),
                  fee(0.0),
                  accrueDates(0),
                  noAccruedOnBreach(false),
                  useAccrueDates(false),
                  intraDayMonitor(false),
                  isBreachedUp(false),
                  isBreachedDown(false){}; 

    CallableEDS(const CallableEDS& rhs);
    CallableEDS& operator=(const CallableEDS& rhs);
   
   // PV of any unsettled cash - i.e. fees "paid" but not yet settled
    double unsettledCash() const {
       static const string method("CallableEDS::unsettledCash");
        try {
            // attention to sign: add fee payments and deduct rebate payments
            double cash = 0.0;
            
            // unsettled rebate payments
            if (!breachDates.empty()) {
                for (int dTrancheKO=0; dTrancheKO<breachDates.size(); dTrancheKO++) {
                    if (instSettle->settles(breachDates[dTrancheKO], asset.get()) > valueDate) {
                        // accrued at breach date
                        double accrued = accruedFee(breachDates[dTrancheKO]);
                        // pv for low date = value date and hi date = breachdate
                        double settlePV = instSettle->pv(valueDate, 
                                                         instSettle->settles(breachDates[dTrancheKO], 
                                                                             asset.get()),
                                                         discount.get(), asset.get());
                        // Up or Down rebate at breach date 
                        double stepUpOrDownRebate = (breachUpOrDown[dTrancheKO] == CallableEDS::UP) ? 
                            upRebate->interpolate(breachDates[dTrancheKO]) :
                            downRebate->interpolate(breachDates[dTrancheKO]);
                        
                        cash -= pcNotionalTranche[dTrancheKO]*(stepUpOrDownRebate - accrued)*settlePV;
                    }
                }
            }
            
            return cash;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    // what is the i-th accrue fee
    DateTime getAccrueDate(int idx) const {
        static const string method("EquityDistressSwap::getAccruedDate");
        if (idx>=feeDates.size()) {
            throw ModelException(method, "there are  "+ Format::toString(feeDates.size()+1) +
                                 "accrue dates and "+ " the accrue fee date of index " + 
                                 Format::toString(idx+1) +" is required");
        }
        else {
            // select the correct accrue dates
            return useAccrueDates?accrueDates[idx]:feeDates[idx];
        }
    }

    // what is the notional for a fee date
    double notionalFee(const DateTime& payDate) const {
        int      feeidx;
        // find the fee date immediately on or after this date
        bool found = false;
        for (int i = 0; i < feeDates.size() && !found; i++) {
            if (feeDates[i] == payDate) {
                found   = true;
                feeidx  = i;
            }
        }
        
        if (!found) {
            // this might actually happen now we allow barrier to go beyond last fee on last day
            // on last day, after the fee has gone, accrued is zero
            return 0.0;
        }
        else {
            double notionalFee = 1.0;
            int numTrancheKO = breachDates.size();  
            int iTrancheKO;
            for (iTrancheKO = 0; iTrancheKO < numTrancheKO && 
                     breachDates[iTrancheKO] < getAccrueDate(feeidx); iTrancheKO++) {
                notionalFee -= pcNotionalTranche[iTrancheKO];
            }
                        
            return notionalFee;
        }
    }

    
    
    // what's the accrued fee at a given date ?
    double accruedFee(const DateTime& hitDate) const {
        static const string method("EquityDistressSwap::accruedFee");
        try {
            // pay the accrued fee only if required
            if (noAccruedOnBreach) {
                return 0.0;
            }
            else {
                int      payidx;
                
                // find the accrue date immediately strictly after this date
                bool found = false;
                for (int i = 0; i < feeDates.size() && !found; i++) {
                    if (getAccrueDate(i) > hitDate) {
                        found   = true;
                        payidx  = i;
                    }
                }
                
                if (!found) {
                    // this might actually happen now we allow barrier to go beyond last fee on last day
                    // on last day, after the fee has gone, accrued is zero
                    return 0.0;
                }
                else {
                    DateTime lo = payidx == 0 ? initialAccrueDate : getAccrueDate(payidx-1);
                    double accrued = fee * daycount->years(lo, hitDate);
                    return accrued;
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // what's the fee to be paid at a given date ?
    double feeToPay(const DateTime& payDate) const {
        static const string method("EquityDistressSwap::feeToPay");
        try {
            int      feeidx;
            // find the fee date immediately on or after this date
            bool found = false;
            for (int i = 0; i < feeDates.size() && !found; i++) {
                if (feeDates[i] == payDate) {
                    found   = true;
                    feeidx  = i;
                }
            }
            
            if (!found) {
                // this might actually happen now we allow barrier to go beyond last fee on last day
                // on last day, after the fee has gone, accrued is zero
                return 0.0;
            }
            else {
                DateTime lo = feeidx == 0 ? initialAccrueDate : getAccrueDate(feeidx-1);
                double accrued = fee * daycount->years(lo, getAccrueDate(feeidx));
                return accrued;
            }
        }
        
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // is there a fee due between notification and call ?
    // return PV of any fee as of notification date
    double dueFee(int notifIdx) const {
        static const string method("CallableEDS::dueFee");
        try {
            const DateTime& notifDate = callSchedule.notification[notifIdx];          
            const DateTime& callDate  = callSchedule.schedule[notifIdx].date;

            int payidx = 0;
            double pv = 0.0;
            // find the fee date immediately after notif date
            bool found = false;
            for (int i = 0; i < feeDates.size() && !found; i++) {
                if (feeDates[i] > notifDate && feeDates[i] < callDate){
                    found = true;
                    pv = discount->pv(notifDate, feeDates[i]);
                    payidx = i;
                }
            }
            
            DateTime low = payidx == 0 ? initialAccrueDate : feeDates[payidx-1];
            DateTime high = feeDates[payidx];

            return (found ? fee*pv*daycount->years(low,high) : 0.0);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }



    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        if (CTree1fLV::TYPE->isInstance(model)) {
            return true; // do pointwise instead
        }
        return false;
    }

    // make a LN vol request - not in real use as prefer LV tree
    CVolRequest* makeVolRequest() const {
        static const string method("CallableEDS::makeVolRequest");
        try {
            double level = downBarrier->lastValue();
            if (!fwdStarting) {
                level *= initialSpot;
            }
            // do it at last barrier level for want of anywhere better
            return new LinearStrikeVolRequest(level,
                                              startDate,
                                              downBarrier->lastDate(),
                                              fwdStarting);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns all strikes the CallableEDS is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method("CallableEDS::getSensitiveStrikes");
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
    DateTime endDate(const Sensitivity* sensControl) const{
        DateTime maturity = downBarrier->lastDate();
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;    
    }

    /** Satisfy LegalTerms::Shift interface */
    bool sensShift(LegalTerms* shift) {
        // Set the barriers for pricing equal to the economic barriers
        if (upEcoBarrier.get()) {
            upBarrier = upEcoBarrier;
        }
        if (downEcoBarrier.get()) {
            downBarrier = downEcoBarrier;
        }
        
        return false;
    }

    bool sensShift(Theta* shift) {
        static const string method = "CallableEDS::sensShift";
        try  {
            // nothing to do here
            // roll the parent (updates value date etc)
            Generic1Factor::sensShift(shift);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
        return true; // our components have theta type sensitivity
    }
       
    // convenient fct for priceDeadInstrument
    double shortHitPrice(int tranche, double rebate, const DateTime& feeDate) const {
        return pcNotionalTranche[tranche] * (rebate - accruedFee(feeDate));
    }
    
    bool priceDeadInstrument(Control* control, Results* results) const{
        static const string method = "CallableEDS::priceDeadInstrument";
        try  {
            bool        deadInstrument = false;
            
            // important convention: 
            // if KO on valueDate, then it must be included in breach schedule
            // => it is not necessary to check whether KO on finalDate or on callDate!
            
            // finalDate = last date on barrier schedule = last date on fee schedule
            // finalDateSettlement = settlement date of finalDate
            DateTime    finalDate           = downBarrier->lastDate();
            DateTime    finalDateSettlement = instSettle->settles(finalDate,asset.get());
            
            double      value = 0.0;
    
            bool    breachScheduleEmpty = breachDates.empty();
            int     numBreachDates = breachDates.size();
            int     numNotionalTranche = pcNotionalTranche.size();

            const   CallSchedule& call  = callSchedule;
            DateTime    callDate; 
            if(call.isCalled) {
                callDate = call.schedule[call.callidx].date;
            }
            
            // if sthg is to be priced via priceDeadInstrument, then all KO are
            // already included in the KO schedule
            double pcNotionalRemaining = 1.0;
            int iBreachDate = 0;
            for (iBreachDate=0; iBreachDate<numBreachDates; iBreachDate++) {
                pcNotionalRemaining -= pcNotionalTranche[iBreachDate];
            }

            
            
            
            if(call.isCalled && valueDate >= instSettle->settles(callDate,asset.get())){
                // everything is zero
                results->storePrice(0.0, discount->getCcy());
                if(control && control->isPricing()) {
                    recordOutputRequests(control, results, 0.0);
                }
                deadInstrument = true; 
            // ************************************************************
            } 
            else if (numBreachDates == numNotionalTranche) {
                // all tranches are used up, however, not everything has been paid
                // => structure is not terminated by calling!
                // 1. remaining payments due to breach dates
                for (int iBreachDate = 1; iBreachDate<=numBreachDates; iBreachDate++) {
                    // for every breach date: check whether there are unsettled payments or not 
                    DateTime breachDate = breachDates[numBreachDates-iBreachDate];
                    DateTime settleBreach = instSettle->settles(breachDate, asset.get());
                    if (settleBreach > valueDate) {
                        // if not yet paid, get fee accrued and rebate
                        double rebate;
                        if (breachUpOrDown[numBreachDates-iBreachDate] == UP) {
                            rebate = upRebate->interpolate(breachDate);
                        } else if (breachUpOrDown[numBreachDates-iBreachDate] == DOWN) {
                            rebate = downRebate->interpolate(breachDate);
                        }
                        value += shortHitPrice(numNotionalTranche-iBreachDate,rebate,breachDate)
                            * notional * discount->pv(settleBreach);
                    }
                }
                // 2. remaining payments due to fee dates
                DateTime settleFee;
                for (int i=0; i<feeDates.size(); i++) {
                    settleFee = feeDates[i];
                    if(settleFee > valueDate && getAccrueDate(i) <= breachDates[breachDates.size()-1]) {
                        // determine corresponding notional for fee
                        double feeNotional = 1.0;
                        for (int j=0; j<breachDates.size(); j++) {
                            if (breachDates[j] < getAccrueDate(i)) {
                                feeNotional -= pcNotionalTranche[j];               
                            }
                        }
                        
                        // pay if strictly after today
                        if (feeDates[i] > valueDate) {
                            value -= feeNotional * feeToPay(feeDates[i])
                                * notional * discount->pv(feeDates[i]);
                        }
                    }
                }
                
                results->storePrice(value, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, value);
                }              
                deadInstrument = true;  
                // ************************************************************
            } else if (call.isCalled &&  
                       callDate <= valueDate && valueDate < instSettle->settles(callDate,asset.get())) {

                double currentTranche = pcNotionalTranche[numNotionalTranche-numBreachDates-1];

                // recall that callDate = finalDate is excluded
                // call Date <= value Date < settlement(callDate)
                double callAmount = call.schedule[call.callidx].amount;
                double accrued = accruedFee(callDate);
                // if callDate = breachDate, then breach is already included in breachSchedule
                if(call.isIssuerCall) {
                    value += callAmount * pcNotionalRemaining;
                } else if (!call.isIssuerCall) {
                    value -= callAmount * pcNotionalRemaining;
                }
                // check breach schedule whether callDate = breachDate
                for (int i=0; i<numBreachDates; i++) {
                    if (breachDates[i] == callDate) {
                        if(breachUpOrDown[i] == UP) {
                            value += currentTranche * upRebate->interpolate(callDate);
                        } else if (breachUpOrDown[i] == DOWN) {
                            value += currentTranche * downRebate->interpolate(callDate);
                        }
                    }
                }
                // now the fee
                value -= accrued * pcNotionalRemaining;
                // finally the settlement
                DateTime callDateSettlement = instSettle->settles(callDate,asset.get());
                value *= discount->pv(callDateSettlement) * notional;
                
                // now, everything on callDate is done
                // check for "old things", i.e. old fees and old rebates 
                
                if (!breachScheduleEmpty) {
                    for (int i=0; i<breachDates.size() && breachDates[i]<callDate; i++) {
                        if (valueDate < breachSettleDates[i]) {
                            // client receives rebate and pays fee accrued
                            DateTime settle = instSettle->settles(breachDates[i],asset.get());
                            double rebate;
                            double accrued = accruedFee(breachDates[i]);
                            if (breachUpOrDown[i]==UP) {
                                rebate = upRebate->interpolate(breachDates[i]);
                                value += pcNotionalTranche[i] * (rebate-accrued) 
                                    *notional*discount->pv(settle);
                            } else if (breachUpOrDown[i]==DOWN) {
                                rebate = downRebate->interpolate(breachDates[i]);
                                value += pcNotionalTranche[i] * (rebate-accrued) 
                                    *notional*discount->pv(settle);
                            }
                        }
                    }
                }
                DateTime payFee;
                for (int iFee=0; iFee<feeDates.size(); iFee++) {
                    payFee = feeDates[iFee];
                    if(payFee > valueDate) {
                        // determine corresponding notional for fee
                        double feeNotional = 1.0;
                        for (int j=0; j<breachDates.size(); j++) {
                            if (breachDates[j] < getAccrueDate(iFee)) {
                                feeNotional -= pcNotionalTranche[j];               
                            }
                        }
                        // pay fee if strictly after today
                        if (payFee > valueDate) {
                            value -= notional * feeNotional * 
                                feeToPay(payFee) * discount->pv(payFee);
                        }
                    }
                }

                results->storePrice(value, discount->getCcy());
                if(control && control->isPricing()) {
                    recordOutputRequests(control, results, 0.0);
                }
                deadInstrument = true;
                
                // ************************************************************
            } else if (valueDate >= finalDate) {
                // structure is definitely NOT called (see previous case)
                // situation:   we are at or past finalDate (= last date of barrier schedule)
                //              and there are unsettled payments
                // I. check whether sthg from finalDate is unsettled
                // II. check whether even "older things" are unsettled
                
                double currentTranche = pcNotionalTranche[numNotionalTranche-numBreachDates-1];

                // I.A) check whether there is a KO on finalDate
                // note: any KO on finalDate has to be entered into breachSchedule! 
                
                if (!breachScheduleEmpty && breachDates.back()==finalDate) {
                    double rebate;
                    if (breachUpOrDown.back() == UP) {
                        rebate = upRebate->lastValue();
                        value += currentTranche * (rebate-accruedFee(finalDate)); 
                    } else if (breachUpOrDown.back() == DOWN) {
                        rebate = downRebate->lastValue();
                        value += currentTranche * (rebate-accruedFee(finalDate));
                    }
                }
                value *=notional*discount->pv(finalDateSettlement);
                
                
                // II.A) check whether previous rebate payments need to be settled
                if (!breachScheduleEmpty) {
                    for (int i=0; i<breachDates.size() && breachDates[i]<finalDate; i++) {
                        if (valueDate < instSettle->settles(breachDates[i],asset.get())) {
                            // client receives rebate and pays fee accrued
                            DateTime settle = instSettle->settles(breachDates[i],asset.get());
                            double rebate;
                            if (breachUpOrDown[i]==UP) {
                                rebate = upRebate->interpolate(breachDates[i]);
                                value += shortHitPrice(i,rebate,breachDates[i])
                                    *notional*discount->pv(settle);
                            } else if (breachUpOrDown[i]==DOWN) {
                                rebate = downRebate->interpolate(breachDates[i]);
                                value += shortHitPrice(i,rebate,breachDates[i])
                                    *notional*discount->pv(settle);
                            }
                        }
                    }
                }
                // II.B) check whether previous fee payments need to be settled
                DateTime payFee;
                for (int i=0; i<feeDates.size(); i++) {
                    payFee = feeDates[i];
                    if(payFee > valueDate) {
                        // determine corresponding notional for fee
                        double feeNotional = 1.0;
                        for (int j=0; j<breachDates.size(); j++) {
                            if (breachDates[j] < getAccrueDate(i)) {
                                feeNotional -= pcNotionalTranche[j];               
                            }
                        }
                        
                        value -= feeNotional * feeToPay(payFee)
                            * notional * discount->pv(payFee);
                    }
                }

                results->storePrice(value, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, value);
                }
                deadInstrument = true;  
            }     
            return deadInstrument;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }            
    }

    /** extra output requests */
    void recordOutputRequests(Control* control, 
                              Results* results, 
                              double   fairValue) const {
        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(
            control,
            results,
            downBarrier->lastDate(),
            valueDate,
            asset.get());
        
        // for shorthand
        const CallSchedule& call  = callSchedule;
        
        bool breached = (breachDates.size() > 0);
        int numBreachDates = breachDates.size();
        OutputRequest* request = NULL;
        
        // I. PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArray paydates;
            
            // if (!breached), then all fee dates (+ settlement) 
            if (!breached) {
                for (int iFee = 0;
                     (iFee < feeDates.size() && (call.isCalled ? getAccrueDate(iFee) <= call.schedule[call.callidx].date : true) ); 
                     iFee++) {
                    paydates.push_back(feeDates[iFee]);
                }
            }
            // if (breached), then all fee dates (+ settlement) up to any breach date,
            // then breach date and so on ...
            // finally: fee dates between last breach date and last fee date 
            // if # breach dates < # tranches
            // recall: we have max(breachDates) < callDate
            if(breached) {
                int iFee = 0;
                for (int jBreachDate=0 ; jBreachDate<numBreachDates; jBreachDate++) {
                    for (; iFee < feeDates.size() && getAccrueDate(iFee) <= breachDates[jBreachDate]; iFee++){
                        paydates.push_back(feeDates[iFee]);
                    }
                    paydates.push_back(instSettle->settles(breachDates[jBreachDate],asset.get()));
                }
                if (breachDates.size() < pcNotionalTranche.size()) {
                    for (; 
                         (iFee < feeDates.size() && (call.isCalled ? getAccrueDate(iFee) <= call.schedule[call.callidx].date : true) );
                         iFee++) {
                        paydates.push_back(feeDates[iFee]);
                    }
                }
            }
            // finally: include call date if isCalled
            if (call.isCalled) {
                paydates.push_back(instSettle->settles(call.schedule[call.callidx].date,asset.get()));
            } 
            
            OutputRequestUtil::recordPaymentDates(control,results,&paydates); 
        }
        
        // II. KNOWN_CASHFLOWS
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArraySP cfl1(new CashFlowArray(0));
            CashFlowArraySP cfl2(new CashFlowArray(0));
            CashFlowArraySP cfl3(new CashFlowArray(0));
            CashFlowArraySP cfl4(new CashFlowArray(0));
            CashFlowArraySP cfl5(new CashFlowArray(0));
            
            // if (!breached), then we know the fees, stop at call date
            if (!breached) {
                for (int i = 0;
                     ( i<feeDates.size() && (call.isCalled ? getAccrueDate(i) <= call.schedule[call.callidx].date : true) );
                     i++) {
                    DateTime pays = feeDates[i];
                    CashFlow cf(pays, -feeToPay(feeDates[i])*notional);
                    cfl1->push_back(cf);
                }
            }
            
            // if breached, then we know the fees and the rebate
            // recall: we have max(breachDates) < callDate
            double pcNotionalRemaining = 1.0;
            if (breached) {
                
                
                double accrued;
                double reb;
                double value;
                
                int iFee = 0;
                int jBreachDate = 0;
                
                for (; jBreachDate<numBreachDates; jBreachDate++) {
                    pcNotionalRemaining -= (jBreachDate != 0) ? pcNotionalTranche[jBreachDate-1] : 0;
                    // cashflow for fee dates
                    for (; iFee < feeDates.size() && getAccrueDate(iFee) <= breachDates[jBreachDate]; iFee++) {
                        DateTime pays = feeDates[iFee];
                        CashFlow cf(pays, -pcNotionalRemaining*feeToPay(feeDates[iFee])*notional);
                        cfl5->push_back(cf);
                    }
                    // cashflow for breach date
                    DateTime pays = instSettle->settles(breachDates[jBreachDate], asset.get());
                    accrued = accruedFee(breachDates[jBreachDate]);
                    if (breachUpOrDown[jBreachDate] == UP) {
                        reb = upRebate->interpolate(breachDates[jBreachDate]);
                    }
                    else {
                        reb = downRebate->interpolate(breachDates[jBreachDate]);
                    }
                    value = pcNotionalTranche[jBreachDate]*(reb - accrued);
                    CashFlow cf(pays, value*notional);
                    cfl2->push_back(cf);
                }
                
                if (breachDates.size() < pcNotionalTranche.size()) {
                    pcNotionalRemaining -= (jBreachDate != 0) ? pcNotionalTranche[jBreachDate-1] : 0;
                    for (; 
                         ( iFee < feeDates.size() && (call.isCalled ? getAccrueDate(iFee) <= call.schedule[call.callidx].date : true) ); 
                         iFee++) {
                        DateTime pays = feeDates[iFee];
                        CashFlow cf(pays, -pcNotionalRemaining*feeToPay(feeDates[iFee])*notional);
                        cfl4->push_back(cf);
                    }
                }
            }
                        
            // if isCalled, pay accrued fee at callDate guaranteed call amount
            if (call.isCalled) {
                DateTime paysCall    = instSettle->settles(call.schedule[call.callidx].date, asset.get());
                double accrued   = accruedFee(call.schedule[call.callidx].date);
                double callAmount= call.schedule[call.callidx].amount;
                double value;
                if (!call.isIssuerCall) {
                    value = -callAmount - accrued;
                }
                else {
                    value = callAmount - accrued;
                }

                CashFlow cf(paysCall, value*pcNotionalRemaining*notional);
                cfl3->push_back(cf);
            }

            // then we merge the 3 sets of cash flows
            cfl1 = CashFlow::merge(cfl1,cfl2);
            cfl1 = CashFlow::merge(cfl1,cfl3);
            cfl1 = CashFlow::merge(cfl1,cfl4);
            cfl1 = CashFlow::merge(cfl1,cfl5);
            
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl1.get());   
        }     
        
        // III. BARRIER_LEVEL
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        bool dead = numBreachDates == pcNotionalTranche.size();
        
        if (request && !fwdStarting && !dead) {
            // report barrier levels over a date range 
            DateTime upperDate;
            
            upperDate = BarrierLevel::barrierWindow(valueDate); 
            if (callSchedule.isCalled){ 
                DateTime callDate = callSchedule.schedule[call.callidx].date;
                if (callDate<BarrierLevel::barrierWindow(valueDate)) {
                    upperDate = callDate;
                }
            }
            
            BarrierLevelArraySP levels(new BarrierLevelArray(0));
            int i;
            
            // lower barrier - use economic barrier (if it exists)
            Schedule* downS = downEcoBarrier.get() ? downEcoBarrier.get(): downBarrier.get();
            CashFlowArraySP downSubset(downS->subset(valueDate, upperDate));
            for (i = 0; i < downSubset->size(); i++) {
                BarrierLevel bl(false,(*downSubset)[i].date,(*downSubset)[i].amount*initialSpot,intraDayMonitor);
                levels->push_back(bl);
            }
            
            // upper barrier - use economic barrier (if it exists)
            Schedule* upS = upEcoBarrier.get() ? upEcoBarrier.get(): upBarrier.get();
            CashFlowArraySP upSubset(upS->subset(valueDate, upperDate));
            for (i = 0; i < upSubset->size(); i++) {
                BarrierLevel bl(true,(*upSubset)[i].date,(*upSubset)[i].amount*initialSpot,intraDayMonitor);
                levels->push_back(bl);
            }
            
            OutputRequestUtil::recordBarrierLevels(control,
                                                   results,
                                                   asset->getTrueName(),
                                                   levels.get());
        }
    }  
    
    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;
    
private:
    // define fee payments
    double        fee;
    string        dcc;
    DateTime      initialAccrueDate;
    DateTimeArray accrueDates;
    DateTimeArray feeDates;
    bool          noAccruedOnBreach;
    bool          useAccrueDates;
    
    CallSchedule  callSchedule;

    // define barriers
    bool            intraDayMonitor;
    bool            isBreachedUp;   // old definition
    bool            isBreachedDown; // old definition
    DateTime        breachDate;     // old definition
    StringArray     breachUpOrDown; // corresponding new definition
    DateTimeArray   breachDates;    // corresponding new definition
    ScheduleSP      upBarrier;
    ScheduleSP      upRebate;
    ScheduleSP      downBarrier;
    ScheduleSP      downRebate;
    ScheduleSP      upEcoBarrier;    // economic barrier-up barrier is for risk
    ScheduleSP      downEcoBarrier;  // economic barrier-down barrier is for risk

    // define the percentage of notional represented by each tranche
    DoubleArray     pcNotionalTranche;

    // internal
    DayCountConventionSP daycount;
    DateTimeArray breachSettleDates;
};  // end of class CallableEDS

/***********************************************************************************************************/
class CEDFDProd: public LatticeProdEDRIns
{
    public:

    CEDFDProd(const CallableEDS* ced, FDModel* model) :
        LatticeProdEDRIns( model, 2, 2 + ced->pcNotionalTranche.size() - ced->breachDates.size() ),
        isCalled(false), inst(ced) 
    {
        if( ! tree1f )
        {
            throw ModelException( "CEDFDProd::CEDFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );
    }
    
    // for mapping step index to notification date index
    typedef map<int, int> NotifDateMap;

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
    {
       return inst->ccyTreatment;
    }
     
    /** calculate at barriers for tree */
    virtual void preCalc(int step); 

    /** product payoff method at maturity */
    void prod_BWD_T(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);
    
    /** premium scaling */
    virtual double scalePremium(double fairValue, 
                                YieldCurveConstSP disc);

    /** extra output requests */    
    void recordOutput(Control* control, 
                      YieldCurveConstSP disc, 
                      Results* results); 

    void AdjustDeltaShift(CControl* control);

    virtual CVolRequestConstSP GetLNRequest() const;

    virtual void update(int& step, FDProduct::UpdateType type);

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing 
    virtual void initProd();    

    virtual void addBarrierCritDates(DateTimeArray  & criticalDates,
                                     const Schedule * barrier) const;

    void compileCriticalDates(DateTimeArray & criticalDates) const;

    void truncateCriticalDates(DateTimeArray & criticalDates) const;

    double getCostOfCall(int step);

    private:
    const CallableEDS*        inst;
    vector<bool>              stepIsKO;
    vector<bool>              stepIsLast;
    vector<bool>              stepIsFee;
    vector<bool>              stepIsNotif;
    vector<double>            stepUpBarrier;
    vector<double>            stepUpRebate;
    vector<double>            stepDnBarrier;
    vector<double>            stepDnRebate;
    double                    refLevel;

    NotifDateMap              notifDateMap;
    DoubleArray               costOfCall;

    bool                      isCalled;
    int                       callidx;
    DateTime                  callDate;

    // barrier at a tree slice after discrete monitor adjustment
    double                    adjUpBarrier[2]; 
    double                    adjDnBarrier[2]; 
};

/** create a fd payoff product */
FDProductSP CallableEDS::createProduct(FDModel* model) const
{
    return FDProductSP( new CEDFDProd(this, model) );
}

/** this sets up the timeline */
void CEDFDProd::init(CControl* control) const 
{
    static const string method = "CEDFDProd::Init";
    try 
    {
        if( tree1f )
        {
            // default to NODE_INSERTION smoothing
            if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) 
            {
                tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
            }

            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode =
                ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? numIns : 0 );

            tree1f->SetDivAmountTreatment(false);

            if ( inst->fwdStarting ) 
                tree1f->controlSameGridFwdStart(inst->ccyTreatment);            
        }        

        // compile list of critical dates
        DateTimeArray critDates;
        compileCriticalDates(critDates);

        // kill off control var - unlikely to be very useful here
        tree1f->DEBUG_UseCtrlVar = false;
           
        // if structure is called, truncate critial dates at call date
        // it does not make sense to simulate until maturity
        if (inst->callSchedule.isCalled) 
        {
            truncateCriticalDates(critDates);
        }

        // add critical dates
        model->addCritDates( critDates );

        // define start and end points of tree
        DateTimeArray segDates;
        segDates.resize(2);

        segDates[0] = getStartDate();

        // if structure is called, simulation is stopped at call date
        if (inst->callSchedule.isCalled) 
        {
            segDates[1] = inst->callSchedule.schedule[inst->callSchedule.callidx].date;
        } 
        else 
        {
            segDates[1] = inst->downBarrier->lastDate();
        }   

        IntArray density( 1, 1 );

        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables */
// this is called per pricing call before tree sweep call (after InitTree)
void CEDFDProd::initProd() 
{
    static const string method = "CEDFDProd::InitProd";
    try 
    {
        int i;
        int lastStep = model->getLastStep();
        
        initSlices( numPrices );
        initInsertNode();

        // set up an array of flags indicating if time step is a KO date
        stepIsKO.resize(lastStep+1);

        // ask tree to decide first about steps that are KO 
        AssetUtil::setStepExercise(stepIsKO,
                                 model->getDates(),
                                 inst->downBarrier, 
                                 true,   // 'american'
                                 inst->asset.getSP());

        // set up an array of flags indicating if the step is the last step of the day
        stepIsLast.resize(lastStep+1);

        for (i = 0; i < lastStep; i++) 
        {
            DateTime stepDate      = model->getDate(i);
            DateTime stepDateAfter = model->getDate(i+1);

            stepIsLast[i] = !stepDateAfter.equals(stepDate,0);
        }

        // set the reference level for the barrier
        refLevel = inst->fwdStarting ? 
            inst->asset->fwdValue(inst->startDate) : inst->initialSpot;

        // set up barrier & rebate for each step
        stepUpBarrier.resize(lastStep+1);
        stepDnBarrier.resize(lastStep+1);
        stepIsFee.resize(lastStep);
        stepIsNotif.resize(lastStep);
        stepUpRebate.resize(lastStep+1);
        stepDnRebate.resize(lastStep+1);

        for (i = 0; i < lastStep; i++) 
        {
            if (stepIsKO[i]) 
            {
                DateTime stepDate = model->getDate(i);

                stepUpRebate[i] =inst->upRebate->interpolate(stepDate);
                stepUpBarrier[i]=inst->upBarrier->interpolate(stepDate)*refLevel;

                stepDnRebate[i] =inst->downRebate->interpolate(stepDate);
                stepDnBarrier[i]=inst->downBarrier->interpolate(stepDate)*refLevel;

                if (stepDnBarrier[i]>stepUpBarrier[i]) 
                {
                    throw ModelException(method,
                        "it must not be the case that the lower barrier is higher than the upper barrier, e.g. " +
                        stepDate.toString() + ".");
                }
            }
            stepIsFee[i]   = false; // for starters
            stepIsNotif[i] = false; // for starters
        }
        // and now the last value (interpolation does not work for i=numSteps)
        if (inst->callSchedule.isCalled) 
        {
            DateTime callDate = inst->callSchedule.schedule[inst->callSchedule.callidx].date;
            stepUpRebate[lastStep] = inst->upRebate->interpolate(callDate);
            stepUpBarrier[lastStep] = inst->upBarrier->interpolate(callDate)*refLevel;
            stepDnRebate[lastStep] = inst->downRebate->interpolate(callDate);
            stepDnBarrier[lastStep] = inst->downBarrier->interpolate(callDate)*refLevel;
        } 
        else if (!inst->callSchedule.isCalled) 
        {
            stepUpRebate[lastStep] = inst->upRebate->lastValue();
            stepUpBarrier[lastStep] = inst->upBarrier->lastValue()*refLevel;
            stepDnRebate[lastStep] = inst->downRebate->lastValue();
            stepDnBarrier[lastStep] = inst->downBarrier->lastValue()*refLevel;
        }

        // and now identify fee dates - skip last one as handled in
        // in maturity payoff, and must occur on last barrier date
        int nextFee = 0;
        int j;
        for (j = 0; j < inst->feeDates.size()-1; j++) 
        {
            for (i = nextFee; i < lastStep; i++) 
            {
                if (inst->feeDates[j] == model->getDate(i)) 
                {
                    stepIsFee[i] = true;
                    nextFee = i+1;
                }
            }
        }    

        // and now call notification - don't care if there's one at maturity
        // as it can't affect the payoff (as you always pay accrued)
        const DateTimeArray& notification = inst->callSchedule.notification;
        int nextNotif = 0;
        for (j = 0; j < notification.size(); j++) 
        {
            for (i = nextNotif; i < lastStep; i++) 
            {
                if (notification[j] == model->getDate(i)) 
                {
                    notifDateMap[i] = j; // map step to notification index
                    stepIsNotif[i] = true;
                    nextNotif = i+1;
                }
            }
        }    

        // store cost of call (i.e. PV(call amount + accrued)) for each
        // notification date
        costOfCall.resize(notification.size());
        const CashFlowArray& schedule = inst->callSchedule.schedule;
        
        for (j = 0; j < notification.size(); j++) 
        {
            const DateTime& callDate = schedule[j].date;
            if (!inst->callSchedule.isIssuerCall) 
            {
                // investor is PUTTING this back
                // pay accrued and put amount
                costOfCall[j] = -schedule[j].amount-inst->accruedFee(callDate) -inst->feeToPay(callDate);
            }
            else 
            {
                // issuer is CALLING it in
                // receive accrued, pay call amount
                costOfCall[j] = schedule[j].amount-inst->accruedFee(callDate)-inst->feeToPay(callDate);;
            }

            costOfCall[j] *= inst->discount->pv(notification[j], callDate);

            // add on PV of any fees due between notif and call
            costOfCall[j] -= inst->dueFee(j);
        }
        if (inst->callSchedule.isCalled) 
        {
            isCalled = true;
            callidx  = inst->callSchedule.callidx;
            callDate = schedule[callidx].date;
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

// here just take care of scaling by notional and any additional 
// discounting for fwd starting case

double CEDFDProd::scalePremium(double fairValue, 
                               YieldCurveConstSP disc)
{
    fairValue -= inst->unsettledCash();
    double fwdStartDF = 1.0;
    if (inst->fwdStarting)
        fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));

    return inst->notional*fwdStartDF*fairValue;
}

/** extra output requests */
void CEDFDProd::recordOutput(Control* control, 
                             YieldCurveConstSP disc, 
                             Results* results)                                 
{
    // get prices at t=0
    // save price
    double price = scalePremium(model->getPrice0( *slices[0] ), disc);
    results->storePrice(price, disc->getCcy());

    // throw this back to the instrument itself
    inst->recordOutputRequests(control, results, price);    
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP CEDFDProd::GetLNRequest() const
{
    CVolRequestConstSP volRequest(inst->makeVolRequest());
    return volRequest;
}


//---------------------------------------------------------------------
// PayoffFD methods
//---------------------------------------------------------------------

/** called before PayoffAtMat, PayoffBeforeMat and tree roll() */
// calculate barrier and place barrier at inserted node if needed
void CEDFDProd::preCalc(int step) 
{
    static const string method("CEDFDProd::preCalc");

    try 
    {    
        int idx = tree1f->getSliceIndex(step);

        adjUpBarrier[idx] = stepUpBarrier[step];
        adjDnBarrier[idx] = stepDnBarrier[step];
        // taking care of discrete monitoring, treat as daily monitor

        if (!inst->intraDayMonitor) 
        {
            vector<double> vol;
            // get vol at barrier
            tree1f->GetStepVol(step, vol, &adjUpBarrier[idx], 0, 0);
            Barrier::BarrierAdjustment(vol[0], true, adjUpBarrier[idx]);
            tree1f->GetStepVol(step, vol, &adjDnBarrier[idx], 0, 0);
            Barrier::BarrierAdjustment(vol[0], false, adjDnBarrier[idx]);
        }

        if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
        { 
            //assuming Lower Barrier is index 0 and Upper uses 1
            // apparently last param = 1 means KO
            tree1f->SetInsertNode(idx, 0, adjDnBarrier[idx], 1); 
            tree1f->SetInsertNode(idx, 1, adjUpBarrier[idx], 1); 
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** product payoff method at maturity */
void CEDFDProd::prod_BWD_T(const TreeSlice & spot,
                                 int step, 
                                 int bot, 
                                 int top, 
                                 int pStart, 
                                 int pEnd,
                                 const vector< TreeSliceSP > & price) 
{
    static const string method("CEDFDProd::prod_BWD_T");
    try 
    {      
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        // if structure is called, then maturity = callDate (simulation stops at callDate)
        // if structure is not called, then maturity = last date on barrier schedule
        DateTime maturity;
        DateTime lastKODate = inst->downBarrier->lastDate();
        
        if (inst->callSchedule.isCalled) 
        {
            maturity = inst->callSchedule.schedule[inst->callSchedule.callidx].date;
        } else 
        {
            maturity = lastKODate;
        }
       
        // settlements
        double settlePV = inst->instSettle->pvAdjust(maturity,inst->discount.get(),inst->asset.get());
        
        // barriers and rebates
        double upKOCM = adjUpBarrier[tree1f->GetSliceIdx()];
        double upKO   = stepUpBarrier[step];   
        double upReb  = stepUpRebate[step];

        double dnKOCM = adjDnBarrier[tree1f->GetSliceIdx()];
        double dnKO   = stepDnBarrier[step];
        double dnReb  = stepDnRebate[step];

        // number of tranches knocked out resp. remaining
        int numTrancheKO = inst->breachDates.size();                         
        
        // remaining notional
        double pcNotionalRemaining = 1.0;
        for (int iTrancheKO = 0; iTrancheKO < numTrancheKO; iTrancheKO++) 
        {
            pcNotionalRemaining -= inst->pcNotionalTranche[iTrancheKO];
        }
        double feeNotional = inst->notionalFee(maturity);

        double pcNotionalRemainingFix = pcNotionalRemaining;
        
        //call amount
        double callAmount = 0.0;
        if (isCalled) 
        {
            callAmount = (inst->callSchedule.isIssuerCall? 1.0:-1.0) *
                inst->callSchedule.schedule[callidx].amount;
        }
        
        // accrued fee on that date
        double accFee = inst->accruedFee(maturity);
        double stepFee = inst->feeToPay(maturity);

        // calculate the value of the fees that have still to be paid
        // if the instrument has not breach
        // i.e. if (fee date[i] >= accrued date[i] > step date)
        int i;
        double extraFees(0.0);
        for (i=0; i<inst->feeDates.size(); i++) 
        {
            if (inst->getAccrueDate(i) > lastKODate && inst->feeDates[i] > lastKODate) 
            {
                extraFees += inst->feeToPay(inst->feeDates[i]) * 
                    inst->discount->pv(lastKODate,inst->feeDates[i]);
            }
        }
        
        double feeNotionalSubstract = 0.0;

        // payoff
        for (int iPrice=pStart; iPrice<=pEnd; iPrice++) 
        {

            // define specific things for specific layer iPrice
            if (iPrice <= pEnd-2 && iPrice > 0 ) 
            {
                pcNotionalRemaining -= (inst->pcNotionalTranche)[iPrice+numTrancheKO-1]; 
                feeNotionalSubstract -= (inst->pcNotionalTranche)[iPrice+numTrancheKO-1];
            }

            // paranoia
            double currentTranche;
            int trancheIdx = iPrice+numTrancheKO;
            if (trancheIdx >= 0 && trancheIdx < inst->pcNotionalTranche.size() && iPrice <=pEnd-2) 
            {
                currentTranche = (inst->pcNotionalTranche)[trancheIdx];
            }
            else 
            {
                currentTranche = pcNotionalRemainingFix;
            }

            // which up and down barriers to use
            double upKObis = (iPrice <= pEnd-1?upKO:upKOCM) - FP_MIN;
            double dnKObis = (iPrice <= pEnd-1?dnKO:dnKOCM) + FP_MIN;

            // notional remaining
            double pcNotRemaining = (iPrice < pEnd-1)? pcNotionalRemaining:pcNotionalRemainingFix;
            feeNotionalSubstract = (iPrice < pEnd-1)?feeNotionalSubstract :0.0;
           

            for (int jNode=bot; jNode<=top; jNode++) 
            {
                // initial step
                (p[iPrice])[jNode] = 0.0;

                // if a tranche is breached
                if (s[jNode]>upKObis || s[jNode]<dnKObis) 
                {
                    // chose up or down rebate
                    double reb = (s[jNode] > upKObis)? upReb : dnReb;
                        
                    // pay the rebate and the accrued fee
                    (p[iPrice])[jNode] += currentTranche*(reb-accFee)*settlePV;
                        
                    // if is called, remove the tranche that has just been breached
                    if (isCalled) 
                    {
                        (p[iPrice])[jNode] -= currentTranche*(callAmount-accFee)*settlePV;
                    }
                }
                // if not breach and if we are on the last KO date
                // pay the fees after the monitoring date
                else if (!inst->callSchedule.isCalled) 
                {
                    (p[iPrice])[jNode] -= extraFees*pcNotRemaining;;  
                }
                
                    
                // client receives (pays) call amount
                if (isCalled) 
                {
                    (p[iPrice])[jNode] += pcNotRemaining*(callAmount-accFee)*settlePV;
                }
                
                // calculate the fees that have to be paid because the accrued period is finished
                // i.e. if (accrued date[i] <= step date <= fee date[i])
                int i;
                double dueFees(0.0);
                for (i=0; i<inst->feeDates.size(); i++) 
                {
                    if (inst->getAccrueDate(i) <= maturity && maturity <= inst->feeDates[i] ) 
                    {
                        dueFees += inst->feeToPay(inst->feeDates[i]) * 
                            inst->discount->pv(maturity,inst->feeDates[i]) * 
                            (inst->notionalFee(inst->feeDates[i])+feeNotionalSubstract);
                    }
                }
                
                // client pay the fees in the future 
                (p[iPrice])[jNode] -= dueFees;
            }
        }
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}
    
/** product payoff method at steps earlier than maturity */
void CEDFDProd::prod_BWD(const TreeSlice & spot,
                               int step, 
                               int bot, 
                               int top, 
                               int pStart, 
                               int pEnd,
                               const vector< TreeSliceSP > & price) 
{
    static const string method("CEDFDProd::prod_BWD");
    try 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        DateTime stepDate = model->getDate(step);

        // define here, assign value after (setp==0 & numTrancheKO>0)
        double settlePV = inst->instSettle->pvAdjust(stepDate,inst->discount.get(),inst->asset.get());
        
        // barriers and rebates
        double upKOCM  = adjUpBarrier[tree1f->GetSliceIdx()];
        double upKO    = stepUpBarrier[step];
        double upReb   = stepUpRebate[step];

        double dnKOCM  = adjDnBarrier[tree1f->GetSliceIdx()];
        double dnKO    = stepDnBarrier[step];
        double dnReb   = stepDnRebate[step]; 

        // number of tranches knocked out
        int numTrancheKO = inst->breachDates.size();                         
                       
        // remaining notional for barrier and for fee to pay 
        // the accrued date may be before a past breach
        // in that case the notional is the one before the breach event
        double pcNotionalRemaining = 1.0;
        for (int iTrancheKO = 0; iTrancheKO < numTrancheKO; iTrancheKO++) 
        {
            pcNotionalRemaining -= inst->pcNotionalTranche[iTrancheKO];
        }
        double feeNotional = inst->notionalFee(stepDate);

        double pcNotionalRemainingFix = pcNotionalRemaining;
        double feeNotionalFix = feeNotional;

        // if entry in "already hit" coincides with "value date", 
        // make sure that hit is not counted twice
        bool alreadyHit = false;
        if (step==0 && numTrancheKO>0 && inst->breachDates[numTrancheKO-1]==inst->valueDate)
        {
            alreadyHit = true; 
        }
        
        // control variate
        if (step==0 && !stepIsKO[step] && !stepIsFee[step] && !stepIsNotif[step]) 
        {
            for (int jNode=bot; jNode<=top; jNode++) 
            {
                (p[0])[jNode] += ((p[pEnd])[jNode] - (p[pEnd-1])[jNode]);
            }
            return;
        }

        // fee and accrual fee
        // note that the full fee is not taking in account on the value date
        double accFee = inst->accruedFee(stepDate);
        double stepFee = inst->valueDate<stepDate?inst->feeToPay(stepDate):0.0;
        
        // payoff
        for (int iPrice=pStart; iPrice<=pEnd; iPrice++) 
        {
            // define specific things for specific layer iPrice
            if (iPrice > pEnd-2) 
            {
                pcNotionalRemaining = pcNotionalRemainingFix;
                feeNotional = feeNotionalFix;
            } 
            else if (iPrice > 0) 
            {
                pcNotionalRemaining -= (inst->pcNotionalTranche)[iPrice+numTrancheKO-1];
                feeNotional -= (inst->pcNotionalTranche)[iPrice+numTrancheKO-1];
            }

            // paranoia
            double currentTranche;
            int trancheIdx = iPrice+numTrancheKO;
            if (trancheIdx >= 0 && trancheIdx < inst->pcNotionalTranche.size() && iPrice <=pEnd-2) {
                currentTranche = (inst->pcNotionalTranche)[trancheIdx];
            }
            else 
            {
                currentTranche = pcNotionalRemainingFix;
            }
          
            // up and down barrier
            double upKObis = (iPrice <= pEnd-1?upKO:upKOCM) - FP_MIN;
            double dnKObis = (iPrice <= pEnd-1?dnKO:dnKOCM) + FP_MIN;

            for (int jNode=bot; jNode<=top; jNode++) 
            {
                
                // for callable
                bool doComparision = true;
          
                double curNotional = pcNotionalRemaining;
                
                // if breach up or down
                if (stepIsKO[step] && (iPrice == pEnd||stepIsLast[step]) && !alreadyHit && 
                    (s[jNode] > upKObis || s[jNode] < dnKObis)) 
                {
                   
                    // down or up rebate rebate
                    double reb = (s[jNode] > upKObis)? upReb : dnReb;
                    
                    // pay the rebate and the accrued fee
                    (p[iPrice])[jNode] = currentTranche*(reb-accFee)*settlePV;
                                                                    
                    curNotional -= currentTranche;
                                            
                    if ((iPrice<pEnd-2 )) 
                    {
                        (p[iPrice])[jNode] += (p[iPrice+1])[jNode];
                    }
                    else 
                    {
                        doComparision  = false;
                    }
                }

                // client pays fee at fee date using the correct notional
                (p[iPrice])[jNode] -= feeNotional*stepFee;
                                       
                // notification
                if (stepIsNotif[step]&&doComparision) 
                {
                    double sign = inst->callSchedule.isIssuerCall? 1:-1;
                
                    (p[iPrice])[jNode] = sign * Maths::min(
                        sign * (p[iPrice])[jNode],
                        sign * getCostOfCall(step) * curNotional); 
                }
            } // end of jNode
        } // end of iPrice
        
        // control variate
        if (step==0) 
        {
            for (int jNode=bot; jNode<=top; jNode++) 
            {
                (p[0])[jNode] += ((p[pEnd])[jNode] - (p[pEnd-1])[jNode]);
            }
            return;
        }            
    } // end of "try"

    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}  

void CEDFDProd::addBarrierCritDates(DateTimeArray  & criticalDates, 
                                    const Schedule * barrier) const 
{
    static const string method("CEDFDProd::addBarrierCritDates");
    try 
    {   
        int i;
        const DateTimeArray& barDates = barrier->getDates();
        if (barrier->getInterp() != Schedule::INTERP_NONE) 
        {
            for (i = 0; i < barDates.size(); i++) 
            {
                criticalDates.push_back(barDates[i]);
            }
        }
        else 
        {
            DateTime maturity = barrier->lastDate();
            if (tree1f->GetStepsPerYear() == 0) 
            {
                tree1f->SetStepsPerYear(-1); // default to daily step
            }
            else if ((tree1f->GetStepsPerYear() < 250 && 
                      tree1f->GetStepsPerYear() > 0))  
            {
                // input is less than one step per day or on Date Case, 
                // add steps around barrier dates
                for (i = 0; i < barrier->length(); i++) 
                {
                    if ((barDates[i] > inst->valueDate) && 
                        (barDates[i] <= maturity)) 
                    {
                        // 2 days gap seem to be best for convergence
                        criticalDates.push_back(barDates[i].rollDate(-2)); 
                        criticalDates.push_back(barDates[i]);
                    }
                }
            }
        } 
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

void CEDFDProd::compileCriticalDates(DateTimeArray & criticalDates) const
{
    static const string method("CEDFDProd::compileCriticalDates");
    try 
    {   
        // compile list of critical dates:
        // fee pay dates
        // barrier dates
        // call notification dates
        // call dates

        criticalDates.resize(inst->feeDates.size());
        int i;
        for (i = 0; i < criticalDates.size() && 
                 inst->feeDates[i]<inst->upBarrier->lastDate(); i++) 
        {
            criticalDates[i] = inst->feeDates[i];
        }
        criticalDates.resize(i);
        
        addBarrierCritDates(criticalDates, inst->downBarrier.get());
        addBarrierCritDates(criticalDates, inst->upBarrier.get());
   
        for (i = 0; i < inst->callSchedule.notification.size(); i++) 
        {
            criticalDates.push_back(inst->callSchedule.notification[i]);
            criticalDates.push_back(inst->callSchedule.schedule[i].date);
        }

        // sort critical dates & make unique
        sort(criticalDates.begin(), criticalDates.end());
        DateTime::removeDuplicates(criticalDates, false);
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

void CEDFDProd::truncateCriticalDates(DateTimeArray & criticalDates) const
{
    static const string method("CEDFDProd::truncateCriticalDates");
    try 
    {
        for (int i=criticalDates.size()-1; i>=0;i--)
        {
            if (criticalDates[i] > 
                inst->callSchedule.schedule[inst->callSchedule.callidx].date) 
            {
                criticalDates.pop_back();
            }
        }        
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

// get the cost of calling for a given timepoint
double CEDFDProd::getCostOfCall(int step) 
{
    NotifDateMap::const_iterator iter = notifDateMap.find(step);
    if (iter == notifDateMap.end()) 
    {
        throw ModelException("CEDFDProd::getCostOfCall",
                             "no notification info found for step " +
                             Format::toString(step));
    }
    return costOfCall[iter->second];
}

void CEDFDProd::update(int & step, 
                       FDProduct::UpdateType type)
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
            prod_BWD_T(*insNodes,
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
        prod_BWD(s,
                 step,
                 bot,
                 top,
                 pStart, 
                  pEnd,
                 price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
        {
            prod_BWD(*insNodes,
                      step,
                      0,
                       tree1f->NumOfInsertNode-1,
                      pStart, 
                      pEnd,
                     *insPrices);
        }
    }   
}

class CallableEDSHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Callable Deposit");
        REGISTER(CallableEDS, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        IMPLEMENTS(LegalTerms::Shift);
        EMPTY_SHELL_METHOD(defaultCallableEDS);
        FIELD(fee, "fee");
        FIELD(dcc, "fee day count");
        FIELD(initialAccrueDate, "initial accrue date");
        FIELD(accrueDates, "accrue period end dates");
        FIELD_MAKE_OPTIONAL(accrueDates);
        FIELD(feeDates, "fee payment dates");
        FIELD(noAccruedOnBreach, "if the accrue fee are not paid on breach");
        FIELD_MAKE_OPTIONAL(noAccruedOnBreach);
        FIELD(useAccrueDates, "use the accrue dates to calculate the fee amount");
        FIELD_MAKE_OPTIONAL(useAccrueDates);
        FIELD(intraDayMonitor, "intra-day barrier");
        FIELD(isBreachedUp, "is up barrier breached?");      // to ensure backwards compatibility
        FIELD_MAKE_OPTIONAL(isBreachedUp);
        FIELD(isBreachedDown, "is down barrier breached?");  // to ensure backwards compatibility
        FIELD_MAKE_OPTIONAL(isBreachedDown);
        FIELD(breachDate, "breach date");                    // to ensure backwards compatibility
        FIELD_MAKE_OPTIONAL(breachDate);                    
        FIELD(breachUpOrDown, "are any barriers breached?");
        FIELD_MAKE_OPTIONAL(breachUpOrDown);
        FIELD(breachDates, "breach dates");
        FIELD_MAKE_OPTIONAL(breachDates);
        FIELD(pcNotionalTranche, "percentage representing each tranche");
        FIELD_MAKE_OPTIONAL(pcNotionalTranche);
        FIELD(upBarrier, "risk up barrier");
        FIELD(upRebate, "up rebate");        
        FIELD(downBarrier, "risk down barrier");
        FIELD(downRebate, "down rebate");        
        FIELD(callSchedule, "call schedule");
        FIELD(upEcoBarrier, "economic up barrier");
        FIELD_MAKE_OPTIONAL(upEcoBarrier);
        FIELD(downEcoBarrier, "economic down barrier");
        FIELD_MAKE_OPTIONAL(downEcoBarrier);
        FIELD_NO_DESC(daycount);
        FIELD_MAKE_TRANSIENT(daycount);
        FIELD_NO_DESC(breachSettleDates);
        FIELD_MAKE_TRANSIENT(breachSettleDates);
    }

    static IObject* defaultCallableEDS(){
        return new CallableEDS();
    }
};

CClassConstSP const CallableEDS::TYPE = CClass::registerClassLoadMethod("CallableEDS", typeid(CallableEDS),
CallableEDSHelper::load);

CClassConstSP const CallableEDS::CallSchedule::TYPE = CClass::registerClassLoadMethod(
    "CallableEDS::CallSchedule", typeid(CallableEDS::CallSchedule), CallableEDS::CallSchedule::load);

const string CallableEDS::UP = "Up";
const string CallableEDS::DOWN = "Down";

/* for class loading */
bool CallableEDSLoad() {
    return (CallableEDS::TYPE != 0);
}


DRLIB_END_NAMESPACE
