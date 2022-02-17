//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquityStabilitySwap.cpp
//
//   Description : Equity Stability Swap - double out version of EDS
//
// Client
//   o Pays a fee at fee dates on the remaining notional of the contract
//   o On knockout, receives rebate and pays fee accrued on the marginal tranche. 
//   o KO is either up or down (i.e double out)
//   o We assume no settlement on pay dates but rebates/accrued are subject to settlement
//   o thus, if we breach on a fee day but before the fee has been paid then the fee and 
//   o the rebate are all paid at settlement   
//  
//   Author      : Andrew J Swain
//
//   Date        : 29 May 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Format.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/VegaParallel.hpp"
#include <map>
#include <algorithm>
#include "edginc/LegalTerms.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class EquityStabilitySwap: public Generic1Factor,
                           virtual public LegalTerms::Shift,
                           virtual public FDModel::IIntoProduct,
                           virtual public LastSensDate {
public:
    static CClassConstSP const TYPE;

    static const string UP;
    static const string DOWN;

    virtual void Validate(){
        static const string method = "EquityStabilitySwap::Validate";
        
        // Call Generic1Factor validate()
        validate();

        if (oneContract) {
            throw ModelException(method, 
                                 "Equity Stability Swap must be "
                                 "notionally booked");
        }

        // sum of tranches is 100%
        double sumNotionalTranche = 0.0;

        for (int iTranche=0; iTranche<pcNotionalTranche.size(); iTranche++) {
            sumNotionalTranche += pcNotionalTranche[iTranche];
        }

        if (!Maths::areEqualWithinTol(sumNotionalTranche, 1.0, 0.000001)) {
            throw ModelException(method, 
                                 "the sum of percentages represented by each tranche "
                                 "has to be 1.0");
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


        // fee dates in order
        DateTime::ensureIncreasing(feeDates, "fee payment dates", true);

        // VALIDATE UP & DOWN BARRIERS
        // last ko/rebate date are co-incident
        DateTime lastKO  = downBarrier->lastDate();
        DateTime lastReb = downRebate->lastDate();

        if (lastReb != lastKO) {
            throw ModelException(method, "last barrier (" + lastKO.toString() + 
                                 ") and last rebate (" + lastReb.toString() + 
                                 ") must co-incide");
        }

        // last Fee must be before last KO and on same day
        DateTime lastFee = feeDates[feeDates.size()-1];
        if (lastFee > lastKO) {
            throw ModelException(method, "last fee (" + lastFee.toString() + 
                                 ") must be on or before last barrier (" + lastKO.toString() + 
                                 ")");
        }
        if (!lastFee.equals(lastKO, false)) {
            throw ModelException(method, "last fee (" + lastFee.toString() + 
                                 ") must be on the same day as the last barrier (" + lastKO.toString() + 
                                 ")");
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

        if (firstKO < initialAccrueDate) {
            throw ModelException(method, 
                                 "first barrier (" + firstKO.toString() + 
                                 ") must be on or after initial "
                                 "accrue date (" +
                                 initialAccrueDate.toString() + ")");          
        }

        if (firstKO != firstReb) {

        throw ModelException(method, 
                                 "first barrier (" + firstKO.toString() + 
                                 ") and first rebate (" +firstReb.toString() + 
                                 ") must co-incide"); 
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


        // is the instrument breached at value date?
        bool breached = (breachDates.size() > 0 || breachUpOrDown.size() > 0);

        if (breached) {
            // breach dates in order
            DateTime::ensureIncreasing(breachDates, "breached dates", true);

            // #breached dates inferior or equal to #tranches of notional
            if (breachDates.size() > pcNotionalTranche.size()) {
                throw ModelException(method, 
                                    "number of breached dates must be inferior "
                                    "or equal to the number of tranches of notional");
            }

            // #breached dates equal to the size of breachUpOrDown
            if (breachDates.size() != breachUpOrDown.size()) {
                throw ModelException(method, 
                                    "number of breached dates must be consistent " 
                                    "with the inputs for up/down");
            }

            // content of breachUpOrDown
            bool val = true;
            for (int breach = 0; breach<breachUpOrDown.size(); breach++) {
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

        if (fwdStarting && startDate > initialAccrueDate) {
            throw ModelException(method, 
                                 "start date (" + startDate.toString() + 
                                 ") is after initial accrue date (" +
                                 initialAccrueDate.toString() + ")"); 
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

        // and finally define remaining variables (not know to the outside)
        daycount = DayCountConventionSP(DayCountConventionFactory::make(dcc));
        if (breachDates.size()>0) {
            breachSettleDates.resize(breachDates.size()-1);
            for (int i=0; i<breachDates.size()-1; i++) {
                breachSettleDates[i] = instSettle->settles(breachDates[i], asset.get());
            }
        }
    }

private:
    friend class EquityStabilitySwapHelper;
    friend class ESSFDProd;

    EquityStabilitySwap():Generic1Factor(TYPE),fee(0.0),intraDayMonitor(false) {};

    EquityStabilitySwap(const EquityStabilitySwap& rhs);
    EquityStabilitySwap& operator=(const EquityStabilitySwap& rhs);

    // what's the PV of any unsettled cash - i.e. fees "paid" but not settled
    // This has been moth-balled now we assume no settlement for fees

//    double unsettledCash() const {
//       static const string method("EquityStabilitySwap::unsettledCash");
//        try {
//            double cash = 0.0;
//            for (int iFee = 0; 
//                 iFee < feeDates.size() && (feeDates[iFee] < valueDate); iFee++) {
//                DateTime pays = instSettle->settles(feeDates[iFee], asset.get());
//                
//                if (pays > valueDate) {
//                    double amount = accruedFee(feeDates[iFee]);
//                    // check whether fee has to be reduced due to breaches
//                    double feeNotional = 1.0;
//                    for (int j=0; j<breachDates.size(); j++) {
//                        if (breachDates[j] <= feeDates[iFee]) {
//                            feeNotional -= pcNotionalTranche[j];
//                        }
//                    }
//                    cash += feeNotional * amount * discount->pv(pays);
//                }
//            }
//
//            return cash;
//        }
//        catch (exception& e) {
//            throw ModelException(e, method);
//        }
//    }
   
    // what's the accrued fee at a given date ?
    double accruedFee(const DateTime& hitDate) const {
        static const string method("EquityStabilitySwap::accruedFee");
        try {
            int      payidx;
            DateTime paydate;
            // find the fee date immediately on or after this date
            bool found = false;
            for (int iFee = 0; iFee < feeDates.size() && !found; iFee++) {
                if (feeDates[iFee] >= hitDate) {
                    found   = true;
                    paydate = feeDates[iFee];
                    payidx  = iFee;
                }
            }
            
            if (!found) {
                // this might actually happen now we allow barrier to go beyond last fee on last day
                // on last day, after the fee has gone, accrued is zero
                return 0.0;
            }

            DateTime lo = payidx == 0 ? initialAccrueDate : feeDates[payidx-1];
            
            double accrued = fee * Maths::max(0.0,daycount->years(lo, hitDate));
             
            return accrued;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Indicates whether VEGA_MATRIX is sensible for this instrument */
    bool avoidVegaMatrix(const IModel*model){
        if (CTree1fLV::TYPE->isInstance(model)) {
            return true; // do pointwise instead
        }
        return false;
    }

    // make a LN vol request - not in real use as prefer LV tree
    CVolRequest* makeVolRequest() const {
        static const string method("EquityStabilitySwap::makeVolRequest");
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

    /** Returns all strikes the EquityStabilitySwap is sensitive to */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method("EquityStabilitySwap::getSensitiveStrikes");
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
        static const string method = "EquityStabilitySwap::sensShift";
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
            
    bool isFeeDate(const DateTime& hitDate) const {
        for (int i = 0; i < feeDates.size(); i++) {
            if (feeDates[i] == hitDate) {
                return true;
            }
        }
        return false;
    }

    bool priceDeadInstrument(Control* control, Results* results) const{
        static const string method = "EquityStabilitySwap::priceDeadInstrument";
        try  {
            bool     deadInstrument = false;
            DateTime end = instSettle->settles(downBarrier->lastDate(),
                                               asset.get());
            DateTime finalDate = downBarrier->lastDate();

            // number of tranches and number of breached dates
            int numBreachDates = breachDates.size();
            int numNotionalTranche = pcNotionalTranche.size();
            
            if (numBreachDates == numNotionalTranche) {
                double value = 0.0;
                for (int iBreachDate = 1; iBreachDate<=numBreachDates; iBreachDate++) {
                    // see if already paid
                    DateTime breachDate = breachDates[numBreachDates-iBreachDate];
                    DateTime pays = instSettle->settles(breachDate, asset.get());
                    if (pays > valueDate) {
                        // get fee accrued - note if breached on fee day but after the fee was paid 
                        // this accrued will be zero which matches fact that only the rebate is due
                        // we just need to handle the nasty case that fee == breach then still no fee 
                        // as it's paid so must kill accrued
                        double accrued = 0.0;
                        if (!isFeeDate(breachDate)) {
                            accrued = accruedFee(breachDate);
                        }
                        // get rebate
                        double reb;
                        if (breachUpOrDown[numBreachDates-iBreachDate] == UP) {
                            reb = upRebate->interpolate(breachDate);
                        }
                        else  {
                            reb = downRebate->interpolate(breachDate);
                        }
                        value += pcNotionalTranche[numNotionalTranche-iBreachDate]
                            *notional*(reb - accrued) * discount->pv(pays);
                    }
                }
                // no old fees to be paid now we have no settlement - so we're done

                results->storePrice(value, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, value);
                }
                              
                deadInstrument = true;
            }
            else if (valueDate >= end) {
                // all over, worth zero
                results->storePrice(0.0, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, 0.0);
                }
                              
                deadInstrument = true;  
            }
            else if (valueDate >= finalDate) {

                // situation: we are at or past last date of barrier schedule (= finalDate),
                // but there are some unsettled payments
                // I. check whether sthg from finalDate is unsettled
                // II. check whether even "older things" are unsettled
                // no fees to consider as they've gone             
                
                double value = 0;
                bool breachScheduleEmpty = breachDates.empty();

                double refLevel = fwdStarting ? asset->fwdValue(startDate) : initialSpot;

                // I. check whether sthg from finalDate is unsettled
                // A) valueDate = finalDate (check whether KO at finalDate) 
                // B) valueDate > finalDate (check whether breach schedule indicates KO on finalDate)
                
                // A) valueDate = finalDate
                if (valueDate == finalDate) {
                    // if KO on finalDate and if last entry in breach dates does not correspond to 
                    // finalDate, then this is an input error
                    if (asset->getSpot() > (upBarrier->lastValue()-FP_MIN)*refLevel) {
                        // check whether already included in breach schedule
                        // if not, then add 1 to numBreachDates
                        if ( breachScheduleEmpty || 
                            (!breachScheduleEmpty && !(breachDates.back()==finalDate)) ) {
                                numBreachDates++;
                        }
                        // client receives rebate only (no fee accrued as it's been paid on the last day)
                        value += pcNotionalTranche[numBreachDates-1]*upRebate->lastValue();
                    } else if (asset->getSpot() < (downBarrier->lastValue()+FP_MIN)*refLevel) {
                        // check whether already included in breach schedule
                        // if not, then add 1 to numBreachDates
                        if ( breachScheduleEmpty || 
                            (!breachScheduleEmpty && !(breachDates.back()==finalDate)) ) {
                                numBreachDates++;
                        }
                        // client receives rebate only (no fee accrued as it's been paid on the last day)
                        value += pcNotionalTranche[numBreachDates-1]*downRebate->lastValue();
                    } else if (!breachScheduleEmpty && (breachDates.back()==finalDate)) {
                        // no fee accrued on a breach as it's been paid already today
                        if (breachUpOrDown.back()==UP) {
                            value += pcNotionalTranche[numBreachDates-1]*upRebate->lastValue();
                        } else if (breachUpOrDown.back()==DOWN) {
                            value += pcNotionalTranche[numBreachDates-1]*downRebate->lastValue();
                        }
                    }
                // B) valueDate > finalDate
                } else if (valueDate>finalDate && !breachScheduleEmpty) {
                    // now, valueDate > finalDate
                    // if sthg happened on finalDate, than it's recorded in the breach schedule
                    if (breachDates.back()==finalDate) {
                        // now find the corresponding rebate and 
                        // client receives rebate (again no accrued as the fee was paid on the final day)
                        if (breachUpOrDown[breachUpOrDown.size()-1]==UP) {
                            value += pcNotionalTranche[numBreachDates-1]*upRebate->lastValue();
                        } else if (breachUpOrDown[breachUpOrDown.size()-1]==DOWN) { 
                            value += pcNotionalTranche[numBreachDates-1]*downRebate->lastValue();
                        }
                    }
                }
                //don't forget notional and discount factor
                value *=notional*discount->pv(end);

                // now, everything concerning things happened on finalDate is done

                // II. check whether even "older things" are unsettled, i.e. 
                // if KO, then client still receives "old" rebates and possibly "old" accrued fees 
                // attention: everything on finalDate is already done
                // no other fees to consider as they've gone             
                if (!breachScheduleEmpty) {
                    for (int i=0; i<breachDates.size() && breachDates[i]<finalDate
                        && valueDate < instSettle->settles(breachDates[i],asset.get()); i++) {
                        // client receives rebate and pays fee accrued if breach date wasn't a fee date
                        DateTime settle = instSettle->settles(breachDates[i],asset.get());
                        if (breachUpOrDown[i]==UP) {
                            value += pcNotionalTranche[i]
                                *upRebate->interpolate(breachDates[i])*notional*discount->pv(settle);
                        } else if (breachUpOrDown[i]==DOWN) {
                            value += pcNotionalTranche[i]
                                *downRebate->interpolate(breachDates[i])*notional*discount->pv(settle);
                        }
                        // only get accrued fee if it wasn't already paid (note if breached same day but after fee accrued = 0)
                        if (!isFeeDate(breachDates[i])) {
                            value -= pcNotionalTranche[i] * accruedFee(breachDates[i])*notional*discount->pv(settle);
                        }
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
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       downBarrier->lastDate(),
                                       valueDate,
                                       asset.get());

        bool breached = (breachDates.size() > 0);
        int numBreachDates = breachDates.size();
        OutputRequest* request = NULL;

        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArray paydates;
            
            if (!breached) {
                // all fee dates
                for (int iFee = 0; iFee < feeDates.size(); iFee++) {
                    paydates.push_back(feeDates[iFee]);
                }
            }
            if (breached) {
                // all fee dates up to any breach date then the breach date and so on
                // fee dates between last breach date and last fee date if #breach dates < #tranches of notional
                int iFee = 0;

                for (int jBreachDate=0 ; jBreachDate<numBreachDates; jBreachDate++) {
                    for (; iFee < feeDates.size() && feeDates[iFee] <= breachDates[jBreachDate]; iFee++) {
                        paydates.push_back(feeDates[iFee]);
                    }
                    paydates.push_back(instSettle->settles(breachDates[jBreachDate],asset.get()));
                }
                if (breachDates.size() < pcNotionalTranche.size()) {
                    for (; iFee < feeDates.size(); iFee++) {
                        paydates.push_back(feeDates[iFee]);
                    }
                }
            }

            OutputRequestUtil::recordPaymentDates(control,results,&paydates); 
        }

        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArray cfl;
            
            if (!breached) {
                // if not KO'd then we know the fees
                for (int iFee = 0; iFee<feeDates.size(); iFee++) {
                    CashFlow cf(feeDates[iFee], -accruedFee(feeDates[iFee])*notional);
                    cfl.push_back(cf);
                }
            }

            if (breached) {
                double pcNotionalRemaining = 1.0;

                double accrued;
                double reb;
                double value;

                int iFee = 0;
                int jBreachDate = 0;

                for (; jBreachDate<numBreachDates; jBreachDate++) {
                    pcNotionalRemaining -= (jBreachDate != 0) ? pcNotionalTranche[jBreachDate-1] : 0;
                    // cash flow for fee dates
                    for (; iFee < feeDates.size() && feeDates[iFee] <= breachDates[jBreachDate]; iFee++) {
                        CashFlow cf(feeDates[iFee], -pcNotionalRemaining*accruedFee(feeDates[iFee])*notional);
                        cfl.push_back(cf);
                    }
                    // cash flow for breach date
                    DateTime pays = instSettle->settles(breachDates[jBreachDate], asset.get());
                    accrued = 0.0;
                    if (!isFeeDate(breachDates[jBreachDate])) {
                        accrued = accruedFee(breachDates[jBreachDate]);
                    }
                    if (breachUpOrDown[jBreachDate] == UP) {
                        reb = upRebate->interpolate(breachDates[jBreachDate]);
                    }
                    else {
                        reb = downRebate->interpolate(breachDates[jBreachDate]);
                    }
                    value = pcNotionalTranche[jBreachDate]*(reb - accrued);
                    CashFlow cf(pays, value*notional);
                    cfl.push_back(cf);
                }
                
                if (breachDates.size() < pcNotionalTranche.size()) {
                    pcNotionalRemaining -= (jBreachDate != 0) ? pcNotionalTranche[jBreachDate-1] : 0;
                    for (; iFee < feeDates.size(); iFee++) {
                        CashFlow cf(feeDates[iFee], -pcNotionalRemaining*accruedFee(feeDates[iFee])*notional);
                        cfl.push_back(cf);
                    }
                }
            }

            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    &cfl);   
        }
        
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        bool dead = numBreachDates == pcNotionalTranche.size();
        if (request && !fwdStarting && !dead) {
            // report barrier levels over a date range
            DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
            BarrierLevelArraySP levels(new BarrierLevelArray(0));
            // use economic barrier (if it exists)
            Schedule* s = downEcoBarrier.get() ?
                downEcoBarrier.get(): downBarrier.get();

            CashFlowArraySP subset(s->subset(valueDate, upperDate));
            if (!subset->empty()) {
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(false,(*subset)[i].date,
                                    (*subset)[i].amount*initialSpot, intraDayMonitor);
                    levels->push_back(bl);
                }
            }

            s = upEcoBarrier.get() ? upEcoBarrier.get(): upBarrier.get();
            subset = CashFlowArraySP(s->subset(valueDate, upperDate));
            if (!subset->empty()) {
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(true,(*subset)[i].date,
                                    (*subset)[i].amount*initialSpot, intraDayMonitor);
                    levels->push_back(bl);
                }
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
    DateTimeArray feeDates;

    // define barriers
    bool          intraDayMonitor;
    DateTimeArray breachDates; // breach dates for risk             
    StringArray   breachUpOrDown; // is breach up or is beach down for risk
    ScheduleSP    upBarrier;
    ScheduleSP    upRebate;
    ScheduleSP    downBarrier;
    ScheduleSP    downRebate;
    ScheduleSP    upEcoBarrier; // economic barrier-up barrier is for risk
    ScheduleSP    downEcoBarrier; // economic barrier-down barrier is for risk

    // define the percentage of notional represented by each tranche
    DoubleArray   pcNotionalTranche; // percentage of notional represented by each tranche

    // internal
    DayCountConventionSP daycount;
    DateTimeArray  breachSettleDates;

};


/**********************************************************************************************************************/

class ESSFDProd: public LatticeProdEDRIns 
{
public:
    ESSFDProd(const EquityStabilitySwap* ess, FDModel *model) :
        LatticeProdEDRIns(model, 2, 2 + ess->pcNotionalTranche.size() - ess->breachDates.size()),
        inst(ess)
    {
        if( ! tree1f )
        {
            throw ModelException( "ESSFDProd::ESSFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );
    }
    void init(CControl* control) const;

    void initProd();

    double scalePremium(const double& fairValue,
                        YieldCurveConstSP disc);

    virtual string getCcyTreatment() const
    {
        return inst->ccyTreatment;
    }

    void recordOutput(Control* control, 
                      YieldCurveConstSP disc, 
                      Results* results);

    /** returns a vol request for log-normal vol */
    virtual CVolRequestConstSP GetLNRequest() const
    {
        CVolRequestConstSP volRequest(inst->makeVolRequest());
        return volRequest;
    }

    void preCalc(int step);

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

    virtual bool Positive() 
    {        
        return false;
    }

    void update(int & step, 
                FDProduct::UpdateType type);

    void addBarrierCritDates(DateTimeArray  &  criticalDates, 
                             const Schedule * barrier) const;

    void compileCriticalDates(DateTimeArray & criticalDates) const;

    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

private:
    const EquityStabilitySwap*  inst;
    vector<bool>              stepIsKO;
    vector<bool>              stepIsLast;
    vector<bool>              stepIsFee;
    vector<double>            stepUpBarrier;
    vector<double>            stepUpRebate;
    vector<double>            stepDnBarrier;
    vector<double>            stepDnRebate;
    double                    refLevel;
    // barrier at a tree slice after discrete monitor adjustment
    double                    adjUpBarrier[2]; 
    double                    adjDnBarrier[2]; 
};

FDProductSP EquityStabilitySwap::createProduct(FDModel* model) const
{
    return FDProductSP( new ESSFDProd(this, model) );
}

/** this sets up the timeline */
void ESSFDProd::init(CControl* control) const
{    
    static const string method = "ESSFDProd::Init";
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
        
        // define start and end points of tree
        DateTimeArray segDates;
        segDates.resize(2);
        
        segDates[0] = getStartDate();
        segDates[1] = inst->downBarrier->lastDate();

        // timeline configuration
        // 'density factor' for timeline
        IntArray density( 1, 1 );

        // compile list of critical dates
        DateTimeArray critDates;
        compileCriticalDates(critDates);

        // kill off control var - unlikely to be very useful here
        tree1f->DEBUG_UseCtrlVar = false;

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

/** initialising and setting product variables */
// this is called per pricing call before tree sweep call (after InitTree)
void ESSFDProd::initProd() 
{
    static const string method = "ESSFDProd::InitProd";
    try 
    {
        int iStep;
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

        for (iStep = 0; iStep < lastStep; iStep++) 
        {
            DateTime stepDate      = model->getDate(iStep);
            DateTime stepDateAfter = model->getDate(iStep+1);

            stepIsLast[iStep] = !stepDateAfter.equals(stepDate,0);
        }
        
        // set the reference level for the barrier
        refLevel = inst->fwdStarting ? 
            inst->asset->fwdValue(inst->startDate) : inst->initialSpot;

        // set up barrier & rebate for each step
        stepUpBarrier.resize(lastStep+1);
        stepDnBarrier.resize(lastStep+1);
        stepIsFee.resize(lastStep+1);
        stepUpRebate.resize(lastStep);
        stepDnRebate.resize(lastStep);

        stepUpBarrier[lastStep] = inst->upBarrier->lastValue()*refLevel;
        stepDnBarrier[lastStep] = inst->downBarrier->lastValue()*refLevel;
        int numFees = inst->feeDates.size();
        stepIsFee[lastStep] = false; // for starters
        if (inst->feeDates[numFees-1] == inst->upBarrier->lastDate()) 
        {
            stepIsFee[lastStep] = true;
        }

        for (iStep = 0; iStep < lastStep; iStep++) 
        {
            if (stepIsKO[iStep]) 
            {
                DateTime stepDate = model->getDate(iStep);

                stepUpRebate[iStep] =inst->upRebate->interpolate(stepDate);
                stepUpBarrier[iStep]=inst->upBarrier->interpolate(stepDate)*refLevel;

                stepDnRebate[iStep] =inst->downRebate->interpolate(stepDate);
                stepDnBarrier[iStep]=inst->downBarrier->interpolate(stepDate)*refLevel;
            }
            stepIsFee[iStep] = false; // for starters
        }
        

        // and now identify fee dates - note last one needn't be last step
        int nextFee = 0;
        for (int jFee = 0; jFee<numFees; jFee++) 
        {
            for (iStep = nextFee; iStep < lastStep; iStep++) 
            {
                if (inst->feeDates[jFee] == model->getDate(iStep)) 
                {
                    stepIsFee[iStep] = true;
                    nextFee = iStep+1;
                }
            }
        }                    
        // just in case we're exactly on a fee
        // stop initial point being a fee date - the fee is deemed to have gone
        stepIsFee[0] = false;
   }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

// here just take care of scaling by notional and any additional 
// discounting for fwd starting case
double ESSFDProd::scalePremium(const double& fairValue,
                               YieldCurveConstSP disc)
{
    double fwdStartDF = 1.0;
    if (inst->fwdStarting)
        fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));

    return inst->notional*fwdStartDF*fairValue;
}   

void ESSFDProd::recordOutput(Control* control, 
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

/** called before PayoffAtMat, PayoffBeforeMat and tree roll() */
// calculate barrier and place barrier at inserted node if needed
void ESSFDProd::preCalc(int step) 
{
    static const string method("ESSFDProd::preCalc");

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
void ESSFDProd::prod_BWD_T(const TreeSlice & spot,
                                 int step, 
                                 int bot, 
                                 int top, 
                                 int pStart, 
                                 int pEnd,
                                 const vector< TreeSliceSP > & price)  
{
    static const string method("ESSFDProd::prod_BWD_T");
    try 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        double settlePV = inst->instSettle->pvAdjust(inst->downBarrier->lastDate(),
                                                     inst->discount.get(), 
                                                     inst->asset.get());

        // what's the accrued fee so far ?
        double accrued = inst->accruedFee(inst->downBarrier->lastDate());

        // barriers for continuous/daily monitoring and rebates
        double upKOCM = adjUpBarrier[tree1f->GetSliceIdx()];
        double upKO   = stepUpBarrier[step];   
        double upRebAtMat  = inst->upRebate->lastValue() * settlePV;

        double dnKOCM = adjDnBarrier[tree1f->GetSliceIdx()];
        double dnKO   = stepDnBarrier[step];
        double dnRebAtMat  = inst->downRebate->lastValue() * settlePV;
        
        // notional remaining at stake (if breached dates)
        int numTrancheKO = inst->breachDates.size();
        double pcNotionalRemaining = 1.0;

        for (int iTrancheKO = 0; iTrancheKO < numTrancheKO; iTrancheKO++) 
        {
            pcNotionalRemaining -= inst->pcNotionalTranche[iTrancheKO];
        }

        double pcNotionalRemainingCV = pcNotionalRemaining;

       // pay-off
        int iPrice = pStart; 
        for (; iPrice<=pEnd; iPrice++) 
        {
            if (iPrice <= pEnd-2) 
            {
                pcNotionalRemaining -= (iPrice != 0) ? (inst->pcNotionalTranche)[iPrice-1+numTrancheKO] : 0;
            }
            for (int jNode=bot; jNode<=top; jNode++) 
            {
                (p[iPrice])[jNode] = 0.0;

                // deal with fee if it's due
                if (stepIsFee[step])
                {
                    if (iPrice <= pEnd-2) 
                    {
                        (p[iPrice])[jNode] = -pcNotionalRemaining*accrued;
                    }
                    else 
                    {
                        (p[iPrice])[jNode] = -pcNotionalRemainingCV*accrued;
                    }
                }

                // pay-off for product with several tranches : daily monitoring without barrier adjustment
                if (iPrice <= pEnd-2) 
                {
                    // if breached up 
                    // -> pay rebate for the tranche
                    // -> don't bother with accrued - fee was paid some time today
                    if (s[jNode] > upKO-FP_MIN) 
                    {
                        (p[iPrice])[jNode] += (inst->pcNotionalTranche)[iPrice+numTrancheKO]*upRebAtMat; 
                    }

                    // if breached down
                    // -> pay rebate for the tranche
                    // -> don't bother with accrued - fee was paid some time today
                    if (s[jNode] < dnKO+FP_MIN) 
                    {
                        (p[iPrice])[jNode] += (inst->pcNotionalTranche)[iPrice+numTrancheKO]*dnRebAtMat; 
                    }
                }

                // pay-off for product with one tranche : daily monitoring without barrier adjustment
                else if (iPrice == pEnd-1) 
                {
                    // don't bother with accrued - fee was paid some time today
                    if (s[jNode] > upKO-FP_MIN) 
                    {
                        (p[iPrice])[jNode] += pcNotionalRemainingCV*upRebAtMat; 
                    }
                    if (s[jNode] < dnKO+FP_MIN) 
                    {
                        (p[iPrice])[jNode] += pcNotionalRemainingCV*dnRebAtMat; 
                    }
                }

                // pay-off for product with one tranche : continuous monitoring with barrier adjustment
                else 
                {
                    // don't bother with accrued - fee was paid some time today
                    if (s[jNode] > upKOCM-FP_MIN) 
                    {
                        (p[iPrice])[jNode] += pcNotionalRemainingCV*upRebAtMat; 
                    }
                    if (s[jNode] < dnKOCM+FP_MIN) 
                    {
                        (p[iPrice])[jNode] += pcNotionalRemainingCV*dnRebAtMat; 
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

/** product payoff method at steps earlier than maturity */
void ESSFDProd::prod_BWD(const TreeSlice & spot,
                               int step, 
                               int bot, 
                               int top, 
                               int pStart, 
                               int pEnd,
                               const vector< TreeSliceSP > & price)  
{
    static const string method("ESSFDProd::prod_BWD");
    try 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        DateTime stepDate = model->getDate(step);

        double settlePV;
        
        // what's the accrued fee so far ?
        double accrued;

        // settlements which might happen after the value date
        int numTrancheKO = (inst->breachDates).size();
        
        // if entry in "already hit" coincides with "value date", 
        // make sure that hit is not counted twice
        bool alreadyHit = false;
        if (step==0 && numTrancheKO>0 && inst->breachDates[numTrancheKO-1]==inst->valueDate)
        {
            alreadyHit = true; 
        }
                    
        if (step == 0 && numTrancheKO > 0) 
        {
            // deal with unsettled rebates/accrued for breaches that have already happened
            for (int dTrancheKO=0; dTrancheKO<numTrancheKO; dTrancheKO++) 
            {
                if (inst->instSettle->settles(inst->breachDates[dTrancheKO], inst->asset.get()) > inst->valueDate) 
                {
                    accrued = 0.0;
                    // if breach was a fee date and breach < fee accrued = 0
                    // if breach was a fee date and breach > want to pay accrued at settlement
                    // if breach was a fee date and breach = fee then fee has already gone - so we kill this case
                    if (!inst->isFeeDate(inst->breachDates[dTrancheKO])) 
                    {
                        accrued = inst->accruedFee(inst->breachDates[dTrancheKO]);
                    }
                    // pv for low date = value date and hi date = breachdate
                    settlePV = inst->instSettle->pv(inst->valueDate, inst->instSettle->settles(inst->breachDates[dTrancheKO], inst->asset.get()),
                        inst->discount.get(), inst->asset.get());
                    // Up or Down rebate at breach date 
                    double stepUpOrDownRebate = (inst->breachUpOrDown[dTrancheKO] == EquityStabilitySwap::UP) ? 
                        inst->upRebate->interpolate(inst->breachDates[dTrancheKO]) :
                        inst->downRebate->interpolate(inst->breachDates[dTrancheKO]);

                    for (int iPrice = pStart; iPrice<=pEnd; iPrice++) 
                    {
                        for (int jNode=bot; jNode<=top; jNode++) 
                        {
                            (p[iPrice])[jNode] += inst->pcNotionalTranche[dTrancheKO]*(stepUpOrDownRebate - accrued)*settlePV;
                        }
                    }
                }
            }
        }

        if (!stepIsKO[step] && !stepIsFee[step]) 
        {
            if (step == 0)
            {
                for (int jNode=bot; jNode<=top; jNode++) 
                {
                    (p[0])[jNode] += ((p[pEnd])[jNode] - (p[pEnd-1])[jNode]);
                }
            }
            return;
        }
           
        // barriers for continuous/daily monitoring
        double upKOCM = adjUpBarrier[tree1f->GetSliceIdx()];
        double upKO   = stepUpBarrier[step];   
        
        double dnKOCM = adjDnBarrier[tree1f->GetSliceIdx()];
        double dnKO   = stepDnBarrier[step];
        
        // notional remaining at stake (if breached dates)
        double pcNotionalRemaining = 1.0;

        for (int dTrancheKO = 0; dTrancheKO<numTrancheKO; dTrancheKO++) 
        {
            pcNotionalRemaining -= inst->pcNotionalTranche[dTrancheKO];
        }

        double pcNotionalRemainingCV = pcNotionalRemaining;

        settlePV = inst->instSettle->pvAdjust(stepDate,
                                              inst->discount.get(), 
                                              inst->asset.get());
        accrued = inst->accruedFee(stepDate);

        // pay-off
        int iPrice = pStart; 
        for (; iPrice<=pEnd; iPrice++) 
        {
            if (iPrice <= pEnd-2) 
            {
                pcNotionalRemaining -= (iPrice != 0) ? (inst->pcNotionalTranche)[iPrice-1+numTrancheKO] : 0;
            }
            for (int jNode=bot; jNode<=top; jNode++) 
            {
                bool weKOHere = false;
                // pay-off for product with several tranches : daily monitoring without barrier adjustment
                if (iPrice <= pEnd-2) {
                    // if it's a KO date and state variable < number of tranches
                    if (stepIsKO[step] && stepIsLast[step] && iPrice != pEnd-2 && !alreadyHit) 
                    {
                        if (s[jNode] > upKO-FP_MIN) 
                        {
                            (p[iPrice])[jNode] = (p[iPrice+1])[jNode] 
                                + inst->pcNotionalTranche[iPrice+numTrancheKO]*stepUpRebate[step]*settlePV;
                            weKOHere = true;
                        }
                        else if (s[jNode] < dnKO+FP_MIN) {
                            (p[iPrice])[jNode] = (p[iPrice+1])[jNode] 
                                + inst->pcNotionalTranche[iPrice+numTrancheKO]*stepDnRebate[step]*settlePV;                                    
                            weKOHere = true;
                        }
                    }
                    // if it's a KO date and state variable = number of tranches
                    if (stepIsKO[step] && stepIsLast[step] && iPrice == pEnd-2 && !alreadyHit) 
                    {
                        if (s[jNode] > upKO-FP_MIN) 
                        {
                            (p[iPrice])[jNode] = 
                                inst->pcNotionalTranche[iPrice+numTrancheKO]*stepUpRebate[step]*settlePV;
                            weKOHere = true;
                        }
                        else if (s[jNode] < dnKO+FP_MIN) 
                        {
                            (p[iPrice])[jNode] = 
                                inst->pcNotionalTranche[iPrice+numTrancheKO]*stepDnRebate[step]*settlePV;
                            weKOHere = true;
                        }
                    }
                     // if it's a fee date, pay the fee regardlinst of KO - no settlement
                    if (stepIsFee[step]) 
                    {
                        (p[iPrice])[jNode] -= pcNotionalRemaining*accrued;
                    }
                    // if we breached only pay accrued on this tranche if it's not a fee date
                    // fee dates are handled below and if we're same day but after a fee accrued will be zero
                    else if (weKOHere) {
                        (p[iPrice])[jNode] -= inst->pcNotionalTranche[iPrice+numTrancheKO]*accrued*settlePV;
                    }
                }

                // pay-off for product with one tranche : daily monitoring without barrier adjustment
                else if (iPrice == pEnd-1) 
                {
                    if (stepIsKO[step] && stepIsLast[step]) 
                    {
                        if (s[jNode] > upKO-FP_MIN) 
                        {
                            (p[iPrice])[jNode] = pcNotionalRemainingCV*stepUpRebate[step]*settlePV;
                            weKOHere = true;
                        }
                        else if (s[jNode] < dnKO+FP_MIN) 
                        {
                            (p[iPrice])[jNode] = pcNotionalRemainingCV*stepDnRebate[step]*settlePV;
                            weKOHere = true;
                        }
                    }
                    if (stepIsFee[step]) 
                    {// no settlement
                        (p[iPrice])[jNode] -= pcNotionalRemainingCV*accrued;
                    }
                    // if we breached only pay accrued on this tranche if it's not a fee date
                    // fee dates are handled below and if we're same day but after a fee accrued will be zero
                    else if (weKOHere) {
                        (p[iPrice])[jNode] -= pcNotionalRemainingCV*accrued*settlePV;
                    }
                }

                // pay-off for product with one tranche : continuous monitoring with barrier adjustment
                else 
                {
                    if (stepIsKO[step]) 
                    {
                        if (s[jNode] > upKOCM-FP_MIN) 
                        {
                            (p[iPrice])[jNode] = pcNotionalRemainingCV*stepUpRebate[step]*settlePV;
                            weKOHere = true;
                        }
                        else if (s[jNode] < dnKOCM+FP_MIN) 
                        {
                            (p[iPrice])[jNode] = pcNotionalRemainingCV*stepDnRebate[step]*settlePV;
                            weKOHere = true;
                        }
                    }
                    if (stepIsFee[step]) 
                    { // no settlement
                        (p[iPrice])[jNode] -= pcNotionalRemainingCV*accrued;
                    }
                    // if we breached only pay accrued on this tranche if it's not a fee date
                    // fee dates are handled below and if we're same day but after a fee accrued will be zero
                    else if (weKOHere) {
                        (p[iPrice])[jNode] -= pcNotionalRemainingCV*accrued*settlePV;
                    }
                }
            }
        }

        // control variate
        if (step == 0)
        {
            for (int jNode=bot; jNode<=top; jNode++) 
            {
                (p[0])[jNode] += ((p[pEnd])[jNode] - (p[pEnd-1])[jNode]);
            }
            return;
        }
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

void ESSFDProd::update(int & step, 
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

void ESSFDProd::addBarrierCritDates(DateTimeArray  &  criticalDates, 
                                    const Schedule * barrier) const 
{
    static const string method("ESSFDProd::addBarrierCritDates");
    try 
    {   
        int iBarDates;
        const DateTimeArray& barDates = barrier->getDates();
        if (barrier->getInterp() != Schedule::INTERP_NONE) 
        {
            for (iBarDates = 0; iBarDates < barDates.size(); iBarDates++) 
            {
                criticalDates.push_back(barDates[iBarDates]);
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
                // input is linst than one step per day or on Date Case, 
                // add steps around barrier dates
                for (iBarDates = 0; iBarDates < barrier->length(); iBarDates++) 
                {
                    if ((barDates[iBarDates] > inst->valueDate) && 
                        (barDates[iBarDates] <= maturity)) 
                    {
                        // 2 days gap seem to be best for convergence
                        criticalDates.push_back(barDates[iBarDates].rollDate(-2)); 
                        criticalDates.push_back(barDates[iBarDates]);
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

void ESSFDProd::compileCriticalDates(DateTimeArray & criticalDates) const 
{
    static const string method("ESSFDProd::compileCriticalDates");
    try 
    {   
        // compile list of critical dates:
        // fee pay dates
        // barrier dates

        criticalDates.resize(inst->feeDates.size());

        int iCritDates;
        for (iCritDates = 0; iCritDates < criticalDates.size(); iCritDates++)
        {
            criticalDates[iCritDates] = inst->feeDates[iCritDates];
        }
        
        addBarrierCritDates(criticalDates, inst->downBarrier.get());
        addBarrierCritDates(criticalDates, inst->upBarrier.get());
   
        // sort critical dates & make unique
        sort(criticalDates.begin(), criticalDates.end());
        DateTime::removeDuplicates(criticalDates, false);
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

class EquityStabilitySwapHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Equity Stability Swap");
        REGISTER(EquityStabilitySwap, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        IMPLEMENTS(LegalTerms::Shift);
        EMPTY_SHELL_METHOD(defaultEquityStabilitySwap);
        FIELD(fee, "fee");
        FIELD(dcc, "fee day count");
        FIELD(initialAccrueDate, "initial accrue date");
        FIELD(feeDates, "fee payment dates");
        FIELD(intraDayMonitor, "intra-day barrier");
        FIELD(breachDates, "breach dates");
        FIELD(breachUpOrDown, "is up barrier breached or is down barrier breached");
        FIELD(upBarrier, "risk up barrier");
        FIELD(upRebate, "up rebate");        
        FIELD(downBarrier, "risk down barrier");
        FIELD(downRebate, "down rebate");        
        FIELD(upEcoBarrier, "economic up barrier");
        FIELD_MAKE_OPTIONAL(upEcoBarrier);
        FIELD(downEcoBarrier, "economic down barrier");
        FIELD_MAKE_OPTIONAL(downEcoBarrier);
        FIELD(pcNotionalTranche, "percentage of notional for each tranche");
        FIELD_NO_DESC(daycount);
        FIELD_MAKE_TRANSIENT(daycount);
        FIELD_NO_DESC(breachSettleDates);
        FIELD_MAKE_TRANSIENT(breachSettleDates);
    }

    static IObject* defaultEquityStabilitySwap(){
        return new EquityStabilitySwap();
    }
};

CClassConstSP const EquityStabilitySwap::TYPE = CClass::registerClassLoadMethod("EquityStabilitySwap", typeid(EquityStabilitySwap),
EquityStabilitySwapHelper::load);

const string EquityStabilitySwap::UP   = "Up";
const string EquityStabilitySwap::DOWN = "Down";

/* for class loading */
bool EquityStabilitySwapLoad() {
    return (EquityStabilitySwap::TYPE != 0);
}


DRLIB_END_NAMESPACE

