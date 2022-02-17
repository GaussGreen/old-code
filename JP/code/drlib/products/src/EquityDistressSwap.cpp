//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquityDistressSwap.cpp
//
//   Description : Equity Distress Swap (generic 1 factor)
//
// Client 
//   o PAYS a fee until earliest of maturity/knockout.
//   o On knockout, RECEIVES rebate, pays fee accrued. 
//   o We assume no settlement on pay dates but rebates/accrued are subject to settlement
//   o thus, if we breach on a fee day but before the fee has been paid then the fee and 
//   o the rebate are all paid at settlement   
// 
//   Author      : Andrew J Swain
//
//   Date        : 18 March 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/Delta.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/AssetUtil.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class EquityDistressSwap;
typedef smartPtr<EquityDistressSwap> EquityDistressSwapSP;

class EquityDistressSwap: public Generic1Factor,
                          virtual public LegalTerms::Shift,
                          virtual public FDModel::IIntoProduct, 
                          virtual public LastSensDate, 
                          virtual public BarrierBreach::IEventHandler,
                          virtual public KnownCashflows::IEventHandler {
public:
    static CClassConstSP const TYPE;  

    virtual void Validate(){
        static const string method = "EquityDistressSwap::Validate";
        
        // Call Generic1Factor validate()
        validate();

        if (oneContract) {
            throw ModelException(method, 
                                 "Equity Distress Swap must be "
                                 "notionally booked");
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

        if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
            throw ModelException(method, "ccy struck not supported");
        }
          
        if (barrier->getInterp() == Schedule::INTERP_NONE) {
            throw ModelException(method, "'dates only' barrier not supported");
        }
        if (rebate->getInterp() == Schedule::INTERP_NONE) {
            throw ModelException(method, "'dates only' rebate not supported");
        }
             

        // fee dates in order
        DateTime::ensureIncreasing(feeDates, "fee payment dates", true);

        // first/last ko/rebate date are co-incident
        DateTime lastKO  = barrier->lastDate();
        DateTime lastReb = rebate->lastDate();
        DateTime firstKO  = barrier->firstDate();
        DateTime firstReb = rebate->firstDate();
       
        if (lastReb != lastKO) {
            throw ModelException(method, "last barrier (" + lastKO.toString() + 
                                 ") and last rebate (" + lastReb.toString() + 
                                 ") must co-incide");
        }

        if (firstKO != firstReb) {
            throw ModelException(method, 
                                 "first barrier (" + firstKO.toString() + 
                                 ") and first rebate (" +firstReb.toString() + 
                                 ") must co-incide"); 
        }

        // no breach in the future
        if (fwdStarting && isBreached) {
            throw ModelException(method, 
                                 "instrument is fwd starting and breached");
        }

        // again no breach in the future
        if (isBreached && breachDate > valueDate) {
            throw ModelException(method, "instrument is knocked out but ko "
                                 "date (" + breachDate.toString() + ") is "
                                 "after today (" + valueDate.toString() + ")");
        }          

        // if there's an economic barrier, check it lines up with risk barrier
        if (ecoBarrier.get()) {
            DateTime firstEco = ecoBarrier->firstDate();
            DateTime lastEco  = ecoBarrier->lastDate();
            if (firstEco != firstKO || lastEco != lastKO) {
                throw ModelException(method, "mis-matched dates between "
                                     "risk barrier ("+firstKO.toString()+
                                     ", " + lastKO.toString() + ") and "
                                     "economic barrier ("+firstEco.toString()+
                                     ", " + lastEco.toString() + ")");
            }
        }
        
        // first fee date must be after the initial accrue date
        if (feeDates[0] <= initialAccrueDate) {
            throw ModelException(method, "first fee date (" + feeDates[0].toString() + 
                                 ") must be after initial accrue date (" +
                                 initialAccrueDate.toString() + ")");
        }


        if (useAccrueDates) {
            // basic validation
            DateTime::ensureIncreasing(accrueDates, "accrue dates", true);

            // accrueDate[0] is the end of the first accrue period
            if (accrueDates[0] <= initialAccrueDate) {
                throw ModelException(method, "first accrue date (" + accrueDates[0].toString() + 
                                     ") must be after initial accrue date (" +
                                     initialAccrueDate.toString() + ")");
            }

            // one fee paiment date for each accrue period
            if (accrueDates.size()!=feeDates.size()) {
                throw ModelException(method, "there are  "+ Format::toString(accrueDates.size()) +
                                     "accrue dates and "+ Format::toString(feeDates.size()) +" fee dates");
            }

            // accrue period period must finish before the fee is paid
            int i;
            for (i=0; i<accrueDates.size(); i++) {
                if (accrueDates[i]>feeDates[i]) {
                    throw ModelException(method, "the accrue date ("
                                         + accrueDates[i].toString() + ") is "
                                         "after fee date (" + feeDates[i].toString() + ")");
                }
            }
        }
        
        // and finally
        daycount = DayCountConventionSP(DayCountConventionFactory::make(dcc));
    }

private:
    friend class EquityDistressSwapHelper;
    friend class EDSFDProd;  

    EquityDistressSwap():Generic1Factor(TYPE),
                         fee(0.0),
                         noAccruedOnBreach(false),
                         useAccrueDates(false),
                         accrueDates(0),
                         intraDayMonitor(false),
                         isUp(false), 
                         isBreached(false) {}; 

    EquityDistressSwap(const EquityDistressSwap& rhs);
    EquityDistressSwap& operator=(const EquityDistressSwap& rhs);
   
    // what is the i-th accrue fee
    DateTime getAccrueDate(int idx) const {
        static const string method("EquityDistressSwap::getAccruedDate");
        if (idx>=feeDates.size()) {
            throw ModelException(method, "there are  "+ Format::toString(feeDates.size()+1) +
                                 "accrue dates and "+ " the accrue fee date of index " + 
                                 Format::toString(idx+1) +" is required");
        }
        else {
            return useAccrueDates?accrueDates[idx]:feeDates[idx];
        }
    }
    
    // what's the accrued fee at a given date ?
    double accruedFee(const DateTime& hitDate) const {
        static const string method("EquityDistressSwap::accruedFee");
        try {
            if (noAccruedOnBreach) {
                return 0.0;
            }
            else {
                int      payidx;
                DateTime paydate;
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
                    // we have to make sure there is no negative accrual, hence the flooring with zero
                    double accrued = fee * Maths::max(daycount->years(lo, hitDate),0.0); 
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


    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        if (CTree1fLV::TYPE->isInstance(model)) {
            return true; // do pointwise instead
        }
        return false;
    }

    // make a LN vol request - not in real use as prefer LV tree
    CVolRequest* makeVolRequest() const {
        static const string method("EquityDistressSwap::makeVolRequest");
        try {
            double level = barrier->lastValue();
            if (!fwdStarting) {
                level *= initialSpot;
            }
            // do it at last barrier level for want of anywhere better
            return new LinearStrikeVolRequest(level,
                                              startDate,
                                              barrier->lastDate(),
                                              fwdStarting);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Returns all strikes the EquityDistressSwap is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method("EquityDistressSwap::getSensitiveStrikes");
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
        DateTime maturity = barrier->lastDate();
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;    
    }

    /** Satisfy LegalTerms::Shift interface */
    bool sensShift(LegalTerms* shift) {
        // Set the barriers for pricing equal to the economic barriers
        if (ecoBarrier.get()) {
            barrier = ecoBarrier;
        }
              
        return false;
    }

    bool sensShift(Theta* shift) {
        static const string method = "EquityDistressSwap::sensShift";
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

    bool priceDeadInstrument(CControl* control, CResults* results) const{
        static const string method = "EquityDistressSwap::priceDeadInstrument";
        try  {
            bool     deadInstrument = false;
            DateTime end = instSettle->settles(barrier->lastDate(),asset.get());

            if (isBreached) {
                double value = 0.0;
                // see if already paid
                DateTime pays = instSettle->settles(breachDate, asset.get());

                // look for the fees to pay in the future
                // in the case accrue date <= value date < fee date
                int i;
                double nextFee(0.0);
                for (i=0; i<feeDates.size(); i++){
                    if (getAccrueDate(i) <= valueDate && valueDate < feeDates[i]) {
                        nextFee += feeToPay(feeDates[i]) * discount->pv(feeDates[i]);
                    }
                }
                value -= nextFee;
                
                if (pays > valueDate) {
                    // get fee accrued - note if breached on fee day but after the fee was paid 
                    // this accrued will be zero which matches fact that only the rebate is due
                    double accrued(accruedFee(breachDate));
                    
                    // if physical settle, then the rebate is 0.0
                    // as a physical delivery event has been created
                    double reb;
                    if (instSettle->isPhysical()) {
                        reb =0.0;
                    }
                    else {
                        reb = rebate->interpolate(breachDate);
                    }

                    value += notional*(reb - accrued) * discount->pv(pays);
                }
                results->storePrice(value, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, value);
                }
                              
                deadInstrument = true;  
            }
            else if (valueDate > barrier->lastDate()) {
                // all over or no breach and last barrier passed
                // so last fees must be paid too
                // calculate the value of the fees that have still to be paid in the future
                // if the instrument has not breach
                double extraFees(0.0);
                int i;
                for (i=0; i<feeDates.size(); i++) {
                    if (feeDates[i] > valueDate) {
                        extraFees += feeToPay(feeDates[i])*discount->pv(feeDates[i]);
                    }
                }
                
                results->storePrice(-extraFees, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, -extraFees);
                }
                              
                deadInstrument = true;  
            }
            return deadInstrument;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }            
    }

    // gets known cashflows
    void getKnownCashflows(CashFlowArray& cfl) const {
        int i;
        // known fees up to breach
        for (i = 0; 
             (isBreached ? feeDates[i] <= breachDate : true) && 
                 i < feeDates.size(); i++) {
            CashFlow cf(feeDates[i], -feeToPay(feeDates[i])*notional);
            cfl.push_back(cf);
        }
        
        // at breach + settlement it's rebate-accrued
        if (isBreached) {
            DateTime pays = instSettle->settles(breachDate, asset.get());
            double accrued = accruedFee(breachDate);
            double reb     = rebate->interpolate(breachDate);
            // on breach
            
            double value;
            if (instSettle->isPhysical()) {
                // if physically settle, the rebate is not cash settle
                value = 0.0 - accrued;
            }
            else {
                value   = (reb - accrued);
            }
            
            CashFlow cf(pays, value*notional);
            cfl.push_back(cf);
            
            // known fees after the breach
            for (i = 0; i < feeDates.size(); i++) {
                if (getAccrueDate(i) <= breachDate&& feeDates[i] > breachDate) {
                    CashFlow cf(feeDates[i], -feeToPay(feeDates[i])*notional);
                    cfl.push_back(cf);
                }
            }
        }
    }

    /** extra output requests */
    void recordOutputRequests(Control* control, 
                              Results* results, 
                              double   fairValue) const {
        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       barrier->lastDate(),
                                       valueDate,
                                       asset.get());

        OutputRequest* request = NULL;
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        int i;

        if (request) {
            DateTimeArray paydates;
            // all fee dates up to any breach date
            // then the breach date (if it exists)
            for (i = 0; i < feeDates.size(); i++) {
                if (!isBreached || getAccrueDate(i) <= breachDate) {
                    // no settlement on fees
                    paydates.push_back(feeDates[i]);
                }
            }
                       
            if (isBreached) {
                paydates.push_back(instSettle->settles(breachDate,asset.get()));
            } 

            OutputRequestUtil::recordPaymentDates(control,results,&paydates); 
        }
        
        request = control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY);
        if (request) {
            if (instSettle->isPhysical()&&isBreached) {
                DateTime pays = instSettle->settles(breachDate, asset.get());
                double reb     = rebate->interpolate(breachDate);
                                
                // use the economic barrier (if it exists)
                double bar = (ecoBarrier.get()? 
                              ecoBarrier->interpolate(breachDate):
                              barrier->interpolate(breachDate));
                
                // the number of units to deliver is 
                // notional*rebate/(barrier level * spot at start)
                double nbToDeliver = notional*reb/(initialSpot*bar);
                
                PhysicalDelivery::recordPhysicalDelivery(
                    nbToDeliver, 
                    0.0, 
                    pays, 
                    asset.get(),
                    control,
                    results);
            }
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArray cfl;

            getKnownCashflows(cfl);
                                
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    &cfl);   
        }   
    
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request && !fwdStarting && !isBreached) {
            // report barrier levels over a date range
            DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
            // use economic barrier (if it exists)
            Schedule* s = ecoBarrier.get() ?ecoBarrier.get():barrier.get();

            CashFlowArraySP subset(s->subset(valueDate, upperDate));
            if (!subset->empty()) {
                BarrierLevelArraySP levels(new BarrierLevelArray(0));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(isUp,(*subset)[i].date,
                                    (*subset)[i].amount*initialSpot, intraDayMonitor);
                    levels->push_back(bl);
                }

                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       levels.get());
            }
        }
    }  
    
    // a private helper that makes a barrier event once one is detected
    void makeBarrierEvent(DateTime breachDate, double level, double barrLvl, 
                          EventResults* events) const{
        StringArraySP assetNames(new StringArray(1, asset->getTrueName()));
        DoubleArraySP assetLevels(new DoubleArray(1, level));
        DoubleArraySP barrLevels(new DoubleArray(1, barrLvl));

        events->addEvent(new BarrierBreach(breachDate, "EDS Barrier",
            intraDayMonitor ? BarrierBreach::CONTINUOUS : BarrierBreach::DAILY,
            BarrierBreach::KNOCK_OUT, isUp, 0,
            assetNames, assetLevels, barrLevels));
    }

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;
 
public:

    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach,
                   IModel* model, 
                   const DateTime& eventDate,
                   EventResults* events) const {
        static const string method = "EquityDistressSwap::getEvents";

        // no barrier events unless the barrier is active
        if (barrier->firstDate() <= eventDate &&
                eventDate <= barrier->lastDate()) {
            // if breached flag is set report event as best we can
            if (isBreached && breachDate < eventDate) {
                makeBarrierEvent(breachDate, 0.0, 
                            barrier->interpolate(eventDate) * initialSpot, 
                            events);
            }
            else {
                // check for breach right now
                // if not doing cts monitoring we can only have event
                // if time of day is right
                if (intraDayMonitor ||
                        eventDate.getTime() == barrier->firstDate().getTime()) {
                    double level = asset->fwdValue(eventDate);
                    double barrLvl = barrier->interpolate(eventDate) * initialSpot;
                    if (isUp ? (level > barrLvl) : (level < barrLvl)) {
                        makeBarrierEvent(eventDate, level, barrLvl, events);
                    }
                }
            }
        }
    }

    // KnownCashflows::IEventHandler interface
    void getEvents(const KnownCashflows* flows,
                   IModel* model, 
                   const DateTime& eventDate,
                   EventResults* events) const {
        static const string method = "EquityDistressSwap::getEvents";
        CashFlowArraySP cfl(new CashFlowArray(0));
        getKnownCashflows(*cfl);
        if (!cfl->empty()) {
            events->addEvent(new KnownCashflows(eventDate, cfl, 
                                                discount->getCcy()));
        }
    }

    // for EDS/CDS spread curves we want a tree without too many steps or we just
    // die in performance terms - centralise it here
    static CTree1f* quickModel(const CTree1f* model) {
        static int QUICK_TREE_STEPS = 100;
        int steps = model->GetStepsPerYear();
        int stepsToUse = QUICK_TREE_STEPS < steps ? QUICK_TREE_STEPS : steps;
        CTree1fSP tree(copy(model));
        tree->SetStepsPerYear(stepsToUse);
        return (tree.release());
    }
        

    // create a simple, started EDS from parameters
    static EquityDistressSwap* makeEDS(
        const DateTime&             valueDate,
        const DateTime&             maturity,
        double                      fee,
        string                      feeInterval,
        string                      dcc,
        bool                        isUp,
        double                      barrier,
        double                      rebate,
        double                      initialSpot,
        string                      ccyTreatment,
        const InstrumentSettlement* instSettle,
        const CAsset*               asset,
        const YieldCurve*           discount) {
        static const string method("EquityDistressSwap::makeEDS");
        try {
            EquityDistressSwapSP eds(new EquityDistressSwap);

            // plug in the easy bits
            eds->valueDate    = valueDate;
            eds->fwdStarting  = false;
            eds->oneContract  = false;
            eds->notional     = 1.0;
            eds->fee          = fee;
            eds->dcc          = dcc;
            eds->initialAccrueDate = valueDate;
            eds->isUp         = isUp;
            eds->initialSpot  = initialSpot;
            eds->ccyTreatment = ccyTreatment;
            eds->instSettle   = InstrumentSettlementSP(copy(instSettle));
            eds->asset        = CAssetWrapper(copy(asset));
            eds->discount     = YieldCurveWrapper(copy(discount));

            // build up schedules
            DateTimeArray dates;
            DoubleArray   ko;
            DoubleArray   reb;

            dates.push_back(valueDate);
            dates.push_back(maturity);
            ko.push_back(barrier);
            ko.push_back(barrier);
            reb.push_back(rebate);
            reb.push_back(rebate);

            eds->barrier = ScheduleSP(new Schedule(dates, 
                                                   ko, 
                                                   Schedule::INTERP_LINEAR));
        
            eds->rebate = ScheduleSP(new Schedule(dates, 
                                                  reb, 
                                                  Schedule::INTERP_LINEAR));
            
            // now the fee
            DateTimeArraySP feeDates(SwapTool::paymentDates(valueDate,
                                                            maturity,
                                                            1,
                                                            feeInterval,
                                                            false));
            eds->feeDates = (*feeDates);
            
            // check it's OK
            eds->Validate();

            return (eds.release());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // create a simple, started EDS from parameters, and price it
    static double quickPricer(
        const DateTime&             valueDate,
        const DateTime&             maturity,
        double                      fee,
        string                      feeInterval,
        string                      dcc,
        bool                        isUp,
        double                      barrier,
        double                      rebate,
        double                      initialSpot,
        string                      ccyTreatment,
        const InstrumentSettlement* instSettle,
        const CAsset*               asset,
        const YieldCurve*           discount,
        CTree1f*                    model) {
        static const string method("EquityDistressSwap::quickPricer");
        try {
            double price;
            // first, build an EDS
            EquityDistressSwapSP eds(makeEDS(valueDate,
                                             maturity,
                                             fee,
                                             feeInterval,
                                             dcc,
                                             isUp,
                                             barrier,
                                             rebate,
                                             initialSpot,
                                             ccyTreatment,
                                             instSettle,
                                             asset,
                                             discount));

            // price it
            IModelSP tree(copy(model));
            CControlSP ctrl(Control::makeFromFlags("", 0.0));

            price = tree->calcPrice(eds.get(), ctrl.get());
            
            return price;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // price and greeks for simple, started EDS
    static Results* quickPricer(
        const DateTime&             valueDate,
        const DateTime&             maturity,
        double                      fee,
        string                      feeInterval,
        string                      dcc,
        bool                        isUp,
        double                      barrier,
        double                      rebate,
        double                      initialSpot,
        string                      ccyTreatment,
        const InstrumentSettlement* instSettle,
        const CAsset*               asset,
        const YieldCurve*           discount,
        Control*                    ctrl,
        CTree1f*                    model) {
        static const string method("EquityDistressSwap::quickPricer");
        try {
            // first, build an EDS
            EquityDistressSwapSP eds(makeEDS(valueDate,
                                             maturity,
                                             fee,
                                             feeInterval,
                                             dcc,
                                             isUp,
                                             barrier,
                                             rebate,
                                             initialSpot,
                                             ccyTreatment,
                                             instSettle,
                                             asset,
                                             discount));

            // price it
            IModelSP tree(copy(model));
            ResultsSP results(new Results);

            ctrl->calculate(tree.get(), eds.get(), results.get());

            return (results.release());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    void edsSpreads(Control* ctrl, CTree1f* model, Results* results) const {
        static const string method("EquityDistressSwap::edsSpreads");
        try {
            static const string tenors[] = {"1Y", "2Y", "3Y", "4Y", "5Y", "7Y", "10Y"};
            int numTenors = 7;
            OutputRequest* ufRequest = ctrl->requestsOutput(OutputRequest::EDS_UPFRONT);
            OutputRequest* rfRequest = ctrl->requestsOutput(OutputRequest::EDS_RUNNING_FEE);
            OutputRequest* rdRequest = ctrl->requestsOutput(OutputRequest::EDS_RISKY_DURATION);
            OutputRequest* dcRequest = ctrl->requestsOutput(OutputRequest::EDS_DELTA_CURVE);
            OutputRequest* vcRequest = ctrl->requestsOutput(OutputRequest::EDS_VEGA_CURVE);

            bool doFV1    = rfRequest || rdRequest;
            bool doGreeks = dcRequest || vcRequest;

            try {
                // only makes sense if we're already started
                if (!fwdStarting && (ufRequest || rfRequest || rdRequest || doGreeks)) {   
                    CControlSP ctrl(Control::makeFromFlags("", 0.0));
                    CTree1fSP tree(quickModel(model));

                    ExpiryResultArraySP upFrontCurve(new ExpiryResultArray(0));
                    ExpiryResultArraySP runningFeeCurve(new ExpiryResultArray(0));
                    ExpiryResultArraySP durationCurve(new ExpiryResultArray(0));
                    ExpiryResultArraySP deltaCurve(new ExpiryResultArray(0));
                    ExpiryResultArraySP vegaCurve(new ExpiryResultArray(0));

                    SensControlPerNameSP delta(new Delta(Delta::DEFAULT_SHIFT));
                    SensitivitySP vega(
                        new VegaParallel(VegaParallel::DEFAULT_SHIFT));
                    
                    if (dcRequest) {
                        ctrl->addSensitivity(delta);
                    }
                    if (vcRequest) {
                        ctrl->addSensitivity(vega);
                    }
                   
                    OutputNameSP deltaName;
                    OutputNameSP vegaName;

                    if (doGreeks) {
                        // get names used to store delta & vega
                        // can't use asset "true" name as EAS don't use same 
                        // string for asset and vol
                        // fail if there are multiple names
                        deltaName = OutputNameSP(new OutputName(asset->getTrueName()));
                    
                        // round-about route for vols as the name is some random client choice
                        ATMVolRequest   vr;
                        IVolProcessedSP vp(asset->getProcessedVol(&vr));
                        vegaName = OutputNameSP(new OutputName(vp->getName()));
                    }

                    for (int i = 0; i < numTenors; i++) {
                        ExpiryConstSP expiry(new MaturityTimePeriod(tenors[i], DateTime::END_OF_DAY_TIME));
                        DateTime maturity = expiry->toDate(valueDate);
                        double   flatBarrier = barrier->firstValue();
                        double   flatRebate  = rebate->firstValue();

                        // want to price with a zero fee, than a 1% fee
                        ResultsSP results(EquityDistressSwap::quickPricer(valueDate,
                                                                          maturity,
                                                                          0.0,
                                                                          "Q",
                                                                          dcc,
                                                                          isUp,
                                                                          flatBarrier,
                                                                          flatRebate,
                                                                          initialSpot,
                                                                          ccyTreatment,
                                                                          instSettle.get(),
                                                                          asset.get(),
                                                                          discount.get(),
                                                                          ctrl.get(),
                                                                          tree.get()));
                        double upFront = results->retrievePrice();

                        upFrontCurve->push_back(ExpiryResult(expiry, upFront));

                        if (doFV1) {
                            double fv1 = EquityDistressSwap::quickPricer(valueDate,
                                                                         maturity,
                                                                         0.01,
                                                                         "Q",
                                                                         dcc,
                                                                         isUp,
                                                                         flatBarrier,
                                                                         flatRebate,
                                                                         initialSpot,
                                                                         ccyTreatment,
                                                                         instSettle.get(),
                                                                         asset.get(),
                                                                         discount.get(),
                                                                         tree.get());

                            double riskyDuration = (upFront - fv1)/0.01;
                            double runningFee    = upFront/riskyDuration;

                            runningFeeCurve->push_back(ExpiryResult(expiry, runningFee));
                            durationCurve->push_back(ExpiryResult(expiry, riskyDuration));
                        }

                        if (doGreeks) {
                            if (dcRequest) {
                                double edsDelta = results->retrieveScalarGreek(delta->getPacketName(),
                                                                               deltaName);
                                deltaCurve->push_back(ExpiryResult(expiry, edsDelta));
                            }
                            if (vcRequest) {
                                double edsVega = results->retrieveScalarGreek(vega->getPacketName(),
                                                                              vegaName);
                                vegaCurve->push_back(ExpiryResult(expiry, edsVega));
                            }
                        }
                    }
                    if (ufRequest) {
                        results->storeRequestResult(ufRequest, upFrontCurve);
                    }
                    if (rfRequest) {
                        results->storeRequestResult(rfRequest, runningFeeCurve);
                    }
                    if (rdRequest) {
                        results->storeRequestResult(rdRequest, durationCurve);
                    }
                    if (doGreeks) {
                        // already checked if there are multiple names?
                        // report both against asset name
                        if (dcRequest) {
                            results->storeRequestResult(dcRequest, deltaCurve, deltaName);
                        }
                        if (vcRequest) {
                            results->storeRequestResult(vcRequest, vegaCurve, deltaName);
                        }
                    }
                }
            }
            catch (exception& e) {
                IObjectSP dead(new Untweakable(e));
                if (ufRequest) {
                    results->storeRequestResult(ufRequest, dead);
                }
                if (rfRequest) {
                    results->storeRequestResult(rfRequest, dead);
                }
                if (rdRequest) {
                    results->storeRequestResult(rdRequest, dead);
                }                    
                if (dcRequest) {
                    results->storeRequestResult(dcRequest, dead);
                }
                if (vcRequest) {
                    results->storeRequestResult(vcRequest, dead);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    // define fee payments
    double        fee;
    string        dcc;
    DateTime      initialAccrueDate;
    DateTimeArray feeDates;
    bool          noAccruedOnBreach;
    bool          useAccrueDates;
    DateTimeArray accrueDates;

    // define barrier
    bool          intraDayMonitor;
    bool          isUp;
    bool          isBreached;
    DateTime      breachDate;
    ScheduleSP    barrier;
    ScheduleSP    rebate;
    ScheduleSP    ecoBarrier;   // economic barrier - not for pricing

    // internal
    DayCountConventionSP daycount;
};

/**************************************************************************************************************/
class EDSFDProd: public LatticeProdEDRIns 
{
public:
    EDSFDProd(const EquityDistressSwap* eds, FDModel* model) : 
        LatticeProdEDRIns(model, 1, 1), inst(eds)
    {
        if( ! tree1f )
        {
            throw ModelException( "EDSFDProd::EDSFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );
    }

    virtual void update(int& step, FDProduct::UpdateType type);

    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(CControl*  control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing 
    virtual void initProd();    

    void compileCriticalDates(DateTimeArray & criticalDates) const;

    /** returns a vol request for log-normal vol */
    virtual CVolRequestConstSP GetLNRequest() const
    {
        CVolRequestConstSP volRequest(inst->makeVolRequest());
        return volRequest;
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
    virtual double scalePremium(const double & fairValue, 
                                YieldCurveConstSP disc);

    /** extra output requests */    
    void recordOutput(Control* control, 
                      YieldCurveConstSP disc, 
                      Results* results); 

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    virtual string getCcyTreatment() const
    {
        return inst->ccyTreatment;
    }

private:
    const EquityDistressSwap* inst;
    vector<bool>              stepIsKO;
    vector<bool>              stepIsFee;
    vector<double>            stepBarrier;
    vector<double>            stepRebate;
    double                    refLevel;
    // barrier at a tree slice after discrete monitor adjustment
    double                    adjBarrier[2]; 
};

/** this sets up the timeline */
void EDSFDProd::init(CControl* control) const 
{    
    static const string method = "EDSFDProd::Init";
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
        segDates[1] = inst->barrier->lastDate();


        // timeline configuration
        // 'density factor' for timeline
        IntArray density( 1, 1 );

        // compile list of critical dates
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
void EDSFDProd::initProd() 
{
    static const string method = "EDSFDProd::InitProd";
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
                                 inst->barrier, 
                                 true,   // 'american'
                                 inst->asset.getSP());

        // set the reference level for the barrier
        refLevel = inst->fwdStarting ? 
            inst->asset->fwdValue(inst->startDate) : inst->initialSpot;

        // set up barrier & rebate for each step
        stepBarrier.resize(lastStep+1);
        stepIsFee.resize(lastStep+1);
        stepRebate.resize(lastStep);

        // handle maturity point
        stepBarrier[lastStep] = inst->barrier->lastValue()*refLevel;
    
        // set up barrier & rebate for each step
        stepBarrier.resize(lastStep+1);
        stepIsFee.resize(lastStep+1);
        stepRebate.resize(lastStep+1);

        // handle maturity point
        stepBarrier[lastStep] = inst->barrier->lastValue()*refLevel;
        int numFees = inst->feeDates.size();
        stepIsFee[lastStep] = false; // for starters
        if (inst->feeDates[numFees-1] == inst->barrier->lastDate()) 
        {
            stepIsFee[lastStep] = true;
        }

        for (i = 0; i <= lastStep; i++) 
        {
            if (stepIsKO[i]) 
            {
                DateTime stepDate = model->getDate(i);

                stepRebate[i] =inst->rebate->interpolate(stepDate);
                stepBarrier[i]=inst->barrier->interpolate(stepDate)*refLevel;
            }
            stepIsFee[i] = false; // for starters
        }

        // and now identify fee dates - note last fee needn't be last step
        int nextFee = 0;
        for (int j = 0; j < numFees; j++) 
        {
            for (i = nextFee; i <= lastStep; i++) 
            {
                if (inst->feeDates[j] == model->getDate(i)) 
                {
                    stepIsFee[i] = true;
                    nextFee = i+1;
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
double EDSFDProd::scalePremium(const double& fairValue,
                               YieldCurveConstSP disc)
{
    double fwdStartDF = 1.0;
    if (inst->fwdStarting)
        fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));

    return inst->notional*fwdStartDF*fairValue;
}   

void EDSFDProd::recordOutput(Control* control, 
                             YieldCurveConstSP disc, 
                             Results* results)
{
    if (control && control->isPricing()) 
    {
    inst->edsSpreads(control, tree1f, results);
    }

    // get prices at t=0
    // save price
    double price = scalePremium(model->getPrice0( *slices[0] ), disc);
    results->storePrice(price, disc->getCcy());

    // throw this back to the instrument itself
    inst->recordOutputRequests(control, results, price);
}

/** called before PayoffAtMat, PayoffBeforeMat and tree roll() */
// calculate barrier and place barrier at inserted node if needed

void EDSFDProd::preCalc(int step) 
{
    static const string method("EDSFDProd::preCalc");

    try 
    {
        int idx = tree1f->getSliceIndex(step);

        adjBarrier[idx] = stepBarrier[step];
        // taking care of discrete monitoring, treat as daily monitor

        if (!inst->intraDayMonitor) 
        {
            vector<double> vol;
            // get vol at barrier
            tree1f->GetStepVol(step, vol, &adjBarrier[idx], 0, 0);
            Barrier::BarrierAdjustment(vol[0], inst->isUp, adjBarrier[idx]);
        }

        if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
        { 
            // apparently last param = 1 means KO
            tree1f->SetInsertNode(idx, 0, adjBarrier[idx], 1); 
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** product payoff method at maturity */
void EDSFDProd::prod_BWD_T(const TreeSlice & spot,
                                 int step, 
                                 int bot, 
                                 int top, 
                                 int pStart, 
                                 int pEnd,
                                 const vector< TreeSliceSP > & price)  
{
    static const string method("EDSFDProd::prod_BWD_T");
    try 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        DateTime stepDate = model->getDate(step);

        double settlePV = inst->instSettle->pvAdjust(inst->barrier->lastDate(),
                                                     inst->discount.get(), 
                                                     inst->asset.get());

       // what's the accrued fee so far ?
        double accrued = inst->accruedFee(inst->barrier->lastDate());
        double ko      = adjBarrier[tree1f->GetSliceIdx()];
    
        int i;            
    
        // calculate the fees that have to be paid because the accrued period is finished
        // i.e. if (accrued date[i] <= step date < fee date[i])
        double dueFees(0.0);
        for (i=0; i<inst->feeDates.size(); i++) 
        {
            if (inst->getAccrueDate(i) <= stepDate && stepDate < inst->feeDates[i] ) 
            {
                dueFees += inst->feeToPay(inst->feeDates[i]) * 
                    inst->discount->pv(stepDate,inst->feeDates[i]);
            }
        }
    
        // calculate the value of the fees that have still to be paid in the future
        // if the instrument has not breach (does not exclude the due fees)
        // i.e. if (fee date[i] >= accrued date[i] > step date)
        double extraFees(0.0);
        for (i=0; i<inst->feeDates.size(); i++) 
        {
            if (inst->getAccrueDate(i) > stepDate && inst->feeDates[i] > stepDate) 
            {
                extraFees += inst->feeToPay(inst->feeDates[i]) * 
                    inst->discount->pv(stepDate,inst->feeDates[i]);
            }
        }
    
        for (i=pStart; i<=pEnd; i++) 
        {
            for (int j=bot; j<=top; j++) 
            {
                // pay the due fees in all the case
                (p[i])[j] = -dueFees;

                // check if we've breached
                if ((inst->isUp && s[j]>ko-FP_MIN)||
                    (!inst->isUp && s[j] < ko+FP_MIN)) 
                {              
                    // pay the rebate and the accrued fee
                    (p[i])[j] += (stepRebate[step]-accrued)*settlePV; 
                }
                else 
                {
                    // if no breach, extra fees have to be oaid 
                    (p[i])[j] -= extraFees;
                }
            
                // pay the fee if we are on a fee date
                if (stepIsFee[step]) 
                {
                    (p[i])[j] -= inst->feeToPay(stepDate);
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
void EDSFDProd::prod_BWD(const TreeSlice & spot,
                               int step, 
                               int bot, 
                               int top, 
                               int pStart, 
                               int pEnd,
                               const vector< TreeSliceSP > & price) 
{
    static const string method("EDSFDProd::prod_BWD");
    try 
    {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        if (!stepIsKO[step] && !stepIsFee[step]) 
        {
            return;
        }

        DateTime stepDate = model->getDate(step);

        double settlePV = inst->instSettle->pvAdjust(stepDate,
                                                     inst->discount.get(), 
                                                     inst->asset.get());

        // what's the accrued fee so far ?
        double accrued = inst->accruedFee(stepDate);
        double ko      = adjBarrier[tree1f->GetSliceIdx()];

        for (int i=pStart; i<=pEnd; i++) 
        {
            for (int j=bot; j<=top; j++) 
            {
       
                // if it's a KO date, see if client pays accrued and
                // receives rebate instead
                if (stepIsKO[step]) 
                {
                    // if breach pay the rebate and accrued fee
                    if ((inst->isUp && s[j] > ko-FP_MIN) || (!inst->isUp && s[j] < ko+FP_MIN)) 
                    {
                        (p[i])[j] = (stepRebate[step]-accrued)*settlePV;
                    } 
                }

                // if it's a fee date, pay the fee regardless of KO
                if (stepIsFee[step]) 
                {
                    (p[i])[j] -= inst->feeToPay(stepDate);
                }
            }
        }                
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

void EDSFDProd::compileCriticalDates(DateTimeArray & criticalDates) const
{
    static const string method("EDSFDProd::compileCriticalDates");
    try 
    {   
        // compile list of critical dates:
        // fee pay dates
        // barrier dates
    
        DateTime lastKO = inst->barrier->lastDate();

        criticalDates.resize(inst->feeDates.size());

        int i;
        for (i = 0; i < criticalDates.size() && inst->feeDates[i]<=lastKO; i++) 
        {
            criticalDates[i] = inst->feeDates[i];
        }

        criticalDates.resize(i);

        DateTime        maturity = inst->barrier->lastDate();
        const Schedule* barrier  = inst->barrier.get();

        if (barrier->getInterp() != Schedule::INTERP_NONE) 
        {
            const DateTimeArray& barDates = barrier->getDates();
            for (i = 0; i < barDates.size(); i++) 
            {
                criticalDates.push_back(barDates[i]);
            }
        }
        else 
        {
            const DateTimeArray& barDates = barrier->getDates();
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

        // sort critical dates & make unique
        sort(criticalDates.begin(), criticalDates.end());
        DateTime::removeDuplicates(criticalDates, false);        
    }
    catch(exception& e) 
    {
        throw ModelException(e, method);
    }
}

void EDSFDProd::update(int & step, 
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

/** create a fd payoff product */
FDProductSP EquityDistressSwap::createProduct(FDModel* model) const
{
    return FDProductSP( new EDSFDProd(this, model) );
}

class EquityDistressSwapHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Equity Default Swap (EDS)");
        REGISTER(EquityDistressSwap, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        IMPLEMENTS(KnownCashflows::IEventHandler);
        EMPTY_SHELL_METHOD(defaultEquityDistressSwap);
        FIELD(fee, "fee");
        FIELD(dcc, "fee day count");
        FIELD(initialAccrueDate, "initial accrue date");
        FIELD(feeDates, "fee payment dates");
        FIELD(intraDayMonitor, "intra-day barrier");
        FIELD(isUp, "is up barrier?");
        FIELD(isBreached, "is it breached?");
        FIELD(breachDate, "breach date");
        FIELD(barrier, "barrier");
        FIELD(rebate, "rebate");        
        FIELD(ecoBarrier, "economic barrier");
        FIELD_MAKE_OPTIONAL(ecoBarrier);
        FIELD_NO_DESC(daycount);
        FIELD_MAKE_TRANSIENT(daycount);
        FIELD(noAccruedOnBreach, "if true, no accrued fee is paid on breach");
        FIELD_MAKE_OPTIONAL(noAccruedOnBreach);
        FIELD(useAccrueDates, "if true, use the accrue"
                     " dates to compute the fees, otherwise use the fee dates");
        FIELD_MAKE_OPTIONAL(useAccrueDates);       
        FIELD(accrueDates, "these dates are used to compute"
                     " the fees id useAccrueDates is true");
        FIELD_MAKE_OPTIONAL(accrueDates);
    }

    static IObject* defaultEquityDistressSwap(){
        return new EquityDistressSwap();
    }
};

CClassConstSP const EquityDistressSwap::TYPE = CClass::registerClassLoadMethod(
    "EquityDistressSwap", typeid(EquityDistressSwap), 
    EquityDistressSwapHelper::load);


/* for class loading */
bool EquityDistressSwapLoad() {
    return (EquityDistressSwap::TYPE != 0);
}

DRLIB_END_NAMESPACE

