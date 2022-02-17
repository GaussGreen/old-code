//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CallableDeposit.cpp
//
//   Description : Callable Deposit Instrument (generic 1 factor)
// 
// Client 
//   o RECEIVES any fixed coupons.
//   o PAYS floating rate plus a spread. 
//   o OWNS a portfolio of European calls/puts of same maturity.
// 
// Issuer
//   o Can CALL (cancel) the structure for preset amounts at preset dates.
//
//   Author      : Andrew J Swain
//
//   Date        : 18 March 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/Average.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class CallableDeposit: public Generic1Factor,
                       virtual public LegalTerms::Shift,
                       virtual public FDModel::IIntoProduct, 
                       virtual public LastSensDate, 
                       virtual public Callability::IEventHandler {
public:
    static CClassConstSP const TYPE;  

    // represents the portolfio of (avg spot) options the client holds
    // strikes are percentages 
    class Portfolio: public CObject {
    public:
        static CClassConstSP const TYPE;  

        // methods
        double intrinsic(double spot, double refLevel) const {
            double iv = 0.0;
            double performance = spot/refLevel;
            for (int i = 0; i < strike.size(); i++) {
                iv += longShort[i]*participation[i]*
                    Maths::max(0.0, callPut[i]*(performance - strike[i]));
            }

            return iv;
        }

        bool isEmpty() const {
            // move this out of validatePop2Obj 
            return longShort.empty();
        }  

        void validatePop2Object() {
            static const string method("CallableDeposit::Portfolio::validatePop2Object");
            try {
                if (longShort.size() != callPut.size()) {
                    throw ModelException(method, "long/short and call/put "
                                         "flags are different lengths");
                }     
                if (longShort.size() != participation.size()) {
                    throw ModelException(method, "long/short and participation "
                                         "are different lengths");
                }     
                if (longShort.size() != strike.size()) {
                    throw ModelException(method, "long/short and strike "
                                         "are different lengths");
                }    

                int i;
                for (i = 0; i < longShort.size(); i++) {
                    if (longShort[i] != 1 && longShort[i] != -1) {
                        throw ModelException(method, "long/short " +
                                             Format::toString(i+1) + " (" + 
                                             Format::toString(longShort[i]) + 
                                             ") must be +1 or -1");
                    }
                       
                    if (callPut[i] != 1 && callPut[i] != -1) {
                        throw ModelException(method, "call/put " +
                                             Format::toString(i+1) + " (" + 
                                             Format::toString(callPut[i]) + 
                                             ") must be +1 (call) or -1 (put)");
                    }
                }

                if ((isPerCallDate && !callDateId.get()) ||
                    (isPerCallDate && callDateId->empty())) {
                        throw ModelException(method, "portfolios per call date "
                                             "selected, but none supplied");
                }
                                             
                if (isPerCallDate) {
                    if (longShort.size() != callDateId->size()) {
                        throw ModelException(method, "long/short and call "
                                             "date id are different lengths");
                    }    
                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

        // make a portfolio for a given call date id
        Portfolio* makePortfolio(int callId, int numCallDates) const {
            static const string method("CallableDeposit::Portfolio::makePortfolio");
            try {
                Portfolio* portfolio = new Portfolio();
                // callDateId uses user indexing 1->N
                // might be missing in which case only maturity
                // portfolio exists
                bool justMaturity = !isPerCallDate;

                for (int i = 0; i < longShort.size(); i++) {
                    int iCall = justMaturity ? numCallDates : (*callDateId)[i]-1;
                    if (callId == iCall) {
                        portfolio->longShort.push_back(longShort[i]);
                        portfolio->callPut.push_back(callPut[i]);
                        portfolio->participation.push_back(participation[i]);
                        portfolio->strike.push_back(strike[i]);
                    }
                }
                return portfolio;
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
         }

         // data
        IntArray    longShort;
        IntArray    callPut;
        DoubleArray participation;
        DoubleArray strike;
        bool        isPerCallDate;  // determines whether or not callDateId is used
        // optional - allows multiple portfolios by specify which call date
        // each of the above applies to. If missing, then we just have a 
        // maturity portfolio, otherwise it's N+1 items, 1 per call date
        // plus maturity
        IntArraySP  callDateId;
        // optional - flag indicating that the portfolio intrinsic value is paid on the
        // call date, not the notification date. The value of the portfolio is still calculated
        // on the notification date.
        bool payPortfolioOnCallDate;

        Portfolio():CObject(TYPE), isPerCallDate(false), payPortfolioOnCallDate(false) {}; 
   
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("Callable Deposit Portfolio");
            REGISTER(CallableDeposit::Portfolio, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultPortfolio);
            FIELD(longShort, "long or short");
            FIELD(callPut, "call or put");
            FIELD(participation, "participation");
            FIELD(strike, "strike");
            FIELD(isPerCallDate, "do we have multiple portfolios");
            FIELD(callDateId, "call date ID for multiple portfolios");
            FIELD_MAKE_OPTIONAL(callDateId);
            FIELD(payPortfolioOnCallDate, "pay intrinsic portfolio value on call date, not notification date");
            FIELD_MAKE_OPTIONAL(payPortfolioOnCallDate);
        }
        
        static IObject* defaultPortfolio(){
            return new CallableDeposit::Portfolio();
        }
    };

    typedef smartPtr<Portfolio> PortfolioSP;

    // useful data holder
    class CallSchedule: public CObject {
    public:
        static CClassConstSP const TYPE;  

        bool          isPuttable;
        bool          includeAccrued;
        DateTimeArray notification;
        CashFlowArray schedule;

        CallSchedule():CObject(TYPE),isPuttable(false),includeAccrued(false) {}; 

        void validatePop2Object() {
            static const string method("CallableDeposit::CallSchedule::validatePop2Object");
            try {
                if (notification.size() != schedule.size()) {
                    throw ModelException(method, "notification dates and call "
                                         "schedule are different lengths");
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
            clazz->setDescription("Callable Deposit Call Schedule");
            REGISTER(CallableDeposit::CallSchedule, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultCallSchedule);
            FIELD(isPuttable, "is puttable");
            FIELD(includeAccrued, "include accrued");
            FIELD(notification, "notification dates");
            FIELD(schedule, "schedule");
        }
        
        static IObject* defaultCallSchedule(){
            return new CallableDeposit::CallSchedule();
        }
    };

    bool hasMonitor() const {
        return !!monitorDates && !monitorDates->empty();
    }

    virtual void Validate(){
        static const string method = "CallableDeposit::Validate";
        int i;
        // Call Generic1Factor validate()
        validate();

        if (oneContract) {
            throw ModelException(method, 
                                 "Callable Deposit must be notionally booked");
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

        if (portfolio.isEmpty() ) {
            throw ModelException(method, "portfolio is empty");
        }
            
        if (hasFloater && !floater.get()) {
            throw ModelException(method, "floating leg is missing");
        }

        if (callSchedule.includeAccrued && !hasFloater) {
            throw ModelException(method, 
                                 "Cannot include accrued payment in call "
                                 "amount if there is no floating leg");
        }

        if (hasFloater && fwdStarting) {
            if (startDate > floater->RefixDates[0]) {
                throw ModelException(method, 
                                     "start date ("+ startDate.toString() + 
                                     ") is after first floating fixing (" + 
                                     floater->RefixDates[0].toString() + ")");
            }
        }

        CashFlow::ensureDatesIncreasing(coupons, "coupons", false);
        CashFlow::ensureDatesIncreasing(callSchedule.schedule,
                                        "call schedule", false);
        DateTime::ensureIncreasing(callSchedule.notification, 
                                   "call notification", false);

        if (!coupons.empty()) {
            if (coupons[coupons.size()-1].date > maturity) {
                throw ModelException(method, 
                                     "last coupon date (" + 
                                     coupons[coupons.size()-1].date.toString()+
                                     ") is after maturity (" + 
                                     maturity.toString() + ")");
            }                              

            if (fwdStarting && coupons[0].date < startDate) {
                throw ModelException(method, 
                                     "first coupon date (" + 
                                     coupons[0].date.toString()+
                                     ") is before start date (" + 
                                     startDate.toString() + ")");
            }
        }
                
        // can't call before notification
        for (i = 0; i < callSchedule.notification.size(); i++) {
            if (callSchedule.notification[i] > callSchedule.schedule[i].date) {
                throw ModelException(method, 
                                     "call notification date (" + 
                                     callSchedule.notification[i].toString()+
                                     ") must come before call date (" + 
                                     callSchedule.schedule[i].date.toString() + 
                                     ")");
            }              
        }

        if (!callSchedule.schedule.empty()) {
            int numCalls = callSchedule.schedule.size();
            DateTime& lastCall = callSchedule.schedule[numCalls-1].date;
            if (lastCall > maturity) {
                throw ModelException(method, 
                                     "last call date (" + lastCall.toString() + 
                                     ") is after maturity ("+maturity.toString()+")");
            }

            DateTime& lastNotify = callSchedule.notification[numCalls-1];
            if (avgOut->getFirstDate() < lastNotify) {
                throw ModelException(method, 
                                     "first avg out date (" + 
                                     avgOut->getFirstDate().toString() + 
                                     ") is before last call notice date ("+
                                     lastNotify.toString()+")");
            }

            // currently get tied in knots if really do have averaging
            // and the first avg out co-incides with last notification
            // just validate against it for now
            if (avgOut->numDates(0) > 1) {
                DateTime firstAvg = avgOut->getFirstDate();
                if (firstAvg.daysDiff(lastNotify) < 2) {
                    throw ModelException(method,
                                         "With averaging, first avg out ("+
                                         firstAvg.toString() + 
                                         ") must be at least 1 day after "
                                         "last call notice date (" + 
                                         lastNotify.toString()+")");
                }
            }
  
            if (fwdStarting && startDate  > callSchedule.notification[0]) {
                throw ModelException(method, 
                                     "start date (" + startDate.toString() + 
                                     ") is after first call notice date ("+
                                     callSchedule.notification[0].toString()+")");                
            }

            if ((lastCall > valueDate) && (avgOut->getLastDate() != maturity)) {
                throw ModelException(method, 
                                     "last avg out date (" + 
                                     avgOut->getLastDate().toString() + 
                                     ") must equal maturity ("+
                                     maturity.toString()+")");
            }           
        }

        if (hasFloater) {
            DateTime lastFloat = floater->PayDates[floater->PayDates.size()-1];
            if (lastFloat > maturity) {
                throw ModelException(method, 
                                     "last floating pay date (" + 
                                     lastFloat.toString() + 
                                     ") is after maturity ("+
                                     maturity.toString()+")");
            }        
        }  

        if (portfolio.isPerCallDate) {
            // each id is the Nth call date - use user style 1->N 
            // rather than C style 0->N-1
            // N+1 is maturity 
            int numCalls = callSchedule.schedule.size();
            bool hasMaturityPortfolio = false;
            bool hasFinalPortfolio = false;
            for (i = 0; i < portfolio.callDateId->size(); i++) {
                int userid = (*portfolio.callDateId)[i];
                if (userid < 1 || userid > numCalls+1) {
                    throw ModelException(method, "portfolio: call date id " + 
                                         Format::toString(i+1) + "("+
                                         Format::toString(userid)+") "
                                         "is out of range [1, " + 
                                         Format::toString(numCalls+1) +
                                         "] inclusive");
                }
            }
        }

        if( isHit ) {
            if( !hasMonitor() )
                throw ModelException(method, "has no monitor dates but isHit=true"); 
            if( hitDate > valueDate )
                throw ModelException(method, "isHit but hitDate is in future");
        }
        if( hasMonitor() != (!!portfolioIfHit && !portfolioIfHit->isEmpty()) ||
            hasMonitor() != (!!levels) ){
                string strA = (hasMonitor()?"":"no"), strB = (hasMonitor()?"no":"");
                throw ModelException(method, "has " + strA + " monitor dates, but has " + strB + 
                    " levels or portfolioIfHit");
        }
        if( hasMonitor() ) {
            // validation for barrier
            if( !levels->isNonNegative() )
                throw ModelException(method, "barrier level must be non negative");
            
            if( !!economicLevels ) {
                if( !economicLevels->isNonNegative() )
                    throw ModelException(method, "economic barrier level must be non negative");
            } else {
                economicLevels = ScheduleSP(levels.clone());
            }

            if (fwdStarting && startDate > monitorDates->front())
                throw ModelException(method, "instrument is fwd starting with monitor dates, but "
                            "start date ("+ startDate.toString() + ") is > first monitor date (" +
                            monitorDates->front().toString());

            // validation for historical sample
            if( !histMonSamples ) histMonSamples = DoubleArraySP(new DoubleArray(0));
            int numPastSampleDates = monitorDates->size() - valueDate.numFutureDates(*monitorDates);
            if (numPastSampleDates > 0) {
                if (histMonSamples->size() < numPastSampleDates)
                    throw ModelException(method,
                                         "There are " + Format::toString(numPastSampleDates) +
                                         " historical monitoring dates but only "+
                                         Format::toString(histMonSamples->size()) +
                                         " historical sample values have been provided.");
                
                // check that all the historical samples are non-zero
                for (i = 0; i < numPastSampleDates; i++)
                {
                    if (!Maths::isPositive((*histMonSamples)[i]))
                        throw ModelException(method,
                                             "Historical sample " + Format::toString(i+1)+
                                             " for monitoring date " + (*monitorDates)[i].toString()+
                                             " must be greater than zero.");
                    
                    // check for historical hitting status
                    // we know that the trade is not fwd starting anymore and can directly use barrier level
                    if( !isHit ) {
                        double level = levels->interpolate((*monitorDates)[i]);
                        if( (isUp && !Maths::isNegative((*histMonSamples)[i]-level) ) ||
                            (!isUp && !Maths::isPositive((*histMonSamples)[i]-level) ) ) {
                            isHit = true;
                        }
                    }
                }

                if( numPastSampleDates < histMonSamples->size() )
                    histMonSamples->resize(numPastSampleDates);
            }

            if (portfolioIfHit->isPerCallDate) {
                int numCalls = callSchedule.schedule.size();
                for (i = 0; i < portfolioIfHit->callDateId->size(); i++) {
                    int userid = (*portfolioIfHit->callDateId)[i];
                    if (userid < 1 || userid > numCalls+1) {
                        throw ModelException(method, "portfolioIfHit: call date id " + 
                                             Format::toString(i+1) + "("+
                                             Format::toString(userid)+") "
                                             "is out of range [1, " + 
                                             Format::toString(numCalls+1) +
                                             "] inclusive");
                    }
                }
            }
        } // end of checking hasMonitor()
    }

    // Callability::IEventHandler interface
    void getEvents(const Callability* call, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const;

private:
    friend class CallableDepositHelper;
    friend class CDFDProd; 

    CallableDeposit():Generic1Factor(TYPE),guarantee(0.0),hasFloater(false), 
        isUp(false), isHit(false), intraDayMonitor(false) {}; 
    CallableDeposit(const CallableDeposit& rhs);
    CallableDeposit& operator=(const CallableDeposit& rhs);
   
    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        if (CTree1fLV::TYPE->isInstance(model)) {
            return true; // do pointwise instead
        }
        return false;
    }

    /** Returns all strikes the CallableDeposit is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method = "CallableDeposit::getSensitiveStrikes";       
        try {
            DoubleArraySP sensStrikes(new DoubleArray(0));
            if (avoidVegaMatrix(model)) {
                throw ModelException(method, 
                                     "VEGA_MATRIX is not valid for this instrument");
            }
            CVolRequestConstSP volRequest(new ATMVolRequest());

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
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;    
    }

    bool sensShift(Theta* shift) {
        static const string method = "CallableDeposit::sensShift";
        try  {
            const DateTime& newDate = shift->rollDate(valueDate);
            // fill in any samples
            avgOut->roll(shift->getUtil(valueDate), 0, asset.get());
            // and fixings
            if (hasFloater) {
                floater->setFixingforThetaShift(valueDate,
                                                discount.get(),
                                                newDate);
            }
            
            // convert barrier from fwd starting and populate hist mon sample
            if ( hasMonitor() ) {
                if( fwdStarting && newDate.isGreaterOrEqual(startDate) &&
                    startDate.isGreaterOrEqual(valueDate))
                {
                    fwdStarting = false;
                    initialSpot = asset->getThetaSpotOnDate(shift, startDate);
                    levels->scale(initialSpot);
                }
                if (valueDate == newDate) {
                    // only set stuff if it was not set before
                    for (int i=0; i<monitorDates->size(); i++) {
                        if (newDate == (*monitorDates)[i] && Maths::isZero((*histMonSamples)[i])) {
                            (*histMonSamples)[i] = asset->getThetaSpotOnDate(shift, (*monitorDates)[i]);
                        }
                    }
                    
                } else {
                    
                    for (int i=0; i<monitorDates->size(); i++) {
                        if (newDate.isGreaterOrEqual((*monitorDates)[i]) && (*monitorDates)[i].isGreater(valueDate)) {
                            (*histMonSamples)[i] = asset->getThetaSpotOnDate(shift, (*monitorDates)[i]);                    
                        }
                    }
                }
            }

            // roll the parent (updates value date etc)
            Generic1Factor::sensShift(shift);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
        return true; // our components have theta type sensitivity
    }

    /** Satisfy LegalTerms::Shift interface */
    virtual bool sensShift(LegalTerms* shift) {
        // Set the barriers for pricing equal to the economic barriers
        if (!!economicLevels) {
            levels = economicLevels;
        }
        return true;
    }

    /**  a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    bool priceDeadInstrument(CControl* control, CResults* results) const{
        static const string method = "CallableDeposit::priceDeadInstrument";
        try  {
            bool deadInstrument = false;

            if (valueDate >= avgOut->getFirstDate()) {
                // if we're into the averaging region, don't even bother with
                // a tree (as we ignore equity dependence after this date)
                // and just value it all closed form and be done with it

                // cash first
                double value = 0.0;
                for (int i = 0; i < coupons.size(); i++) {
                    if (coupons[i].date > valueDate) {
                        value += notional*coupons[i].amount * discount->pv(coupons[i].date);
                    }
                }
                if (hasFloater) {
                    value -= notional*floater->getPV(valueDate, discount.get());
                }

                // and then maturity portfolio
                int numCalls = callSchedule.notification.size();
                PortfolioSP matPrtf((hasMonitor()&&isHit?*portfolioIfHit:portfolio).makePortfolio(numCalls, numCalls));

                value += portfolioValue(matPrtf.get(), 
                                        initialSpot, 
                                        valueDate, 
                                        notional, 
                                        initialSpot);

                // the above all handle maturity < today < settle
                // and finally the guarantee
                if (valueDate < maturity) {
                    value += notional*guarantee * discount->pv(maturity);
                }

                results->storePrice(value, discount->getCcy());
                if (control && control->isPricing()) {
                    recordOutputRequests(control, results, value);
                }
                              
                deadInstrument = true;  // as far as the tree is concerned anyway
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
                                       maturity,
                                       valueDate,
                                       asset.get());

        OutputRequest* request = NULL;
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTime paymentDate = instSettle->settles(maturity, asset.get());
            DateTimeArray date(1, paymentDate);

            int i;
            for (i = 0; i < coupons.size(); i++) {
                date.push_back(coupons[i].date);
            }
            if (hasFloater) {
                for (i = 0; i < floater->PayDates.size(); i++) {
                    date.push_back(floater->PayDates[i]);
                }        
            }        

            OutputRequestUtil::recordPaymentDates(control,results,&date); 
        }
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            // these need to be in increasing date order to aggregate
            // hence the merge rather than push_back which handles this
            CashFlowArraySP coupCFL(new CashFlowArray(0));
            int i;
            // first coupons
            for (i = 0; i < coupons.size(); i++) {
                coupCFL->push_back(coupons[i]);
                (*coupCFL)[i].amount*= notional;
            }
            // then guarantee
            CashFlow cf(maturity, guarantee*notional);
            CashFlowArraySP gteeCFL(new CashFlowArray(1, cf));

            CashFlowArraySP fltCFL(new CashFlowArray(0));
            // grab any floaters
            if (hasFloater) {
                PayStreamSP fltpay(floater->makePayStream(discount.get()));
                
                CashFlowArrayConstSP floatcf(fltpay->knownCashflows(valueDate,
                                                                    0,
                                                                    false,
                                                                    valueDate,
                                                                    0.0));
                
                for (i = 0; i < floatcf->size(); i++) {
                    fltCFL->push_back((*floatcf)[i]);
                    (*fltCFL)[i].amount *= -notional;
                }              
            }

            CashFlowArraySP prtfCFL(new CashFlowArray(0));;                              
            // if we're past the last sample then handle the maturity portfolio
            if (valueDate.isGreater(avgOut->getLastDate())) {
                int numCalls = callSchedule.notification.size();
                PortfolioSP matPrtf((hasMonitor()&&isHit?*portfolioIfHit:portfolio).makePortfolio(numCalls, numCalls));

                if (!matPrtf.get()) {
                    throw ModelException("CallableDeposit::recordOutputRequests", 
                                         "can't build maturity portfolio");
                }

                // can't reliabaly call portfolioValue as we are  
                // after option maturity so this values to zero

                double prtf = 0.0;
                double avg  = avgOut->averageToDate(valueDate);
                double perf = avg/initialSpot;

                for (i = 0; i < matPrtf->strike.size(); i++) {
                    prtf += matPrtf->longShort[i]*matPrtf->participation[i]*
                        Maths::max(0.0, 
                                   matPrtf->callPut[i]*(perf-matPrtf->strike[i]));
                }

                prtf *= notional;
                CashFlow prtfcf(instSettle->settles(maturity,asset.get()),prtf);
                prtfCFL->push_back(prtfcf);
            }

            // now glue them all together
            CashFlowArraySP merge1(CashFlow::merge(coupCFL, gteeCFL));
            CashFlowArraySP merge2(CashFlow::merge(merge1, fltCFL));
            CashFlowArraySP cfl(CashFlow::merge(merge2, prtfCFL));

            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get());   
        }        

        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request && hasMonitor() && !isHit && !fwdStarting) {
            // report barrier levels over a date range
            DateTime upperDate = BarrierLevel::barrierWindow(valueDate);

            BarrierLevelArraySP lvls(new BarrierLevelArray(0));

            // use economic barrier (if it exists)
            Schedule* s = economicLevels.get() ? economicLevels.get(): levels.get();
            CashFlowArraySP subset(s->subset(valueDate, upperDate));
            for (int i = 0; i < subset->size(); i++) {
                BarrierLevel bl(isUp,(*subset)[i].date,(*subset)[i].amount, intraDayMonitor);
                lvls->push_back(bl);
            }
            
            if (!lvls->empty()) {
                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       lvls.get());
            }
        }            

    }  

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const; 

    // value the portfolio of avg options
    double portfolioValue(const Portfolio*   prtf,
                          double             spot, 
                          const    DateTime& when, 
                          double             refLevel) const {
        static const string method = "CallableDeposit::portfolioValue";
        try {   
            // survive spot = 0
            spot = Maths::max(spot, 0.0001);
            return portfolioValue(prtf, refLevel/spot, when, spot/refLevel, spot);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    double portfolioValue(const Portfolio* prtf,
                          double           strikeScale, 
                          const DateTime&  startDate,
                          double           notl,
                          double           refLevel) const {
        static const string method = "CallableDeposit::portfolioValue";
        try {     
            double        value = 0.0;

            if (!prtf) {
                // null really means no portfolio
                return 0.0;
            }

            CClosedFormLN cfln("VolPreferred");

            for (int i = 0; i < prtf->strike.size(); i++) {
                InstrumentSP avg(Average::makeAvgSpot(prtf->callPut[i] == 1,
                                                      maturity,
                                                      prtf->strike[i]*strikeScale,
                                                      avgOut.get(),
                                                      instSettle.get(),
                                                      premiumSettle.get(),
                                                      asset.get(),
                                                      ccyTreatment,
                                                      discount.get(),
                                                      valueDate,
                                                      startDate > valueDate,
                                                      startDate,
                                                      false, // always notional
                                                      notl,
                                                      refLevel));
                                                   
                value += cfln.calcPrice(avg.get(), 0) * 
                    prtf->longShort[i] * prtf->participation[i];
            }
            return value;                  
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

private:
   
    DateTime      maturity;
    double        guarantee;
    CashFlowArray coupons;
    LiborLegSP    floater;
    bool          hasFloater;
    SampleListSP  avgOut;
    CallSchedule  callSchedule;
    Portfolio     portfolio;

    // optional portfolio that can be KI
    // so if barrier is hit, this REPLACEs the portfolio
    bool          isUp;
    ScheduleSP    levels;            // the schedule of dates and levels
    ScheduleSP    economicLevels;
    DateTimeArraySP     monitorDates; 
    DoubleArraySP       histMonSamples;  // spot price on hist sample dates
    smartPtr<Portfolio> portfolioIfHit;
    bool          isHit; // indicate if barrier was already hit
    DateTime      hitDate;
    bool          intraDayMonitor;
                                  
};

typedef smartPtr<CallableDeposit::Portfolio> CDPortfolioSP;

/*/////////////////////////////////////////////////////////////////////////////////////////////////////*/
/*/////////////////////////////////////////////////////////////////////////////////////////////////////*/

class CDFDProd: public LatticeProdEDRIns
{
public:    
    CDFDProd(const CallableDeposit* cd, FDModel* mdl) : 
        LatticeProdEDRIns(mdl, 1, cd->hasMonitor()&&(!cd->isHit)?2:1),    
        inst( cd )
    {         
        if( ! tree1f )
        {
            throw ModelException( "CDFDProd::CDFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

        // compile list of critical dates
        compileCriticalDates();
    }

    // hold key data for each critical date
    class ExerInfo 
    {
    public:
        bool     isNotification; // is now a call notification date
        DateTime callDate;       // if so, this is the corresponding call date
        bool     isFirstAvgOut;  // is it the first avg out date ?
        // Cost of calling today: includes call amount + coupons paid after
        // today upto and including call date - floating payments received 
        // during the same period all pv'd to today 

        double   pv;             // Discount factor for notification to call date
        double   matPV;          // Discount factor for maturity to settlement
        double   costOfCall;     // issuer's cost of calling
        double   callAmount;     // as percent of notional 
        double   issuerPayment;  // fixed payments - floating payments
        double   monLevel;      // negative if not monitor date
        vector<CDPortfolioSP> notifPortfolios; // array of portfolio for call date
        vector<CDPortfolioSP> matPortfolios; // array of portfolio for maturity

        ExerInfo(): isNotification(false), isFirstAvgOut(false), costOfCall(0.0), 
            callAmount(0.0), issuerPayment(0.0), pv(0.0), matPV(0.0) {
                notifPortfolios.clear();
                matPortfolios.clear();
        }
    };

    //---------------------------------------------------------------------
    // Product FD methods
    //---------------------------------------------------------------------    

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }
   
    /** this sets up the timeline */
    virtual void init(CControl* control) const
    {    
        static const string method = "CDFDProd::Init";
        try 
        {    
            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode =
                ( tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ? numIns : 0 );

            if ( inst->fwdStarting ) 
                tree1f->controlSameGridFwdStart(inst->ccyTreatment);            
              
            // don't bother with dollar divs
            tree1f->SetDivAmountTreatment(false);

            // kill off control var - unlikely to be very useful here
            tree1f->DEBUG_UseCtrlVar = false;

            // add critical dates
            model->addCritDates( criticalDates );

            // define start and end points of tree
            DateTimeArray segDates(2);
            segDates[0] = getStartDate();
            segDates[1] = inst->maturity;

            // timeline configuration
            // 'density factor' for timeline
            IntArray density( 1, 1 );

            // prepare model set up
            model->initSegments( segDates, density );
        }
        catch (exception& e) 
        {
            throw ModelException(e, method);
        }
    }
           
    /** initialising and setting product variables */
    // this is called per pricing call before tree sweep call (after InitTree)
    virtual void initProd() 
    {
        static const string method = "CDFDProd::InitProd";
        try 
        {
            // if recycling timeline, tree doesn't call Init() again
            // but I need critical dates to exist for the mappings below
            if (criticalDates.empty()) 
            {
                compileCriticalDates();
            }

            int i;          
            int lastStep = model->getLastStep();
            
            initSlices( numPrices );
            initInsertNode();

            // set up an array of flags indicating if time step is a critical date or not
            stepIsCritical.resize(lastStep+1);

            // always set maturity point as a critical date
            stepIsCritical[lastStep] = true;
            // and set everything else to false for starters
            for (i = 0; i < lastStep-1; i++) 
            {
                stepIsCritical[i] = false;
            }

            // find first critical date in timeline
            int loIdx = 0;
            int hiIdx = criticalDates.size() - 1;
            while (loIdx < hiIdx && 
                   (criticalDates[loIdx] < model->getDate(0))) 
            {
                loIdx++;
            }

            // maturity is 'magic' i.e. not on timeline per se
            critDateMap[lastStep] = hiIdx;

            for (i = 0; i < lastStep && (loIdx < hiIdx); i++) 
            {
                if (model->getDate(i).equals(criticalDates[loIdx], true)) 
                {
                    critDateMap[i] = loIdx;
                    stepIsCritical[i] = true;
                    loIdx++;
                }
            }   

            // store key product data at each critical date
            exerInfo.resize(criticalDates.size());
            int callIdx  = inst->callSchedule.notification.size()-1;
            int numCalls = inst->callSchedule.notification.size();
            int monitorIdx = (inst->hasMonitor()?inst->monitorDates->size():0)-1;

            for (i = exerInfo.size() - 1; i >= 0; i--) {
                if ((callIdx >= 0) &&  
                    criticalDates[i] == inst->callSchedule.notification[callIdx]) 
                {
                    // Process Notification dates
                    exerInfo[i].isNotification = true;
                    exerInfo[i].callDate   = inst->callSchedule.schedule[callIdx].date;
                    exerInfo[i].callAmount = inst->callSchedule.schedule[callIdx].amount;

                    // calc pv of all payments made after notice up to call date
                    // - coupons paid
                    // - floating payments received
                    // - accrued received on call

                    double pv = inst->discount->pv(criticalDates[i], exerInfo[i].callDate);

                    // Store the pv factor for discounting the payment of the portfolio
                    if( inst->portfolio.payPortfolioOnCallDate ) {
                        exerInfo[i].pv = pv;
                    } else {
                        // Portfolio pays on notification date, so no discount
                        exerInfo[i].pv = 1.0;
                    }                        

                    exerInfo[i].costOfCall = (inst->guarantee+exerInfo[i].callAmount) * pv;
                        
                    exerInfo[i].costOfCall += pvCoupons(criticalDates[i], 
                                                        exerInfo[i].callDate);
                    if (inst->hasFloater) 
                    {
                        exerInfo[i].costOfCall -= inst->floater->getPV(inst->valueDate,
                                                                     criticalDates[i], 
                                                                     exerInfo[i].callDate, 
                                                                     inst->discount.get());
                        
                        if (inst->callSchedule.includeAccrued) 
                        {
                            // subtract accrued
                            exerInfo[i].costOfCall -= inst->floater->getAccrued(inst->valueDate,
                                                                              exerInfo[i].callDate,
                                                                              inst->discount.get())*pv;
                        }
                    }

                    exerInfo[i].notifPortfolios.clear();
                    if( !inst->hasMonitor() || !inst->isHit )
                        exerInfo[i].notifPortfolios.push_back(CDPortfolioSP(inst->portfolio.makePortfolio(callIdx, 
                                                                                      numCalls)));
                    if( inst->hasMonitor() )
                        exerInfo[i].notifPortfolios.push_back(CDPortfolioSP(inst->portfolioIfHit->makePortfolio(callIdx, 
                                                                                      numCalls)));
                    callIdx--;
                }
                
                if (criticalDates[i].equals(avgStart)) 
                {
                    // Process portfolio at maturity
                    exerInfo[i].isFirstAvgOut = true;

                    exerInfo[i].matPV = inst->instSettle->pvAdjust(inst->maturity,
                                                                   inst->discount.get(), 
                                                                   inst->asset.get());

                    // get maturity portfolio
                    exerInfo[i].matPortfolios.clear();
                    if( !inst->hasMonitor() || !inst->isHit )
                        exerInfo[i].matPortfolios.push_back(CDPortfolioSP(inst->portfolio.makePortfolio(numCalls, 
                                                                                      numCalls)));
                    if( inst->hasMonitor() )
                        exerInfo[i].matPortfolios.push_back(CDPortfolioSP(inst->portfolioIfHit->makePortfolio(numCalls, 
                                                                                      numCalls)));
                }

                // and now the issuer payment
                // any fixed or floating coupons today
                // and the guarantee if at maturity
                int j;
                if (inst->hasFloater) 
                {
                    CashFlowArrayConstSP cf(inst->floater->getCashFlowArray(inst->valueDate,
                                                                          inst->discount.get()));

                    for (j = 0; j < cf->size(); j++) 
                    {
                        if ((*cf)[j].date == criticalDates[i]) 
                        {
                            exerInfo[i].issuerPayment -= (*cf)[j].amount;
                        }
                    }
                }

                for (j = 0; j < inst->coupons.size(); j++) 
                {
                    if (inst->coupons[j].date == criticalDates[i]) 
                    {
                        exerInfo[i].issuerPayment += inst->coupons[j].amount;
                    }
                }

                if (criticalDates[i] == inst->maturity) 
                {
                    exerInfo[i].issuerPayment += inst->guarantee;
                }

                // for monitoring, get barrier level
                if (monitorIdx >= 0 &&
                    criticalDates[i] == (*inst->monitorDates)[monitorIdx])
                {
                    exerInfo[i].monLevel = inst->levels->interpolate(criticalDates[i]);
                    if( inst->fwdStarting )
                        exerInfo[i].monLevel *= inst->asset->fwdValue(inst->startDate);
                    monitorIdx--;
                } else 
                {
                    exerInfo[i].monLevel = -1;
                }
            }         

            // and finally the reference level for the portfolio
            refLevel = inst->fwdStarting ? 
                inst->asset->fwdValue(inst->startDate) : inst->initialSpot;
        }
        catch (exception& e) 
        {
            throw ModelException(e, method);
        }
    }

    // here just take care of scaling by notional and any additional 
    // discounting for fwd starting case    
    
    double scalePremium(const double& fairValue,
                        YieldCurveConstSP disc)
    {
        double fwdStartDF = 1.0;
        if (inst->fwdStarting)
            fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));

        return inst->notional*fwdStartDF*fairValue;
    }   
    
    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
    {
        return inst->ccyTreatment;
    }

    /** extra output requests */
    virtual void recordOutput(Control* control, 
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
    virtual CVolRequestConstSP GetLNRequest() const 
    {
        // go in ATM as we're really meant to be LV
        CVolRequestConstSP volRequest(new ATMVolRequest());
        return volRequest;
    }

    //---------------------------------------------------------------------
    // Payoff FD methods
    //---------------------------------------------------------------------
 
    /** called before PayoffAtMat, PayoffBeforeMat and tree roll() */
    virtual void preCalc(int step) 
    {
        int idx = tree1f->getSliceIndex(step);

        if (!stepIsCritical[step])
            return;

        ExerInfo& exer = exerLookup(step);
        adjBarrier = exer.monLevel;
        if( !Maths::isNegative(adjBarrier) ) 
        {
        
            // taking care of discrete monitoring, treat as daily monitor
            if (!inst->intraDayMonitor) 
            {
                vector<double> vol;
                // get vol at barrier
                tree1f->GetStepVol(step, vol, &adjBarrier, 0, 0);
                Barrier::BarrierAdjustment(vol[0], inst->isUp, adjBarrier);
            }
        
            if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION)
                tree1f->SetInsertNode(idx, 0, adjBarrier, 0);
        }
    }

    void mergeMonPriceArray(const double* s, 
                                  int bot, 
                                  int top, 
                                  double level, 
                                  double *price0, 
                                  double *price1)
    {
        if( inst->isUp ) 
        {
            int j=top; 
            while(j>=bot && s[j]>=level) {
                price0[j] = price1[j];
                j--;
            }
        } 
        else 
        {
            int j=bot; 
            while(j<=top && s[j]<=level) {
                price0[j] = price1[j];
                j++;
            }
        }
    }

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type)
    {   
        // we assume just need one und level for spot here
        const TreeSlice & s = payoffIndex->getValue( step );
        int bot, top;
        s.getCalcRange( bot, top );

        const vector< TreeSliceSP > & price = slices;
        int pStart = 0, pEnd = price.size() - 1;

        if (type == FDProduct::BWD_T){
            prod_BWD_T(   s,
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
        else if(type == FDProduct::BWD){
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

    /** product payoff method at maturity */
    virtual void prod_BWD_T(const TreeSlice & spot,
                            int step, 
                            int bot, 
                            int top, 
                            int pStart, 
                            int pEnd,
                            const vector< TreeSliceSP > & price) 
    {
        static const string method("CDFDProd::prod_BWD_T");
        try 
        {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            double prtf = 0.0;

            // get the exer info for this step
            ExerInfo& exer = exerLookup(step);

            // payoff at mat depends on the average.  Since we are using a tree, 
            // we will add the equity participation on the first ave day. 
            // So, value today is cash payment. 
            // If today is first avg date, add equity participation here.
            for (int i=pStart; i<=pEnd; i++) 
            {
                for (int j=bot; j<=top; j++) 
                {
                    // This is the maturity step, so the nodes start with 0 value.
                    (p[i])[j] = 0.0;

                    // if first avg out date add portfolio value to node value
                    if (exer.isFirstAvgOut) 
                    {
                        // XXX This must be a portfolio with no averaging out period,
                        // so we just discount by the precalculated pv according to the
                        // instrument settlement.
                        // We don't use the average spot calculation, since there is only a
                        // single averaging out date.
                        (p[i])[j] += exer.matPortfolios[i]->intrinsic(s[j], refLevel)*exer.matPV;
                    }

                    // adjust node value by cash payments, which are made 
                    // regardless of call.  Value is value to note holder,
                    // so we adjust the value down by the amount received by the holder
                    (p[i])[j] += exer.issuerPayment;                    

                    // if isNotification, issuer will call if it is cheaper.  
                    // That means that if payoff > costOfCall issuer will call 
                    // and pay the equivalent of costOfCall to note holder. 
                    // Note: this decision is made after all other payments are 
                    // made today.  Since we are going backwards in the tree, that
                    // means we should test for call before "undoing" the issuerPayments 
                    if (exer.isNotification) 
                    {
                        if (exer.notifPortfolios[i].get()) 
                        {
                            // Calculate portfolio value and discount to call date (today is notification)
                            prtf = exer.notifPortfolios[i]->intrinsic(s[j], refLevel)*exer.pv;
                        }

                        if (inst->callSchedule.isPuttable) 
                        {
                            (p[i])[j] = Maths::max((p[i])[j], exer.costOfCall + prtf);
                        }
                        else 
                        {
                            (p[i])[j] = Maths::min((p[i])[j], exer.costOfCall + prtf);
                        }                            
                    }
                }
            }

            // if has monitor, check for monitor condition
            if( inst->hasMonitor() && !inst->isHit && !Maths::isNegative(adjBarrier) ) 
            {
                mergeMonPriceArray(s, bot, top, adjBarrier, (p[0]), (p[1]));
            }
        }
        catch(exception& e) 
        {
            throw ModelException(e, method);
        }
    }
    
    /** product payoff method at steps earlier than maturity */
    virtual void prod_BWD(const TreeSlice & spot,
                          int step, 
                          int bot, 
                          int top, 
                          int pStart, 
                          int pEnd,
                          const vector< TreeSliceSP > & price) 
    {
        static const string method("CDFDProd::prod_BWD");
        try 
        {
            double * s = spot.getValues();
            const vector< double * > & p = getValues( price );

            if (!stepIsCritical[step]) 
            {
                return;
            }
            // get the exer info for this step
            ExerInfo& exer = exerLookup(step);

            double prtf = 0.0;

            /* "Early exercise" procedures: Because the instrument supports averaging
               which inherently goes against the use of a tree for pricing, we have to
               make some fudges.  What we do is this: we assume there is no dependence
               on equity prices after the first avg out date.  So, the value after
               the first avg out date is just the pv of cashflows.  Consequently,
               after first avg out date, the only "exercises" are cashflows. On
               these dates, add the value of the cashflow to the node value. On the
               first avg out date, add all the equity dependence.  To do this, 
               increase the node value with the sum of fwd starting averages with 
               strikes equal to strike = (actual strike for the note)/UnderlyingValue
               All cashflows prior to the first avg out date are "exercises", where
               the node value is adjusted by the cashflow.  All call notice dates
               (which must be prior to first avg out date) are real exercises.  For
               these, min node value with pv of call amount, which may include accrued
               interest.
               
               Recap.  Early exercises are:
               - all cashflows
               - value of portfolio on first AveSample date
               - call notice dates (which must be prior to first AveSample date)
               
            */
            double   pvf = 1.0;
            if (exer.isFirstAvgOut) 
            {
                // XXX Calculate the discount of cashflow as calculated from today...
                pvf = inst->discount->pv(avgStart);
            }

            for (int i=pStart; i<=pEnd; i++) 
            {
                for (int j=bot; j<=top; j++) 
                {
                    // if first avg out date add portfolio value to node value
                    if (exer.isFirstAvgOut) 
                    {
                        // XXX ...to cancel it out from the portfolio value calculation
                        (p[i])[j] += inst->portfolioValue(exer.matPortfolios[i].get(), s[j], avgStart, refLevel)/pvf;
                    }

                    // adjust node value by cash payments, which are made 
                    // regardless of call.  Value is value to note holder,
                    // so we adjust the value down by the amount received by the holder
                    (p[i])[j] += exer.issuerPayment;                    

                    // if isNotification, issuer will call if it is cheaper.  
                    // That means that if payoff > costOfCall issuer will call 
                    // and pay the equivalent of costOfCall to note holder. 
                    // Note: this decision is made after all other payments are 
                    // made today.  Since we are going backwards in the tree, that
                    // means we should test for call before "undoing" the issuerPayments 
                    if (exer.isNotification) 
                    {
                        if (exer.notifPortfolios[i].get()) 
                        {
                            // Calculate portfolio value and discount to call date (today is notification)
                            prtf = exer.notifPortfolios[i]->intrinsic(s[j], refLevel)*exer.pv;
                        }

                        if (inst->callSchedule.isPuttable) 
                        {
                            (p[i])[j] = Maths::max((p[i])[j], exer.costOfCall + prtf);
                        }
                        else 
                        {
                            (p[i])[j] = Maths::min((p[i])[j], exer.costOfCall + prtf);
                        }                            
                    }
                }
            }
            
            // if has monitor, check for monitor condition
            if( inst->hasMonitor() && !inst->isHit && !Maths::isNegative(adjBarrier) ) 
            {
                mergeMonPriceArray(s, bot, top, adjBarrier, (p[0]), (p[1]));
            }
        }
        catch(exception& e) 
        {
            throw ModelException(e, method);
        }
    }
    

    //---------------------------------------------------------------------
    // TreeFD::IProduct methods
    //---------------------------------------------------------------------

    // run off default implementations


    // own methods
    double pvCoupons(const DateTime& fromDate, 
                     const DateTime& toDate) const 
    {
        // get the pv of coupons between fromDate & toDate as of fromDate
    
        double pv = 0.0;
        for (int i = 0; i < inst->coupons.size(); i++) 
        {
            const DateTime& paydate = inst->coupons[i].date;
            if (  (paydate >= fromDate && paydate <= toDate) &&
                (paydate > inst->valueDate)  ) 
            {
                pv += inst->coupons[i].amount * inst->discount->pv(fromDate,paydate);
            }
        }
        return pv;
    }

    // get the ExerInfo for a given timepoint
    ExerInfo& exerLookup(int step) 
    {
        CritDateMap::const_iterator iter = critDateMap.find(step);
        if (iter == critDateMap.end()) 
        {
            throw ModelException("CDFDProd::exerLookup",
                                 "no exercise info found for step " +
                                 Format::toString(step));
        }
        return exerInfo[iter->second];
    }

    void compileCriticalDates()
    {
        static const string method("CDFDProd::compileCriticalDates");
        try 
        {       
            // compile list of critical dates
            // already have maturity plus
            // any coupons 
            // any floating pay dates
            // call notification dates
            // and the first avg out date
            // actually want just before first avg out date as we approx the equity
            // dependence using fwd starting options,so if they start on first 
            // avg out dates we get tied in knots as the first fixing is 'past'
        
            criticalDates.resize(inst->coupons.size());
            int i;
            for (i = 0; i < criticalDates.size(); i++) 
            {
                criticalDates[i] = inst->coupons[i].date;
            }
            
            if (inst->hasFloater) 
            {
                for (i = 0; i < inst->floater->getSize(); i++) 
                {
                    criticalDates.push_back(inst->floater->PayDates[i]);
                }
            }

            for (i = 0; i < inst->callSchedule.notification.size(); i++) 
            {
                criticalDates.push_back(inst->callSchedule.notification[i]);
            }
           
            // now if there's just 1 avg out date (which must be maturity)
            // we don't need to worry about pricing fwd starting avg options
            // later on as we just do the intrinsic
            
            DateTime firstAvgOut = inst->avgOut->getFirstDate();
            if (inst->avgOut->numDates(0) == 1) 
            {
                avgStart = firstAvgOut;
            }
            else 
            {
                // get the end of the day before avg out starts
                avgStart = DateTime(firstAvgOut.getDate()-1, DateTime::END_OF_DAY_TIME);
            }

            criticalDates.push_back(avgStart);
            criticalDates.push_back(inst->maturity);

            // add monitor dates if has barrier to be monitored
            if ( inst->hasMonitor() && !inst->isHit )
            {
                for (i = 0; i < inst->monitorDates->size(); i++) 
                {
                    criticalDates.push_back((*inst->monitorDates)[i]);
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

private:
    // for mapping step index to critical date index
    typedef map<int, int> CritDateMap;
      
    const CallableDeposit* inst;
    vector<bool>           stepIsCritical;
    DateTimeArray          criticalDates;
    vector<ExerInfo>       exerInfo;
    CritDateMap            critDateMap;
    double                 refLevel;
    DateTime               avgStart;
    // barrier at a tree slice after discrete monitor adjustment
    double  adjBarrier;
};

/* create a fd payoff tree - new for all fd/tree state variable interface */
 
FDProductSP CallableDeposit::createProduct(FDModel* model) const
{
    return FDProductSP( new CDFDProd(this, model) );
} 

// Callability::IEventHandler interface
void CallableDeposit::getEvents(const Callability* call, 
                                IModel*            model, 
                                const DateTime&    eventDate,
                                EventResults*      events) const{
    static const string method = "CallableDeposit::getEvents";

    // returns ALL callability dates
    // no comment is made on the desirability of calling
    for (int callIdx = 0; callIdx < callSchedule.notification.size(); callIdx++) {
        events->addEvent(new Callability(callSchedule.notification[callIdx], 
                            callSchedule.schedule[callIdx].date,
                            callSchedule.isPuttable ? Callability::PUTTABLE :
                                                      Callability::CALLABLE,
                            callSchedule.schedule[callIdx].amount));
    }
}

class CallableDepositHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Callable Deposit");
        REGISTER(CallableDeposit, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(Callability::IEventHandler);
        EMPTY_SHELL_METHOD(defaultCallableDeposit);
        FIELD(maturity, "maturity");
        FIELD(guarantee, "guarantee");
        FIELD(coupons, "coupons");
        FIELD(floater, "floater");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(hasFloater, "hasFloater");
        FIELD(avgOut, "average-out samples");
        FIELD(callSchedule, "call schedule");
        FIELD(portfolio, "portfolio");
        FIELD(isUp, "isUp");
        FIELD_MAKE_OPTIONAL(isUp);
        FIELD(monitorDates, "monitorDates");
        FIELD_MAKE_OPTIONAL(monitorDates);
        FIELD(levels, "levels");
        FIELD_MAKE_OPTIONAL(levels);
        FIELD(economicLevels, "economicLevels");
        FIELD_MAKE_OPTIONAL(economicLevels);
        FIELD(portfolioIfHit, "portfolioIfHit that replaces portfolio if barrier is hit");
        FIELD_MAKE_OPTIONAL(portfolioIfHit);
        FIELD(isHit, "isHit");
        FIELD_MAKE_OPTIONAL(isHit);
        FIELD(hitDate, "hitDate");
        FIELD_MAKE_OPTIONAL(hitDate);
        FIELD(histMonSamples, "histMonSamples");
        FIELD_MAKE_OPTIONAL(histMonSamples);
        FIELD(intraDayMonitor, "intraDayMonitor");
        FIELD_MAKE_OPTIONAL(intraDayMonitor);
    }

    static IObject* defaultCallableDeposit(){
        return new CallableDeposit();
    }
};

CClassConstSP const CallableDeposit::TYPE = CClass::registerClassLoadMethod(
    "CallableDeposit", typeid(CallableDeposit), CallableDepositHelper::load);

CClassConstSP const CallableDeposit::CallSchedule::TYPE = CClass::registerClassLoadMethod(
    "CallableDeposit::CallSchedule", typeid(CallableDeposit::CallSchedule), CallableDeposit::CallSchedule::load);

CClassConstSP const CallableDeposit::Portfolio::TYPE = CClass::registerClassLoadMethod(
    "CallableDeposit::Portfolio", typeid(CallableDeposit::Portfolio), CallableDeposit::Portfolio::load);

/* for class loading */
bool CallableDepositLoad() {
    return (CallableDeposit::TYPE != 0);
}


DRLIB_END_NAMESPACE

