//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRebate.cpp
//
//   Description : 
//
//   Date        : Oct 2003
//
//
//   $Log: IRebate.cpp,v $
//   Revision 1.13  2005/05/03 11:16:38  qhou
//   init rebatePayAtHit with true
//
//   Revision 1.12  2005/05/03 09:15:01  aswain
//   fix UMR in FlatRebateMaker
//
//   Revision 1.11  2005/04/29 00:02:47  qhou
//   relax validation on rebatePayPeriod dates
//
//   Revision 1.10  2005/04/19 05:29:46  qhou
//   move payAtHit/payDates to flat/schedule rebateMaker respectively. add trivial classes to
//   expose these in IMS without change existing rebate interface
//
//   Revision 1.9  2005/04/15 04:21:42  qhou
//   initialize stub rule in streamrebate. fix core if floatStreamRebate is NULL in getMarket
//
//   Revision 1.8  2004/09/10 18:58:17  snesbitt
//   Adding StreamRebate to IRebate
//
//   Revision 1.7  2004/08/11 09:14:23  Kkitazaw
//   change getKOPV.
//
//   Revision 1.6  2004/08/03 02:20:33  Kkitazaw
//   Added streamRebateMaker into RebateMakerWrapper.
//
//   Revision 1.5  2004/07/29 14:12:20  jmusset
//   Added discount field to the IRebate::getRebate function
//
//   Revision 1.4  2004/02/18 17:25:47  mrobson
//   Use new macros for interfaces
//
//   Revision 1.3  2004/02/16 13:37:43  snesbitt
//   ScheduleRebate now takes dates/values/interp instead of Schedule
//   to better fit in IMS. Also cleaned up validation.
//
//   Revision 1.2  2003/12/31 12:02:09  snesbitt
//   add ScheduleRebate class
//
//   Revision 1.1  2003/10/10 17:25:12  snesbitt
//   Generalised rebate
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/IRebate.hpp"
#include "edginc/BarrierUtil.hpp"
#include "edginc/MCPathGenerator.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/KOStubRule.hpp"
#include "edginc/KOFixedLeg.hpp"
#include "edginc/KOLiborLeg.hpp"

DRLIB_BEGIN_NAMESPACE

void IRebateMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IRebateMaker, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IRebateMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "IRebateMaker", typeid(IRebateMaker), load);

class FlatRebate : virtual public IRebate {
public:
    FlatRebate(FlatRebateMaker* maker, const DateTimeArray& simDates):
        maker(FlatRebateMakerSP(copy(maker))) {
    };

    virtual double getLevel(int simDateIdx) const {
        // a constant
        return maker->rebateAmount;
    };

    virtual DateTime getPayDate(const DateTime& hitDate) const {
        return hitDate;
    };

private:
    FlatRebateMakerSP maker;
};

typedef refCountPtr<FlatRebate> FlatRebateSP;

// The published face
FlatRebateMaker::FlatRebateMaker(double rebateAmount): 
    CObject(TYPE), rebateAmount(rebateAmount) {}

FlatRebateMaker::FlatRebateMaker(CClassConstSP clazz):
    CObject(clazz), rebateAmount(0.0) {}

IRebateSP FlatRebateMaker::getRebate(const DateTimeArray& simDates,
                                     const YieldCurve*    discount)
{
    // independent of the sim dates
    return FlatRebateSP(new FlatRebate(this, simDates));
}

void FlatRebateMaker::getMarket(const IModel*     model, 
                                const MarketData* market) {}

double FlatRebateMaker::getRebateAmount() const{
    return rebateAmount;
}

class FlatRebateMakerHelper{
public:
    static IObject* defaultFlatRebateMaker(){
        return new FlatRebateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FlatRebateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRebateMaker);
        FIELD(rebateAmount,          "rebate amount");
        FIELD_MAKE_OPTIONAL(rebateAmount);
        EMPTY_SHELL_METHOD(defaultFlatRebateMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const FlatRebateMaker::TYPE =
CClass::registerClassLoadMethod("FlatRebateMaker", 
                                typeid(FlatRebateMaker), FlatRebateMakerHelper::load);


/*********************************************************************/

class ScheduleRebate : virtual public IRebate {
public:
    ScheduleRebate(ScheduleRebateMaker* maker,
                   const DateTimeArray& simDates) {
        static const string method = "ScheduleRebate::ScheduleRebate";

        Schedule sched(maker->rebateSchedDates,
                       maker->rebateSchedValues,
                       maker->rebateSchedInterp);

#if 0
        // Check the range of the rebate schedule is sufficient
        // Tried to use sched.coversDateRange but it's kind of meaningless
        // Instead do : if N type then make sure schedule has dates corresponding
        // to all simDates; if L/S make sure range of rebateSchedDates covers
        // simDates
        if (maker->rebateSchedInterp == INTERP_NONE) {
            // This requires that the simDates are in order
            // Assume rebateSchedDates are in order
            if (maker->rebateSchedDates.size() < simDates.size()) {
                throw ModelException(method,
                                     "For interp " + INTERP_NONE +
                                     " require at least as many dates in schedule as will need a value : " +
                                     Format::toString(maker->rebateSchedDates.size()) + " given but need " +
                                     Format::toString(simDates.size()));
            }
            int s=0;
            int r=0;
            while (s<simDates.size()) {
                while (r<maker->rebateSchedDates.size()) {
                    if (simDates[s] < maker->rebateSchedDates[r]) {
                        // missed one
                        throw ModelException(method,
                                             simDates[s].toString + " is missing from schedule");
                    }
                    if (simDates[s] == maker->rebateSchedDates[r]) {
                        // match - so move on both
                        s++;
                        r++;
                    } else if (simDates[s] > maker->rebateSchedDates[r]) {
                        // extra rebate sched date is ok - keep looking
                        r++;
                    }
                }
            }
        } else {
            if (maker->rebateSchedDates.front() > simDates.front() ||
                maker->rebateSchedDates.back() < simDates.back()) {
                throw ModelException("ScheduleRebate::ScheduleRebate",
                                     "Schedule does not cover required date range between " +
                                     simDates.front().toString() + " and " +
                                     simDates.back().toString());
            }
        }
#endif
        // To Allow the simDates are not corresponding to rebate schedule,
        // intoroducing bool array.
        levels = DoubleArray(simDates.size());
        isRebateStep.resize(simDates.size(),false);
        for (int i=0; i<levels.size(); i++) {
            try {
                levels[i] = sched.interpolate(simDates[i]);
                isRebateStep[i] = true;
            }
            catch (exception& e){
                isRebateStep[i] = false;
            }            
        }

        if( !maker->rebatePayPeriodDates.empty() && maker->rebatePayPeriodDates.back() < simDates.back() )
            throw ModelException("ScheduleRebate::ScheduleRebate",
                                 "If has pay period, last pay period date must be on or after last monitor date");
        payDates = maker->rebatePayPeriodDates;
    };

    virtual double getLevel(int simDateIdx) const {
        if (isRebateStep[simDateIdx])
            return levels[simDateIdx];
        else
            throw ModelException("ScheduleRebate", "There is no interpolated level or rebate found");
    };

    virtual DateTime getPayDate(const DateTime& hitDate) const {
        if( payDates.empty() ) return hitDate;
        
        for( int i=0; i<payDates.size(); i++) {
            if( hitDate <= payDates[i] ) return payDates[i]; 
        }

        return DateTime(0,0); // should never reach this step actually
    };

private:
    DoubleArray levels; // need only this
    DateTimeArray payDates;
    vector<bool> isRebateStep;
};

typedef refCountPtr<ScheduleRebate> ScheduleRebateSP;

ScheduleRebateMaker::ScheduleRebateMaker(CClassConstSP clazz):
    CObject(clazz){}

void ScheduleRebateMaker::validatePop2Object() {
    if( !rebatePayPeriodDates.empty() ) {
        DateTime::ensureIncreasing(rebatePayPeriodDates, "ScheduleRebateMaker::validatePop2Object."
            "pay period dates must in strict ascending order", true);   
    }
}

IRebateSP ScheduleRebateMaker::getRebate(const DateTimeArray& simDates,
                                         const YieldCurve*    discount)
{
    return ScheduleRebateSP(new ScheduleRebate(this, simDates));
}

void ScheduleRebateMaker::getMarket(const IModel*     model, 
                                    const MarketData* market) {}

ScheduleSP ScheduleRebateMaker::makeSchedule() const{
    ScheduleSP mySchedule = ScheduleSP(new Schedule(rebateSchedDates,rebateSchedValues,rebateSchedInterp));
    return mySchedule;
}

class ScheduleRebateMakerHelper{
public:
    static IObject* defaultScheduleRebateMaker(){
        return new ScheduleRebateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ScheduleRebateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRebateMaker);
        FIELD(rebateSchedDates,    "Schedule Dates");
        FIELD(rebateSchedValues,   "Schedule Values");
        FIELD(rebateSchedInterp,   "Interpolation Type");
        FIELD(rebatePayPeriodDates,"Pay period end dates");
        FIELD_MAKE_OPTIONAL(rebateSchedDates);
        FIELD_MAKE_OPTIONAL(rebateSchedValues);
        FIELD_MAKE_OPTIONAL(rebateSchedInterp);
        FIELD_MAKE_OPTIONAL(rebatePayPeriodDates);
        EMPTY_SHELL_METHOD(defaultScheduleRebateMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const ScheduleRebateMaker::TYPE =
CClass::registerClassLoadMethod("ScheduleRebateMaker", 
                                typeid(ScheduleRebateMaker), ScheduleRebateMakerHelper::load);

/*********************************************************************/
/** Rebate is "FloatStream + FixedStream + Fixed Amount"             */
/*********************************************************************/

class StreamRebate : virtual public IRebate {
public:
    StreamRebate(StreamRebateMaker*   maker,
                 const DateTimeArray& simDates,
                 const YieldCurve*    discount) {
        static const string method = "StreamRebate::StreamRebate";
        
        // No need for validation code since trying to interpolate will catch any problems.
        Schedule sched(maker->rebateSchedDates,
                       maker->rebateSchedValues,
                       maker->rebateSchedInterp);
        levels = DoubleArray(simDates.size());


        KOStubRuleSP koFltRule = KOStubRuleSP(new KOStubRule(maker->floatKoStubRule,false,"AsSchedule"));    //koStubRule,isAccrueUpToSettle,payTimingRule
        KOLiborLegSP koFlt = KOLiborLegSP(new KOLiborLeg(maker->floatStreamRebate,koFltRule,simDates[0],discount));
        
        KOStubRuleSP koFixRule = KOStubRuleSP(new KOStubRule(maker->fixedKoStubRule,false,"AsSchedule"));    //koStubRule,isAccrueUpToSettle,payTimingRule        
        KOFixedLegSP koFix = KOFixedLegSP(new KOFixedLeg(maker->fixedStreamRebate,koFixRule));

        for (int i=0; i<levels.size(); i++) {
            double fltVal = (!maker->floatStreamRebate) ? 0.0 :
                  koFlt->getKOPV(simDates[i]);
            double fixVal = (!maker->fixedStreamRebate)? 0.0 :
                  koFix->getKOPV(simDates[i],discount);

            // These are the values at the hit date
            levels[i] = fltVal + fixVal + 
                sched.interpolate(simDates[i]);
        }
    };

    virtual double getLevel(int simDateIdx) const {
        return levels[simDateIdx];
    };

    virtual DateTime getPayDate(const DateTime& hitDate) const {
        throw ModelException("StreamRebate::getPayDate", "Not supported");
    };

private:

    DoubleArray levels; // need only this
};

typedef refCountPtr<StreamRebate> StreamRebateSP;

StreamRebateMaker::StreamRebateMaker()
: CObject(TYPE), floatKoStubRule("N"), fixedKoStubRule("N")
{} 

IRebateSP StreamRebateMaker::getRebate(const DateTimeArray& simDates,
                                       const YieldCurve*    discount) 
{
    return StreamRebateSP(new StreamRebate(this, simDates, discount));
}

// populate from market cache 
void StreamRebateMaker::getMarket(const IModel*     model, 
                                  const MarketData* market) {
    if( !!floatStreamRebate )
        floatStreamRebate->getMarket(model, market);
}

class StreamRebateMakerHelper{
public:
    static IObject* defaultStreamRebateMaker(){
        return new StreamRebateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(StreamRebateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRebateMaker);
        FIELD(floatStreamRebate,          "Floating rebate stream");
        FIELD(floatKoStubRule,     "KO Stub rule for floating stream");
        FIELD(fixedStreamRebate,          "Fixed rebate stream");
        FIELD(fixedKoStubRule,     "KO Stub rule for fixed stream");
        FIELD(rebateSchedDates,    "Schedule Dates");
        FIELD(rebateSchedValues,   "Schedule Values");
        FIELD(rebateSchedInterp,   "Interpolation Type");
        // IMS imposes certain restrictions...
        FIELD_MAKE_OPTIONAL(floatStreamRebate);
        FIELD_MAKE_OPTIONAL(floatKoStubRule);
        FIELD_MAKE_OPTIONAL(fixedStreamRebate);
        FIELD_MAKE_OPTIONAL(fixedKoStubRule);
        FIELD_MAKE_OPTIONAL(rebateSchedDates);
        FIELD_MAKE_OPTIONAL(rebateSchedValues);
        FIELD_MAKE_OPTIONAL(rebateSchedInterp);
        EMPTY_SHELL_METHOD(defaultStreamRebateMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const StreamRebateMaker::TYPE =
CClass::registerClassLoadMethod("StreamRebateMaker", 
                                typeid(StreamRebateMaker), StreamRebateMakerHelper::load);


/******************************************************************/
// class to Flat/ScheduleRebate to expose payAtHit/payPeriodDates in IMS
/******************************************************************/

class FlatRebateMaker2 : public FlatRebateMaker {    
public:
    static CClassConstSP const TYPE;

    FlatRebateMaker2() : FlatRebateMaker(TYPE){};
    ~FlatRebateMaker2(){};

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(FlatRebateMaker2, clazz);
        SUPERCLASS(FlatRebateMaker);
        EMPTY_SHELL_METHOD(defaultFlatRebateMaker2);
    }

    static IObject* defaultFlatRebateMaker2(){
        return new FlatRebateMaker2();
    }
};

CClassConstSP const FlatRebateMaker2::TYPE = CClass::registerClassLoadMethod(
    "FlatRebateMaker2", typeid(FlatRebateMaker2), FlatRebateMaker2::load);


class ScheduleRebateMaker2 : public ScheduleRebateMaker {    
public:
    static CClassConstSP const TYPE;

    ScheduleRebateMaker2() : ScheduleRebateMaker(TYPE){};
    ~ScheduleRebateMaker2(){};

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(ScheduleRebateMaker2, clazz);
        SUPERCLASS(ScheduleRebateMaker);
        EMPTY_SHELL_METHOD(defaultScheduleRebateMaker2);
    }

    static IObject* defaultScheduleRebateMaker2(){
        return new ScheduleRebateMaker2();
    }
};

CClassConstSP const ScheduleRebateMaker2::TYPE = CClass::registerClassLoadMethod(
    "ScheduleRebateMaker2", typeid(ScheduleRebateMaker2), ScheduleRebateMaker2::load);


/*********************************************************************/

/* The ultimate wrapping of IRebateMaker's mainly for use in Pyramid 
 */
#define REBATE_TYPE_FLAT     "Flat"
#define REBATE_TYPE_SCHEDULE "Schedule"
#define REBATE_TYPE_STREAM   "Stream"

class RebateMakerWrapper : public CObject,
                           virtual public IRebateMaker,
                           virtual public Theta::Shift {
public: // how can I have this protected or private?
    string                   rebateMakerType;
    FlatRebateMakerSP        flatMaker;
    ScheduleRebateMakerSP    scheduleRebateMaker;
    StreamRebateMakerSP      streamRebateMaker;

private:
    IRebateMakerSP  realMaker; // $unregistered

public:
    static CClassConstSP const TYPE;

    virtual IRebateSP getRebate(const DateTimeArray& simDates,
                                const YieldCurve*    discount) {
        return realMaker->getRebate(simDates, discount);
    }

    virtual void getMarket(const IModel*     model, 
                           const MarketData* market) {
        realMaker->getMarket(model, market);
    }

    // Implementation of the Theta shift interface
    virtual bool sensShift(Theta* shift) {
        // This class needs to filter all
        // info through the selector field
        if (rebateMakerType==REBATE_TYPE_STREAM) {
            return true; // continue tweaking
        } 
        return false; // nothing to do so avoid trying with "dummy" streamRebateMaker
    }
    
    // validation
    void validatePop2Object(){
        static const string routine = "RebateMakerWrapper::validatePop2Object";

        if (rebateMakerType.empty()){
            throw ModelException(routine, "Blank Rebate Maker specified!");
        }
        if (rebateMakerType==REBATE_TYPE_FLAT) {
            if (flatMaker.get()) {
                realMaker = flatMaker;
            } else {
                throw ModelException(routine, "Expected flat rebateMaker but none supplied!");
            }
        } else if (rebateMakerType==REBATE_TYPE_SCHEDULE) {
            if (scheduleRebateMaker.get()) {
                realMaker = scheduleRebateMaker;
            } else {
                throw ModelException(routine, "Expected schedule rebateMaker but none supplied!");
            }
        } else if (rebateMakerType==REBATE_TYPE_STREAM) {
            if (streamRebateMaker.get()) {
                realMaker = streamRebateMaker;
            } else {
                throw ModelException(routine, "Expected stream rebateMaker but none supplied!");
            }
        } else {
            throw ModelException(routine, "Unrecognised Rebate Maker " + rebateMakerType);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RebateMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRebateMaker);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultRebateMakerWrapper);
        FIELD(rebateMakerType, "Flat, Schedule or Stream");
        FIELD(flatMaker,  "Flat Rebate Maker");
        FIELD_MAKE_OPTIONAL(flatMaker);
        FIELD(scheduleRebateMaker,  "Schedule Rebate Maker");
        FIELD_MAKE_OPTIONAL(scheduleRebateMaker);
        FIELD(streamRebateMaker,  "Stream Rebate Maker");
        FIELD_MAKE_OPTIONAL(streamRebateMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    RebateMakerWrapper(): CObject(TYPE){}

    RebateMakerWrapper(CClassConstSP clazz): CObject(clazz){}

    static IObject* defaultRebateMakerWrapper(){
        return new RebateMakerWrapper();
    }

    /** Override clone method to make sure realMaker references
        the maker inside the clone and not the original */
    IObject* clone() const{
        // first clone all the registered fields
        IObject*  copy = CObject::clone();
        RebateMakerWrapper* rmw = dynamic_cast<RebateMakerWrapper*>(copy);
        if (!rmw){
            throw ModelException("RebateMakerWrapper::clone"); // shouldn't happen
        }
        rmw->validatePop2Object();
        return copy;
    }

};

typedef smartPtr<RebateMakerWrapper> RebateMakerWrapperSP;

CClassConstSP const RebateMakerWrapper::TYPE =
CClass::registerClassLoadMethod("RebateMakerWrapper", 
                                typeid(RebateMakerWrapper), load);

/******************************************************************/
// Wrapper to be able for alternative view of rebatemakerwrapper
/******************************************************************/

class RebateMakerWrapper2 : public RebateMakerWrapper {    
public:
    static CClassConstSP const TYPE;

    RebateMakerWrapper2() : RebateMakerWrapper(TYPE){};
    ~RebateMakerWrapper2(){};

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(RebateMakerWrapper2, clazz);
        SUPERCLASS(RebateMakerWrapper);
        EMPTY_SHELL_METHOD(defaultRebateMakerWrapper2);
    }

    static IObject* defaultRebateMakerWrapper2(){
        return new RebateMakerWrapper2();
    }
    
};

CClassConstSP const RebateMakerWrapper2::TYPE = CClass::registerClassLoadMethod(
    "RebateMakerWrapper2", typeid(RebateMakerWrapper2), RebateMakerWrapper2::load);



/******************************************************************/
// Function for BarrierUtil/BarrierPay
/******************************************************************/


// rebate payment is once hit for KO or never hit for KI
class RebateBarrierPay : public IBarrierPay
{
public:
    // this is the standard rebate for KO barrier. if isOut, pay once when hit,
    // for KI barrier, pay rebate if not KI at end
    // in both case, only pay if hitValue = 0, ie payIfHit=false
    RebateBarrierPay(const IRebateMaker* rebateMaker, bool isOut, const InstrumentSettlement* instSettle)
        : rebateMaker(copy(rebateMaker)), isOut(isOut), settle(copy(instSettle))
    {
        if( !rebateMaker )
            throw ModelException("RebateBarrierPay", "empty rebate maker");
    }
    ~RebateBarrierPay(){};

    void preprocess(DateTime valueDt, DateTimeArrayConstSP simDates, const DateTimeArray& monitorDates, const YieldCurve* discount)
    {
        relevantMonDts = isOut?monitorDates:DateTimeArray(1, monitorDates.back());
        IRebateSP rebate(rebateMaker->getRebate(relevantMonDts, discount));

        payDtPerStep.resize(simDates->size());
        dfs.resize(simDates->size());
        levels.resize(simDates->size());
        for(int i=0, j=0; i<simDates->size(); i++)
        {
            if( j < relevantMonDts.size() && (*simDates)[i] == relevantMonDts[j] )
            {
                payDtPerStep[i] = settle->settles(rebate->getPayDate((*simDates)[i]), 0);
                dfs[i] = (payDtPerStep[i]<=valueDt)?0.0:discount->pv(valueDt, payDtPerStep[i]);
                levels[i] = rebate->getLevel(j);
                j++;
            }
            else
            {
                payDtPerStep[i] = DateTime(0, 0);
                dfs[i] = 0.0;
                levels[i] = 0.0;
            }
        }
    }

    bool payOnce() const 
    { return isOut; } // payOnce if isOut, otherwise pay if not KI

    bool payIfHitOne() const 
    { return false; }

    // we don't need to observe any spot price for rebate calculation
    void payObsDates(DateTimeArray &obsDts, bool &hasContinuous) const
    {
        obsDts.clear();
        hasContinuous = false; 
    }

    // all rebate dates for now
    DateTimeArray criticalMonDates(const DateTimeArray& monitorDates, 
                                    bool isContinuous) const
    { return isOut?monitorDates:DateTimeArray(1, monitorDates.back()); }

    const DateTimeArray& relevantMonDates() const
    { return relevantMonDts; }

    double value(int step, DateTime valueDate, const IMCPathGenerator* pathGen)
    { 
        // add to value only if payment in future. fv=0 if paydate in past
        return dfs[step] * levels[step];
    }

    void collect(int step, const IMCPathGenerator* pathGenIn, double scalingFactor, 
                 CashFlowArray& knownCFs, PhysicalDeliveryArray& phyDs)
    {
        knownCFs.push_back(CashFlow(payDtPerStep[step], scalingFactor*levels[step]));
    }

private:
    IRebateMakerSP      rebateMaker;
    bool                isOut;
    DateTimeArraySP     obsDates;
    InstrumentSettlementSP settle;

    // *** transient fields below ***
    // pay dates with same size as sim dates
    DoubleArray         levels;
    DoubleArray         dfs;
    DateTimeArray       relevantMonDts;
    DateTimeArray       payDtPerStep;
};

refCountPtr<IBarrierPay> IRebateMaker::createBarrierPay(const IRebateMaker* rebateMaker, 
                                                     bool isOut, 
                                                     const InstrumentSettlement* instSettle)
{
    return refCountPtr<IBarrierPay>(new RebateBarrierPay(rebateMaker, isOut, instSettle));
}

DRLIB_END_NAMESPACE

    

