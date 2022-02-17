//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ITimeAggregate.cpp
//
//   Description : Captures collapse of TIME dimension from array to single value in various ways
//
//   Date        : Mar 04
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/ITimeAggregate.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

void ITimeAggregateMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(ITimeAggregateMaker, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const ITimeAggregateMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "ITimeAggregateMaker", typeid(ITimeAggregateMaker), load);

/*********************************************************************/
// Class to help support PAYMENT_DATES and KNOWN_CASHFLOWS request
class TimeAggKnownFlows {
public:
    TimeAggKnownFlows(const DateTimeArray& paymentDates) : 
        uniqPayDates(new DateTimeArray(0)), 
        knownFlows(new CashFlowArray(0)), 
        mapToUniq(new IntArray(paymentDates.size(), 0)) {
        if (paymentDates.size()>0) {
            // cope also with out-of-order pay dates
            DateTimeArray dl = paymentDates;
            sort(dl.begin(), dl.end());
            uniqPayDates->push_back(dl[0]);
    
            for (int i = 1; i < dl.size(); i++) {
                if (dl[i] > dl[i-1]) {
                    uniqPayDates->push_back(dl[i]);
                }
            } 
            // inefficient but done only "once"
            for(int k=0; k<paymentDates.size(); k++) {
                for(int j=0; j<uniqPayDates->size(); j++) {
                    if (paymentDates[k] == (*uniqPayDates)[j]) {
                        (*mapToUniq)[k] = j;
                        break;
                    }
                }
            }
        }
    }
    
    void addFlow(int             idx,       // trusting caller
                 double          howMuch) {
        CashFlow cf((*uniqPayDates)[(*mapToUniq)[idx]], howMuch);
        CashFlowArraySP cfa(new CashFlowArray(0));
        cfa->push_back(cf);
        if (knownFlows.get()) {
            knownFlows = CashFlow::merge(cfa, knownFlows);
        } else {
            knownFlows = cfa;
        }
    }
    
    const DateTimeArray* getPayDates() const {
        return uniqPayDates.get();
    }

    const CashFlowArray* getFlows() const {
        return knownFlows.get();
    }
    
private:
    DateTimeArraySP   uniqPayDates;
    CashFlowArraySP   knownFlows;
    IntArraySP        mapToUniq; // orig -> uniq
};

/*********************************************************************/
// XXX Do we need to worry about settlement? -> yes and done
// XXX What about PAYMENT_DATE requests? -> should be sorted ... TEST!
// XXX If date is in the past? -> caching introduced ... TEST!
// XXX Pay at hit?
// XXX PV to today, or FV to mat? 
// XXX If some past, what about pre-computation? Current IAggregate leaves such caching to product. -> sorted ... TEST!
///////////
// So here we have :
// explicit list of payment dates, but value to final one 
class SumPVAggregate : virtual public ITimeAggregate {
public:
    SumPVAggregate(SumPVAggregateMaker* maker,
                   const DoubleArray&   fvFactors,
                   IDoubleArray*        components,
                   const DateTime&      today,
                   const DateTimeArray& realPayDates):
        maker(maker),
        fvs(fvFactors), 
        components(components),
        knownFlows(realPayDates),
        donePast(false),
        valPast(0.0) {

        for(iFuture=0; iFuture<maker->aggDates.size(); iFuture++) {
            if (maker->aggDates[iFuture]>today) {
                break;
            }
        }
    };
    
    virtual double aggregate() {
        if (!donePast) {
            for(int i=0; i<iFuture; i++) {
                valPast += fvs[i] * (*components)[i];
                // Fill in knownFlows for known flows!
                // Note no use of fvs[]
                knownFlows.addFlow(i, (*components)[i]);
            }
            donePast = true;
        }
        double val = valPast;
        for(int i=iFuture; i<components->size(); i++) {
            val += fvs[i] * (*components)[i];
        }
        return val;
    };

    // To support the PAYMENT_DATES request
    virtual const DateTimeArray* getPaymentDates() const {
        return knownFlows.getPayDates();
    };

    // To support the KNOWN_CASHFLOWS request
    virtual const CashFlowArray* getKnownFlows() const {
        return knownFlows.getFlows();
    };

private:
    SumPVAggregateMaker*     maker;
    DoubleArray              fvs;
    IDoubleArray*            components;
    TimeAggKnownFlows        knownFlows;

    int                      iFuture; // first future index
    bool                     donePast;
    double                   valPast;
};

void SumPVAggregateMaker::validatePop2Object(){
    static const string routine = "SumPVAggregateMaker::validatePop2Object";
    if (aggDates.size()==0) {
        throw ModelException(routine,
                             "Must have at least one aggregation date (none supplied)");
    }
    if (aggDates.size() != aggPayDates.size()) {
        throw ModelException(routine,
                             "Must have equal numbers of aggregate dates (#= " +
                             Format::toString(aggDates.size()) + 
                             ") and pay dates (#= " +
                             Format::toString(aggPayDates.size()) + ")");
    }
    for (int i=0; i<aggDates.size(); i++) {
        if (aggDates[i] > aggPayDates[i]) {
            throw ModelException(routine,
                                 "Pay dates must not be earlier than aggregate dates. Check date #" +
                                 Format::toString(i+1) + " = " +
                                 aggDates[i].toString() + " but pays at " +
                                 aggPayDates[i].toString());
        }
        if (i>0 &&
            aggDates[i-1] >= aggDates[i]) {
            throw ModelException(routine,
                                 "Aggregate dates must be strictly increasing. Check date #" +
                                 Format::toString(i) + " = " +
                                 aggDates[i-1].toString() + " versus date #" +
                                 Format::toString(i+1) + " = " +
                                 aggDates[i].toString());
        }
    }
}

const DateTimeArray& SumPVAggregateMaker::getDates() const {
    return aggDates;
}

// It's possible a lot of this could be held within the product and reduce
// the need for prodView - there is perhaps too much movement of data here. XXX
ITimeAggregate* SumPVAggregateMaker::getAggregate(const DateTime&            valueAtDate,
                                                  const ITimeAggProductView* prodView,
                                                  IDoubleArray*              comps) {
    static const string routine = "SumPVAggregateMaker::getAggregate";
    int len = comps->size();

    if (aggDates.size() != len) {
        throw ModelException(routine, 
                             "Number of dates (" + Format::toString(aggDates.size()) +
                             ") must equal number of components (" + 
                             Format::toString(comps->size()) + ")");
    }

    DoubleArray fvFactors(len, 1.0);
    // This is when it is all over
    DateTime finalSettleDate = prodView->settles(valueAtDate);
    DateTimeArray myPayDates; // possibly with settlement included
    if (adjustPayDatesForInstSettlement) {
        // This is how it was done historically - painful and
        // confusing. Left here for backwards compatability
        // allowing a gradual switch to the non-adjustment approach.
        DateTime finalAggSettleDate = prodView->settles(aggPayDates.back());
        if (finalAggSettleDate > finalSettleDate) {
            throw ModelException(routine, 
                                 "Aggregate pay dates cannot be after final settlement date but " +
                                 aggPayDates.back().toString() + " settles on " +
                                 finalAggSettleDate.toString() + " which is after " +
                                 finalSettleDate.toString());
        }
        myPayDates.resize(len);
        for (int i=0; i<len; i++) {
            myPayDates[i] = prodView->settles(aggPayDates[i]);
        }
    } else {
        // Prefer no adjustment but allow a choice to ease the 
        // upgrade in prod. New instruments should not adjust.
        const DateTime& finalAggSettleDate = aggPayDates.back();
        if (finalAggSettleDate > finalSettleDate) {
            throw ModelException(routine, 
                                 "Aggregate pay dates cannot be after final settlement date but " +
                                 aggDates.back().toString() + " settles on " +
                                 finalAggSettleDate.toString() + " which is after " +
                                 finalSettleDate.toString());
        }
        myPayDates = aggPayDates;
    }
    const YieldCurve* yc = prodView->getYieldCurve();
    const DateTime& today = prodView->getValueDate();
    for (int i=0; i<len; i++) {
        if (myPayDates[i] <= today) {
            fvFactors[i] = 0.0;
        } else {
            fvFactors[i] = yc->pv(finalSettleDate, myPayDates[i]); // note the parameter order
        }
    }
    
    return new SumPVAggregate(this, fvFactors, comps,
                              today, myPayDates);
}

class SumPVAggregateMakerHelper{
public:
    static IObject* defaultSumPVAggregateMaker(){
        return new SumPVAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SumPVAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITimeAggregateMaker);
        EMPTY_SHELL_METHOD(defaultSumPVAggregateMaker);
        FIELD(aggDates,          "Dates for which values will be provided");
        FIELD(aggPayDates,       "Actual payment dates (poss settlement adjust)");
        FIELD(adjustPayDatesForInstSettlement, "[Optional]Default = true");
        FIELD_MAKE_OPTIONAL(adjustPayDatesForInstSettlement);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const SumPVAggregateMaker::TYPE =
CClass::registerClassLoadMethod("SumPVAggregateMaker", 
                                typeid(SumPVAggregateMaker), SumPVAggregateMakerHelper::load);

/*********************************************************************/
/*********************************************************************/


/* The ultimate wrapping of IDoubleArrayModifierMaker's mainly for use in Pyramid 
 */
#define TIME_AGGREGATE_TYPE_SUMPV    "SumPV"

class TimeAggregateMakerWrapper : public CObject,
                                  virtual public ITimeAggregateMaker {
public: 
    string                           timeAggregateMakerType;
    SumPVAggregateMakerSP            sumPVMaker;

private:
    ITimeAggregateMakerSP  realMaker;

public:
    static CClassConstSP const TYPE;

    virtual const DateTimeArray& getDates() const {
        return realMaker->getDates();
    }

    virtual ITimeAggregate* getAggregate(const DateTime&            valueAtDate,
                                         const ITimeAggProductView* prodView,
                                         IDoubleArray*              comps){
        return realMaker->getAggregate(valueAtDate, prodView, comps);
    }

    // validation
    void validatePop2Object(){
        static const string routine = "TimeAggregateMakerWrapper::validatePop2Object";

        if (timeAggregateMakerType.empty()){
            throw ModelException(routine, "Blank Time Aggregate Maker specified!");
        }
        if (timeAggregateMakerType==TIME_AGGREGATE_TYPE_SUMPV) {
            if (sumPVMaker.get()) {
                realMaker = sumPVMaker;
            } else {
                throw ModelException(routine, "Expected simple SumPV but none supplied!");
            }
        } else {
            throw ModelException(routine, "Unrecognised Time Aggregate Maker " + timeAggregateMakerType);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TimeAggregateMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITimeAggregateMaker);
        EMPTY_SHELL_METHOD(defaultTimeAggregateMakerWrapper);
        FIELD(timeAggregateMakerType, "Name of time agg maker type");
        FIELD(sumPVMaker,  "sumPVMaker");
        FIELD_MAKE_OPTIONAL(sumPVMaker);
        FIELD(realMaker, "realMaker");
        FIELD_MAKE_TRANSIENT(realMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    TimeAggregateMakerWrapper(): CObject(TYPE){}

    static IObject* defaultTimeAggregateMakerWrapper(){
        return new TimeAggregateMakerWrapper();
    }
};

typedef smartPtr<TimeAggregateMakerWrapper> TimeAggregateMakerWrapperSP;

CClassConstSP const TimeAggregateMakerWrapper::TYPE =
CClass::registerClassLoadMethod("TimeAggregateMakerWrapper", 
                                typeid(TimeAggregateMakerWrapper), load);

DRLIB_END_NAMESPACE

    

