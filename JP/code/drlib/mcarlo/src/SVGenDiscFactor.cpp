//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenDiscFactor.cpp
//
//   Description : A Generator of MC Discount Factor State Variables
//
//   Date        : 20 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor - from  an array of dates */
SVGenDiscFactor::SVGenDiscFactor(const DateTime&      today, // when to discount to
                           YieldCurveConstSP    yieldCurve, 
                           const DateTimeArray& dates):
    today(today), isMargin(false), 
    yieldCurve(yieldCurve), originalDates(dates), adjustedDates(dates) {
    static const string routine = "SVGenDiscFactor::SVGenDiscFactor";
    ASSERT(!dates.empty());

    try {
        validate();
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Constructor - from a single date */
SVGenDiscFactor::SVGenDiscFactor(
    const DateTime&   today, // when to discount to
    YieldCurveConstSP yieldCurve,
    const DateTime&   date):
    today(today), isMargin(false), 
yieldCurve(yieldCurve), originalDates(1, date), adjustedDates(1, date) {
    static const string routine = "SVGenDiscFactor::SVGenDiscFactor";
    try {
        validate();
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** work out dates from settlement object */
void SVGenDiscFactor::initDatesFromInstSettlement(
    InstrumentSettlementSP instSettle,
    const DateTimeArray&   theDates,
    CAsset*                asset) {
    
    static const string routine = "SVGenDiscFactor::initDatesFromInstSettlement";
    
    try {
        isMargin = instSettle->isMargin();
        if (!isMargin){
            adjustedDates = DateTimeArray(theDates.size());
            for (int i = 0; i < adjustedDates.size(); i++){
                adjustedDates[i] = instSettle->settles(theDates[i], asset);
            }
        } else {
            adjustedDates = theDates;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Constructor - from  an array of dates and a settlement */
SVGenDiscFactor::SVGenDiscFactor(const DateTime&        today, // when to discount to
                           YieldCurveConstSP      yieldCurve, 
                           InstrumentSettlementSP instSettle,
                           const DateTimeArray&   theDates):
    today(today), yieldCurve(yieldCurve), originalDates(theDates) {
    static const string routine = "SVGenDiscFactor::SVGenDiscFactor";
    ASSERT(!theDates.empty());

    try {
        initDatesFromInstSettlement(instSettle, theDates, 0);
        validate();        
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Constructor - from a single date and a settlement. Does not support
    physical settlement */
SVGenDiscFactor::SVGenDiscFactor(const DateTime&        today, // when to discount to
                           YieldCurveConstSP      yieldCurve,
                           InstrumentSettlementSP instSettle,
                           const DateTime&        date):
    today(today), yieldCurve(yieldCurve), originalDates(1, date) {
    static const string routine = "SVGenDiscFactor::SVGenDiscFactor";
    try {
        initDatesFromInstSettlement(instSettle, originalDates, 0);
        validate();
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}
    
/** Constructor - from a single date and a settlement. Does support
    physical settlement */
SVGenDiscFactor::SVGenDiscFactor(const DateTime&        today, // when to discount to
                           YieldCurveConstSP      yieldCurve,
                           InstrumentSettlementSP instSettle,
                           const DateTime&        date,
                           CAssetSP               asset):
    today(today), yieldCurve(yieldCurve), originalDates(1, date) {
    static const string routine = "SVGenDiscFactor::SVGenDiscFactor";
    try {
        initDatesFromInstSettlement(instSettle, originalDates, asset.get());
        validate();
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Retrieve the yield curve associated with this SVGenDiscFactor */
YieldCurveConstSP SVGenDiscFactor::getYieldCurve() const{
    return yieldCurve;
}
/** Retrieve the dates for which a discount factor is required */
const DateTimeArray& SVGenDiscFactor::getDates() const{
    return adjustedDates;
}

/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenDiscFactor::create(IStateVariableSP             oldStateVar,
                                   IStateVariableGen::IStateGen* pathGen) const{
    return getSVDiscFactor(pathGen);
}


/** Returns a MC Discount Factor state variable which then
    provides access to the path etc. This is the method that
    products should call to get an SVQmcDiscFactorSP. */
SVDiscFactorSP SVGenDiscFactor::getSVDiscFactor(IStateVariableGen::IStateGen* pathGen) const
{
    SVDiscFactorSP discFactorSV(&dynamic_cast<SVDiscFactor&>(*pathGen->create(this)));
    return discFactorSV;
}

void SVGenDiscFactor::validate() {
    static const string routine = "SVGenDiscFactor::validate";

    try {
        DateTime::ensureIncreasing(originalDates, "Discount dates", 0);
        DateTime::ensureIncreasing(adjustedDates, "Settlement adjusted discount dates", 0);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


class SVGenDiscFactor::DeterminsticSV: public SVQmcDiscFactor {
public:

    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const{
        return doingThePast;
    }
    /** All elements are inside array */
    virtual double element(int idx) const {
        return (*this)[idx];
    }

    virtual double getDF(int i)const { return element(i);}

    DeterminsticSV(YieldCurveConstSP    yieldCurve, 
                   const DateTime&      today,
                   const DateTimeArray& originalDates,
                   const DateTimeArray& adjustedDates,
                   bool                 doingThePast,
                   bool                 isMargin):
    SVQmcDiscFactor(NULL, adjustedDates), // TODO: derive DeterminsticSV from a correct base class
    firstDiscFactor(0.0), 
    doingThePast(doingThePast)
    {
        static const string routine = "SVGenDiscFactor::DeterminsticSV::DeterminsticSV";
        
        try {
            // Figure out the start and end of loop
            int numOrigDatesInPast = today.numPastDates(originalDates);
            int end = doingThePast? numOrigDatesInPast: originalDates.size();
            // but for non-SRM models there is a need for discount factors after 
            // non-IR simulation dates, yet while some payment dates remain in the
            // future. This is cheap and while not strictly correct seems to cover
            // the case well enough.
            int endAlloc = originalDates.size();

            // Only compute discount factors for future adjustedDates
            discountFactors = vector<double>(endAlloc, 0.0);
            int numAdjDatesInPast = today.numPastDates(adjustedDates);
            for (int i = numAdjDatesInPast; i < endAlloc; i++){
                discountFactors[i] = isMargin ? 1.0 : yieldCurve->pv(adjustedDates[i]);
            }
        
            // Create path and first discount factor
            if (!discountFactors.empty()){
                firstDiscFactor = discountFactors.front();
            }
            //Path::initialize(....)
            SVPath::initialize(&discountFactors[0],
                             doingThePast? 0: numOrigDatesInPast,
                             end);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }
private:
    double         firstDiscFactor;
    vector<double> discountFactors;
    bool           doingThePast;
};

/** For use by Path Generators that want to use determinstic rates. */
SVDiscFactor* SVGenDiscFactor::determinsticSV(bool doingPast) const{
    return new DeterminsticSV(yieldCurve, today, originalDates, adjustedDates, doingPast, isMargin);
}

/** implementing 'visitor' model */
void SVGenDiscFactor::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE
