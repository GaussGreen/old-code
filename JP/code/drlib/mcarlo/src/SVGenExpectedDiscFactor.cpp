//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedDiscFactor.cpp
//
//   Description : A Generator of MC Discount Factor State Variables
//
//   Date        : 6 Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"
#include "edginc/SVQmcImplemented.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor - for computing expected discount values on calcDate between
    pvDate and each date in dates */
SVGenExpectedDiscFactor::SVGenExpectedDiscFactor(
    const DateTime&      calcDate, // when to compute
    const DateTime&      pvDate, // when to discount to
    YieldCurveConstSP    yieldCurve,
    const DateTimeArray& dates,
    bool                 computeLog):
    VirtualDestructorBase(), //refCount(0),
    calcDate(calcDate), pvDate(pvDate),
    dates(dates), yieldCurve(yieldCurve), computeLog(computeLog)
{
    static const string method("SVGenExpectedDiscFactor::validate");
    ASSERT(!dates.empty());
    if (calcDate.isGreater(pvDate)){
        throw ModelException(method, "calcDate "+
                             calcDate.toString()+" is after pvDate "+
                             pvDate.toString());
    }
    if (!dates.empty() && pvDate.isGreater(dates.front())){
        throw ModelException(method, "pvDate "+pvDate.toString()+
                             " is after first date to discount from");
    }
}

/** Retrieve the yield curve associated with this SVGenExpectedDiscFactor */
YieldCurveConstSP SVGenExpectedDiscFactor::getYieldCurve() const{
    return yieldCurve;
}
/** Retrieve the dates for which a discount factor is required */
const DateTimeArray& SVGenExpectedDiscFactor::getDates() const{
    return dates;
}
/** Returns the date on which to compute expected value */
const DateTime& SVGenExpectedDiscFactor::getCalcDate() const{
    return calcDate;
}
/** Returns the date to which we pv to */
const DateTime& SVGenExpectedDiscFactor::getPVDate() const{
    return pvDate;
}

/** is the log of the expected discount factor wanted */
bool SVGenExpectedDiscFactor::logRequired() const{
    return computeLog;
}

/** Returns a MC Discount Factor state variable which then
    provides access to the path etc. This is the method that
    products should call to get an SVExpectedDiscFactorSP. */

 SVExpectedDiscFactorSP SVGenExpectedDiscFactor::getSVExpectedDiscFactor(
    IStateVariableGen::IStateGen* pathGen) const
{
     SVExpectedDiscFactorSP discFactorSV(
        & dynamic_cast<SVExpectedDiscFactor&>(* pathGen->create(this)));
    return discFactorSV;
}

class SVGenExpectedDiscFactor::DeterminsticSV: public SVQmcExpectedDiscFactor
{
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
    /** asset specific name of one of the generic functions */
    virtual double getFwdDF( int idx ) const {
        return element(idx);
    }

    /** asset specific name of one of the generic functions */
    virtual double firstDF() const {
        return element(0);
    }

    /** Once the calcDate is in the past we report a value of zero for
        each of the discount factors (currently we're not supporting
        historic expected values) */
    DeterminsticSV(const SVGenExpectedDiscFactor* expDiscFactor,
                   const DateTime&             today,
                   bool                        doingThePast):
        SVQmcExpectedDiscFactor(NULL, today, expDiscFactor->dates,
                                expDiscFactor->computeLog,
                                expDiscFactor->getYieldCurve()),
        firstDiscFactor(0.0), doingThePast(doingThePast){
        /*
          Formerly 
          "bool calcDateInPast = today.isGreaterOrEqual(expDiscFactor->calcDate);"
          but relaxed somewhat. Seems to me the condition is better expressed as
          whether the "pvFrom" date is in the past or not, but this will do me for now.
        */
        bool calcDateInPast = today.isGreater(expDiscFactor->calcDate); 
        int start, end;
        int numDates = expDiscFactor->dates.size();
        if (doingThePast){
            start = 0;
            end = calcDateInPast? numDates: 0;
            discountFactors.resize(end);
        } else {
            bool computeLog = expDiscFactor->computeLog;
            start = calcDateInPast? numDates: 0;
            end = numDates;
            discountFactors.resize(end);
            if (!calcDateInPast){
                for (int i = start; i < end; i++){
                    discountFactors[i] = expDiscFactor->yieldCurve->
                        pv(expDiscFactor->pvDate, expDiscFactor->dates[i]);
                    if (computeLog){
                        discountFactors[i] = log(discountFactors[i]);
                    }
                }
            }
        }
        // Create path and first discount factor
        // in case there is no DFs we add one element, so that firstDF() works as expected (returns 0).

        if (discountFactors.empty())
            discountFactors.push_back(firstDiscFactor); // note that in this case we still keep start=end=0
        firstDiscFactor = discountFactors.front();

        SVPath::initialize(&discountFactors[0],start, end); // initialize path with direct access constructor

    }
private:
    double         firstDiscFactor;
    vector<double> discountFactors;
    bool           doingThePast;
};

/** For use by Path Generators that want to use determinstic rates. */
SVExpectedDiscFactor * SVGenExpectedDiscFactor::determinsticSV(
    const DateTime& today,
    bool            doingPast) const{
    return new DeterminsticSV(this, today, doingPast);
}

/** implementing 'visitor' model */
void SVGenExpectedDiscFactor::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE
