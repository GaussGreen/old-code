//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedSurvivalDiscFactor.cpp
//
//   Description : A Generator of MC EXPECTED Energy Futures prices State Variables
//
//
//   Date        : 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE


SVGenExpectedEnergyFuture::SVGenExpectedEnergyFuture(
    const DateTime&				_calcDate,		// when to compute
    EnergyFuturesCurveConstSP   _futureCurve,
    const DateTimeArray&		_dates,
    bool						_computeLog)		// true: do log
	:
    calcDate(_calcDate), /*pvDate(pvDate),*/ dates(_dates), 
    futureCurve(_futureCurve), computeLog(_computeLog)
{    
    static const string method("SVGenExpectedEnergyFuture::SVGenExpectedEnergyFuture");
    ASSERT(!dates.empty());

    /*if (calcDate.isGreater(pvDate)){
        throw ModelException(method, "calcDate "+
                             calcDate.toString()+" is after pvDate "+
                             pvDate.toString());
    }
    if (!dates.empty() && pvDate.isGreater(dates.front())){
        throw ModelException(method, "pvDate "+pvDate.toString()+
                             " is after first date to discount from");
    }*/
	// need to check if calcDate.isGreater(dates.front()) ??  --> measurement date is after 1st future maturity date
}


/** Retrieve the CDS par spread curve associated with this SVGenExpectedSurvivalDiscFactor */
EnergyFuturesCurveConstSP SVGenExpectedEnergyFuture::getEnergyFutureCurve() const 
{
    return futureCurve;
}

/** Retrieve the dates for which a discount factor is required */
const DateTimeArray& SVGenExpectedEnergyFuture::getDates() const {
    return dates;
}

/** Returns the date on which the expected value should be computed */
const DateTime& SVGenExpectedEnergyFuture::getCalcDate() const {
    return calcDate;
}

/** is the log of the expected discount factor wanted */
bool SVGenExpectedEnergyFuture::logRequired() const {
    return computeLog;
}

/** name of underlying credit */
const string SVGenExpectedEnergyFuture::getName() const 
{
    return futureCurve->getName();
}

/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in. 
    The return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenExpectedEnergyFuture::create(
    IStateVariableSP oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const
{
    return getSVExpEnergyFuture(pathGen);
}

/** Returns a MC Expected Survival Discount Factor state variable which then
    provides access to the values etc. This is the method that
    products should call to get an SVGenExpectedSurvivalDiscFactor::IStateVar. */
SVExpEnergyFutureSP SVGenExpectedEnergyFuture::getSVExpEnergyFuture(
    IStateVariableGen::IStateGen* pathGen) const
{
    SVExpEnergyFutureSP expEnergyFutureSV( &dynamic_cast<SVExpEnergyFuture &>(*pathGen->create(this)) );
    return expEnergyFutureSV;
}

class SVGenExpectedEnergyFuture::PastSV: public SVQmcExpEnergyFuture 
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

    /** Once the calcDate is in the past we report a value of zero for
    each of the future prices (currently we're not supporting
    historic expected values) */
    PastSV(
        const SVGenExpectedEnergyFuture* expFutureGen,
        const DateTime& today,
        bool doingThePast)
        :
        SVQmcExpEnergyFuture(
            NULL, 
            expFutureGen->calcDate, 
            expFutureGen->dates,
            expFutureGen->computeLog ),
        firstFuturePrice(0.0), 
        doingThePast(doingThePast)
    {
        bool calcDateInPast = today.isGreaterOrEqual(expFutureGen->calcDate);
        int start, end;
        int numDates = expFutureGen->dates.size();
        if (doingThePast){
            start = 0;
            end = calcDateInPast? numDates: 0;
            futurePrices.resize(end);
        } else {
            bool computeLog = expFutureGen->computeLog;
            start = calcDateInPast? numDates: 0;
            end = numDates;
            futurePrices.resize(end);
            if (!calcDateInPast){
                for (int i = start; i < end; i++){
                    futurePrices[i] = expFutureGen->futureCurve->getFwdRate(i); // FIX THIS!?
                    if (computeLog){
                        futurePrices[i] = log(futurePrices[i]);
                    }
                }
            }
        }
        // Create path and first future price
        // in case there is no DFs we add one element, so that firstDF() works as expected (returns 0).

        if (futurePrices.empty())
            futurePrices.push_back(firstFuturePrice); // note that in this case we still keep start=end=0
        firstFuturePrice = futurePrices.front();

        SVPath::initialize(&futurePrices[0],start, end); // initialize path with direct access constructor
        prepared = true; //FIXME: debug only
    }

private:
    double         firstFuturePrice;
    vector<double> futurePrices;
    bool           doingThePast;
};

/** For use by Path Generators that want to use determinstic future prices. */
SVExpEnergyFuture * SVGenExpectedEnergyFuture::getPastSV(
    const DateTime& today,
    bool doingPast) const
{
    return new PastSV(this, today, doingPast);
}


/** implementing 'visitor' model */
void SVGenExpectedEnergyFuture::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
	sv->processSVGen(this);
}


DRLIB_END_NAMESPACE

