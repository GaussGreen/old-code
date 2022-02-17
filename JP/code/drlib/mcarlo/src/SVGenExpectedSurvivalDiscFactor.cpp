//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedSurvivalDiscFactor.cpp
//
//   Description : A Generator of MC EXPECTED Survival Discount Factor State Variables
//
//   Author      : Eva X Strasser
//
//   Date        : 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"
//#include "edginc/MCPathConfigSRMGenSV.hpp"

#include <iostream>
using namespace std;

DRLIB_BEGIN_NAMESPACE

/** Constructor - for computing expected survival discount values on calcDate between
        pvDate and each date in dates. Note that currently past values for this
        are not supported. In particular a value of 0 is returned when calcDate
        is in the past (even if computeLog is true) */
SVGenExpectedSurvivalDiscFactor::SVGenExpectedSurvivalDiscFactor(
    const DateTime&         calcDate,       // when to compute
    const DateTime&         pvDate,         // when to discount to
    ICDSParSpreadsConstSP   cdsParSpreadCurve,
    const DateTimeArray&    dates,
    bool                    computeLog):    // true: do log
    calcDate(calcDate), pvDate(pvDate), dates(dates),
    cdsParSpreadCurve(cdsParSpreadCurve), computeLog(computeLog)
{
    
    static const string method("SVGenExpectedSurvivalDiscFactor::validate");
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
    if (0) {
        cerr << "SVGenExpectedSurvivalDiscFactor: " << pvDate.toString();
        for(int i=0; i < dates.size(); ++i)
            cerr << " " << dates[i].toString();
        cerr << endl;
    }
}


/** Retrieve the CDS par spread curve associated with this SVGenExpectedSurvivalDiscFactor */
ICDSParSpreadsConstSP SVGenExpectedSurvivalDiscFactor::getCDSParSpreadCurve() const {
    return cdsParSpreadCurve;
}

/** Retrieve the dates for which a discount factor is required */
const DateTimeArray& SVGenExpectedSurvivalDiscFactor::getDates() const {
    return dates;
}
/** Returns the date to which we pv to */
const DateTime& SVGenExpectedSurvivalDiscFactor::getPVDate() const{
    return pvDate;
}
/** Returns the date on which the expected value should be computed */
const DateTime& SVGenExpectedSurvivalDiscFactor::getCalcDate() const {
    return calcDate;
}

/** is the log of the expected discount factor wanted */
bool SVGenExpectedSurvivalDiscFactor::logRequired() const {
    return computeLog;
}

/** name of underlying credit */
const string SVGenExpectedSurvivalDiscFactor::getName() const 
{
    return cdsParSpreadCurve->getName();
}

/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in. 
    The return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenExpectedSurvivalDiscFactor::create(
    IStateVariableSP              oldStateVar,
    IStateVariableGen::IStateGen*  pathGen) const
{
    return getSVExpectedSurvivalDiscFactor(pathGen);
}

/** Returns a MC Expected Survival Discount Factor state variable which then
    provides access to the values etc. This is the method that
    products should call to get an SVGenExpectedSurvivalDiscFactor::IStateVar. */

SVExpSurvDiscFactorSP SVGenExpectedSurvivalDiscFactor::getSVExpectedSurvivalDiscFactor(
    IStateVariableGen::IStateGen* pathGen) const
{
    SVExpSurvDiscFactorSP survivalDiscFactorSV(
        &dynamic_cast<SVExpSurvDiscFactor&>(
            *pathGen->create(this)));
    return survivalDiscFactorSV;
}

class SVGenExpectedSurvivalDiscFactor::DeterminsticSV: public SVQmcExpSurvDiscFactor{
public:
    /** All elements are inside array */
    virtual double element(int idx) const {
        return (*this)[idx];
    }

    /** Returns true if this state variable is being used to 'simulate' the past. 
        This is a useful method for users of state variables - need to see how hard it is to implement */
    virtual bool doingPast() const{
        return doingThePast;
    }
    /** Once the calcDate is in the past we report a value of zero for each of the discount factors 
        (currently we're not supporting historic expected values) */
    DeterminsticSV(const SVGenExpectedSurvivalDiscFactor*  expSurvivalDiscFactor,
                   const DateTime&                      today,
                   bool                                 doingThePast):
        SVQmcExpSurvDiscFactor(NULL,
                                expSurvivalDiscFactor->calcDate,
                                expSurvivalDiscFactor->getDates(),
                                expSurvivalDiscFactor->logRequired()),
        doingThePast(doingThePast)
        {
        bool calcDateInPast = today.isGreaterOrEqual(expSurvivalDiscFactor->calcDate);
        int start, end;
        int numDates = expSurvivalDiscFactor->dates.size();
        if (doingThePast){
            start = 0;
            end = calcDateInPast? numDates: 0;
            survivalDiscountFactors.resize(end);
        } else {
            bool computeLog = expSurvivalDiscFactor->computeLog;
            start = calcDateInPast? numDates: 0;
            end = numDates;
            survivalDiscountFactors.resize(end);
            if (!calcDateInPast){
                DefaultRatesSP defaultRates = expSurvivalDiscFactor->cdsParSpreadCurve->defaultRates();
                for (int i = start; i < end; i++){
                    survivalDiscountFactors[i] = 
                        defaultRates->calcDefaultPV(expSurvivalDiscFactor->pvDate, expSurvivalDiscFactor->dates[i]);
                    if (computeLog){
                        survivalDiscountFactors[i] = log(survivalDiscountFactors[i]);
                    }
                }
            }
        }
        
        // Create path and first discount factor
        if (survivalDiscountFactors.empty())
            survivalDiscountFactors.push_back(0.0); // so firstSDF() returns 0
        
        SVPath::initialize(&survivalDiscountFactors[0], start, end);
    }
private:
    vector<double> survivalDiscountFactors;
    bool           doingThePast;
};

/** For use by Path Generators that want to use determinstic rates. */
SVExpSurvDiscFactor* SVGenExpectedSurvivalDiscFactor::determinsticSV(
    const DateTime& today,
    bool            doingPast) const{
    return new DeterminsticSV(this, today, doingPast);
}
/** implementing 'visitor' model */
void SVGenExpectedSurvivalDiscFactor::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE

