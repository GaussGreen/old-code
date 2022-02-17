//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedBasisForward.cpp
//
//   Description : A Generator of MC EXPECTED Basis Forward Spread State Variables
//
//   Author      : Anatoly V Morosov
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

#include <iostream>
using namespace std;

DRLIB_BEGIN_NAMESPACE

/** Constructor - for computing expected survival discount values on calcDate between
        pvDate and each date in dates. Note that currently past values for this
        are not supported. In particular a value of 0 is returned when calcDate
        is in the past (even if computeLog is true) */
SVGenExpectedBasisFwdSpread::SVGenExpectedBasisFwdSpread(
    IBasisIndexCurveConstSP   basisCurve,
    const DateTime&         calcDate,       // when to compute
    const DateTimeArray&    resetDates):    
    calcDate(calcDate), dates(resetDates),
    basisCurve(basisCurve) {
    
    static const string method("SVGenExpectedBasisFwdSpread::SVGenExpectedBasisFwdSpread");
    ASSERT(!resetDates.empty());

    if (0) {
        cerr << "SVGenExpectedBasisFwdSpread: " << calcDate.toString();
        for(int i=0; i < dates.size(); ++i)
            cerr << " " << dates[i].toString();
        cerr << endl;
    }
}


/** Retrieve the CDS par spread curve associated with this SVGenExpectedSurvivalDiscFactor */
IBasisIndexCurveConstSP SVGenExpectedBasisFwdSpread::getBasisCurve() const {
    return basisCurve;
}

/** Retrieve the dates for which a discount factor is required */
const DateTimeArray& SVGenExpectedBasisFwdSpread::getDates() const {
    return dates;
}
/** Returns the date to which we pv to */
const DateTime& SVGenExpectedBasisFwdSpread::getPVDate() const{
    return calcDate;
}
/** Returns the date on which the expected value should be computed */
const DateTime& SVGenExpectedBasisFwdSpread::getCalcDate() const {
    return calcDate;
}

/** name of underlying credit */
const string SVGenExpectedBasisFwdSpread::getName() const 
{
    return basisCurve->getName();
}

/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in. 
    The return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenExpectedBasisFwdSpread::create(
    IStateVariableSP              oldStateVar,
    IStateVariableGen::IStateGen*  pathGen) const
{
    return getSVExpectedBasisFwdSpread(pathGen);
}

SVExpectedBasisFwdSpreadSP SVGenExpectedBasisFwdSpread::getSVExpectedBasisFwdSpread(
    IStateVariableGen::IStateGen* pathGen) const
{
    SVExpectedBasisFwdSpreadSP basisSpreadSV(
        &dynamic_cast<SVExpectedBasisFwdSpread&>(
            *pathGen->create(this)));
    return basisSpreadSV;
}

class SVGenExpectedBasisFwdSpread::DeterminsticSV: public SVQmcExpectedBasisFwdSpread {
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
    DeterminsticSV(const SVGenExpectedBasisFwdSpread*  expBasisSpread,
                   const DateTime&                      today,
                   bool                                 doingThePast):
        SVQmcExpectedBasisFwdSpread(NULL,
                                    expBasisSpread->getCalcDate(),
                                    expBasisSpread->getDates()),
        doingThePast(doingThePast)
    {
        // not implemented, populates with 0
        bool calcDateInPast = today.isGreaterOrEqual(expBasisSpread->getCalcDate());
        int numDates = expBasisSpread->getDates().size();
        int end = calcDateInPast? numDates: 0;
        storedBasisSpreads.resize(end);

        // Create path and first discount factor
        if (storedBasisSpreads.empty())
            storedBasisSpreads.push_back(0.0); // so first() returns 0
        
        SVPath::initialize(&storedBasisSpreads[0], 0, end);
    }
private:
    vector<double> storedBasisSpreads;
    bool           doingThePast;
};

/** For use by Path Generators that want to use determinstic rates. */
SVExpectedBasisFwdSpread* SVGenExpectedBasisFwdSpread::determinsticSV(
    const DateTime& today,
    bool            doingPast) const{
    return new DeterminsticSV(this, today, doingPast);
}
/** implementing 'visitor' model */
void SVGenExpectedBasisFwdSpread::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE

