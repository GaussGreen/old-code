//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenAggregatedSurvivalDiscFactor.cpp
//
//   Description : Generator that captures all dates from a set of SDFs and ESDFs
//
//   Author      : Vladimir Grebinskiy
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE

SVGenAggregatedSurvivalDiscFactor::SVGenAggregatedSurvivalDiscFactor(
            DateTimeArraySP   _sdfDates,            // set of dates for discount factors
            SpotIdxArraySP    _sdfIdxSP,
            DateTimeArraySP   _esdfRequestedDates,  //  union of {t_i} for all expected factors
            FwdIdxArraySP     _esdfReqIdxSP,
            DateTimeArraySP   _esdfForwardDates,    // union of all {T_j} for all expected factors
            FwdIdxArraySP     _esdfForIdxSP,
            const DateTime&        _maxDiffDate,        // asset specific max diffusion
            const DateTime&        _maxCurveDate,       // asset specific max curve
            ICDSParSpreadsConstSP  _cdsParSpreadCurve,
            bool                   _computeLog)    :
            sdfDates(_sdfDates),
            esdfRequestedDates(_esdfRequestedDates),
            esdfForwardDates(_esdfForwardDates),
            maxMaturityDate(_maxDiffDate),
            maxCurveDate(_maxCurveDate),
            cdsParSpreadCurve(_cdsParSpreadCurve),
            computeLog(_computeLog),
            sdfIdxSP(_sdfIdxSP),
            esdfReqIdxSP(_esdfReqIdxSP),
            esdfForIdxSP(_esdfForIdxSP)
{
    // Check that requested dates are not fully in the future
    check(sdfDates, maxMaturityDate);
    check(esdfRequestedDates, maxMaturityDate);
    check(esdfForwardDates, maxCurveDate);
}
            
void SVGenAggregatedSurvivalDiscFactor::check(DateTimeArraySP dates, const DateTime& bound)
{
    const string method = "SVGenAggregatedSurvivalDiscFactor::check(DateTimeArraySP dates, const DateTime& bound)";
    if (!dates && !(*dates).empty() && (*dates)[0] > bound)
        throw ModelException(method, "All dates are beyond max date! Dates[0]= " + (*dates)[0].toString() + " maxDate= " + bound.toString());
}
        
/** Returns a Expected Survival Discount Factor state variable which then
provides access to the values etc. This is the method that
products should call to get an SVExpSurvDiscFactor. */

SVAggregatedSurvDiscFactorSP SVGenAggregatedSurvivalDiscFactor::getIQSVGenAggregatedSurvivalDiscFactorSV(
                            IStateVariableGen::IStateGen* pathGen) const
{
    SVAggregatedSurvDiscFactorSP aggregatedDiscFactorSV(
            &dynamic_cast<SVAggregatedSurvDiscFactor&>(
            *pathGen->create(this)));
    return aggregatedDiscFactorSV;
}
IStateVariableSP SVGenAggregatedSurvivalDiscFactor::create(
        IStateVariableSP /*oldStateVar*/,
        IStateGen*     stateGen) const
{
    return getIQSVGenAggregatedSurvivalDiscFactorSV(stateGen);
}

/** Retrieve the CDS par spread curve associated with this SVGenAggregatedSurvivalDiscFactor */
ICDSParSpreadsConstSP SVGenAggregatedSurvivalDiscFactor::getCDSParSpreadCurve() const
{
    return cdsParSpreadCurve;
}

/** Retrieve the dates for which a discount factor is required */
const DateTimeArray& SVGenAggregatedSurvivalDiscFactor::getSDFDates() const
{
    return *sdfDates;
}
/** Retrieve measurement the dates for which a expected discount factor is required */
const DateTimeArray& SVGenAggregatedSurvivalDiscFactor::getESDFRequestedDates() const
{
    return *esdfRequestedDates;
}

/** Retrieve forward dates for which a expected discount factor is required */
const DateTimeArray& SVGenAggregatedSurvivalDiscFactor::getESDFForwardDates() const
{
    return *esdfForwardDates;
}

/** is the log of the expected discount factor wanted */
bool SVGenAggregatedSurvivalDiscFactor::logRequired() const
{
    return computeLog;
}

/** name of underlying credit */
const string SVGenAggregatedSurvivalDiscFactor::getName() const
{
    return cdsParSpreadCurve->getName();
}

void SVGenAggregatedSurvivalDiscFactor::attachSVGen(IElemStateVariableGenVisitor* sv) const    // support visitor pattern
{
    sv->processSVGen(this);
}

// Return last date which is <= "m". If no such date exist, return the smallest possible one
static DateTime lastNonGreater(const DateTime & m, const DateTimeArray& dates)
{
    if (dates.empty() ||  dates.front() > m)
        return DateTime(); // everything is bigger than m -- return the smallest dates
    else
        return dates[m.findLower(dates)];
}

// if the dates are d1< d2< d3 <= maxMaturityDate < d4 < ...
//    return d3
DateTime SVGenAggregatedSurvivalDiscFactor::getMaxMatDate() const
{
    DateTime date= max(lastNonGreater(maxMaturityDate, *sdfDates), lastNonGreater(maxMaturityDate, *esdfRequestedDates));
    return date;
}

// return the last forward date which is not in the future vs maxCurveDate
DateTime SVGenAggregatedSurvivalDiscFactor::getMaxCurveDate() const
{
    DateTime date= lastNonGreater(maxCurveDate, *esdfForwardDates);
    return date;
}

DRLIB_END_NAMESPACE
