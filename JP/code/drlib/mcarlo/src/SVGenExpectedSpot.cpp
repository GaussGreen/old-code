//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedSpot.cpp
//
//   Description : A Generator of MC EXPECTED Asset Spot State Variables
//                 aka Forward Prices
//
//   Date        : 19 Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/SVGenExpectedSpot.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor - for computing on calcDate the expected spot values on
    dates */
SVGenExpectedSpot::SVGenExpectedSpot(
    int                   assetIdx, // into MultiMarketFactors
    const DateTime&       calcDate, // When to compute the
    const DateTimeArray&  dates): // expected spot on these dates
    assetIdx(assetIdx), calcDate(calcDate), dates(dates)
{
    ASSERT(!dates.empty());
}

/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenExpectedSpot::create(IStateVariableSP              oldStateVar,
                                     IStateVariableGen::IStateGen*  pathGen) const{
    return getExpSpotSV(pathGen);
}

/** Returns a MC Expected Discount Factor state variable which then
    provides access to the values etc. This is the method that
    products should call to get an SVGenExpectedSpot::IStateVar. */
SVPathSP SVGenExpectedSpot::getExpSpotSV(
    IStateVariableGen::IStateGen* pathGen) const{
    IStateVariableSP  sv = pathGen->create(this);
    SVPathSP expSpotSV(
        &dynamic_cast<SVPath&>(*sv));
    return expSpotSV;
}
SVExpectedFXSP SVGenExpectedSpot::getExpFXSV(
    IStateVariableGen::IStateGen* pathGen) const{
    IStateVariableSP  sv = pathGen->create(this);
    SVExpectedFXSP expSpotSV(
        &dynamic_cast<SVExpectedFX&>(*sv));
    return expSpotSV;
}
SVExpectedEQSP SVGenExpectedSpot::getExpEQSV(
    IStateVariableGen::IStateGen* pathGen) const{
    IStateVariableSP  sv = pathGen->create(this);
    SVExpectedEQSP expSpotSV(
        &dynamic_cast<SVExpectedEQ&>(*sv));
    return expSpotSV;
}

/** Returns the date on which the expected value should be computed */
const DateTime&   SVGenExpectedSpot::getCalcDate() const{
    return calcDate;
}

/** Returns the dates for when the expected spot should be computed */
const DateTimeArray& SVGenExpectedSpot::getDates() const{
    return dates;
}

/** Returns the number of dates */
int SVGenExpectedSpot::numDates() const{
    return dates.size();
}
 
/** Returns the index of the associated asset (into MultiMarketFactors) */
int SVGenExpectedSpot::getAssetIdx() const{
    return assetIdx;
}
   
class SVGenExpectedSpot::ZeroSV: public SVQmcExpectedFX
//virtual IStateVar
{
public:
    virtual ~ZeroSV(){}

    virtual bool doingPast() const{
        return true;
    }
    
    /** Gives a zero path for historic dates */
    virtual const SVPath& path() const{
        return thePath;
    }

    /** Returns the expected spot value at date idx given 
    the current state of the world, ie. current time point 
    on current path.  Instruments set up the state variable
    so they should know the mapping between idx and actual dates
    (also, getDates can be used to look these dates back up 
    if needed). */
    virtual double getFwd( int idx ) const {
        return thePath[idx];
    }

    ZeroSV(int start, int end): SVQmcExpectedFX(NULL, DateTime(), DateTimeArray(), false),
                                zeroValues(end),
                                thePath((end > 0) ? &zeroValues[0] : NULL, start, end)
                                {}
 
private:
    vector<double> zeroValues;
    SVPath   thePath;
};
    
/** Create an IStateVar which has zero for all dates when the calcDate is in the
    past */
SVExpectedFXSP SVGenExpectedSpot::createHistoricZeroSV(
    const DateTime& today) const{
    bool calcDateInPast = today.isGreaterOrEqual(calcDate);
    int end = calcDateInPast? dates.size(): 0;
    return SVExpectedFXSP(new ZeroSV(0, end));
}

/** implementing 'visitor' model */
void SVGenExpectedSpot::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

// temporary placement for SVGenSpot... 
// TODO: sort out the directory dependencies, move SVGenSpot files to mcarlo

DRLIB_END_NAMESPACE


