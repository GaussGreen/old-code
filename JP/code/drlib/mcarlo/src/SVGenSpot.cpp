//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenSpot.cpp
//
//   Description : A Generator of MC Spot State Variables
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_MCSPOT_CPP
#include "edginc/SVGenSpot.hpp"
//#include "edginc/MCPathConfigSRMGenSV.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE


MCPath::IStateVar::IStateVar(){}
MCPath::IStateVar::~IStateVar(){}



int MCPath::PathByAsset::numAssets() const{
    return paths.size();
}
        
bool MCPath::PathByAsset::doingPast() const{
    return isDoingPast;
}

const SVPath& MCPath::PathByAsset::path(int iAsset) const { 
    return (*(paths[iAsset])); 
}  

SVPathSP MCPath::PathByAsset::getSV(int iAsset)  { 
    return paths[iAsset]; 
} 

double MCPath::PathByAsset::maxDriftEstimate(int iAsset) const{
    return maxDrifts[iAsset];
}

double MCPath::PathByAsset::getSpotPrice(int iAsset) const{
    return paths[iAsset]->getValue();
}

void MCPath::PathByAsset::getSpotPrices(DoubleArray& prices) const{
    int nbAssets = numAssets();
    for (int i = 0; i < nbAssets; ++i) {
        SVPath* cp = paths[i].get();
        if (cp)
            prices[i] = getSpotPrice(i);
        else
            prices[i] = 0;
        // in some cases, e.g. CcyProt.xml cp can be 0
        // i am not sure why
    }
}

void MCPath::PathByAsset::advance() {
    int nbAssets = numAssets();
    for (int i = 0; i < nbAssets; ++i) {
        SVPath* cp = paths[i].get();
        if (cp)
            cp->advance();
        // in some cases, e.g. CcyProt.xml cp can be 0
        // i am not sure why
    }
}

void MCPath::PathByAsset::reset() {
    int nbAssets = numAssets();
    for (int i = 0; i < nbAssets; ++i) {
        SVPath* cp = paths[i].get(); 
        if (cp)
            cp->reset(); 
        // in some cases, e.g. CcyProt.xml cp can be 0
        // i am not sure why
    }
}

MCPath::PathByAsset::PathByAsset(int numAssets,
                                 const vector<double>& maxDrifts, 
                                 bool isDoingPast):
paths(numAssets), maxDrifts(maxDrifts), isDoingPast(isDoingPast) {}

MCPath::PathByAsset::PathByAsset(const vector<SVPathSP>& paths, 
                                 const vector<double>& maxDrifts,
                                 bool                  isDoingPast): 
    paths(paths), maxDrifts(maxDrifts), isDoingPast(isDoingPast){}


MCPath::~MCPath(){}

/** Returns the sim series for this MCPath */
SimSeriesConstSP MCPath::getSimSeries() const{
    return simSeries;
}

/** Returns the number of assets */
int MCPath::numAssets() const{
    return simSeries->getNumAssets();
}

/** Returns the number of dates (including all assets) */
int MCPath::numDates() const{
    return simSeries->getAllDates().size();
}

/** Returns the number of dates for the specified asset */
int MCPath::numDates(int iAsset) const{
    return simSeries->numDates(iAsset);
}

/** Constructor - takes reference to supplied simSeries */
MCPath::MCPath(SimSeriesConstSP simSeries): simSeries(simSeries){}

/** Constructor - from numAssets and an array of dates (same dates per
    asset) */
MCPath::MCPath(int numAssets, const DateTimeArray& dates){
    SimSeriesSP theSimSeries(new SimSeries(numAssets));
    theSimSeries->addDates(dates);
    simSeries = theSimSeries;
}

/** Constructor - from numAssets and a single date */
MCPath::MCPath(int numAssets, const DateTime& date){
    SimSeriesSP theSimSeries(new SimSeries(numAssets));
    theSimSeries->addDates(DateTimeArray(1, date));
    simSeries = theSimSeries;
}


//////////////////////////////////////////////////////////////////////////


/** Constructor - takes reference to supplied simSeries */
SVGenSpot::SVGenSpot(SimSeriesConstSP simSeries): 
MCPath(simSeries) {}

/** Constructor - from numAssets and an array of dates (same dates per
    asset) */
SVGenSpot::SVGenSpot(int numAssets, const DateTimeArray& dates): 
MCPath(numAssets, dates) {}

/** Constructor - from numAssets and a single date */
SVGenSpot::SVGenSpot(int numAssets, const DateTime& date):
MCPath(numAssets, date) {}

/** Create the corresponding State Variable for this State Variable
    Generator (from IStateVariableGen interface). */
IStateVariableSP SVGenSpot::create(IStateVariableSP             oldStateVar,
                             IStateVariableGen::IStateGen* pathGen) const{
    return getSpotSV(pathGen);
}

/** Returns a MC Spot state variable which then provides access to the
    path etc. This is the method that products should call to get an
    MCPath::IStateVar. */
MCPath::IStateVarSP SVGenSpot::getSpotSV(IStateVariableGen::IStateGen* pathGen) const{
    IStateVarSP spotPathSV(&dynamic_cast<IStateVar&>(*pathGen->create(this)));
    return spotPathSV;
}

/** implementing 'visitor' pattern */
void SVGenSpot::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

/** Constructor - takes reference to supplied simSeries */
MCQuadVar::MCQuadVar(SimSeriesConstSP simSeries): 
MCPath(simSeries) {}

/** Constructor - from numAssets and an array of dates (same dates per
    asset) */
MCQuadVar::MCQuadVar(int numAssets, const DateTimeArray& dates): 
MCPath(numAssets, dates) {}

/** Constructor - from numAssets and a single date */
MCQuadVar::MCQuadVar(int numAssets, const DateTime& date):
MCPath(numAssets, date) {}

/** Create the corresponding State Variable for this State Variable
    Generator (from IStateVariableGen interface). */
IStateVariableSP MCQuadVar::create(IStateVariableSP             oldStateVar,
                                IStateVariableGen::IStateGen* pathGen) const{
    return getQuadVarSV(pathGen);
}

/** Returns a MC Spot state variable which then provides access to the
    path etc. This is the method that products should call to get an
    MCPath::IStateVar. */
MCPath::IStateVarSP MCQuadVar::getQuadVarSV(IStateVariableGen::IStateGen* pathGen) const{
    IStateVarSP quadVarPathSV(&dynamic_cast<IStateVar&>(*pathGen->create(this)));
    return quadVarPathSV;
}

void MCQuadVar::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

/** Constructor - takes reference to supplied simSeries */
MCSqrtAnnualQuadVar::MCSqrtAnnualQuadVar(SimSeriesConstSP simSeries): 
MCPath(simSeries) {}

/** Constructor - from numAssets and an array of dates (same dates per
    asset) */
MCSqrtAnnualQuadVar::MCSqrtAnnualQuadVar(int numAssets, const DateTimeArray& dates): 
MCPath(numAssets, dates) {}

/** Constructor - from numAssets and a single date */
MCSqrtAnnualQuadVar::MCSqrtAnnualQuadVar(int numAssets, const DateTime& date):
MCPath(numAssets, date) {}

/** Create the corresponding State Variable for this State Variable
    Generator (from IStateVariableGen interface). */
IStateVariableSP MCSqrtAnnualQuadVar::create(IStateVariableSP             oldStateVar,
                                          IStateVariableGen::IStateGen* pathGen) const{
    return getSqrtAnnualQuadVarSV(pathGen);
}

/** Returns a MC Spot state variable which then provides access to the
    path etc. This is the method that products should call to get an
    MCPath::IStateVar. */
MCPath::IStateVarSP MCSqrtAnnualQuadVar::getSqrtAnnualQuadVarSV(IStateVariableGen::IStateGen* pathGen) const{
    IStateVarSP sqrtAnnualQuadVarPathSV(&dynamic_cast<IStateVar&>(*pathGen->create(this)));
    return sqrtAnnualQuadVarPathSV;
}

void MCSqrtAnnualQuadVar::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE
