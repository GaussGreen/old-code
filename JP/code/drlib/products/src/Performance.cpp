//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Performance.cpp
//
//   Description : Measures performance of assets
//
//   Date        : 24 Oct 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Performance.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"

#ifdef STATE_VARIABLES
#endif

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////

#ifdef STATE_VARIABLES

bool PerformanceSV::StateVar::doingPast() const {
    return isPast;
}


int PerformanceSV::StateVar::numPerformances() const {
    return nbAssets;
}


void PerformanceSV::StateVar::calcPerf(DoubleArray& perf) const {
    for(int iAsset = 0; iAsset < nbAssets; iAsset++) {
        perf[iAsset] = calcPerf(iAsset);
    }
}
        

CVolRequestLN* PerformanceSV::StateVar::getVolInterp(const MCProductClient*  mcProduct,
                                                     int                     iAsset) const {
    double strike = performance->strikePerAsset[iAsset];
    const DateTime& startDate = refLevelGen->getDates(iAsset).front();
    const DateTime& today = mcProduct->getToday();
    bool fwdStarting = startDate.isGreater(today);
    const DateTimeArray& assetDates = spotPathGen->getSimSeries()->getDates(iAsset);

    double interpLevel;
    if (fwdStarting){
        interpLevel = strike;
    } else {
        /* not forward starting - some samples have fixed already
           (this includes averaging in) */
        int numInSample = spotPathGen->numDates(iAsset);
        int numRemaining = 
            today.numFutureDates(assetDates);
        // here we're just using iPath = 0 since at this point
        // all the paths are the same
        double ref = refLevel->refLevel(iAsset);
        interpLevel = (numInSample * ref * strike - 
                       sumOut[iAsset]) / numRemaining;
    }

    return new LinearStrikeTSVolRequest(interpLevel,
                                        startDate,
                                        assetDates.back(),
                                        fwdStarting);
}

double PerformanceSV::StateVar::maxDeltaScalingFactor(int iAsset) const {
    double ref = refLevel->refLevel(iAsset);
    return (performance->participationPerAsset[iAsset] / ref);
}

void PerformanceSV::StateVar::update(IStateVariableGen::IStateGen* pathGen) {
    static const string routine = "SVGenBarrierHVEur::StateVar::update";

    try { 
        // Pull the state variables out of the PathGen:
        isPast   = pathGen->doingPast();
        spotPath = spotPathGen->getSpotSV(pathGen);
        refLevel = refLevelGen->getRefLevelSV(refLevel, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

PerformanceSV::StateVar::StateVar(const PerformanceSV*             performance, 
                                  const IRefLevel::IStateVarGenSP& refLevelGen,
                                  const SVGenSpotSP&                  spotPathGen):
performance(performance),
nbAssets(spotPathGen->numAssets()),
sumOut(spotPathGen->numAssets()),
spotPathGen(spotPathGen),
refLevelGen(refLevelGen) {}


double PerformanceSV::StateVar::calcPerf(int iAsset) const {
    static const string method = "PerformanceSV::calcPerf";
    try {
        int    numSamples = spotPathGen->numDates(iAsset);
        double thisPerf = sumOut[iAsset]; // for performance

        const SVPath& path = spotPath->path(iAsset);
        int endIdx = path.end();
        int iStep  = path.begin();

        for (; iStep < endIdx; iStep++) {
            thisPerf += path[iStep];
        }

        if(isPast) { // preserve values
            sumOut[iAsset] = thisPerf;
        }

        // get simple percentage based performance
        double ref = refLevel->refLevel(iAsset);
        thisPerf = (thisPerf / numSamples) / ref - 
            performance->strikePerAsset[iAsset];
        // then cope with performance type 
        double perf = performance->participationPerAsset[iAsset] * 
            IPerformance::Util::calcPerf(performance->perfTypePerAsset[iAsset], thisPerf);

        return perf;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

    
PerformanceSV::PerformanceSV(const StringArray&        perfTypePerAsset,
                             const DoubleArray&        strikePerAsset,
                             const DoubleArray&        participationPerAsset,
                             SVGenSpotSP                  spotPathGen,
                             IRefLevel::IStateVarGenSP refLevelGen):
perfTypePerAsset(perfTypePerAsset), strikePerAsset(strikePerAsset), 
participationPerAsset(participationPerAsset), spotPathGen(spotPathGen), 
refLevelGen(refLevelGen) {
    static const string routine = "PerformanceSV::PerformanceSV";

    try {
        // Validate dimensions
        int numAssets = spotPathGen->numAssets();

        if (numAssets != perfTypePerAsset.size() ||
            numAssets != strikePerAsset.size() ||
            numAssets != participationPerAsset.size()) {
            string message("Mismatch between number of assets (" +
                           Format::toString(numAssets) + 
                           "), number of performance types (" +
                           Format::toString(perfTypePerAsset.size()) +
                           "), number of strikes"+
                           Format::toString(strikePerAsset.size()) +
                           "), and number of participations (" +
                           Format::toString(participationPerAsset.size())+")");
            throw ModelException(routine, message);
        }

        // Validate on a per asset basis
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            const DateTimeArray& assetDates = spotPathGen->getSimSeries()->getDates(iAsset);
            IPerformance::Util::validatePerfFlags(
                perfTypePerAsset[iAsset], strikePerAsset[iAsset], assetDates);
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void PerformanceSV::collectStateVars(IStateVariableCollectorSP svCollector) const {
    static const string routine = "PerformanceSV::collectStateVars";

    try {
        svCollector->append(spotPathGen.get());
        svCollector->append(refLevelGen.get());
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


IStateVariableSP PerformanceSV::create(IStateVariableSP             oldStateVar,
                                    IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "PerformanceSV::create";

    try {
        IStateVarSP oldPerfStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
        return getPerfStateVar(oldPerfStateVar, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


IPerformanceSV::IStateVarSP PerformanceSV::getPerfStateVar(
    IStateVarSP               oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "PerformanceSV::getPerfStateVar";

    try {
        if(oldStateVar.get()) {
            // Update european barrier with new state variables etc.
            StateVar* spotSV = dynamic_cast<StateVar*>(oldStateVar.get());
            if(!spotSV) {
                throw ModelException(
                    "Expected class of type PerformanceSV::StateVar");
            }
            spotSV->update(pathGen);
            return oldStateVar;
        } else {
            // Create a new one
            StateVarSP newSpotSV(new StateVar(this, refLevelGen, spotPathGen));
            newSpotSV->update(pathGen);
            return newSpotSV;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

void PerformanceSV::recordDates(SimSeries* simSeries, int iAsset) const {
    simSeries->addDates(iAsset, spotPathGen->getSimSeries()->getDates(iAsset));
}

PerformanceSVSP PerformanceSV::mergePerfs(const PerformanceSVArray& perfs,
                                          IRefLevel::IStateVarGenSP refLevelGen) {
    int nbAssets = perfs.size();

    StringArray perfTypePerAsset(nbAssets);
    DoubleArray strikePerAsset(nbAssets);
    DoubleArray participationPerAsset(nbAssets);
    SimSeriesSP series(new SimSeries(nbAssets));

    for(int iAsset = 0; iAsset < nbAssets; iAsset++) {
        const PerformanceSV& src = *perfs[iAsset];
        perfTypePerAsset[iAsset] = src.perfTypePerAsset[iAsset];
        strikePerAsset[iAsset]   = src.strikePerAsset[iAsset];
        participationPerAsset[iAsset] = src.participationPerAsset[iAsset];
        series->addDates(iAsset, src.spotPathGen->getSimSeries()->getDates(iAsset));
    }
    
    SVGenSpotSP spotPathGen(new SVGenSpot(series));

    // Create state var generator
    PerformanceSVSP perfGen(new PerformanceSV(
        perfTypePerAsset,
        strikePerAsset,
        participationPerAsset,
        spotPathGen,
        refLevelGen));

    return perfGen;
}
                                      

//////////////////////////////////////////////////////////////////

#endif

DEFINE_TEMPLATE_TYPE_WITH_NAME("AssetPerformanceArray", IAssetPerformanceArray);

void IAssetPerformance::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IAssetPerformance, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IAssetPerformance::TYPE = CClass::registerInterfaceLoadMethod(
    "AssetPerformance", typeid(IAssetPerformance), load);


void IPerformance::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IPerformance, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IPerformance::TYPE = CClass::registerInterfaceLoadMethod(
    "Performance", typeid(IPerformance), load);

/* in header file (for performance) 
   double IPerformance::Util::calcPerf(const string&        perfType,
   double               originalPerf);
*/

void IPerformance::Util::validatePerfFlags(
    const string&        perfType,
    double               strike,
    const DateTimeArray& sampleDates)
{
    static const string routine("IPerformance::Util::validatePerfFlags");
    // check performance type
    if (perfType.empty()){
        throw ModelException(routine, "Blank performance type specified");
    }
    switch (perfType[0]){
    case PERF_TYPE_FORWARD:
    case PERF_TYPE_CALL:
    case PERF_TYPE_PUT:
    case PERF_TYPE_STRADDLE:
        break;
    default:
        throw ModelException(routine,
                             "Unrecognised performance type "+perfType+
                             ". Must be F, C, P or S");
    }
    if (Maths::isNegative(strike)){
        throw ModelException(routine, "Strike ("+
                             Format::toString(strike)+") is negative");
    }
    // validate dates are not empty and are in order 
    DateTime::ensureIncreasing(sampleDates, "sample dates", true);
}

/** Simple performance where perfType, strike, participation, and 
    sampleDates are the same across all assets (a better name anyone?) */
class SimplePerformance: public CObject,
                         virtual public IPerformance,
                         virtual public IAssetPerformance{
private:
    // all fields apply to all assets 
    string             perfType;         
    double             strike;           
    double             participation;   
    DateTimeArray      sampleDates;

public:
    static CClassConstSP const TYPE;

    class MCPerf: virtual public IPerformance::IMCPerf,
                  virtual public IAssetPerformance::IMCPerf {
    private:
        /// fields ///
        const SimplePerformance*  performance;
        mutable DoubleArray       sumOut;
        int                       nbAssets;
        IntArray                  dateMap;
        bool                      isTrivialMap;
    public:
        /** Returns the number of performances that this object will be
            calculating ie the size of the perf parameter in calcPerf() */
        virtual int numPerformances() const{
            return nbAssets;
        }

        /** Calculates the performance of each asset returning the results
            in perf (which must of the required size) */
        virtual void calcPerf(const IMCPathGenerator*  pathGen,
                              int                      iPath,
                              DoubleArray&             perf) const{
            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                perf[iAsset] = calcPerf(pathGen, iPath, iAsset);
            }
        }

        /** Calculates the performance of ith asset */
        virtual double calcPerf(const IMCPathGenerator*  pathGen,
                                int                      iPath,
                                int                      iAsset) const {
            static const string method = "SimplePerformance::calcPerf";
            try {
                double perf;
                int    numSamples = performance->sampleDates.size();
                double thisPerf = sumOut[iAsset]; // for performance
                int    endIdx = pathGen->end(iAsset);
                int    iStep = pathGen->begin(iAsset);
                const double* path = pathGen->Path(iAsset, iPath);
                if (isTrivialMap){ // for performance
                    for (; iStep < endIdx; iStep++){
                        thisPerf += path[iStep];
                    }
                } else {
                    // loop over our dates (a subset of path)
                    for (iStep += dateMap[iStep]; iStep < endIdx; 
                         iStep++, iStep += dateMap[iStep]){
                        thisPerf +=path[iStep];
                    }
                }
                if (pathGen->doingPast()){ // preserve values
                    sumOut[iAsset] = thisPerf;
                }
                // get simple percentage based performance
                double refLevel = pathGen->refLevel(iAsset, iPath);
                thisPerf = (thisPerf / numSamples) /refLevel - 
                    performance->strike;
                // then cope with performance type 
                perf = performance->participation *
                    IPerformance::Util::calcPerf(performance->perfType, 
                                                 thisPerf);
                return perf;
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

        // for the LogNormal path generator
        CVolRequestLN* getVolInterp(const IMCProduct*        mcProduct,
                                    const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
            double strike = performance->strike;
            const IRefLevel* refLevel = mcProduct->getRefLevel();
            const DateTime& startDate = refLevel->getAllDates().front();
            const DateTime& today = mcProduct->getToday();
            bool fwdStarting = startDate.isGreater(today);
            double interpLevel;
            if (fwdStarting){
                interpLevel = strike;
            } else {
                /* not forward starting - some samples have fixed already
                   (this includes averaging in) */
                int numInSample = performance->sampleDates.size();
                int numRemaining = 
                    today.numFutureDates(performance->sampleDates);
                // here we're just using iPath = 0 since at this point
                // all the paths are the same
                interpLevel = (numInSample * pathGen->refLevel(iAsset, 0) * 
                               strike - sumOut[iAsset])/ numRemaining;
            }
            const SimSeries* simSeries = mcProduct->getSimSeries();
            return new LinearStrikeTSVolRequest(interpLevel,
                                                startDate,
                                                simSeries->getLastDate(),
                                                fwdStarting);
        }
        
        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level  */
        virtual double maxDeltaScalingFactor(
            const IMCPathGenerator* futurePathGen, 
            int                     iAsset) const{
            double refLevel = futurePathGen->refLevel(iAsset,
                                                      0 /* path irrelevant*/ );
            return (performance->participation/refLevel);
        }

        MCPerf(const SimplePerformance* performance, 
               const DateTime&          today,
               const SimSeries*         allDates,
               int                      iAsset): 
            performance(performance), 
            sumOut(allDates->getNumAssets()),
            nbAssets(allDates->getNumAssets()) {
            static const string routine("SimplePerformance::MCPerf");
            try{
                /* get mapping for dates (tells us our how dates map onto all
                   simulation dates) */
                dateMap = allDates->createMap(iAsset,
                                              performance->sampleDates, 
                                              isTrivialMap);
            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }
    };
    
    friend class MCPerf;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SimplePerformance, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerformance);
        IMPLEMENTS(IAssetPerformance);
        EMPTY_SHELL_METHOD(defaultSimplePerformance);
        FIELD(perfType, "Performance Type (for all assets)");
        FIELD(strike,   "Strike (for all assets)");
        FIELD(participation, "Participation (for all assets)");
        FIELD(sampleDates,  "Samples (for all assets)");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    SimplePerformance(): CObject(TYPE), strike(0.0), participation(0.0){}

    static IObject* defaultSimplePerformance(){
        return new SimplePerformance();
    }

public:
    IPerformance::IMCPerf* getMCPerf(const DateTime&    today,
                       const SimSeries*   allDates) const{
        return new MCPerf(this, today, allDates, 0); //same dates for all assets
    }

    IAssetPerformance::IMCPerf* getMCPerf(const DateTime&    today,
                       const SimSeries*   allDates,
                       int                iAsset) const{
        return new MCPerf(this, today, allDates, iAsset);
    }

    virtual PerformanceSVSP getStateVarGen(
        IRefLevel::IStateVarGenSP refLevelGen, int nbAssets) {
        // Vectorize everything
        StringArray perfTypePerAsset(nbAssets, perfType);
        DoubleArray strikePerAsset(nbAssets, strike);
        DoubleArray participationPerAsset(nbAssets, participation);
        SVGenSpotSP spotPathGen(new SVGenSpot(nbAssets, sampleDates));

        // Create state var generator
        PerformanceSVSP perfGen(new PerformanceSV(
            perfTypePerAsset,
            strikePerAsset,
            participationPerAsset,
            spotPathGen,
            refLevelGen));

        return perfGen;
    }
    
    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries) const{
        simSeries->addDates(sampleDates);
    }
    virtual void recordDates(SimSeries* simSeries, int iAsset) const{
        simSeries->addDates(iAsset, sampleDates);
    }
    
       
    // validation
    void validatePop2Object(){
        IPerformance::Util::validatePerfFlags(perfType, strike, sampleDates);
    }
};

typedef smartPtr<SimplePerformance> SimplePerformanceSP;
typedef array<SimplePerformanceSP, SimplePerformance> SimplePerformanceArray;
typedef smartPtr<SimplePerformanceArray> SimplePerformanceArraySP;

DEFINE_TEMPLATE_TYPE(SimplePerformanceArray);

CClassConstSP const SimplePerformance::TYPE =
CClass::registerClassLoadMethod("SimplePerformance", 
                                typeid(SimplePerformance), load);


/** Performance where perfType, strike, participation, and 
    sampleDates can be different per asset (a better name anyone?) */
class PerformancePerAsset: public CObject,
                           virtual public IPerformance{
private:
    // all fields apply on a per asset basis
    StringArray         perfTypePerAsset;         
    DoubleArray         strikePerAsset;           
    DoubleArray         participationPerAsset;   
    DateTimeCluster     sampleDatesPerAsset;    
    
public:
    static CClassConstSP const TYPE;

    class MCPerf: virtual public IPerformance::IMCPerf{
    private:
        /// fields ///
        const PerformancePerAsset *performance;
        int                        nbAssets;
        mutable DoubleArray        sumOut;
        vector<IntArray>           dateMapPerAsset;
        bool                       isTrivialMap;
    public:
        /** Returns the number of performances that this object will be
            calculating ie the size of the perf parameter in calcPerf() */
        virtual int numPerformances() const{
            return nbAssets;
        }

        /** Calculates the performance of each asset returning the results
            in perf (which must of the required size) */
        virtual void calcPerf(const IMCPathGenerator*  pathGen,
                              int                      iPath,
                              DoubleArray&             perf) const {
            static const string method = "PerformancePerAsset::calcPerf";
            try {
                for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                    double thisPerf = sumOut[iAsset]; // for performance
                    int    endIdx = pathGen->end(iAsset);
                    int    iStep = pathGen->begin(iAsset);
                    // loop over our dates (a subset of path)
                    const double* path = pathGen->Path(iAsset, iPath);
                    if (isTrivialMap){ // for performance
                        for (; iStep < endIdx; iStep++){
                            thisPerf += path[iStep];
                        }
                    } else {
                        const IntArray& dateMap = dateMapPerAsset[iAsset];
                        for (iStep += dateMap[iStep]; iStep < endIdx; 
                             iStep++, iStep += dateMap[iStep]){
                            thisPerf += path[iStep];
                        }
                    }
                    if (pathGen->doingPast()){ // preserve values
                        sumOut[iAsset] = thisPerf;
                    }
                    // get simple percentage based performance
                    double refLevel = pathGen->refLevel(iAsset, iPath);
                    int numSamples =
                        performance->sampleDatesPerAsset[iAsset].size();
                    thisPerf = (thisPerf / numSamples) /
                        refLevel - performance->strikePerAsset[iAsset];
                    // then cope with performance type on a per asset basis
                    perf[iAsset] = performance->participationPerAsset[iAsset] *
                        IPerformance::Util::calcPerf(
                            performance->perfTypePerAsset[iAsset], thisPerf);

                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }            
        }

        // for the LogNormal path generator
        CVolRequestLN* getVolInterp(const IMCProduct*        mcProduct,
                                    const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
            double strike = performance->strikePerAsset[iAsset];
            const IRefLevel* refLevel = mcProduct->getRefLevel();
            const DateTime& startDate = refLevel->getAllDates().front();
            const DateTime& today = mcProduct->getToday();
            bool fwdStarting = startDate.isGreater(today);
            double interpLevel;
            if (fwdStarting){
                interpLevel = strike;
            } else {
                /* not forward starting - some samples have fixed already
                   (this includes averaging in) */
                int numInSample = 
                    performance->sampleDatesPerAsset[iAsset].size();
                int numRemaining = 
                    today.numFutureDates(performance->
                                         sampleDatesPerAsset[iAsset]);
                // here we're just using iPath = 0 since at this point
                // all the paths are the same
                interpLevel = (numInSample * pathGen->refLevel(iAsset, 0) * 
                               strike - sumOut[iAsset])/ numRemaining;
            }
            const SimSeries* simSeries = mcProduct->getSimSeries();
            return new LinearStrikeTSVolRequest(interpLevel,
                                                startDate,
                                                simSeries->getLastDate(),
                                                fwdStarting);
        }
        
        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level  */
        virtual double maxDeltaScalingFactor(
            const IMCPathGenerator* futurePathGen, 
            int                     iAsset) const{
            double refLevel = futurePathGen->refLevel(iAsset,
                                                      0 /* path irrelevant*/ );
            return (performance->participationPerAsset[iAsset]/refLevel);
        }

        MCPerf(const PerformancePerAsset* performance, 
               const DateTime&            today,
               const SimSeries*           allDates): 
            performance(performance), 
            nbAssets(allDates->getNumAssets()),
            sumOut(nbAssets),
            dateMapPerAsset(nbAssets){
            static const string routine("PerformancePerAsset::MCPerf");
            try{
                /* get mapping for dates (tells us our how dates map onto all
                   simulation dates) */
                isTrivialMap = true;
                for (int i = 0; i < nbAssets; i++){
                    bool assetIsTrivialMap = false;
                    dateMapPerAsset[i] = allDates->createMap(
                        i, performance->sampleDatesPerAsset[i], 
                        assetIsTrivialMap);
                    if (isTrivialMap && !assetIsTrivialMap){
                        isTrivialMap = false;
                    }
               }
            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }
    };
    
    friend class MCPerf;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerformancePerAsset, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerformance);
        EMPTY_SHELL_METHOD(defaultPerformancePerAsset);
        FIELD(perfTypePerAsset, "Performance Type (for all assets)");
        FIELD(strikePerAsset,   "Strike (for all assets)");
        FIELD(participationPerAsset, "Participation (for all assets)");
        FIELD(sampleDatesPerAsset,  "Samples (for all assets)");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    PerformancePerAsset(): CObject(TYPE) {}

    static IObject* defaultPerformancePerAsset(){
        return new PerformancePerAsset();
    }

public:
    IMCPerf* getMCPerf(const DateTime&    today,
                       const SimSeries*   allDates) const{
        return new MCPerf(this, today, allDates);
    }

    virtual PerformanceSVSP getStateVarGen(
        IRefLevel::IStateVarGenSP refLevelGen, int nbAssets) {
        
        // Create a simseries object and stuff in the dates per asset
        SimSeriesSP series(new SimSeries(nbAssets));
        series->addDates(sampleDatesPerAsset);
        SVGenSpotSP spotPathGen(new SVGenSpot(series));

        // Create state var generator
        PerformanceSVSP perfGen(new PerformanceSV(
            perfTypePerAsset,
            strikePerAsset,
            participationPerAsset,
            spotPathGen,
            refLevelGen));

        return perfGen;        
    }
    
    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries) const{
        simSeries->addDates(sampleDatesPerAsset);
    }

    // validation
    void validatePop2Object(){
        static const string routine("PerformancePerAsset::validatePop2Object");
        if (perfTypePerAsset.size() != strikePerAsset.size() ||
            perfTypePerAsset.size() != sampleDatesPerAsset.size()){
            string m("Mismatch between number of performance types ("+
                     Format::toString(perfTypePerAsset.size())+"), number "
                     "of strikes"+
                     Format::toString(strikePerAsset.size())+"), and number "
                     "of arrays of dates ("+
                     Format::toString(sampleDatesPerAsset.size())+")");
            throw ModelException(routine, m);
        }
        for (int i = 0; i < perfTypePerAsset.size(); i++){
            try{
                IPerformance::Util::validatePerfFlags(perfTypePerAsset[i], 
                                                      strikePerAsset[i],
                                                      sampleDatesPerAsset[i]);
            } catch (exception& e){
                string m("Failed for asset number "+Format::toString(i));
                throw ModelException(e, routine, m);
            }
        }
    }
};

CClassConstSP const PerformancePerAsset::TYPE =
CClass::registerClassLoadMethod("PerformancePerAsset", 
                                typeid(PerformancePerAsset), load);


/** Performance where perfType, strike, participation, and 
    sampleDates can be different per asset (a better name anyone?) */
class PerformanceByAsset: public CObject,
                          virtual public IPerformance{
private:
    IAssetPerformanceArray perfPerAsset;
    
public:
    static CClassConstSP const TYPE;

    PerformanceByAsset(IAssetPerformanceArray& perfs) : CObject(TYPE),
                                                        perfPerAsset(perfs){}

    class MCPerf: virtual public IPerformance::IMCPerf{
    private:
        /// fields ///
        vector<IAssetPerformance::IMCPerfSP> perfByAsset;
        int                                  nbAssets;
    public:
        /** Returns the number of performances that this object will be
            calculating ie the size of the perf parameter in calcPerf() */
        virtual int numPerformances() const{
            return nbAssets;
        }

        /** Calculates the performance of each asset returning the results
            in perf (which must of the required size) */
        virtual void calcPerf(const IMCPathGenerator*  pathGen,
                              int                      iPath,
                              DoubleArray&             perf) const{
            static const string method = "PerformanceByAsset::calcPerf";
            try {
                for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                    perf[iAsset] = 
                        perfByAsset[iAsset]->calcPerf(pathGen, iPath, iAsset);
                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

        // for the LogNormal path generator
        CVolRequestLN* getVolInterp(const IMCProduct*        mcProduct,
                                    const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
            return perfByAsset[iAsset]->getVolInterp(mcProduct, pathGen, 
                                                     iAsset);
        }
        
        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level  */
        virtual double maxDeltaScalingFactor(
            const IMCPathGenerator* futurePathGen, 
            int                     iAsset) const{
            return perfByAsset[iAsset]->maxDeltaScalingFactor(futurePathGen,
                                                              iAsset);
        }

        MCPerf(const PerformanceByAsset* performance, 
               const DateTime&            today,
               const SimSeries*           allDates): 
            nbAssets(allDates->getNumAssets()){
            static const string routine("PerformanceByAsset::MCPerf");
            try {
                for (int i = 0; i < nbAssets; i++){
                    IAssetPerformance::IMCPerfSP mcp(
                        performance->perfPerAsset[i]->getMCPerf(today,
                                                                allDates,
                                                                i));
                    perfByAsset.push_back(mcp);
                }
            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }
    };
    
    friend class MCPerf;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerformanceByAsset, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerformance);
        EMPTY_SHELL_METHOD(defaultPerformanceByAsset);
        FIELD(perfPerAsset, "Performance (for all assets)");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        registerObjectFromArrayMethod(IAssetPerformanceArray::TYPE,
                                      TYPE,
                                      &fromArray);
        registerObjectFromArrayMethod(SimplePerformanceArray::TYPE,
                                      TYPE,
                                      &fromSimpleArray);
    }
     
    // for reflection
    PerformanceByAsset(): CObject(TYPE) {}

    static IObject* defaultPerformanceByAsset(){
        return new PerformanceByAsset();
    }

public:
    IMCPerf* getMCPerf(const DateTime&    today,
                       const SimSeries*   allDates) const{
        return new MCPerf(this, today, allDates);
    }

    virtual PerformanceSVSP getStateVarGen(
        IRefLevel::IStateVarGenSP refLevelGen, int nbAssets) {
        
        PerformanceSVArray perfs(nbAssets);
        for(int iAsset = 0; iAsset < nbAssets; iAsset++) {
            perfs[iAsset] = perfPerAsset[iAsset]->getStateVarGen(refLevelGen, nbAssets);
        }

        PerformanceSVSP perfGen = PerformanceSV::mergePerfs(perfs, refLevelGen);
        
        return perfGen;
    }
    
    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries) const{
        DateTimeArrayArray(perfPerAsset.size());
        for (int i = 0; i < perfPerAsset.size(); i++){
            perfPerAsset[i]->recordDates(simSeries, i);
        }
    }

    static IObjectSP fromArray(const IObjectSP& object, 
                               CClassConstSP requiredType) {
        IAssetPerformanceArray& perfs = 
            dynamic_cast<IAssetPerformanceArray&>(*object);
        return IObjectSP(new PerformanceByAsset(perfs));
    }

    static IObjectSP fromSimpleArray(const IObjectSP&  object,
                                     CClassConstSP     requiredType) {
        SimplePerformanceArray& sp = 
            dynamic_cast<SimplePerformanceArray&>(*object);
        IAssetPerformanceArray perfs(sp.size());
        for (int i = 0; i < sp.size(); i++) {
            perfs[i] = IAssetPerformanceSP::dynamicCast((IObjectSP)sp[i]);
        }
        return IObjectSP(new PerformanceByAsset(perfs));
    }  
};

CClassConstSP const PerformanceByAsset::TYPE =
CClass::registerClassLoadMethod("PerformanceByAsset", 
                                typeid(PerformanceByAsset), load);


DEFINE_TEMPLATE_TYPE_WITH_NAME("PerformanceArray", IPerformanceArray);


DRLIB_END_NAMESPACE

    
