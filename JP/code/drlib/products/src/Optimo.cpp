/*
 Product description

The Optimo is an option on the best of several baskets.
The baskets have the same underlyings but different weight associated to each underlying.
Each basket performance is the maximum of the average basket level on different observation dates.

*/

#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

class Optimo: public GenericNFBase, 
              virtual public IMCIntoProduct {
protected:
    /// fields ////////
    DoubleArrayArray            weights;           // weights for each strategy
    DateTimeArray               averagingDates;    // dates used to calculate the level of the underlyings
    DateTimeArray               monitoringDates;   // dates when the level of the baskets is observed basket
    IAggregateMakerSP           basketsAggregate;  // on each monitoring date, agregate the performance of the different baskets
    IAggregateMakerSP           timeAgregate;      // agregate the performances calculated on each monitoring date
    IDoubleArrayModifierMakerSP overallOption;     // option at maturity
    
public:
    static CClassConstSP const TYPE;
    friend class OptimoMC;
    friend class OptimoSVMC;

    // validation
    void validatePop2Object(){
        static const string method = "Optimo::validatePop2Object";
        GenericNFBase::validatePop2Object();

        // validate dates are not empty - order is handled by SimSeries
        if (monitoringDates.empty()) {
            throw ModelException(method, "No monitoring dates supplied!");
        }
        
        // check monitoring dates are amongst the averaging dates
        if (!DateTime::isSubset(averagingDates, monitoringDates)) {
            throw ModelException(method, "Monitoring dates should be a subset of averaging dates");
        }
        
        if (weights.size()!=assets->numFactors()) {
            throw ModelException(method, "There number of assets ("+Format::toString(assets->numFactors())+
                                 ") should be equal to the number of array of weights ("+Format::toString(weights.size())+")");
        }
        else if (weights[0].size()==0) {
            throw ModelException(method, "The asset 0 should have at least one weight");
        }
        else {
            int nbStrategies = weights[0].size();
            for (int i=1;i<weights.size();i++) {
                if (weights[i].size()!= nbStrategies) {
                    throw ModelException(method, "There number of weights ("+Format::toString(nbStrategies)+
                                         ") for the asset 0 should be equal to the number of weights ("+
                                         Format::toString(weights[i].size())+") for the asset "+
                                         Format::toString(i));
                }
            }
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return averagingDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    Optimo(): GenericNFBase(TYPE) {} // for reflection
    Optimo(const Optimo& rhs);     // not implemented
    Optimo& operator=(const Optimo& rhs); // not implemented

    static IObject* defaultOptimo(){
        return new Optimo();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Optimo, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultOptimo);
        FIELD(weights, "Weights defining each strategies, given for each asset");
        FIELD(averagingDates, "Dates used to calculate the averaged performance of each asset");
        FIELD(monitoringDates, "Dates when the performance of each strategy is observed");
        FIELD(basketsAggregate, "Defines the way the different strategy performances are aggregated on a given monitoring date");
        FIELD(timeAgregate, "Defines the way the performance on each monitoring dates are aggregated");
        FIELD(overallOption, "Overall option on the final performance at maturity");
    }
};

/* MC product class for Optimo */
class OptimoMC: public IMCProduct,
                virtual public IMCProductLN,
                virtual public IMCProductImplied {
private:
    const Optimo*      inst;                   // Instrument
    int                nbAssets;               // convenient
    int                nbBaskets;              // convenient
    int                nbMonitorings;          // convenient
    DoubleArray        assetsPerf;             // [nbAssets] current perf per asset
    DoubleArray        sum;                    // [nbAssets] sum of the performances of each asset
    DoubleArray        refLevel;               // [nbAssets], saves alloc 
    DoubleArrayArray   weightsPerBaskets;      // [nbBaskets] X [nbAssets]
    SimpleDoubleArray  monitoringPerformances; // [nbMonitorings]

    // maintain state of the instrument in the past
    DoubleArray        sumSoFar;                    // [nbAssets]
    DoubleArray        refLevelSoFar;               // [nbAssets], saves alloc 
    SimpleDoubleArray  monitoringPerformancesSoFar; // [nbMonitorings]
    int                iMonitoringSoFar;
    
    // Operational aggregation and performance calcs
    IAggregateSP           basketsAggregate;
    IAggregateSP           timeAgregate;
 
    IDoubleArrayModifierSP overallOption;
    SimpleDoubleArray  basketsPerformance;     // [nbBaskets] performance for each basket
    TrivialDoubleArray finalPerformance;

    IntArray               nbAvgPerMonitoringDate; // [nbMonitorings]
    IntArray               monitoringMap;          // [nbAvgDates+1] - convenient way to track coupon dates
    
public:
    
    OptimoMC(const Optimo*      inst,
             const SimSeriesSP&   simSeries):
        IMCProduct(inst->assets.get(),
                   inst->valueDate,
                   inst->discount.get(),
                   inst->refLevel,
                   simSeries,
                   inst->pastValues,
                   inst->instSettle.get(),
                   simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        nbBaskets(inst->weights[0].size()),
        nbMonitorings(inst->monitoringDates.size()),
        assetsPerf(nbAssets,0.0),
        sum(nbAssets,0.0),
        refLevel(nbAssets, 0.0),
        weightsPerBaskets(nbBaskets,DoubleArray(nbAssets,0.0)),
        sumSoFar(nbAssets,0.0),
        refLevelSoFar(nbAssets, 0.0),
        monitoringPerformances(inst->monitoringDates.size(),0.0),
        monitoringPerformancesSoFar(inst->monitoringDates.size(),0.0),
        iMonitoringSoFar(0),
        basketsPerformance(nbBaskets, 0.0), 
        finalPerformance(0.0),
        nbAvgPerMonitoringDate(nbMonitorings,0)
        {
            // transpose the weigth matrix - more convenient
            for(int i=0;i<nbAssets;i++){
                for(int j=0;j<nbBaskets;j++) {
                    weightsPerBaskets[j][i] = inst->weights[i][j];
                }
            } 

            // cumulative number of Averaging Dates
            basketsAggregate = IAggregateSP(inst->basketsAggregate->getAggregate(&basketsPerformance));
            timeAgregate     = IAggregateSP(inst->timeAgregate->getAggregate(&monitoringPerformances));
            overallOption    = IDoubleArrayModifierSP(inst->overallOption->getModifier(&finalPerformance));
            
            bool isTrivial;
            monitoringMap = DateTime::createMapping(inst->averagingDates,
                                                    inst->monitoringDates,
                                                    isTrivial);
        }
    
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&             prices) {
        
        // update with preserved values from doing past
        
        int beginIdx = pathGen->begin(0); // same for all assets
        int endIdx   = pathGen->end(0);
        int iAsset, iStep;

        DoubleArray avg(nbAssets,0.0);
        
        sum = sumSoFar;
        int iMonitoring = iMonitoringSoFar;
        monitoringPerformances = monitoringPerformancesSoFar;

        for(iAsset=0; iAsset<nbAssets; iAsset++) {
            refLevel[iAsset] =  pathGen->refLevel(iAsset, 0/*iPath*/);
        }
        
        // Form the asset basket at each monitoring date first
        for (iStep=beginIdx; iStep<endIdx; iStep++) { 
            bool isMonitoringDate = (monitoringMap[iStep]==0);   // true iff a monitoring date
            
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];
                if (isMonitoringDate) {
                    assetsPerf[iAsset] = sum[iAsset] / ((iStep+1)* refLevel[iAsset]);
                }
            }
            
            if (isMonitoringDate) {
                for (int iBasket=0;iBasket<nbBaskets;iBasket++) {
                    basketsPerformance[iBasket] = Maths::ArrayScalarProduct(weightsPerBaskets[iBasket],assetsPerf);
                }
                // calculate baskets agregate
                monitoringPerformances[iMonitoring] = basketsAggregate->aggregate();
                iMonitoring ++;
            }
        }

        // preserve values for past
        if (doingPast()){ 
            sumSoFar = sum;
            refLevelSoFar = refLevel;
            iMonitoringSoFar = iMonitoring;
            monitoringPerformancesSoFar = monitoringPerformances;
        }

        if (!doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            
            finalPerformance() = timeAgregate->aggregate();
            overallOption->apply();
            
            // now scale by notional
            prices.add(inst->notional * finalPerformance()); 
        }
    }
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        static const string routine = "OptimoMC::initialiseLN";
        throw ModelException(routine, "Methodology not supported");
    }

    // any old level so that MC implied works
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "OptimoMC::getVolInterp";

        try {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = getRefLevel()->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);

            double interpLevel  = 1.0;

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));

            return reqarr;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

};

/* MC product class for OptimoSV */
class OptimoSVMC: public MCProductClient,
                  virtual public IMCProductLN,
                  virtual public IMCProductImplied {
private:
    const Optimo*      inst;                   // Instrument
    int                nbAssets;               // convenient
    int                nbBaskets;              // convenient
    int                nbMonitorings;          // convenient
    DoubleArray        assetsPerf;             // [nbAssets] current perf per asset
    DoubleArray        sum;                    // [nbAssets] sum of the performances of each asset
    DoubleArray        refLevel;               // [nbAssets], saves alloc 
    DoubleArrayArray   weightsPerBaskets;      // [nbBaskets] X [nbAssets]
    SimpleDoubleArray  monitoringPerformances; // [nbMonitorings]

    // maintain state of the instrument in the past
    DoubleArray        sumSoFar;                    // [nbAssets]
    DoubleArray        refLevelSoFar;               // [nbAssets], saves alloc 
    SimpleDoubleArray  monitoringPerformancesSoFar; // [nbMonitorings]
    int                iMonitoringSoFar;
    
    // Operational aggregation and performance calcs
    IAggregateSP           basketsAggregate;
    IAggregateSP           timeAgregate;
 
    IDoubleArrayModifierSP overallOption;
    SimpleDoubleArray  basketsPerformance;     // [nbBaskets] performance for each basket
    TrivialDoubleArray finalPerformance;

    IntArray               nbAvgPerMonitoringDate; // [nbMonitorings]
    IntArray               monitoringMap;          // [nbAvgDates+1] - convenient way to track coupon dates

    // State variables and generators
    SVGenSpotSP                  spotGen;      // Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  // Generator for ref level
    SVGenDiscFactorSP            dfGen;        // Generator for discount factors
    SVGenSpot::IStateVarSP       spotSV;       // Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   // Ref level state variable
    SVDiscFactorSP         dfSV;         // Df state variable

public:

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        static const string routine = "OptimoSVMC::collectStateVars";
        try{
            svCollector->append(spotGen.get());             // spot level
            svCollector->append(refLevelGen.get());         // reference level
            svCollector->append(dfGen.get());               // and a DiscFactor one
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) {
        static const string routine = "OptimoSVMC::pathGenUpdated";
        try{
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    OptimoSVMC(const Optimo*      inst,
               const SimSeriesSP& simSeries):
        MCProductClient(inst->assets.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        inst->refLevel,
                        simSeries,
                        inst->pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        nbBaskets(inst->weights[0].size()),
        nbMonitorings(inst->monitoringDates.size()),
        assetsPerf(nbAssets,0.0),
        sum(nbAssets,0.0),
        refLevel(nbAssets, 0.0),
        weightsPerBaskets(nbBaskets,DoubleArray(nbAssets,0.0)),
        sumSoFar(nbAssets,0.0),
        refLevelSoFar(nbAssets, 0.0),
        monitoringPerformancesSoFar(inst->monitoringDates.size(),0.0),
        monitoringPerformances(inst->monitoringDates.size(),0.0),
        iMonitoringSoFar(0),
        basketsPerformance(nbBaskets, 0.0), 
        finalPerformance(0.0),
        nbAvgPerMonitoringDate(nbMonitorings,0),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())) {
        // transpose the weigth matrix - more convenient
        for(int i=0;i<nbAssets;i++){
            for(int j=0;j<nbBaskets;j++) {
                weightsPerBaskets[j][i] = inst->weights[i][j];
            }
        } 

        // cumulative number of Averaging Dates
        basketsAggregate = IAggregateSP(inst->basketsAggregate->getAggregate(&basketsPerformance));
        timeAgregate     = IAggregateSP(inst->timeAgregate->getAggregate(&monitoringPerformances));
        overallOption    = IDoubleArrayModifierSP(inst->overallOption->getModifier(&finalPerformance));
        
        bool isTrivial;
        monitoringMap = DateTime::createMapping(inst->averagingDates,
                                                inst->monitoringDates,
                                                isTrivial);
    }
    
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&             prices) {
        
        DoubleArray avg(nbAssets,0.0);
        
        sum = sumSoFar;
        int iMonitoring = iMonitoringSoFar;
        int iAsset,iStep;
        monitoringPerformances = monitoringPerformancesSoFar;

        // Same begin & end for all assets, so read from the first
        const SVPath& path = spotSV->path(0);
        int    beginIdx = path.begin();
        int    endIdx   = path.end();

        for(iAsset=0; iAsset<nbAssets; iAsset++) {
            refLevel[iAsset] =  refLevelSV->refLevel(iAsset);
        }

        // Form the asset basket at each monitoring date first
        for (iStep=beginIdx; iStep<endIdx; iStep++) { 
            bool isMonitoringDate = (monitoringMap[iStep]==0);   // true iff a monitoring date
            
            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                const SVPath& path = spotSV->path(iAsset);
                sum[iAsset] += path[iStep];
                if (isMonitoringDate) {
                    assetsPerf[iAsset] = sum[iAsset] / ((iStep+1)* refLevel[iAsset]);
                }
            }

            if (isMonitoringDate) {
                for (int iBasket=0;iBasket<nbBaskets;iBasket++) {
                    basketsPerformance[iBasket] = Maths::ArrayScalarProduct(weightsPerBaskets[iBasket],assetsPerf);
                }
                // calculate baskets aggregate
                monitoringPerformances[iMonitoring] = basketsAggregate->aggregate();
                iMonitoring ++;
            }
        }
        
        // preserve values for past
        if (doingPast()){ 
            sumSoFar = sum;
            refLevelSoFar = refLevel;
            iMonitoringSoFar = iMonitoring;
            monitoringPerformancesSoFar = monitoringPerformances;
        }

        if (!doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            
            finalPerformance() = timeAgregate->aggregate();
            overallOption->apply();
            
            // now scale by notional
            prices.add(inst->notional * finalPerformance() * dfSV->firstDF()); 
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseImplied(const IMCPathGenerator* pathGen)const{}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const IMCPathGenerator* pathGen)const{
        static const string routine = "OptimoSVMC::initialiseLN";
        throw ModelException(routine, "Methodology not supported");
    }
    
    // any old level so that MC implied works
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "OptimoSVMC::getVolInterp";

        try {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = refLevelGen->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);
            
            double interpLevel  = 1.0;
            
            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
            
            return reqarr;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

//////////////////////////////////////////////////////////////////////////

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* Optimo::createProduct(const MonteCarlo* model) const {
    static const string method = "Optimo::createProduct";

    try {
        int nbAssets = assets->NbAssets();
        
        // Create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(nbAssets));

        // do we crop these dates to only be future ones??
        simSeries->addDates(averagingDates);
    
        if(model->stateVarUsed()) {
            return new OptimoSVMC(this, simSeries);
        } else {
            // Otherwise, use old methodology
            return new OptimoMC(this, simSeries);
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


CClassConstSP const Optimo::TYPE = CClass::registerClassLoadMethod(
    "Optimo", typeid(Optimo), Optimo::load);

// * for class loading (avoid having header file) */
bool OptimoLoad() {
    return (Optimo::TYPE != 0);
}

DRLIB_END_NAMESPACE
