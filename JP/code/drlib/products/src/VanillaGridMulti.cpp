//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : VanillaGridMulti.cpp
//
//   Description : Multi Factor Version of VanillaGrid
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/VanillaGridMulti.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE_WITH_NAME("IAggregateMakerArray", IAggregateMakerArray);
DEFINE_TEMPLATE_TYPE_WITH_NAME("IDoubleArrayModifierMakerArray", IDoubleArrayModifierMakerArray);

/** validation */
void VanillaGridMulti::validatePop2Object(){
    static const string method("VanillaGridMulti::validatePop2Object");
    try {
        /** call parent method */
        GenericNFBase::validatePop2Object();
        
        /** check dimensions */
        int nbBaskets = baskets->size();
        int nbPerfTypes = perfTypes->size();
        if ((nbBaskets>1) && (nbPerfTypes>1)) {
            if (!(nbBaskets==nbPerfTypes)) {
                throw ModelException(method, "If more than one basket type and more than one performance type "
                    "are supplied, then the two arrays have to have the same length.");
            }
        }

        nbPerfPerMat = max(nbBaskets,nbPerfTypes);

        /** check that average out dates are increasing */
        DateTime::ensureIncreasing(avgOutDates, "Average Out Dates", false /*failIfEmpty*/);        
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}
   
/** mandatory, since derivation from GenericNFBase */ 
const DateTimeArray VanillaGridMulti::samplingDates() const {
    return avgOutDates;
}

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* VanillaGridMulti::createProduct(const MonteCarlo* model) const {
    static const string method("VanillaGridMulti::createProduct");
    try {
        // we need to create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
        simSeries->addDates(avgOutDates);
        return new VanillaGridMultiMC(this, simSeries);
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** additional constructor */
VanillaGridMulti::VanillaGridMulti(VanillaGridMulti*                grid,
                                   const DateTimeArray&             avgOutDate,
                                   IAggregateMakerArraySP           basket,
                                   IDoubleArrayModifierMakerArraySP perfType) : GenericNFBase(TYPE){
    this->valueDate = grid->valueDate;
    this->notional = grid->notional;
    this->instSettle = grid->instSettle;
    this->assets = grid->assets;
    this->discount = grid->discount;
    this->refLevel = grid->refLevel;
    this->pastValues = grid->pastValues;
    this->isdaDateAdjust = grid->isdaDateAdjust;
    this->avgOutDates = avgOutDate;
    this->baskets = basket;
    this->perfTypes = perfType;
}
   
VanillaGridMulti::VanillaGridMulti(): GenericNFBase(TYPE) {}

IObject* VanillaGridMulti::defaultVanillaGridMulti(){
    return new VanillaGridMulti();
}

/** Invoked when Class is 'loaded' */
void VanillaGridMulti::load(CClassSP& clazz){
    REGISTER(VanillaGridMulti, clazz);
    SUPERCLASS(GenericNFBase);
    IMPLEMENTS(IMCIntoProduct);
    EMPTY_SHELL_METHOD(defaultVanillaGridMulti);                
    FIELD(avgOutDates, "avg out dates");
    FIELD(baskets, "avg basket or rainbow basket"); 
    FIELD(perfTypes, "performances");        
    FIELD(nbPerfPerMat, "nb performances per maturity");
    FIELD_MAKE_TRANSIENT(nbPerfPerMat);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const VanillaGridMulti::TYPE = CClass::registerClassLoadMethod(
    "VanillaGridMulti", typeid(VanillaGridMulti), VanillaGridMulti::load);

bool VanillaGridMultiLoad() {
    return (VanillaGridMulti::TYPE != 0);
}

/** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
void VanillaGridMultiMC::initialiseLN(const  MCPathGenerator*  pathGen) const {}

/** equivalent to InstIntoMCProduct. Need to call parent's constructor */
VanillaGridMultiMC::VanillaGridMultiMC(const VanillaGridMulti*  inst,
                                       const SimSeriesSP&       simSeries) :
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
        basketHelper(nbAssets, 0.0),
        performanceHelper(0.0),
        refLevel(nbAssets, 0.0),
        refLevelSoFar(nbAssets, 0.0),
        futureBeginIdx(0), 
        countPaths(0.0) {         
    static const string method("VanillaGridMultiMC::VanillaGridMultiMC");
    try {
        int nbAvgOutDates = inst->avgOutDates.size();
        int nbPerfTypes = inst->perfTypes->size();
        int nbBaskets = inst->baskets->size();

        int nbPerformances = max(nbPerfTypes, nbBaskets);                
        
        /** allocate space for value and value_ste */
        payoffMatrix = CDoubleMatrixSP(new DoubleMatrix(nbPerformances, nbAvgOutDates));
        payoffMatrixSq = CDoubleMatrixSP(new DoubleMatrix(nbPerformances, nbAvgOutDates));
        payoffMatrixSqHelper = CDoubleMatrixSP(new DoubleMatrix(nbPerformances, nbAvgOutDates));
    
        basket.resize(nbPerformances);
        performance.resize(nbPerformances);

        for (int iPerf=0; iPerf<nbPerformances; iPerf++) {
            basket[iPerf] = (nbBaskets>1) ? 
                IAggregateSP((*inst->baskets)[iPerf]->getAggregate(&basketHelper)) : 
                IAggregateSP((*inst->baskets)[0]->getAggregate(&basketHelper));
            performance[iPerf] = (nbPerfTypes>1) ? 
                IDoubleArrayModifierSP((*inst->perfTypes)[iPerf]->getModifier(&performanceHelper)) :
                IDoubleArrayModifierSP((*inst->perfTypes)[0]->getModifier(&performanceHelper));
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Called within the simulation loop */
void VanillaGridMultiMC::payoff(const MCPathGenerator*  pathGen,
                                IMCPrices&                prices) {
    static const string method("VanillaGridMultiMC::payoff");
    try {
        int     beginIdx = pathGen->begin(0); // same for all assets
        int     endIdx   = pathGen->end(0);                        

        refLevel = refLevelSoFar;

        /** determination of reference level */
        for (int iAsset = 0; iAsset < nbAssets; iAsset ++) {
            refLevel[iAsset] = pathGen->refLevel(iAsset, 0);
        }                    
        
        SimpleDoubleArray tmpBasketHelper(nbAssets,0.0);
        for (int iStep = beginIdx; iStep < endIdx; iStep++) {
            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                basketHelper[iAsset] = pathGen->Path(iAsset,0)[iStep] / refLevel[iAsset];				
            }
            // copy 
            tmpBasketHelper = basketHelper;
            /** aggregation */
            int nbPerformances = payoffMatrix->numCols();
            for (int iPerf = 0; iPerf < nbPerformances; iPerf++) {
                performanceHelper() = basket[iPerf]->aggregate();
                performance[iPerf]->apply();
                (*payoffMatrix)[iPerf][iStep] += performanceHelper();
                if ((int)(countPaths) % 2) {
                    (*payoffMatrixSqHelper)[iPerf][iStep] += performanceHelper();
                    (*payoffMatrixSqHelper)[iPerf][iStep] /= 2.0; 
                    (*payoffMatrixSq)[iPerf][iStep] 
                        += (*payoffMatrixSqHelper)[iPerf][iStep]*(*payoffMatrixSqHelper)[iPerf][iStep];
                } else {
                    (*payoffMatrixSqHelper)[iPerf][iStep] = performanceHelper(); 
                }
                // original
                basketHelper = tmpBasketHelper;
            }                
        }
        
        if (pathGen->doingPast()){ 
            countPaths = 0.0;
            refLevelSoFar = refLevel;
        } else {
            futureBeginIdx = beginIdx;
            countPaths += 1.0;
        }

        performanceHelper() = basket[0]->aggregate();
        performance[0]->apply();
        prices.add(performanceHelper()); 

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** invoked after final simulated path is run. */
void VanillaGridMultiMC::recordExtraOutput(CControl*        control,
                                           Results*         results,
                                           const IMCPrices& prices) const {
    static const string method("VanillaGridMultiMC::recordExtraOutput");
    try {
        bool hasFuture = inst->avgOutDates.back().isGreater(inst->valueDate);
        if (control) {
            if (hasFuture) {
                /** discounting has to be manullay */
                int nbSteps = payoffMatrix->numRows();
                int nbPerfs = payoffMatrix->numCols();
                CDoubleMatrixSP stdErrMatrix(new CDoubleMatrix(nbPerfs, nbSteps));
                for (int iStep=futureBeginIdx; iStep<nbSteps && Maths::isPositive(countPaths); iStep++) {
                    double thisDiscount = inst->discount->pv(inst->valueDate, inst->avgOutDates[iStep]);
                    for (int iPerf=0; iPerf<nbPerfs; iPerf++) {
                        (*payoffMatrix)[iPerf][iStep] /= countPaths;
                        (*payoffMatrix)[iPerf][iStep] *= thisDiscount;
                        (*payoffMatrixSq)[iPerf][iStep] *= thisDiscount*thisDiscount;
                        (*payoffMatrixSq)[iPerf][iStep] /= (countPaths/2.0); 
                        (*payoffMatrixSq)[iPerf][iStep] -= (*payoffMatrix)[iPerf][iStep]*(*payoffMatrix)[iPerf][iStep];
                        (*payoffMatrixSq)[iPerf][iStep] /= (countPaths/2.0-1.0);
                        if (Maths::isPositive((*payoffMatrixSq)[iPerf][iStep])) {
                            (*payoffMatrixSq)[iPerf][iStep] = sqrt((*payoffMatrixSq)[iPerf][iStep]);
                        } else {
                            (*payoffMatrixSq)[iPerf][iStep] = 0.0;
                        }
                    }
                } 
            } else {
                payoffMatrix->fill(0.0); 
                payoffMatrixSq->fill(0.0);
            }
            OutputRequest* request = control->requestsOutput(OutputRequest::OPTION_PRICE);
            if (request) {
                VanillaGrid::OutputSP pricesGrid;
                pricesGrid = VanillaGrid::OutputSP(new 
                    VanillaGrid::Output(payoffMatrix->numCols(), payoffMatrix->numRows()));
                
                for (int iStep = 0 ; iStep < payoffMatrix->numRows(); iStep++) {                    
                    for (int iPerf = 0 ; iPerf < payoffMatrix->numCols() ; iPerf++) {                        
                        pricesGrid->setValue(iPerf, iStep, (*payoffMatrix)[iPerf][iStep]);
                    }
                }
                results->storeRequestResult(request, pricesGrid);                
            }
            request = control->requestsOutput(OutputRequest::VALUE_STE);
            if (request) {
                VanillaGrid::OutputSP pricesGrid;
                pricesGrid = VanillaGrid::OutputSP(new 
                    VanillaGrid::Output(payoffMatrixSq->numCols(), payoffMatrixSq->numRows()));
                
                for (int iStep = 0 ; iStep < payoffMatrixSq->numRows(); iStep++) {                    
                    for (int iPerf = 0 ; iPerf < payoffMatrixSq->numCols() ; iPerf++) {                        
                        pricesGrid->setValue(iPerf, iStep, (*payoffMatrixSq)[iPerf][iStep]);
                    }
                }
                results->storeRequestResult(request, pricesGrid);                    
            }
        }
    } catch (exception &e) {
        throw ModelException(e, method);
    }
}

/** need interpolation level for the LogNormal path generator */ 
CVolRequestLNArray VanillaGridMultiMC::getVolInterp(const MCPathGenerator* pathGen,
                                                    int                     iAsset) const {
    const static string method("VanillaGridMultiMC::getVolInterp");
    try {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here
        
        double interpLevel = (*inst->perfTypes)[0]->getInterpLevel(iAsset);

        const DateTime& today = getToday();
        const DateTime& startDate = getRefLevel()->getAllDates().front();
        bool fwdStarting = startDate.isGreater(today);
        if (!fwdStarting) {        
            interpLevel *= pathGen->refLevel(iAsset, 0);
        }

        if (!fwdStarting){
            interpLevel = interpLevel * pathGen->refLevel(iAsset, 0);
        }
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                    startDate, 
                                                                    inst->avgOutDates.back(),
                                                                    fwdStarting));
        return reqarr;
    } catch (exception &e) {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE