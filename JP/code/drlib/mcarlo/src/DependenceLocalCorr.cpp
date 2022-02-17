//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceLocalCorr.hpp
//
//   Description : Holds LocalCorr Dependence (Maker)
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_DEPENDENCE_LOCAL_CORR_CPP
#include "edginc/DependenceLocalCorr.hpp"

#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/MDFUtil.hpp"

#include "edginc/MCRandom.hpp"

#include <fstream>

DRLIB_BEGIN_NAMESPACE

/*************************/
/* LOCAL CORR DEPENDENCE */
/*************************/

/** Destructor */
LocalCorr::~LocalCorr(){}

/** Default Constructor */
LocalCorr::LocalCorr(){}

/** Proper Constructor */
LocalCorr::LocalCorr(DateTimeSP             valueDate,
                     int                    nbIF,
                     IMCRandomSP            randomGenIF,
                     DateTimeArraySP        timelineIF,
                     int                    nbMF,
                     IMCRandomSP            randomGenMF, 
                     DateTimeArraySP        timelineMF,
                     IntArraySP             fwdCorrIndexArray,
                     DoubleMatrixArraySP    fwdCorrBetaFactorArray,
                     IntArraySP             indexArrayMF, 
                     LocalCorrSqueezeArray  localCorrSqueezeArray,
                     IntArray               localCorrSqueezeIndexArray) :
    nbIF(nbIF), nbTimeStepsIF(timelineIF->size()),
    randomGenIF(randomGenIF),
    nbMF(nbMF), nbTimeStepsMF(timelineMF->size()),    
    randomGenMF(randomGenMF), 
    fwdCorrIndexArray(fwdCorrIndexArray), 
    indexArrayMF(indexArrayMF),
    fwdCorrBetaFactorArray(fwdCorrBetaFactorArray),
    localCorrSqueezeArray(localCorrSqueezeArray),
    localCorrSqueezeIndexArray(localCorrSqueezeIndexArray) {
    static const string method("LocalCorr::LocalCorr");
    try {
        /** some validation */                
        if (nbTimeStepsMF < nbTimeStepsIF) {
            throw ModelException(method, "Mismatch between sim timeline and market factor timeline");
        }        
        
        /** save deltaT (need to create dummy time metric) */
        string name = "Holiday";
        const HolidaySP holidayTmp(new Holiday(name, DateTimeArray(0), true));
        TimeMetricSP dummyTimeMetric(new TimeMetric(0.0, holidayTmp.get()));
        
        deltaTimeIF = DoubleArraySP(new DoubleArray(nbTimeStepsIF));
        sqrtDeltaTimeIF = DoubleArraySP(new DoubleArray(nbTimeStepsIF));
        deltaTimeMF = DoubleArraySP(new DoubleArray(nbTimeStepsMF));
        sqrtDeltaTimeMF = DoubleArraySP(new DoubleArray(nbTimeStepsMF));

        /** first deltaT from valueDate to first market factor sim date */
        (*deltaTimeMF)[0] = dummyTimeMetric->yearFrac(*valueDate, timelineMF->front());
        (*sqrtDeltaTimeMF)[0] = sqrt((*deltaTimeMF)[0]);
        for (int iStep=1; iStep < nbTimeStepsMF; iStep++) {
            (*deltaTimeMF)[iStep] = 
                dummyTimeMetric->yearFrac((*timelineMF)[iStep-1], // start of interval
                                          (*timelineMF)[iStep]);  // end of interval
            (*sqrtDeltaTimeMF)[iStep] = sqrt((*deltaTimeMF)[iStep]);
        }

        (*deltaTimeIF)[0] = dummyTimeMetric->yearFrac(*valueDate, timelineIF->front());
        (*sqrtDeltaTimeIF)[0] = sqrt((*deltaTimeIF)[0]); 
        for (int iStep=1; iStep < nbTimeStepsIF; iStep++) {
            (*deltaTimeIF)[iStep] = 
                dummyTimeMetric->yearFrac((*timelineIF)[iStep-1], // start of interval
                                          (*timelineIF)[iStep]);  // end of interval
            (*sqrtDeltaTimeIF)[iStep] = sqrt((*deltaTimeIF)[iStep]);
        }        

        /** Y_t^1 and Y_t^2 ... */
        storedMF = CDoubleMatrixSP(new DoubleMatrix(nbMF, nbTimeStepsIF)); 
        /** int sq(Y_u) dY_u, int sq(Y_u) du, int sq^2(Y_u) du */
        nbRegions = localCorrSqueezeArray.size();
        storedIntegralsMF = DoubleMatrixArraySP(new DoubleMatrixArray(nbRegions));
        for (int iRegion=0; iRegion<nbRegions; iRegion++) {
            (*storedIntegralsMF)[iRegion] = CDoubleMatrixSP(new DoubleMatrix(3, nbTimeStepsIF));
        }                
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/** validation */
void LocalCorr::validatePop2Object() {}

void LocalCorr::correlateSeries( DoubleMatrix& noise, int pathIdx ) {
    /** noise[iAsset][jStep] is only placeholder (ie full of zeros)
        pathIdx does not need to be multiplied by 2, since we are NOT doing antithetics */
    static const string method("LocalCorr::correlateSeries");
    try{        
         if (!(pathIdx % 2)) { // we are on an odd path and ahve to update storedIntegralsMF
            
            randomGenMF->generate(pathIdx/2); // generate
            const DoubleMatrix& randomsMF = randomGenMF->getRandomNumbers(); // retrieve
            
            /** reset stored values */
            storedMF->fill(0.0); // DoubleMatrix(nbMF, nbTimeStepsIF)   
            for (int iRegion=0; iRegion<nbRegions; iRegion++) {
                (*storedIntegralsMF)[iRegion]->fill(0.0); // DoubleMatrix(3, nbTimeStepsIF)
            }

            double totalTimeSofar = 0.0;            
            double leftLimitMainMF = 0.0;
            int marketFactorIndex = 0;

            /** outer loop sim dates - inner loop market factor dates */
            for (int iSimDate=0; iSimDate<nbTimeStepsIF; iSimDate++) {
                
                while ( marketFactorIndex <= (*indexArrayMF)[iSimDate] ) {                    
                    
                    double thisDeltaT = (*deltaTimeMF)[marketFactorIndex]; // increment
                    double sqrtThisDeltaT = (*sqrtDeltaTimeMF)[marketFactorIndex];

                    /** first thing todo: determine squeeze(Y(t-)) ... ! */
                    double thisSqueeze = 0.0;
                    if (Maths::isPositive(totalTimeSofar)) {
                        double tmp = leftLimitMainMF / sqrt(totalTimeSofar);
                        
                        for (int iRegion=1; iRegion<nbRegions; iRegion++) {                            
                            thisSqueeze = 
                                localCorrSqueezeArray[iRegion]->computeSqueeze(tmp, totalTimeSofar);

                            (*(*storedIntegralsMF)[iRegion])[0][iSimDate] 
                                += thisSqueeze * randomsMF[0][marketFactorIndex] * sqrtThisDeltaT;
                            (*(*storedIntegralsMF)[iRegion])[1][iSimDate]
                                += thisSqueeze * thisDeltaT;
                            (*(*storedIntegralsMF)[iRegion])[2][iSimDate]
                                += thisSqueeze * thisSqueeze * thisDeltaT;     
                        }
                    } 

                    for (int iMF=0; iMF<nbMF; iMF++) {
                        (*storedMF)[iMF][iSimDate] += randomsMF[iMF][marketFactorIndex] * sqrtThisDeltaT; 
                    }
                    
                    leftLimitMainMF += randomsMF[0][marketFactorIndex] * sqrtThisDeltaT; 
                    marketFactorIndex += 1;                    
                    totalTimeSofar += thisDeltaT;
                }                
            }
        }

        /** outer look fwd corr dates - inner lool sim dates */
        int fwdCorrIndex = -1, marketFactorIndex = 0; // helper indices                

        randomGenIF->generate(pathIdx);  // generate 
        const DoubleMatrix& randomsIF = randomGenIF->getRandomNumbers(); // retrieve


        /** loop through fwd corr dates */
        for (int fwdCorrDate=0; fwdCorrDate<fwdCorrIndexArray->size(); fwdCorrDate++) {
            CDoubleMatrixConstSP thisBetaMatrix = (*fwdCorrBetaFactorArray)[fwdCorrDate];

            /** loop through sim dates (superset of fwd corr dates) */
            for (int iSimDate=fwdCorrIndex+1; iSimDate<=(*fwdCorrIndexArray)[fwdCorrDate]; iSimDate++) { 
                double thisDeltaT = (*deltaTimeIF)[iSimDate];
                double sqrtThisDeltaT = (*sqrtDeltaTimeIF)[iSimDate];

                if (Maths::isPositive(thisDeltaT)) {

                    const double* thisMainBetaArray = (*thisBetaMatrix)[0]; 
                                   
                    for (int iAsset=0; iAsset<nbIF; iAsset++) {

                        int thisIdx = localCorrSqueezeIndexArray[iAsset];

                        double intSqYdY = (*(*storedIntegralsMF)[thisIdx])[0][iSimDate];
                        double intSqYdt = (*(*storedIntegralsMF)[thisIdx])[1][iSimDate];
                        double intSq2Ydt = (*(*storedIntegralsMF)[thisIdx])[2][iSimDate];
                        
                        double mainBeta = thisMainBetaArray[iAsset];

                        double tmp = 1.0 * thisDeltaT 
                            - 2.0 * mainBeta * (1.0-mainBeta) * intSqYdt
                            - (1.0-mainBeta)*(1.0-mainBeta) * intSq2Ydt;

                        noise[iAsset][iSimDate] = (1.0-mainBeta) * intSqYdY;                    
                        
                        for (int iMF=0; iMF<nbMF; iMF++) {
                            double thisBeta = (*thisBetaMatrix)[iMF][iAsset];                        
                            tmp -= thisBeta*thisBeta * thisDeltaT;
                            noise[iAsset][iSimDate] += thisBeta * (*storedMF)[iMF][iSimDate]; 
                        }
                        
                        noise[iAsset][iSimDate] += randomsIF[iAsset][iSimDate] * sqrt(tmp);
                        noise[iAsset][iSimDate] /= sqrtThisDeltaT;
                        noise[iAsset][iSimDate] *= (-1.0);
                    }
                } else {
                    if (iSimDate > 0) {
                        for (int iAsset=0; iAsset<nbIF; iAsset++) {
                            noise[iAsset][iSimDate] = noise[iAsset][iSimDate-1];
                        }
                    } else {
                        for (int iAsset=0; iAsset<nbIF; iAsset++) {
                            noise[iAsset][iSimDate] = 0.0;
                        }
                    }
                }
                fwdCorrIndex += 1;
            }
        }

        /*if (pathIdx<5000) {
            for (int iAsset=0; iAsset<nbIF; iAsset++) {
                (*storeCorrRN)[iAsset][pathIdx] = noise[iAsset][nbTimeStepsIF-1];
            }
            (*storeCorrRN)[nbIF][pathIdx] = (*storedMF)[0][nbTimeStepsIF-1] / sqrt((*deltaTimeIF)[nbTimeStepsIF-1]);
        }*/
        
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/**  see description in pure virtual declaration */
DoubleMatrix LocalCorr::getCorrelations(int index){
    return DoubleMatrix(0,0); // TODO    
}
double LocalCorr::getCorrelationIJ(int i, int j, int index){
    return 0.0; // TODO	
}

/*******************************/
/* LOCAL CORR DEPENDENCE MAKER */
/*******************************/

const string DependenceMakerLocalCorr::DEFAULT_MARKET_FACTOR_FREQ = "1W";
const int DependenceMakerLocalCorr::DEFAULT_NB_FACTORS = 1;
const double DependenceMakerLocalCorr::DEFAULT_PRECISION = 1E-04;
const int DependenceMakerLocalCorr::DEFAULT_MAX_ITER = 10;

void DependenceMakerLocalCorr::validatePop2Object() {
    static const string method("DependenceMakerLocalCorr::validatePop2Object");
    try {
        if (nbFactors==0) {
            throw ModelException(method, "Nb of factors used for beta factorization must at least be one, "
                "but is zero.");
        }
        if (Maths::isNegative(precision)) {
            throw ModelException(method, "Precision input for beta factorization must not be negative, " 
                "but equals " + Format::toString(precision) + ".");
        }
        if (maxIter < 10) {
            throw ModelException(method, "Max nb of iterations for beta factoriazion mut be at least 10, "
                "but equals " + Format::toString(maxIter) + ".");
        }
    } catch (exception& e){
        throw ModelException(e,method);
    } 
}

/** modifies MDF if necessary */
void DependenceMakerLocalCorr::modifyMarketDataFetcher(MarketDataFetcherSP mdf) {    
    MDFUtil::setCorrSkewMode(*mdf, true);
}

DependenceSP DependenceMakerLocalCorr::createDependence(const DependenceMaker::ISupport* support) const {
    static const string method("DependenceMakerLocalCorr::createDependence");
    try {
        /** support used for validation and to get specific data */
        const DependenceMakerLocalCorr::Support* supportDmlc = 
            dynamic_cast<const DependenceMakerLocalCorr::Support*>(support);
        if (!supportDmlc) {
            throw ModelException(method, "DependenceMakerLocalCorr is not a valid dependence maker.");
        }

        DependenceMakerSP dependenceMaker = supportDmlc->getDependenceMaker();
      
        /** retrieve mAsset and check whether only one asset */
        const IMultiFactors* mAsset = supportDmlc->getMultiFactors();
        int nbIF = mAsset->NbAssets();
        if (nbIF == 1) { // no CorrTS necessary, no CorrSkew necessary
            DoubleMatrix correlation(1,1);
            correlation[0][0] = 1.0;
            return DependenceSP(new Gauss(correlation)); // Gauss as good as anything else 
        }

        /** retrieve timeline & valuedate & simdates */
        DateTimeArray timeline = supportDmlc->getSimDates(); 
        DateTime valueDate = timeline.front();
        DateTimeArray timelineIF(timeline.begin()+1, timeline.end()); // sim dates
        
        /** GAUSS or GAUSSTERM */
        IntArray            fwdCorrIndexArray;        
        DoubleMatrixArray   fwdCorrBetaFactorArray;

        /** squeeze for input correlation matrix */
        LocalCorrSqueezeArray localCorrSqueezeArray = mAsset->getLocalCorrSqueezeArray();        
        
        /** GaussTerm ? */ 
        const DependenceMakerGaussTerm* dpmGT = 
            dynamic_cast<const DependenceMakerGaussTerm*>(dependenceMaker.get());
        if (!dpmGT) {
            /** Gauss ? */
            const DependenceMakerGauss* dpmG =
                dynamic_cast<const DependenceMakerGauss*>(dependenceMaker.get());
            if (!dpmG) {
                throw ModelException(method, "Wrong DependenceMaker supplied to DependenceMakerLocalCorr");
            }

            /** GAUSS  */           
           CDoubleMatrixSP correlationMatrix = 
               CDoubleMatrixSP(new DoubleMatrix(*mAsset->factorsCorrelationMatrix()));            
            
            BetaFactorizationSP betaFactorization = 
                correlationMatrix->computeBetaFactorization(nbFactors, precision, maxIter);

            double allowedMax = betaFactorization->squeezeLimitHigh;
            double allowedMin = betaFactorization->squeezeLimitLow;
            for (int iRegion=0; iRegion<localCorrSqueezeArray.size(); iRegion++) {
                double usedMax = localCorrSqueezeArray[iRegion]->getMaxSqueeze();
                if ( Maths::isPositive(usedMax - allowedMax) ) {
                    throw ModelException(method, "The max allowed squeeze equals " + 
                        Format::toString(allowedMax) + 
                        " but the max used squeeze equals " +
                        Format::toString(usedMax) + ".\n This is the case for region" + 
                            localCorrSqueezeArray[iRegion]->getName() + ".");
                }
                double usedMin = localCorrSqueezeArray[iRegion]->getMinSqueeze(); 
                if ( Maths::isPositive(allowedMin - usedMin) ) {
                    throw ModelException(method, "The min allowed squeeze equals " + 
                        Format::toString(allowedMin) + 
                        " but the min used squeeze equals " +
                        Format::toString(usedMin) + ".\n This is the case for region" + 
                            localCorrSqueezeArray[iRegion]->getName() + ".");
                }
            }
            /** a bit clumsy ... */
            fwdCorrIndexArray.resize(timelineIF.size());
            fwdCorrBetaFactorArray.resize(timelineIF.size());
            for (int iStep=0; iStep<timelineIF.size(); iStep++) {
                fwdCorrIndexArray[iStep] = iStep;                    
                fwdCorrBetaFactorArray[iStep] = betaFactorization->betaFactors;
            }            
        } else {
            dpmGT->checkInputs(); // a bit clumsy ...            

            /** retrieve correlation market data */
            CorrTermDataSP corrTermData = 
                CorrelationTerm::getCorrelationTermSqueezesAndExpiries (
                    valueDate, 
                    nbIF,
                    mAsset->getCorrObjArray(), 
                    mAsset->getCorrTermObjArray());                 

            DateTimeArray fwdCorrDatesArray;
            /** create arrays of dates and index array */
            DependenceMakerGaussTerm::createFwdCorrIdxAndDatesArray(
                timeline, // valuedate + simdates
                dpmGT->getFwdCorrModeUsed(),
                mAsset->getTimeMetricArray(),
                fwdCorrIndexArray,     // output
                fwdCorrDatesArray);        
            
            /** compute fwd correls at fwd dates */
            DoubleMatrixArraySP fwdCorrTermMatrixArray = 
                dpmGT->computeFwdCorrMatrixArray(
                    fwdCorrIndexArray, 
                    fwdCorrDatesArray,
                    corrTermData,
                    mAsset,
                    supportDmlc->getFwdVarAtDates(dpmGT->getFwdCorrModeUsed()->doInterpolateAtmFwd()), 
                    timeline); // timeline (valuedate + simdates)

            fwdCorrBetaFactorArray.resize(fwdCorrIndexArray.size());
                                 
            for (int i=0; i<fwdCorrIndexArray.size(); i++) {                
                BetaFactorizationSP betaFactorization = 
                    (*fwdCorrTermMatrixArray)[i]->
                        computeBetaFactorization(nbFactors, precision, maxIter);
                
                double allowedMax = betaFactorization->squeezeLimitHigh;
                double allowedMin = betaFactorization->squeezeLimitLow;
                for (int iRegion=0; iRegion<localCorrSqueezeArray.size(); iRegion++) {
                    double usedMax = localCorrSqueezeArray[iRegion]->getMaxSqueeze();
                    if ( Maths::isPositive(usedMax - allowedMax) ) {
                        throw ModelException(method, "The max allowed squeeze equals " + 
                            Format::toString(allowedMax) + 
                            " but the max used squeeze equals " +
                            Format::toString(usedMax) + ".\n This is the case for region" + 
                            localCorrSqueezeArray[iRegion]->getName() + ".");
                    }
                    double usedMin = localCorrSqueezeArray[iRegion]->getMinSqueeze(); 
                    if ( Maths::isPositive(allowedMin - usedMin) ) {
                        throw ModelException(method, "The min allowed squeeze equals " + 
                            Format::toString(allowedMin) + 
                            " but the min used squeeze equals " +
                            Format::toString(usedMin) + ".\n This is the case for region" + 
                            localCorrSqueezeArray[iRegion]->getName() + ".");
                    }
                }
                fwdCorrBetaFactorArray[i] = betaFactorization->betaFactors;
            }

        } /** we now have fwdCorrIndexArray and corresponding fwdCorrBetaFactorArray */
        
        /** need to check whether we need to insert additional dates for market factor */
        MaturityPeriodSP marketFactorInterval(new MaturityPeriod(samplingFrequency));
        
        int basicCount;
        string basicInterval;
        marketFactorInterval->decompose(basicCount, basicInterval);        

        /** timeline for the MARKET FACTOR */
        DateTimeArraySP timelineMF(new DateTimeArray(1, timelineIF.back()));         
        while(timelineMF->front() > valueDate) {
            DateTime dateToInsert = 
                marketFactorInterval->toDate(-basicCount, basicInterval, timelineMF->front());
            int indexLow = dateToInsert.findLower(timelineIF) + 1; // idx of first date greater than this one
            if (timelineIF[indexLow] < timelineMF->front()) {
                // want to insert this simdate and all inbetween
                int indexUp = timelineMF->front().findLower(timelineIF);
                if (timelineIF[indexUp].equals(timelineMF->front())) {
                    indexUp--; // to avoid duplicates
                }
                for (int iIdx=indexUp; iIdx>=indexLow; iIdx--) {
                    timelineMF->insert(timelineMF->begin(), timelineIF[iIdx]); 
                }               
            } else {
                timelineMF->insert(timelineMF->begin(), dateToInsert);
            }
        }
        /** by construction: timelineMF[0] <= valueDate, hence, remove first element again */
        timelineMF->erase(timelineMF->begin());

        vector<int> indexArrayTmp = DateTime::getIndexes(*timelineMF, timelineIF);
        IntArraySP indexArrayMF(new IntArray(indexArrayTmp.begin(),indexArrayTmp.end()));

        /** need to create 2 MCRandoms, one for the idiosyncratic / market factor */
        CDoubleMatrix helperMatrixIF(nbIF,nbIF);
        for (int iIF=0; iIF<nbIF; iIF++) {
            helperMatrixIF[iIF][iIF] = 1.0;
        }
        IMCRandomSP randomGenIF = IMCRandomSP(new
            MCRandom(0, // path generator
                     DependenceSP(new Gauss(helperMatrixIF)),
                     supportDmlc->getRandomGenerator(),         // IRandomSP
                     supportDmlc->getIdiosynFactorRandomCache(),// FactorRandomCache
                     supportDmlc->carefulRandoms(),             // isCarefulRandoms
                     timelineIF.size(),                         // numDates
                     nbIF,                                      // numFactors
                     supportDmlc->getNbPastDates()));           // numPastDates

        CDoubleMatrix helperMatrixMF(nbFactors,nbFactors);
        for (int iMF=0; iMF<nbFactors; iMF++) {
            helperMatrixMF[iMF][iMF] = 1.0;
        }
        IMCRandomSP randomGenMF = IMCRandomSP(new 
            MCRandom(0, // path generator 
                     DependenceSP(new Gauss(helperMatrixMF)),
                     supportDmlc->getRandomGenerator(),         // IRandomSP
                     supportDmlc->getMarketFactorRandomCache(), // FactorRandomCache
                     false,                                     // isCarefulRandoms
                     timelineMF->size(),                        // numDates
                     nbFactors,                                 // numFactors
                     0,                                         // numPastDates
                     false));                                   // doAntithetics

        return DependenceSP(new LocalCorr(
            DateTimeSP(new DateTime(valueDate)),
            nbIF,
            randomGenIF,
            DateTimeArraySP(new DateTimeArray(timelineIF)), 
            nbFactors,
            randomGenMF,
            timelineMF,
            IntArraySP(new IntArray(fwdCorrIndexArray)),      // some info for correl TS
            DoubleMatrixArraySP(new DoubleMatrixArray(fwdCorrBetaFactorArray)),
            indexArrayMF, // some info for market factor
            localCorrSqueezeArray,
            mAsset->getLocalCorrSqueezeIndexArray()));            
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// Invoked when Class is 'loaded'
void DependenceMakerLocalCorr::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DependenceMakerLocalCorr, clazz);        
    SUPERCLASS(SkewMaker);    
    EMPTY_SHELL_METHOD(defaultDependenceMakerLocalCorr);
    FIELD(samplingFrequency, "sampling frequency of market factor");        
    FIELD(nbFactors, "nb of factors for beta factorization");
    FIELD(precision, "precision for beta factorization");
    FIELD(maxIter, "maxIter for beta factorization");
}

IObject* DependenceMakerLocalCorr::defaultDependenceMakerLocalCorr(){
    return new DependenceMakerLocalCorr();
}

DependenceMakerLocalCorr::DependenceMakerLocalCorr() : 
SkewMaker(TYPE) {}

DependenceMakerLocalCorr::DependenceMakerLocalCorr(const CClassConstSP& clazz) :
SkewMaker(clazz) {}

CClassConstSP const DependenceMakerLocalCorr::TYPE = CClass::registerClassLoadMethod(
	"DependenceMakerLocalCorr", typeid(DependenceMakerLocalCorr), DependenceMakerLocalCorr::load);

DRLIB_END_NAMESPACE

