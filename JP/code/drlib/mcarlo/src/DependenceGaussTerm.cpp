//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceGaussTerm.hpp
//
//   Description : Holds GaussTerm (SRM) Dependence (Maker)
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_DEPENDENCE_GAUSS_TERM_CPP
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
#include "edginc/MDFUtil.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////////////////////
// HELPER METHOD -- TO BE RETIRED ONE DAY ... 
///////////////////////////////////////////////////////////////////////////////////////////////

DoubleArray decomposeFwdCorrModeString(const string& fwdCorrMode) {
    static const string method("CorrelationTerm::decomposeFwdCorrModeString");
    try {
    DoubleArray fwdCorrModeBarriers;
    const char* position = fwdCorrMode.c_str();
    while (*position != '\0') { // that indicated the end of the string
        char   c;
        // skip whitespace
        while ((c = *position) == ' ' || c == '\t'){
            position++;
        }
        if (c == '-' || c == '.' || isdigit (c)) {
            char* tmpPos;
            double value = strtod(position, &tmpPos);
            position = tmpPos;
            // skip whitespace
            char cc;
            while ((cc = *position) == ' ' || cc == '\t'){
                position++;
            }
            fwdCorrModeBarriers.push_back(value);
        } else if (c == ','){
            // skip over our chosen separator and continue
            position++;
        } else {
            // not recognised
            throw ModelException(method, "Failed to parse fwdCorrMode " + fwdCorrMode 
                + " as comma separated string with two doubles");
        }
    }
    if (!(fwdCorrModeBarriers.size() == 2)) {
        throw ModelException(method, "Wrong input for fwdCorrMode in DependenceMaker.\n" 
            "Two strictly increasing doubles between zero and one, separated by a comma, are needed");
    }
    if (Maths::isNegative(fwdCorrModeBarriers[0]) || 
        Maths::isNegative(fwdCorrModeBarriers[1]) ||
        Maths::isNegative(fwdCorrModeBarriers[1] - fwdCorrModeBarriers[0]) ||
        Maths::isNegative(1.0 - fwdCorrModeBarriers[1]) )  {
        throw ModelException(method, "Wrong input for fwdCorrMode in DependenceMaker.\n"
              "Both values need to be between zero and one and \n"
              "the second value must be strictly greater than the first value.\n"
              "However, the inputs are "
            + Format::toString(fwdCorrModeBarriers[0]) + " and " 
            + Format::toString(fwdCorrModeBarriers[1]));
    }
    return fwdCorrModeBarriers;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////
// GAUSS TERM
///////////////////////////////////////////////////////////////////////////////////////////////

/** Destructor */
GaussTerm::~GaussTerm(){ /*empty*/ }

/** Validation */
void GaussTerm::validatePop2Object(){  /*empty*/ }

/** Default Constructor */
GaussTerm::GaussTerm(){}

/** constructor from correlation matrices */
GaussTerm::GaussTerm(const DoubleMatrixArraySP& correlationTerm,
                     const IntArray&            fwdCorrIndexArray): 
    correlationTerm(correlationTerm), fwdCorrIndexArray(fwdCorrIndexArray) {
    static const string method("GaussTerm::GaussTerm");
    try{
        // dimensions...
        nbFwdCorrSteps  = correlationTerm->size();
        nbAssets        = (*correlationTerm)[0]->numCols();
        int iStep;
        for( iStep=0; iStep<nbFwdCorrSteps; iStep++ ) {
            if( (*correlationTerm)[0]->numCols() != nbAssets ) {
                throw ModelException(method, 
                                     "In step " + Format::toString(iStep) + 
                                     "have numCols " + Format::toString((*correlationTerm)[0]->numCols()) + 
                                     ", expected " + Format::toString(nbAssets) );
            }
            else if( (*correlationTerm)[0]->numRows() != nbAssets ) {
                throw ModelException(method, 
                                     "In step " + Format::toString(iStep) + 
                                     "have numRows " + Format::toString((*correlationTerm)[0]->numRows()) + 
                                     ", expected " + Format::toString(nbAssets) );
            }
        }

        // allocations for speed in correlateSeries
        correlCoeffsTerm.resize(nbAssets);
        for (int a = 0; a < nbAssets; ++a) {
            correlCoeffsTerm[a].resize(nbAssets);
            for (int b = 0; b < nbAssets; ++b) {
                correlCoeffsTerm[a][b].resize(nbFwdCorrSteps);
            }
        }
        try {
            for( iStep=0; iStep<nbFwdCorrSteps; iStep++ ) {
                DoubleMatrix corrTerm = (*correlationTerm)[iStep]->computeSquareRoot();
                for (int j = 0; j < nbAssets; ++j) {
                    vector<double>* correlCoeffsTermJ = &correlCoeffsTerm[j][0];
                    for (int i = 0; i < nbAssets; ++i) {
                        correlCoeffsTermJ[i][iStep] = corrTerm[j][i];
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e,method, "Failed to compute square root of correlation matrix");
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }

}

void GaussTerm::correlateSeries( DoubleMatrix& noise, int pathIdx ) {
    static const string method("Gauss::correlateSeries");
    try{
        // pathIdx not needed here
        int index, iIndex, iAsset, jAsset, iDate;
        for (iAsset = nbAssets - 1; iAsset > 0; iAsset --) {
            double* mainSeries = noise[iAsset];
            const vector<double>* correlCoeffsI = &correlCoeffsTerm[iAsset][0];
            const double* correlCoeffsII = &correlCoeffsI[iAsset][0];
            index = -1;
            for (iDate = 0; iDate < nbFwdCorrSteps; iDate ++){
                int tmpInt = fwdCorrIndexArray[iDate];
                double tmpCorrelCoeffsII = correlCoeffsII[iDate];
                for (iIndex = index+1; iIndex <= tmpInt; iIndex++) {
                    mainSeries[iIndex] *= tmpCorrelCoeffsII;
                }
                index = fwdCorrIndexArray[iDate];
            }
            for (jAsset = iAsset - 1; jAsset > -1; jAsset --) {
                double* otherSeries = noise[jAsset];
                const double* correlCoeffsIJ = &correlCoeffsI[jAsset][0];
                index = fwdCorrIndexArray.back()+1;
                for (iDate = nbFwdCorrSteps - 1; iDate > 0; --iDate) {
                    int tmpInt = fwdCorrIndexArray[iDate-1];
                    double tmpCorrelCoeffsIJ = correlCoeffsIJ[iDate];
                    for (iIndex = index-1; iIndex > tmpInt; iIndex--) {
                        mainSeries[iIndex] += otherSeries[iIndex] * tmpCorrelCoeffsIJ;
                    }
                    index = iIndex+1;
                }
                // special for iDate = 0
                for (iIndex = index-1; iIndex >=0; iIndex--) {
                    mainSeries[iIndex] += otherSeries[iIndex] * correlCoeffsIJ[0];
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

//// see desription in pure virtual declaration
DoubleMatrix GaussTerm::getCorrelations(int index){
    return *(*correlationTerm)[index];
}
double GaussTerm::getCorrelationIJ(int i, int j, int index){
	return (*(*correlationTerm)[index])[i][j];
}

/*****************************************/
/* EXTENDED SPARSE GAUSS TERM DEPENDENCE */
/*****************************************/

/** Same as Gauss but can cope for case where the 'square root' of the
correlation matrix is supplied (and is not in lower triangular form) */

/** constructor */
ExtendedSparseGaussTerm::ExtendedSparseGaussTerm(
    const vector<SparseDoubleMatrixSP>& sqrtCorrs,
    IntArray                            fwdCorrIndexArrayInput) : 
    correlCoeffsSparse(sqrtCorrs) {
        /** the fields correlationTerm and correlCoeffsTerm are not populated */
	nbAssets = sqrtCorrs[0]->numRows();
    nbFwdCorrSteps = sqrtCorrs.size();
    fwdCorrIndexArray = fwdCorrIndexArrayInput;
}

DoubleMatrix ExtendedSparseGaussTerm::getCorrelations(int index) {
    static const string method("ExtendedSparseGaussTerm::getCorrelations");
	throw ModelException(method,
        "ExtendedSparseGaussTerm, does not store the whole correlation matrix. " 
		"Use element-wise getCorrelationIJ instead.");
}

double ExtendedSparseGaussTerm::getCorrelationIJ(int i, int j, int index){
	return correlCoeffsSparse[(unsigned int)index]->GetMMTElem(i,j);
}

void ExtendedSparseGaussTerm::correlateSeries( DoubleMatrix& noise, int pathIdx) {        
    // pathIdx not needed here
    int index = -1;
    for (int idx = 0; idx < fwdCorrIndexArray.size(); idx++) {
        for (int iStep = index + 1; iStep <= fwdCorrIndexArray[idx]; iStep++) {
            correlCoeffsSparse[idx]->leftMultiplyToRow(noise, iStep); 
        }
        index = fwdCorrIndexArray[idx];
    }
}

/*******************************/
/* GAUSS TERM DEPENDENCE MAKER */
/*******************************/

void DependenceMakerGaussTerm::validatePop2Object() {
    static const string method("DependenceMakerGaussTerm::validatePop2Object");
    try {
        init(); 
    } catch (exception& e){
        throw ModelException(e,method);
    } 
}

/** modifies MDF if necessary */
void DependenceMakerGaussTerm::modifyMarketDataFetcher(MarketDataFetcherSP mdf) {
    MDFUtil::setCorrTermMode(*mdf, true);
}


void DependenceMakerGaussTerm::init() {
    /** populate fwdCorrModeUsed */
    /** either fwdCorrMode = NOTUSED or it has to be comma separated string */
    if (!CString::equalsIgnoreCase(fwdCorrMode, "NOTUSED")) {
        DoubleArray barriers = decomposeFwdCorrModeString(fwdCorrMode);
        fwdCorrModeUsed = FwdCorrModeSP(new FwdCorrModeSimdates(barriers[0], barriers[1]));
    } else {
        fwdCorrModeUsed = fwdCorrModeDetails->getFwdCorrMode();
    }
    /** validation - both obj might be supplied, only one is used -> validate only the used one */
    fwdCorrModeUsed->checkInputs();
}

void DependenceMakerGaussTerm::checkInputs() const {
    static const string method("DependenceMakerGaussTerm::checkInputs");
    try {
        /** some validation */
        if (!Maths::isPositive(eigenValueFloor)) {
            throw ModelException(method, "EigenValueFloor has to be positive, but is " + 
                Format::toString(eigenValueFloor));
        }
        if (!Maths::isPositive(maxSqError) || Maths::isPositive(maxSqError - 1.0)) {
            throw ModelException(method, "MaxSqError has to be between zero and one, but is " +
                Format::toString(maxSqError));
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** a static helper method */
void DependenceMakerGaussTerm::createFwdCorrIdxAndDatesArray(
    const DateTimeArray&           timeline, // valuedate + simdates
    const FwdCorrModeSP            fwdCorrMode,    
    const TimeMetricArray&         timeMetricArray,
    IntArray&                      fwdCorrIndexArray,   // output
    DateTimeArray&                 fwdCorrDatesArray) { // output
    static const string method("DependenceMakerGaussTerm::createAllFwdCorrData");
    try {
        /** timeline & valuedate & simdates */
        DateTime valueDate = timeline.front();
        DateTimeArray simDates(timeline.begin()+1, timeline.end()); // first date after value date        
                
        const FwdCorrModeStandard* fwdCorrModeStandard = 
            dynamic_cast<const FwdCorrModeStandard*>(fwdCorrMode.get());
        if(!fwdCorrModeStandard) { // must be fwdCorrModeSimdates now ... 
            const FwdCorrModeSimdates* fwdCorrModeSimdates = 
                dynamic_cast<const FwdCorrModeSimdates*>(fwdCorrMode.get());
            if (!fwdCorrModeSimdates) {
                throw ModelException(method, "Wrong input for supplied" );
            }        
            /** CASE A: fwd corr dates equal sim dates */            
            fwdCorrIndexArray.resize(simDates.size());
            for (int iStep=0; iStep<simDates.size(); iStep++) {
                fwdCorrIndexArray[iStep] = iStep;
            }
            fwdCorrDatesArray = simDates;            
        } else {
            /** CASE B: fwd corr dates according to tenor */                        

            /** decompose fwdCorrInterval */
            int basicCount;
            string basicInterval;
            MaturityPeriodSP fwdCorrInterval(new MaturityPeriod(fwdCorrModeStandard->getFwdCorrInterval()));
            fwdCorrInterval->decompose(basicCount, basicInterval);
        
            /** fwdCorrDatesArray & fwdCorrIndexArray */
            fwdCorrDatesArray = DateTimeArray(1, simDates.back()); 
            fwdCorrIndexArray = IntArray(1, simDates.size()-1);

            int count = -1, index, nbBusDaysDiff;
            while (fwdCorrDatesArray.front() > valueDate) {
                DateTime insertDate = 
                    fwdCorrInterval->toDate(count * basicCount, basicInterval, simDates.back());
                /** simdates does NOT incl valuedate, hence, we are on the safe side */
                index = insertDate.findLower(simDates); 
                nbBusDaysDiff = fwdCorrModeStandard->getNbStubDays();
                /* calc nbBusDaysDiff between insertDate & fwdCorrDatesArray.front() */                
                for (int iAsset=0; iAsset<timeMetricArray.size(); iAsset++) {
                    int tmp = timeMetricArray[iAsset]->getHolidays()->
                        businessDaysDiff(insertDate, fwdCorrDatesArray.front());
                    nbBusDaysDiff = Maths::min(nbBusDaysDiff,tmp);
                }				
                if (index < 0) {
                    break; 
                } else if ( (index < fwdCorrIndexArray[0]) && 
				    (nbBusDaysDiff >=fwdCorrModeStandard->getNbStubDays()) ) {
                    fwdCorrDatesArray.insert(fwdCorrDatesArray.begin(), simDates[index]);
                    fwdCorrIndexArray.insert(fwdCorrIndexArray.begin(), index);
                }
                    count --;
            }
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

CDoubleMatrixSP DependenceMakerGaussTerm::fwdVarAtFwdCorrDates(
    const DateTimeArray&            timeline,
    const IntArray&                 fwdCorrIndexArray,
    const DoubleMatrix&             fwdVarAtSimDates) {
    static const string method("DependenceMakerGaussTerm::fwdVarAtFwdCorrDates");
    try {
        /** check whether we need to do sthg about fwd vars */        
        if (fwdCorrIndexArray.size() == timeline.size()-1) {
            return CDoubleMatrixSP(new DoubleMatrix(fwdVarAtSimDates)); // nthg todo
        } else {
            /** compute fwdVarAtFwdCorrDates from fwdVarAtSimdates */
            CDoubleMatrixSP fwdVarAtFwdCorrDates(new 
                DoubleMatrix(fwdVarAtSimDates.numCols(), fwdCorrIndexArray.size()));
            int index;
            for (int iAsset=0; iAsset<fwdVarAtSimDates.numCols(); iAsset++) {
                index = -1;
                for (int iDate=0; iDate<fwdCorrIndexArray.size(); iDate++) {
                    double fwdVar = 0.0;
                    for (int iIndex = index+1; iIndex <= fwdCorrIndexArray[iDate]; iIndex++) {
                        fwdVar += fwdVarAtSimDates[iAsset][iIndex];
                    }
                    (*fwdVarAtFwdCorrDates)[iAsset][iDate] = fwdVar;
                    index = fwdCorrIndexArray[iDate];
                }
            }    
            return fwdVarAtFwdCorrDates;
        }      
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** a non-static helper method */
DoubleMatrixArraySP DependenceMakerGaussTerm::computeFwdCorrMatrixArray (
    IntArray&               fwdCorrIndexArray,
    DateTimeArray&          fwdCorrDatesArray,    
    CorrTermDataSP          corrTermData,
    const IMultiFactors*    mAsset,
    const DoubleMatrix&     fwdVarAtSimDates,
    const DateTimeArray&    timeline) const {
    static const string method("DependenceMakerGaussTerm::computeFwdCorrMatrixArray");
    try { 
        CDoubleMatrixSP fwdVarAtFwdCorrDates = 
            DependenceMakerGaussTerm::fwdVarAtFwdCorrDates(timeline,
                                                           fwdCorrIndexArray,
                                                           fwdVarAtSimDates);
          
        /** compute fwdCorrelation matrices */   
        return CorrelationTerm::CorrelationTermMatrix(
            timeline.front(), 
            fwdCorrDatesArray, 
            *fwdVarAtFwdCorrDates, 
            mAsset->getTimeMetricArray(), 
            corrTermData,
            CorrelationTerm::forward, 
            eigenValueFloor, 
            maxSqError, 
            fwdCorrModeUsed->getBarriers(), 
            mAsset->getSkipFwdCorrArray());        
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

DependenceSP DependenceMakerGaussTerm::createDependence(const DependenceMaker::ISupport* support) const {
    static const string method("DependenceMakerGaussTerm::createDependence");
    try {    
        /** some validation ... **/
        checkInputs(); // a bit late, should be called from pathConfig -> todo
        const DependenceMakerGaussTerm::Support* supportDmgt = 
            dynamic_cast<const DependenceMakerGaussTerm::Support*>(support);
        if (!supportDmgt) {
            throw ModelException(method, "DependenceMakerGaussTerm is not a valid dependence maker.");
        }

    /** if only one asset -> create GAUSS dependence */        
        int nbAssets = supportDmgt->getMultiFactors()->NbAssets(); 
        if (nbAssets == 1) { // nthg todo really
            DoubleMatrix correlation(1,1);
            correlation[0][0] = 1.0;
            return DependenceSP(new Gauss(correlation)); // Gauss as good as GaussTerm
        }

        DateTimeArray timeline = supportDmgt->getSimDates();
        const IMultiFactors* mAsset = supportDmgt->getMultiFactors();

        /** retrieve correlation corrTermData */
        CorrTermDataSP corrTermData = 
            CorrelationTerm::getCorrelationTermSqueezesAndExpiries (
                timeline.front(), 
                nbAssets,
                mAsset->getCorrObjArray(), 
                mAsset->getCorrTermObjArray());                 

        IntArray            fwdCorrIndexArray;        
        DateTimeArray       fwdCorrDatesArray;      

        DependenceMakerGaussTerm::createFwdCorrIdxAndDatesArray(
            timeline, // valuedate + simdates
            fwdCorrModeUsed,
            mAsset->getTimeMetricArray(),
            fwdCorrIndexArray,     // output
            fwdCorrDatesArray);         
        
        DoubleMatrixArraySP fwdCorrTermMatrixArray = 
            computeFwdCorrMatrixArray(fwdCorrIndexArray, 
                                      fwdCorrDatesArray,
                                      corrTermData,
                                      mAsset,
                                      supportDmgt->getFwdVarAtDates(fwdCorrModeUsed->doInterpolateAtmFwd()), 
                                      timeline); // timeline (valuedate + simdates)


        return DependenceSP(new GaussTerm(fwdCorrTermMatrixArray, fwdCorrIndexArray));
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// Invoked when Class is 'loaded'
void DependenceMakerGaussTerm::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DependenceMakerGaussTerm, clazz);
    SUPERCLASS(DependenceMaker);
    IMPLEMENTS(IRevertTypeConvert);
    EMPTY_SHELL_METHOD(defaultDependenceMakerGaussTerm);
    FIELD(eigenValueFloor, "eigen value floor for pos definite modification");
    FIELD_MAKE_OPTIONAL(eigenValueFloor);
    FIELD(maxSqError, "max sq error for deviation from pos definite modification");
    FIELD_MAKE_OPTIONAL(maxSqError);
    FIELD(fwdCorrMode, "mode for checking fwd corr matrices");
    FIELD_MAKE_OPTIONAL(fwdCorrMode);
    FIELD(fwdCorrModeDetails, "used fwd corr mode");
    FIELD_MAKE_OPTIONAL(fwdCorrModeDetails);
    FIELD(fwdCorrModeUsed, "used fwd corr mode");
    FIELD_MAKE_TRANSIENT(fwdCorrModeUsed);
}

IObject* DependenceMakerGaussTerm::defaultDependenceMakerGaussTerm(){
    return new DependenceMakerGaussTerm();
}


// for the IMS interface (see IRevertTypeConvert)
IObjectSP DependenceMakerGaussTerm::revert(const string& interfaceType) const {
    static const string method = "DependenceMakerGaussTerm::revert";
    try {
        if (interfaceType != IRevertTypeConvert::PYRAMID) {
            throw ModelException( method, 
                                  "Cannot convert a DependenceMakerGaussTerm for"
                                  " the interface " + interfaceType);               
        }
        
        IObjectSP wrapper(new DependenceMakerWrapper(DependenceMakerGaussTermSP(copy(this))));
        return wrapper;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

FwdCorrModeSP DependenceMakerGaussTerm::getFwdCorrModeUsed() const {
    static const string method = "DependenceMakerGaussTerm::getFwdCorrModeUsed";
    if (!fwdCorrModeUsed) {
        throw ModelException(method, "Internal error");
    }
    return fwdCorrModeUsed;
}

DependenceMakerGaussTerm::DependenceMakerGaussTerm(): DependenceMaker(TYPE), 
eigenValueFloor(CorrelationTerm::EIGEN_VALUE_FLOOR), 
maxSqError(CorrelationTerm::MAX_SQ_ERROR),
fwdCorrMode("NOTUSED"), 
fwdCorrModeDetails(new FwdCorrModeWrapper()) {
    init(); 
}

CClassConstSP const DependenceMakerGaussTerm::TYPE = CClass::registerClassLoadMethod(
	"DependenceMakerGaussTerm", typeid(DependenceMakerGaussTerm), DependenceMakerGaussTerm::load);

/***********************************/
/* GAUSS TERM SRM DEPENDENCE MAKER */
/***********************************/

DependenceMakerGaussTermSrm::DependenceMakerGaussTermSrm() : DependenceMakerGaussSrm(TYPE) {
    corrTermStructure = true;
}

void DependenceMakerGaussTermSrm::load(CClassSP& clazz){
    clazz->setPublic(); 
    REGISTER(DependenceMakerGaussTermSrm, clazz);
    SUPERCLASS(DependenceMakerGaussSrm);
    EMPTY_SHELL_METHOD(defaultDependenceMakerGaussTermSrm);    
}

IObject* DependenceMakerGaussTermSrm::defaultDependenceMakerGaussTermSrm(){
    return new DependenceMakerGaussTermSrm();
}

CClassConstSP const DependenceMakerGaussTermSrm::TYPE = CClass::registerClassLoadMethod(
    "DependenceMakerGaussTermSrm", typeid(DependenceMakerGaussTermSrm), DependenceMakerGaussTermSrm::load);


DRLIB_END_NAMESPACE

