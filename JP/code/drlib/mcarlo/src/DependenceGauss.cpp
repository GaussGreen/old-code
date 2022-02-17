//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceGauss.hpp
//
//   Description : Holds GaussTerm (SRM) Dependence (Maker)
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_DEPENDENCE_GAUSS_CPP
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
#include "edginc/MDFUtil.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************************************/
/* MODES FOR FWDCORRELATION						            */
/* FwdCorrMode -> FwdCorrModeStandard & FwdCorrModeSimdates */
/************************************************************/

void FwdCorrMode::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FwdCorrMode, clazz);
    SUPERCLASS(CObject);
}

/** for reflection */
FwdCorrMode::FwdCorrMode() : CObject(TYPE) {}
FwdCorrMode::FwdCorrMode(const CClassConstSP& clazz) : CObject(clazz) {}

CClassConstSP const FwdCorrMode::TYPE = CClass::registerClassLoadMethod(
    "FwdCorrMode", typeid(FwdCorrMode), FwdCorrMode::load);

void FwdCorrModeStandard::checkInputs() const {
    const static string method = "FwdCorrModeStandard::validatePop2Object";
    try {
        if (nbStubDays < 0) {
            throw ModelException(method, "Nb of stub days for FwdCorrModeStandard must not be negative, "
                "but is " + Format::toString(nbStubDays));
        }
        if (Maths::isNegative(fwdMaxSqError) || Maths::isNegative(1.0-fwdMaxSqError)) {
            throw ModelException(method, "Max sq error for FwdCorrModeStandard must be between 0 and 1, "
                "but is " + Format::toString(fwdMaxSqError));
        }
        // some further validation (using a dummy date, but this doesnt matter)
        MaturityPeriodSP usedFwdCorrInterval(new MaturityPeriod(fwdCorrInterval));
        DateTime testDate1("01-Jan-2000",DateTime::START_OF_DAY);
        DateTime testDate2 = usedFwdCorrInterval->toDate(testDate1);
        string name = "Holiday";
        const HolidaySP holidayTmp(new Holiday(name, DateTimeArray(0), true));
        TimeMetricSP timeMetricTmp(new TimeMetric(0.0, holidayTmp.get()));
        int nbBusDaysDiff = timeMetricTmp->getHolidays()->
                    businessDaysDiff(testDate1, testDate2);
        if (nbBusDaysDiff < nbStubDays) {
            throw ModelException(method, "Mismatch in DependenceMakerGaussTerm for input for nbStubDays = "
                + Format::toString(nbStubDays) + " and fwdCorrInterval " 
                + fwdCorrInterval);
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

DoubleArray FwdCorrModeStandard::getBarriers() const {
    const static string method = "FwdCorrModeStandard::getBarriers";
    try {
        DoubleArray barriers(2);
        barriers[0] = fwdMaxSqError;
        barriers[1] = fwdMaxSqError; // same for lower and upper barrier
        return barriers;
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** invoked when Class is 'loaded' */
void FwdCorrModeStandard::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FwdCorrModeStandard, clazz);
    SUPERCLASS(FwdCorrMode);
    EMPTY_SHELL_METHOD(defaultFwdCorrModeStandard);
    FIELD(fwdCorrInterval, "interval when to compute fwd corrs");
    FIELD_MAKE_OPTIONAL(fwdCorrInterval);
    FIELD(nbStubDays, "nb of business days for stub period for fwd corrs");
    FIELD_MAKE_OPTIONAL(nbStubDays);
    FIELD(fwdMaxSqError, "barrier for fwd corr massaging");
    FIELD_MAKE_OPTIONAL(fwdMaxSqError);
    FIELD(interpolateAtmFwd, "vol interp atm fwd or atm spot");
    FIELD_MAKE_OPTIONAL(interpolateAtmFwd);
}

IObject* FwdCorrModeStandard::defaultFwdCorrModeStandard(){
    return new FwdCorrModeStandard();
}

/** for reflection */
FwdCorrModeStandard::FwdCorrModeStandard() : FwdCorrMode(TYPE), 
fwdCorrInterval("1W"), 
nbStubDays(3),
fwdMaxSqError(0.0025), // corresponds to 5% correlation
interpolateAtmFwd(true) {} 

string FwdCorrModeStandard::getFwdCorrInterval() const {
	return fwdCorrInterval;
}

int FwdCorrModeStandard::getNbStubDays() const {
	return nbStubDays;
}

bool FwdCorrModeStandard::doInterpolateAtmFwd() const {
	return interpolateAtmFwd;
}

CClassConstSP const FwdCorrModeStandard::TYPE = CClass::registerClassLoadMethod(
    "FwdCorrModeStandard", typeid(FwdCorrModeStandard), load);

void FwdCorrModeSimdates::checkInputs() const {
    const static string method = "FwdCorrModeSimdates::validatePop2Object";
    try {
        if (Maths::isNegative(fwdMaxSqErrorLow) || Maths::isNegative(1.0-fwdMaxSqErrorLow) ||
            Maths::isNegative(fwdMaxSqErrorHigh) || Maths::isNegative(1.0-fwdMaxSqErrorHigh) || 
            Maths::isNegative(fwdMaxSqErrorHigh - fwdMaxSqErrorLow)) {
            throw ModelException(method, "The two barriers must be between zero and one, and "
                "fwdMaxSqErrorLow must be less than fwdMaxSqErrorHigh. However, the values are "
                + Format::toString(fwdMaxSqErrorLow) + " and " + Format::toString(fwdMaxSqErrorHigh));
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

DoubleArray FwdCorrModeSimdates::getBarriers() const {
    const static string method = "FwdCorrModeStandard::getBarriers";
    try {
        DoubleArray barriers(2);
        barriers[0] = fwdMaxSqErrorLow;
        barriers[1] = fwdMaxSqErrorHigh;
        return barriers;
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/** invoked when Class is 'loaded' */
void FwdCorrModeSimdates::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FwdCorrModeSimdates, clazz);
    SUPERCLASS(FwdCorrMode);
    EMPTY_SHELL_METHOD(defaultFwdCorrModeSimdates);
    FIELD(fwdMaxSqErrorLow, "lower barrier for fwd corr massaging");
    FIELD_MAKE_OPTIONAL(fwdMaxSqErrorLow);
    FIELD(fwdMaxSqErrorHigh, "lower barrier for fwd corr massaging");
    FIELD_MAKE_OPTIONAL(fwdMaxSqErrorHigh);
    FIELD(interpolateAtmFwd, "vol interp atm fwd or atm spot");
    FIELD_MAKE_OPTIONAL(interpolateAtmFwd);

}

IObject* FwdCorrModeSimdates::defaultFwdCorrModeSimdates(){
    return new FwdCorrModeSimdates();
}

/** for reflection */
FwdCorrModeSimdates::FwdCorrModeSimdates() : FwdCorrMode(TYPE), 
fwdMaxSqErrorLow(0.0025), fwdMaxSqErrorHigh(0.01), interpolateAtmFwd(true) {} // corresponds to 5% and 10% correlation

/** constructor for backwards compatibility */
FwdCorrModeSimdates::FwdCorrModeSimdates(double barrierLow, double barrierHigh) : FwdCorrMode(TYPE),
fwdMaxSqErrorLow(barrierLow), fwdMaxSqErrorHigh(barrierHigh), interpolateAtmFwd(true) {}

bool FwdCorrModeSimdates::doInterpolateAtmFwd() const {
	return interpolateAtmFwd;
}

CClassConstSP const FwdCorrModeSimdates::TYPE = CClass::registerClassLoadMethod(
    "FwdCorrModeSimdates", typeid(FwdCorrModeSimdates), load);

/** Wrapper class for FwdCorrMode */

#define FWD_CORR_MODE_STANDARD      "STANDARD"
#define FWD_CORR_MODE_SIMDATES      "SIMDATES"

/** additional helper method for constructor & validatePop2Object */
void FwdCorrModeWrapper::init() {
    static const string method = "FwdCorrModeWrapper::init";
    try {
        if (fwdCorrModeType.empty()){
            fwdCorrModeType = FWD_CORR_MODE_STANDARD;
            fwdCorrModeStandard = FwdCorrModeStandardSP(new FwdCorrModeStandard());            
            realFwdCorrMode = fwdCorrModeStandard;
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** validatePop2Object is called first */
void FwdCorrModeWrapper ::validatePop2Object(){
    static const string routine = "FwdCorrModeWrapper::validatePop2Object";
    try{
        init();
        if (CString::equalsIgnoreCase(fwdCorrModeType, FWD_CORR_MODE_STANDARD)) {
            if (fwdCorrModeStandard.get()) { 
                realFwdCorrMode = fwdCorrModeStandard;
            } else {
                throw ModelException(routine, "Expected FwdCorrMode STANDARD " 
                                        "but none supplied!");
            }
        } else if (CString::equalsIgnoreCase(fwdCorrModeType, FWD_CORR_MODE_SIMDATES)) {
            if (fwdCorrModeSimdates.get()) { 
                realFwdCorrMode = fwdCorrModeSimdates;
            } else {
                throw ModelException(routine, "Expected FwdCorrMode SIMDATES "
                                        "but none supplied!");
            }
        } else {
            throw ModelException(routine, "Unrecognised FwdCorrMode "
                                 + fwdCorrModeType + ". Expected " 
                                 + FWD_CORR_MODE_STANDARD + ", " 
                                 + FWD_CORR_MODE_SIMDATES);
        }
    }  catch (exception& e){
        throw ModelException(e, routine);
    }
}

FwdCorrModeSP FwdCorrModeWrapper::getFwdCorrMode() const {
    return realFwdCorrMode;
}

/** invoked when Class is 'loaded' */
void FwdCorrModeWrapper::load(CClassSP& clazz){
    REGISTER(FwdCorrModeWrapper, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultFwdCorrModeWrapper);
    FIELD(fwdCorrModeType, "STANDARD or SIMDATES");
    FIELD(fwdCorrModeStandard,  "fwdCorrMode STANDARD");
    FIELD_MAKE_OPTIONAL(fwdCorrModeStandard);
    FIELD(fwdCorrModeSimdates,  "fwdCorrMode SIMDATES");
    FIELD_MAKE_OPTIONAL(fwdCorrModeSimdates);
    FIELD(realFwdCorrMode, "real fwdCorrMode");
    FIELD_MAKE_TRANSIENT(realFwdCorrMode);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
     
/** for reflection */
FwdCorrModeWrapper::FwdCorrModeWrapper() : CObject(TYPE) {
    init();
}

IObject* FwdCorrModeWrapper::defaultFwdCorrModeWrapper(){
    return new FwdCorrModeWrapper();
}

CClassConstSP const FwdCorrModeWrapper::TYPE =
CClass::registerClassLoadMethod("FwdCorrModeWrapper", 
                                typeid(FwdCorrModeWrapper), load);

/***************************************************************/
/* GAUSS DEPENDENCE                                            */
/* independent steps, assets with step independent correlation */
/***************************************************************/

/** Destructor */
Gauss::~Gauss(){ /*empty*/ }

/** Validation */
void Gauss::validatePop2Object(){  /*empty*/ }

/** Default Constructor */
Gauss::Gauss() {}

/** constructor from correlation matrix */
Gauss::Gauss(const DoubleMatrix& correlation ) : correlation(correlation) {
    static const string method("Gauss::Gauss");
    try {
        // dimension checking ...
        nbAssets = correlation.numCols();
        if( correlation.numRows() != nbAssets ) {
            throw ModelException(method, 
                                 "Have numRows " + Format::toString(correlation.numRows()) + 
                                 " and numCols " + Format::toString(nbAssets) );
        }
        isTrivial = false;
        if (correlation.isIdentity()) {
            isTrivial = true; 
            return;
        }
        try {  // allocations for speed in correlateSeries 
            correlCoeffs = correlation.computeSquareRoot(); 
        } catch (exception& e){
            throw ModelException(e,method, "Failed to compute square root of correlation matrix");
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

void Gauss::correlateSeries( DoubleMatrix& noise, int pathIdx ) {
    // pathIdx not needed here
    if (isTrivial) {
        return;
    }
    // Start from last dimension as we are changing them, first one unchanged
    int nbSteps = noise.numRows();
    for (int iAsset = nbAssets - 1; iAsset > 0; iAsset --) {
        double* mainSeries = noise[iAsset];
        double correlCoeffsDiag = correlCoeffs[iAsset][iAsset];
        for (int iStep = 0; iStep < nbSteps; iStep ++){
            mainSeries[iStep] *= correlCoeffsDiag;
        }
        for (int jAsset = iAsset - 1; jAsset > -1; jAsset --) {
            double* otherSeries = noise[jAsset];
            double correlCoeffsOffDiag = correlCoeffs[iAsset][jAsset];
            for (int iStep = 0; iStep < nbSteps; iStep ++) {
                mainSeries[iStep] += otherSeries[iStep] * correlCoeffsOffDiag;
            }
        }
    }
}
/** see description in pure virtual declaration */
DoubleMatrix Gauss::getCorrelations(int index){
    return correlation;
}
double Gauss::getCorrelationIJ(int i, int j, int index){
	return correlation[i][j];
}

///////////////////////////////////////////////////////////////////////////////////////////////
// EXTENDED SPARSE GAUSS
///////////////////////////////////////////////////////////////////////////////////////////////

/** Same as Gauss but can cope for case where the 'square root' of the
correlation matrix is supplied (and is not in lower triangular form) */
class ExtendedSparseGauss: public Gauss {
public:
    /** constructor */
    ExtendedSparseGauss(const SparseDoubleMatrix& sqrtCorrs) : 
      correlCoeffsSparse(sqrtCorrs) {
          /** the fields correlation and correlCoeffs are not populated */
		nbAssets = sqrtCorrs.numRows();
	}

	DoubleMatrix getCorrelations(int index) {
        static const string method("ExtendedSparseGauss::getCorrelations");
		throw ModelException(method,
            "ExtendedSparseGauss does not store the whole correlation matrix. " 
			"Use element-wise getCorrelationIJ instead.");
	}
	
	// setting up the Beta matrix structure, so that the complete correlation matrix
	//    R = B * T * T^ * B^   (where ^ is used for "transposed")
	//
	// beta correlations can be set up such that the "idiosyncratic" part is never too small
	// to that purpose -- 'maxcorr' parameter is used, so that B[a][a] >= sqrt(1-maxcorr^2)

	void setBetaCorrelations(const SparseDoubleMatrix& betas, const DoubleArray& maxcorr) {
        // if called with empty sparse matrix - clean up the member beta corr storage 
        // and do nothing else!
        if (!betas.size()) // it is empty
        {
            betaCoeffsSparse = betas; // make member empty too
            return;
        }


		// verify that betas has the same number of rows and cols as corrs
		if (betas.numCols() != correlCoeffsSparse.numCols() || 
			betas.numRows() != correlCoeffsSparse.numRows())
			throw ModelException("ExtendedSparseGauss::SetBetaCorrelations",
				"Beta Correlation Matrix must have the same size as regular correlations");

		if (betas.numCols() != (unsigned int) maxcorr.size())
			throw ModelException("ExtendedSparseGauss::SetBetaCorrelations",
				"Beta Correlation Matrix must have the same size as maximum corr array");

		if (betas.checkNoDiags())
			throw ModelException("ExtendedSparseGauss::SetBetaCorrelations",
				"Beta Correlation Matrix parameter must not have any diagonal elements."
				" They are dynamically computed instead.");

		// calculating and adding diag elements B[a][a] (and also 1 for B[i][i])
        // Important: the maxcorr is taken into account -- i.e. the values in beta matrix 
        // should be scaled if the calculated value of B[a][a] is less than sqrt(1-maxcorr[a]^2)
		SparseCollection sc;

        // TODO: unit test the calculations, see that all the indices are used properly
        //    (complication -- inverted indices of DoubleMatrix)

        // for each tier 1 or 2 object
        for(size_t a=0; a<betas.numCols(); ++a)
        {

            if (maxcorr[a]>1 || maxcorr[a]<0)
                throw ModelException("ExtendedSparseGauss::SetBetaCorrelations",
                    Format::toString("Maximum corr array has an element maxcorr[%d] = %g outside of range [0,1]", a, maxcorr[a]));

            // getting a single line of Bak
            DoubleArray v(betas.numCols(), 0.0);
            v[a] = 1.0;
            betas.leftMultiplyTo(v);

            double total = 0.0;

            for(size_t kv = 0; kv<betas.numCols(); ++kv)
                for(size_t kw = 0; kw<betas.numCols(); ++kw)
                {
                    double vw = v[kv]*v[kw];
                    total += (vw!=0.0) ? vw*correlCoeffsSparse.GetMMTElem(kv,kw) : 0.0;
                }

            // need rescaling of all betas of this line if total > maxcorr^2
            double corrscale = (total <= maxcorr[a]*maxcorr[a]) ? 1.0 : maxcorr[a]/sqrt(total);
            for(size_t i = 0; i<betas.numCols(); ++i)
                if (v[i] != 0.0)
                    sc.push_back(i, a, v[i]*corrscale);

            sc.push_back(a, a, sqrt(1.0 - corrscale*corrscale*total));
        }

		betaCoeffsSparse = sc;

    }


    void setBetaCorrelations(const SparseDoubleMatrix& betas, double maxcorr = 1.0) {
		setBetaCorrelations(betas, DoubleArray(betas.numCols(), maxcorr));
	}

	// comp correlation is rho = Bia Can Cbn Bjb = (Bi)a (Can Cbn) (Bj)b 
	inline double getCompCorrelation(int i, int j, const SparseDoubleMatrix& betas){

		if (!betas.size())
			return correlCoeffsSparse.GetMMTElem(i,j);

		DoubleArray v(betas.numCols(), 0.0);
		v[i] = 1.0;
		betas.leftMultiplyTo(v);

		DoubleArray w(betas.numCols(), 0.0);
		w[j] = 1.0;
		betas.leftMultiplyTo(w);

        double total = 0.0;

        // the intended number of betas per each tier 2 asset is only 1, 
        // so these cycles below are likely to lead to just one call to expensive GetMMTElem()
		for(size_t kv = 0; kv<betas.numCols(); ++kv)
            for(size_t kw = 0; kw<betas.numCols(); ++kw)
		    {
                double vw = v[kv]*w[kw];
                total += (vw!=0.0) ? vw*correlCoeffsSparse.GetMMTElem(kv,kw) : 0.0;
		    }

		return total;
	}

	double getCorrelationIJ(int i, int j, int index){
		return getCompCorrelation(i,j,betaCoeffsSparse);
	}

	void correlateSeries( DoubleMatrix& noise, int pathIdx) {
		correlCoeffsSparse.leftMultiplyTo(noise);

		// if Tier2 correlations are prescribed as well
		if (betaCoeffsSparse.size())
		 betaCoeffsSparse.leftMultiplyTo(noise);
	}


    // DEBUG
     void print(ostream& out, int mode = 0) // 0-compact, 1-full
     {
         // in the inversed col/row notation of qlib the lower triangular matrix prints out
         // transposed, i.e. upper triangular.

         out<<"\n\nSqrt of corr matrix:\n";
         correlCoeffsSparse.print(out, mode);
         out<<"\n\nBeta matrix:\n";
         betaCoeffsSparse.print(out, mode);
      
     }



protected:
    // only a decomposed (i.e. sqrt of) corr matrix is stored
	SparseDoubleMatrix  correlCoeffsSparse;

    // this is to implement Tier2-Tier1 structure of correlations
    SparseDoubleMatrix  betaCoeffsSparse;
};




///////////////////////////////////////////////////////////////////////////////////////////////
// GAUSS
///////////////////////////////////////////////////////////////////////////////////////////////

DependenceSP DependenceMakerGauss::createDependence(const DependenceMaker::ISupport* support) const {
    static const string method("DependenceMakerGauss::createDependence");
    try {
        /** support used for validation to get specific data */
        const DependenceMakerGauss::Support* supportDmg = 
            dynamic_cast<const DependenceMakerGauss::Support*>(support);
        if (!supportDmg) {
            throw ModelException(method, "DependenceMakerGauss is not a valid dependence maker.");
        }        
        CDoubleMatrix corrMatrix(*supportDmg->getGaussData()); 
        return DependenceSP(new Gauss(corrMatrix));
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// Invoked when Class is 'loaded'
void DependenceMakerGauss::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DependenceMakerGauss, clazz);
    IMPLEMENTS(IRevertTypeConvert);
    SUPERCLASS(DependenceMaker);
    EMPTY_SHELL_METHOD(defaultDependenceMakerGauss);
}

IObject* DependenceMakerGauss::defaultDependenceMakerGauss(){
    return new DependenceMakerGauss();
}

// for the IMS interface (see IRevertTypeConvert)
IObjectSP DependenceMakerGauss::revert(const string& interfaceType) const {
    static const string method = "DependenceMakerGauss::revert";
    try {
        if (interfaceType != IRevertTypeConvert::PYRAMID) {
            throw ModelException( method, 
                                  "Cannot convert a DependenceMakerGauss for"
                                  " the interface " + interfaceType);               
        }
        
        IObjectSP wrapper(new DependenceMakerWrapper(DependenceMakerGaussSP(copy(this))));
        return wrapper;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DependenceMakerGauss::DependenceMakerGauss() : DependenceMaker(TYPE) {};

CClassConstSP const DependenceMakerGauss::TYPE = CClass::registerClassLoadMethod(
	"DependenceMakerGauss", typeid(DependenceMakerGauss), DependenceMakerGauss::load);


/******************************/
/* GAUSS SRM DEPENDENCE MAKER */
/******************************/

void DependenceMakerGaussSrm::validatePop2Object() {
    static const string method("DependenceMakerGaussSrm::validatePop2Object");
    try {
        init();
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

void DependenceMakerGaussSrm::init() {
    static const string method("DependenceMakerGaussSrm::validatePop2Object");
    try {
        fwdCorrModeUsed = fwdCorrModeDetails->getFwdCorrMode();
        fwdCorrModeUsed->checkInputs();
        checkInputs();
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

void DependenceMakerGaussSrm::checkInputs() const {
    static const string method("DependenceMakerGaussSrm::checkInputs");
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

/** modifies MDF if necessary */
void DependenceMakerGaussSrm::modifyMarketDataFetcher(MarketDataFetcherSP mdf) {
    if (corrTermStructure) {
        MDFUtil::setCorrTermMode(*mdf, true);
    }
}

void DependenceMakerGaussSrm::setCorrTermStructureMode(bool doCorrTermStructure) {
    corrTermStructure = doCorrTermStructure; 
}

bool DependenceMakerGaussSrm::getCorrTermStructureMode() const {
    return corrTermStructure;
}

DependenceSP DependenceMakerGaussSrm::createDependence(const DependenceMaker::ISupport* support) const {
    static const string method("DependenceMakerGaussSrm::createDependence");
    try {
        /** support used for validation to get specific data */
        const DependenceMakerGaussSrm::Support* supportDmgSrm = 
            dynamic_cast<const DependenceMakerGaussSrm::Support*>(support);
        if (!supportDmgSrm) {
            throw ModelException(method, "DependenceMakerGaussSrm is not a valid dependence maker.");
        }

        /** case 1 -- NO corr term structure and NO corr mapping */
        if ( (!corrTermStructure && !corrMapping) || (supportDmgSrm->nbEqEqAssets()==1)) {
            /** retrieve sparse corr matrix from pathConfig */
            vector<SparseDoubleMatrixSP> sparseGaussMatrixArray = 
                supportDmgSrm->createSparseGaussMatrixArray(this);
            ExtendedSparseGauss* pGauss = new ExtendedSparseGauss(*sparseGaussMatrixArray[0]);
            SparseDoubleMatrixSP betaCorrelations = supportDmgSrm->getBetaCorrelations(this);
            if (betaCorrelations.get()) // if it is non-empty
                pGauss->setBetaCorrelations(*betaCorrelations, supportDmgSrm->getMaxBetaCorr(this));            
            //pGauss->print(cerr, 1); // DEBUG printout
            return DependenceSP(pGauss);
        }

        /** note: no notion of time metric in srm => create dummy array */
        string name = "Holiday";
        const HolidaySP holiday(new Holiday(name, DateTimeArray(0), false));
        TimeMetricSP timeMetric(new TimeMetric(0.05, holiday.get()));
        TimeMetricArray timeMetricArrayDummy(1, timeMetric);

        IntArray fwdCorrIndexArray;
        DateTimeArray fwdCorrDatesArray;            
        DependenceMakerGaussTerm::createFwdCorrIdxAndDatesArray(
            supportDmgSrm->getSimDates(), // valuedate + simdates
            fwdCorrModeUsed,    // standard (default) or simdates
            timeMetricArrayDummy,  // just a dummy here
            fwdCorrIndexArray,  // output
            fwdCorrDatesArray); // output

        /** case 2 -- NO corr term structure and YES corr mapping */
        if (!corrTermStructure && corrMapping) {
            /** retrieve sparse corr matrix array from pathConfig */
            vector<SparseDoubleMatrixSP> sparseGaussMatrixArray = 
                supportDmgSrm->createSparseGaussMatrixArray(this, 
                                                            fwdCorrIndexArray, 
                                                            fwdCorrDatesArray);
            return DependenceSP(new ExtendedSparseGaussTerm(sparseGaussMatrixArray,
                                                            fwdCorrIndexArray));
        }        

        /** case 3+4 -- YES corr term structure and NO/YES corr mapping 
            distinction is done in pathconfig */

        /** timeline */
        DateTimeArray simDates = supportDmgSrm->getSimDates();
        DateTime valueDate = simDates.front();        

        /** asset names */
        CorrelationCommonArray  corrObjArray;
        CorrelationTermArray    corrTermObjArray; 
        supportDmgSrm->getCorrelationData(corrObjArray,
                                          corrTermObjArray); 
            
        // get further corrTermData for GaussTerm via static method in CorrelationTerm
        CorrTermDataSP corrTermData = 
            CorrelationTerm::getCorrelationTermSqueezesAndExpiries (
            valueDate, supportDmgSrm->nbEqEqAssets(), 
            corrObjArray, corrTermObjArray); 

        /** convert fwdVar at simDates to fwdVar at fwdCorrDates */
        DoubleMatrix fwdVarAtSimDates = supportDmgSrm->getFwdVarAtDates(); // nbAssets x nbDates
        CDoubleMatrixSP fwdVarAtFwdCorrDates = 
            DependenceMakerGaussTerm::fwdVarAtFwdCorrDates(simDates,
                                                           fwdCorrIndexArray,
                                                           fwdVarAtSimDates);
        
        TimeMetricArray timeMetricArray(fwdVarAtFwdCorrDates->numCols(), timeMetric);
        BoolArray skipFwdCorrelation(corrTermData->benchmarkTermExpiries->size(), false);
            
        /** compute fwdCorrelation matrices */   
        DoubleMatrixArraySP fwdCorrelation = 
            CorrelationTerm::CorrelationTermMatrix(valueDate, 
                                                   fwdCorrDatesArray, 
                                                   *fwdVarAtFwdCorrDates, 
                                                   timeMetricArray, 
                                                   corrTermData,
                                                   CorrelationTerm::forward, 
                                                   eigenValueFloor, 
                                                   maxSqError, 
                                                   fwdCorrModeUsed->getBarriers(), 
                                                   skipFwdCorrelation);
           
        vector<SparseDoubleMatrixSP> sparseGaussMatrixArray = 
            supportDmgSrm->createSparseGaussTermMatrixArray(this, 
                                                            fwdCorrelation, 
                                                            fwdCorrIndexArray,
                                                            fwdCorrDatesArray);
        
        return DependenceSP(new ExtendedSparseGaussTerm(sparseGaussMatrixArray,
                                                        fwdCorrIndexArray));    
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// Invoked when Class is 'loaded'
void DependenceMakerGaussSrm::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DependenceMakerGaussSrm, clazz);
    SUPERCLASS(DependenceMaker);
    EMPTY_SHELL_METHOD(defaultDependenceMakerGaussSrm);
    FIELD(corrMapping, "whether or not adj eqeq corrs");
    FIELD_MAKE_OPTIONAL(corrMapping);
    FIELD(corrTermStructure, "whehter or not corr term structure");
    FIELD_MAKE_OPTIONAL(corrTermStructure);
    FIELD(eigenValueFloor, "");
    FIELD_MAKE_OPTIONAL(eigenValueFloor);
    FIELD(maxSqError, "");
    FIELD_MAKE_OPTIONAL(maxSqError);
    FIELD(fwdCorrModeDetails, "");
    FIELD_MAKE_OPTIONAL(fwdCorrModeDetails);
    FIELD(fwdCorrModeUsed, "");
    FIELD_MAKE_TRANSIENT(fwdCorrModeUsed);
}

IObject* DependenceMakerGaussSrm::defaultDependenceMakerGaussSrm(){
    return new DependenceMakerGaussSrm();
}

DependenceMakerGaussSrm::DependenceMakerGaussSrm(): DependenceMaker(TYPE), 
corrMapping(false), 
corrTermStructure(false),
eigenValueFloor(CorrelationTerm::EIGEN_VALUE_FLOOR), 
maxSqError(CorrelationTerm::MAX_SQ_ERROR), 
fwdCorrModeDetails(new FwdCorrModeWrapper()) {
    init(); 
}

DependenceMakerGaussSrm::DependenceMakerGaussSrm(CClassConstSP clazz) : DependenceMaker(clazz),
corrMapping(false), 
corrTermStructure(false),
eigenValueFloor(CorrelationTerm::EIGEN_VALUE_FLOOR), 
maxSqError(CorrelationTerm::MAX_SQ_ERROR), 
fwdCorrModeDetails(new FwdCorrModeWrapper()) {
    init(); 
}

bool DependenceMakerGaussSrm::doCorrMapping() const {
	return corrMapping;
}

double DependenceMakerGaussSrm::getEigenValueFloor() const {
	return eigenValueFloor;
}

double DependenceMakerGaussSrm::getMaxSqError() const {
	return maxSqError;
} 

FwdCorrModeSP DependenceMakerGaussSrm::getFwdCorrModeUsed() const {
    static const string method("DependenceMakerGaussSrm::getFwdCorrModeUsed");
    try {
        if (!fwdCorrModeUsed) {
            throw ModelException("FwdCorrModeUsed requested, but not allocated");
        } else {
            return fwdCorrModeUsed;
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const DependenceMakerGaussSrm::TYPE = CClass::registerClassLoadMethod(
	"DependenceMakerGaussSrm", typeid(DependenceMakerGaussSrm), DependenceMakerGaussSrm::load);

DRLIB_END_NAMESPACE

