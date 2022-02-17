//----------------------------------------------------------------------------
//
//   Group       : QR Equities 
//
//   Filename    : DependenceGauss.cpp
//
//   Description : Holds Basic Dependence (Maker) and POISSON and CLAYTON
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_DEPENDENCE_CPP
#include "edginc/Dependence.hpp"
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/DependenceLocalCorr.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/***************************/
/* BASIC DEPENDENCE OBJECT */
/***************************/

/* Invoked when Class is 'loaded' */
void Dependence::load(CClassSP& clazz) {
	REGISTER(Dependence, clazz);
	SUPERCLASS(CObject);
}

/* constructor */
Dependence::Dependence(): CObject(TYPE) {}
Dependence::Dependence(const CClassConstSP& clazz): CObject(clazz) {}

/* definition of TYPE */
CClassConstSP const Dependence::TYPE = CClass::registerClassLoadMethod(
	"Dependence", typeid(Dependence), load);


/**********************/
/* POISSON DEPENDENCE */
/**********************/

/** Destructor */
Poisson::~Poisson(){ /*empty*/ }

/** Validation */
void Poisson::validatePop2Object(){  /*empty*/ }

/** Default Constructor */
Poisson::Poisson(){}

/** Constructor */
Poisson::Poisson(double crashRate) : crashRate(crashRate) {
    static const string method("Poisson::Poisson");
    try{
        // checking ...
        if( Maths::isNegative(crashRate) ) {
            throw ModelException(method, 
                                 "Have crashRate " +
                                 Format::toString(crashRate) + 
                                 ", should be >= 0." );
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/** Returns the arrival times of Poisson jumps: assumes a matrix with one column */
void Poisson::correlateSeries( DoubleMatrix& noise, int pathIdx ) {
    static const string method("Poisson::correlateSeries");
    try{
        // pathIdx not needed here
        int nbCols = noise.numCols();
        int nbRows = noise.numRows();
        int iRow, iCol;
        double uniform;

        for( iCol=0; iCol<nbCols; iCol++ ) {
            double* mainSeries = noise[iCol];
            uniform = Maths::max( .0001, Maths::min( .9999, N1(mainSeries[0]) ) );
            mainSeries[0] = -log(uniform)/crashRate;
            for ( iRow=1; iRow<nbRows; iRow++ ) {
                uniform = Maths::max( .0001, Maths::min( .9999, N1(mainSeries[iRow]) ) );
                mainSeries[iRow] 
                    = -log(uniform)/crashRate + mainSeries[iRow-1];
            }
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

//// see desription in pure virtual declaration
DoubleMatrix Poisson::getCorrelations(int /*index*/){
    throw ModelException("Poisson::getCorrelations", "Not implemented yet");
}
double Poisson::getCorrelationIJ(int /*i*/, int /*j*/, int /*index*/){
	throw ModelException("Poisson::getCorrelations", "Not implemented yet");
}

/** Calculates the cumulative Poisson distribution with intensity crashRate*time */
void Poisson::numJumpsPreprocess( DoubleArray& timeInYears, int maxJumps ) {
    cumulative = DoubleMatrix(timeInYears.size(),maxJumps);
    int iDate, iStep;
    for( iDate=0; iDate<timeInYears.size(); iDate++ ) {
        double* thisCumul   = cumulative[iDate];
        double proba        = exp(-timeInYears[iDate]*crashRate);
        double dt           = timeInYears[iDate];
        double oneCumul     = 0;
        for( iStep=0; iStep<maxJumps; iStep++ )
        {
            oneCumul += proba;
            thisCumul[iStep] = oneCumul;
            proba *= crashRate * dt / (double)(iStep+1);
        }
    }
}

/** Returns the number of jumps between dates 
    - independent assets, 
    - time intervals as set out in numJumpsPreprocess */
void Poisson::numJumps( DoubleMatrix& noise ) {

    static const string method("Poisson::numJumps");

    try{
        int nbCols = noise.numCols();
        int nbRows = noise.numRows();
        int iCol;

        for( iCol=0; iCol<nbCols; iCol++ ) {
            numJumps(noise[iCol],nbRows);
        }

    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/** Returns the number of jumps between dates 
    - independent assets, 
    - time intervals as set out in numJumpsPreprocess */
void Poisson::numJumps( double* noise, int nbRows ) {
    static const string method("Poisson::numJumps");
    try{
        int iRow;
        double uniform;

        for( iRow=0; iRow<nbRows; iRow++ ) {
            double* thisCumul  = cumulative[iRow];
            uniform = N1(noise[iRow]);
            int iJump = 0;
            while( uniform>thisCumul[iJump] && iJump<cumulative.numRows() ) {
                iJump++;
            }
            noise[iRow] = (double)iJump;
        }

    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/**********************/
/* CLAYTON DEPENDENCE */
/**********************/

/** Destructor */
Clayton::~Clayton(){ /*empty*/ }

/** Validation */
void Clayton::validatePop2Object() { 
    static const string method("Clayton::validatePop2Object");

    try{
        //if( (alpha.numCols() != nbAssets) || (alpha.numRows() != nbAssets) ) {
        //    throw ModelException(method, "dimensions of alpha should be (nbAssets x nbAssets)");
        //}
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

/** Default Constructor */
Clayton::Clayton(){}

/** Constructor */
Clayton::Clayton(
    const double        theta,
    const DoubleMatrix& alpha ):
    theta(theta),
    alpha(alpha) {
    static const string method("Clayton::Clayton");

    // dimension checking ...
    nbAssets = alpha.numCols();

    // allocations for speed in correlateSeries
    tempA = CDoubleArray(nbAssets);
    tempU = CDoubleArray(nbAssets);
    tempX = CDoubleArray(nbAssets);

    // theta in [0,1] assumed
    if( theta>0 && theta<1 )
    {
        thetaUsed = 2.*theta / (1.-theta);
    }
    else if( theta <= 0 )
    {
        thetaUsed = 0.0000001;
    }
    else
    {
        thetaUsed = 50.0;
    }
}


/** Help function for DependenceMatrix: parametrisation of term structure */
void Clayton::correlateSeries( DoubleMatrix& noise, int pathIdx ) {

    static const string method("Clayton::correlateSeries");

    try{
        // pathIdx not needed here
        // do the same at each time step, might need time dependent parms...
        int iAsset, jAsset, iStep;
        int nbSteps = noise.numCols(); // ???

        for( iStep=0; iStep<nbSteps; iStep++ ) {

            // vector of uniform rv's
            for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
		        tempU[iAsset] = N1(noise[iAsset][iStep]);
	        }

            // vector of gamma rv's
	        for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
                gammaRandom = imsl_d_random_gamma(1,1./thetaUsed,0);
                if( !gammaRandom ) {
                    throw ModelException(method, "failed to generate gamma variate");
                }
                tempA[iAsset] = gammaRandom[0];
            }

            // the algorithm
	        for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
                double oneGammaRandom = 0;
	            for( jAsset=0; jAsset<nbAssets; jAsset++ ) {
                    oneGammaRandom += alpha[iAsset][jAsset] * tempA[jAsset];
                }
                double logUoverGamma = log(tempU[iAsset]) / oneGammaRandom;
                tempX[iAsset] = 1;
	            for( jAsset=0; jAsset<nbAssets; jAsset++ ) {
                    tempX[iAsset] 
                        *= pow( 1-alpha[iAsset][jAsset]*logUoverGamma, -1/thetaUsed );
                }

                // vector of normal rv's
    		    noise[iAsset][iStep] = N1Inverse( tempX[iAsset] );
            }
        }

    } catch (exception& e){
        throw ModelException(e,method);
    }
}
//// see desription in pure virtual declaration
DoubleMatrix Clayton::getCorrelations(int /*index*/){
    throw ModelException("Clayton::getCorrelations", "Not implemented yet");
}
double Clayton::getCorrelationIJ(int /*i*/, int /*j*/, int /*index*/){
	throw ModelException("Clayton::getCorrelations", "Not implemented yet");
}

/**************************/
/* BASIC DEPENDENCE MAKER */
/**************************/

/** modifies MDF if necessary, default implementation empty */
void DependenceMaker::modifyMarketDataFetcher(MarketDataFetcherSP mdf) {}

/* definition of TYPE */
CClassConstSP const DependenceMaker::TYPE = CClass::registerClassLoadMethod(
	"DependenceMaker", typeid(DependenceMaker), DependenceMaker::load);

// Invoked when Class is 'loaded'
void DependenceMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DependenceMaker, clazz);
    SUPERCLASS(CObject);
}

DependenceMaker::DependenceMaker(const CClassConstSP& clazz): CObject(clazz) {}

/********************/
/* BASIC SKEW MAKER */
/********************/

CClassConstSP const SkewMaker::TYPE = CClass::registerClassLoadMethod(
	"SkewMaker", typeid(SkewMaker), SkewMaker::load); 

void SkewMaker::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SkewMaker, clazz);
    SUPERCLASS(DependenceMaker);    
}
    
SkewMaker::SkewMaker(const CClassConstSP& clazz) : DependenceMaker(clazz) {}

/***********/
/* for IMS */
/***********/

#define DEPENDENCEMAKER_TYPE_GAUSS          "GAUSS"
#define DEPENDENCEMAKER_TYPE_GAUSS_SRM      "GAUSSSRM"
#define DEPENDENCEMAKER_TYPE_GAUSS_TERM     "GAUSSTERM"
#define DEPENDENCEMAKER_TYPE_GAUSS_TERM_SRM "GAUSSTERMSRM"
#define DEPENDENCEMAKER_TYPE_LOCAL_CORR     "LOCALCORR"

DependenceMakerWrapper::~DependenceMakerWrapper(){}

/** validatePop2Object is called first */
void DependenceMakerWrapper::validatePop2Object(){
	static const string routine = "DependenceMakerWrapper::validatePop2Object";
	try{
		if (dependenceMakerType.empty()){
			throw ModelException(routine, 
				"Blank Dependendence Maker Type specified!");
		}
		if (CString::equalsIgnoreCase(dependenceMakerType, DEPENDENCEMAKER_TYPE_GAUSS)) {
			if (dependenceMakerGauss.get()) { 
				realDependenceMaker = dependenceMakerGauss;
			} else {
				throw ModelException(routine, "Expected Dependence Maker GAUSS "
					"but none supplied!");
			}
		} else if (CString::equalsIgnoreCase(dependenceMakerType, DEPENDENCEMAKER_TYPE_GAUSS_SRM)) {
			if (dependenceMakerGaussSrm.get()) { 
				realDependenceMaker = dependenceMakerGaussSrm;
			} else {
				throw ModelException(routine, "Expected Dependence Maker GAUSSSRM "
					"but none supplied!");
			}
		} else if (CString::equalsIgnoreCase(dependenceMakerType, DEPENDENCEMAKER_TYPE_GAUSS_TERM)) {
			if (dependenceMakerGaussTerm.get()) {
				realDependenceMaker = dependenceMakerGaussTerm;
			} else {
				throw ModelException(routine, "Expected Dependence Maker GAUSSTERM "
					"but none supplied!");
			}
		} else if (CString::equalsIgnoreCase(dependenceMakerType, DEPENDENCEMAKER_TYPE_GAUSS_TERM_SRM)) {
			if (dependenceMakerGaussTermSrm.get()) {
				realDependenceMaker = dependenceMakerGaussTermSrm;
			} else {
				throw ModelException(routine, "Expected Dependence Maker GAUSSTERMSRM "
					"but none supplied!");
			}
        } else if (CString::equalsIgnoreCase(dependenceMakerType, DEPENDENCEMAKER_TYPE_LOCAL_CORR)) {
			if (dependenceMakerLocalCorr.get()) {
				realDependenceMaker = dependenceMakerLocalCorr;
			} else {
				throw ModelException(routine, "Expected Dependence Maker LOCALCORR "
					"but none supplied!");
			}
		} else {
			throw ModelException(routine, "Unrecognised Dependence Maker "
				+ dependenceMakerType + ". Expected " 
				+ DEPENDENCEMAKER_TYPE_GAUSS + ", " 
				+ DEPENDENCEMAKER_TYPE_GAUSS_SRM + ", " 
				+ DEPENDENCEMAKER_TYPE_GAUSS_TERM + ", " 
				+ DEPENDENCEMAKER_TYPE_GAUSS_TERM_SRM + ", "
                + DEPENDENCEMAKER_TYPE_LOCAL_CORR);
		}
	}  catch (exception& e){
		throw ModelException(e, routine);
	}
}

/** convert is called second */
void DependenceMakerWrapper::convert(IObjectSP&    object,
								     CClassConstSP requiredType) const {
	static const string method = "DependenceMakerWrapper::convert";
	try {
		if (requiredType != DependenceMaker::TYPE) {
			throw ModelException(method, 
				"Cannot convert a DependenceMakerWrapper into "
				"object of type "+ requiredType->getName());               
		}
		object = realDependenceMaker;
	}
	catch (exception& e) {
		throw ModelException(e, method);
	}
}

/** Invoked when Class is 'loaded' */
void DependenceMakerWrapper::load(CClassSP& clazz){
	REGISTER(DependenceMakerWrapper, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(ITypeConvert);
	EMPTY_SHELL_METHOD(defaultDependenceMakerWrapper);
	FIELD(dependenceMakerType, "GAUSS or GAUSSTERM");
	FIELD(dependenceMakerGauss,  "dependence maker GAUSS");
	FIELD_MAKE_OPTIONAL(dependenceMakerGauss);
	FIELD(dependenceMakerGaussSrm,  "dependence maker GAUSSSRM");
	FIELD_MAKE_OPTIONAL(dependenceMakerGaussSrm);
	FIELD(dependenceMakerGaussTerm,  "dependence maker GAUSSTERM");
	FIELD_MAKE_OPTIONAL(dependenceMakerGaussTerm);
	FIELD(dependenceMakerGaussTermSrm,  "dependence maker GAUSSTERMSRM");
    FIELD_MAKE_OPTIONAL(dependenceMakerGaussTermSrm);
	FIELD(dependenceMakerLocalCorr,  "dependence maker LOCALCORR");
    FIELD_MAKE_OPTIONAL(dependenceMakerLocalCorr);
	FIELD(realDependenceMaker, "real dependence maker");
	FIELD_MAKE_TRANSIENT(realDependenceMaker);
	clazz->setPublic(); // make visible to EAS/spreadsheet
}

// for reflection
DependenceMakerWrapper::DependenceMakerWrapper(): CObject(TYPE) {}

IObject* DependenceMakerWrapper::defaultDependenceMakerWrapper(){
	return new DependenceMakerWrapper();
}

DependenceMakerWrapper::DependenceMakerWrapper(DependenceMakerGaussSP dependenceMaker): 
CObject(TYPE),
dependenceMakerType(DEPENDENCEMAKER_TYPE_GAUSS),
dependenceMakerGauss(dependenceMaker) {
	realDependenceMaker = dependenceMakerGauss;
}

DependenceMakerWrapper::DependenceMakerWrapper(DependenceMakerGaussTermSP dependenceMaker):
CObject(TYPE),
dependenceMakerType(DEPENDENCEMAKER_TYPE_GAUSS_TERM),
dependenceMakerGaussTerm(dependenceMaker) {
    realDependenceMaker = dependenceMakerGaussTerm;
}

DependenceMakerWrapper::DependenceMakerWrapper(DependenceMakerLocalCorrSP dependenceMaker):
CObject(TYPE),
dependenceMakerType(DEPENDENCEMAKER_TYPE_LOCAL_CORR),
dependenceMakerLocalCorr(dependenceMaker) {
    realDependenceMaker = dependenceMakerLocalCorr;
}

CClassConstSP const DependenceMakerWrapper::TYPE =
CClass::registerClassLoadMethod("DependenceMakerWrapper", 
								typeid(DependenceMakerWrapper), load);

DRLIB_END_NAMESPACE

