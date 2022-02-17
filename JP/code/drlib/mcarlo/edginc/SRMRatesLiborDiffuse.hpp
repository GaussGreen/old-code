//----------------------------------------------------------------------------
//
//	 Group		 : QR
//
//	 Filename	 : SRMRatesLiborDiffuse.hpp 
//
//	 Description : forward libor model path generation
//
//
//----------------------------------------------------------------------------

#ifndef	SRMRATESLIBORDIFF_HPP
#define	SRMRATESLIBORDIFF_HPP

#include <cstdio>
#include <cassert>
#include <set>

#include "edginc/DECLARE.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MCRandomLite.hpp"

#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include "edginc/SRMSwap.hpp"
#include "edginc/SRMRatesLiborUtil.hpp"
#include "edginc/QMCFXBaseDiffuse.hpp"

DRLIB_BEGIN_NAMESPACE

class IQMCHelperTimeLogic;

const double TOO_SMALL = 1.e-10;
const double BGM_MIN_VOL = 1.e-4; 

typedef enum {
	no_forward,
	smallStep,
	predictorCorrector
} simulationType;

typedef enum {
    no_calib,
	diagonal,
	column,
	CMSspread,
} calibRate;

// base class for low level IR path generator
class SRMRatesLiborDiffuse: public QMCRatesDiffuse
{

public:

	virtual	~SRMRatesLiborDiffuse();

	SRMRatesLiborDiffuse();

	virtual	void setSRMRatesLiborDiffuse(
		int					   randomIndex,
		const DateTime&		   today,
		SRMRatesLiborUtilSP	   srmRatesLiborSP,	
		bool				   saveDomLnMoney,
		bool				   saveSigmaR,
		double				   NbSigmasMax,
		double				   NbSigmasMin,
		const vector<double>&  df, // for historic dates
		const vector<double>&  irFxCorr,  // TODO : not used at the moment
		bool				   calibrateAgainstSwaptionVols,
		const string&          calibStyle,
		const string&          calibMaturity,
		const string&          calibMaturityCMS);

	void setLiborModel(
		long nbFactors,
		simulationType simType);

	 // two-step semi-analytical calibration routine :
	 // 1. caplet vols
	 // 2. gradient descent minimization algorithm for swaption vols 
	 void analyticCalib(double     calibAccuracy);
  
    size_t getDiscYCIdx(void) 
    {
        if (discYCIdx < 0)
        {
            ASSERT(liborUtilSP.get()); 
            discYCIdx =  registerYCFlavor(liborUtilSP->getDiscYC());
        }
        return discYCIdx;

    }
    size_t getDiffYCIdx(void) 
    {
        if (diffYCIdx < 0)
        {
            ASSERT(liborUtilSP.get()); 
            diffYCIdx =  registerYCFlavor(liborUtilSP->getDiffYC());
        }
        return diffYCIdx;
    }

    QMCRatesUtilSP getRatesUtil(void) { ASSERT(liborUtilSP.get()); return liborUtilSP;}

	/**	finalize the timeline, allocate	memory */
	// addtional initialization	when allDates is know. no more new dates after this	point
	void finalize(DateTimeArrayConstSP allDates);

    void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

	void generatePathSmallStep(IQMCRNGManagerSP rngMgr);

	void generatePathPredCorr(IQMCRNGManagerSP rngMgr);

	// for capturing state of	diffusion
	struct DiffusedState
	{
		int freezeIdx; // index of last "frozen" libor
		double obsTime; // offset from today
		double dt; // time to first reset, offset from obsTime
		double shortRate;
		vector<double> stateDf; // engine discount factors 
	};
	
	// for capturing diffused state on each date on which calculations performed
	typedef	vector<DiffusedState>::const_iterator DiffusedStateIter;

	// performance critical	functions
	double getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx	j) 
	{
		return log(this->getExpectedDiscFactor(idx, i, j));
	}

	double getExpectedDiscFactor(size_t	idx, FwdIdx	i, FwdIdx j);

    // called after origSVol externally recalibrated via ICE 
	// not defined for BGM
	virtual void recalibrate(QMCRatesUtilSP thisSRMRatesUtilSP) { return; } 

    void printModel(void) const;

protected:

    SRMRatesLiborUtilSP liborUtilSP;

	// setParameters: for given parameter vectors a,b,d,lam it sets up the 
	// matrix members m_A, m_B, m_D and the vector member m_lambda
	void setParameters(const vector<double> &a, 
		               const vector<double> &b, 
		               const vector<double> &d, 
					   const double &lam,
					   const vector<double> &h );

    void choleskyFactor(double timeFirst, double timeLast);

	void calibrateCaplets(void);

	// since there may be more than one libor per caplet, we supplement
	// the previous analytical caplet calibration with a numerical one  
	void numericCalibCaplets(void);

	void fudgeCalib(
		const vector<double> &a, 
		const vector<double> &b,
		const vector<double> &d, 
		const double &lambda,
		const vector<double> &h);

	void bumpFudge(
		double startDate, 
		double endDate,
		double bumpA,
		double bumpD);

	void setCalibDetails(
		const string& calibType,
	    double calibMat,
		double calibMatCMS);

	double getTotalError(void);

    double getCorError(void);

	void setModelCorrelationMatrix(double time1, double time2);

    void setModelSwaptionGrid();

	// getLiborVols: returns the instantaneous factor 
	// loadings for every libor
	void getLiborVols(
		vector< vector<double> > &instLiborVol, 
		double time,
		long   firstLibor,
		long   lastLibor) const;
	
	// getFactorVols: provides building blocks used to compute integrated 
	// Libor covariances and effective correlations. It is used both in model
	// calibration to calculate the model swaption volatility and in the engine
	// to implement the predictor-corrector method. The output is the two-dimensional
	// array facLdArray with num. of Libors rows and 3*num. of factor columns
	// The dimensions of the matrix are assumed to be known
	void getFactorVols(
		vector< vector<double> > &facLdArray, 
		long   firstLibor, 
		long   lastLibor,
		double time1, 
		double time2);

	// returns skew parameters for libor rates
	const vector<double>& getLiborSkew(void) const;

	// sets skew paarameter for each LIBOR rate
	void setLiborSkew(vector<double> &skew);

	// getVolError: as getTotalError, but for a subset of swaptions with
	// column and rows from 0 to row
	double getVolError(long calibRow, long calibColumn);

	// getModelBlackVolatility: 
	// analytical swaption Black volatility using
	// Rebonato's approximation
	double getModelBlackVolatility(const SRMSwapClass &swap);

	double getSwapRateCorr(
		SRMSwapClass &swap1,
		SRMSwapClass &swap2,
		double start,
		double end);

	DateTime today;

	const double m_lastSimTime; // in years
    const int m_liborInterval; // in months

	vector<SRMSwapClassSP> m_gridSwap; 

	double          qLeft;
	double          qRight;
	vector<double>  QPivot;   // [Nall]
	bool            zeroQ;    // true if qLeft = qRight = 0   
	double          pivotRatio;

	vector<DiffusedState> VectorDiffusedStates; //	for	computing expected df [Nedf]

	int 
		m_lastLiborUsed,
		m_nbLibors, // number of state variables
		m_nbFactors; // number of stochastic driving factors

	double 
		m_calibMat,
		m_calibMatCMS,
		m_predCorrStep, // size of long step (in years)
		m_corrLower,
		m_corrUpper,
		m_lambda,
		m_lambdaInit;

	vector< const double * >
		m_gaussianVector;  // Brownian motion increments

	// ASSUMPTION: m_targetCorr[i] is correlation between swaps with 
	// indices m_corrSwapIndices[i] and m_corrSwapIndices[i+1]
	vector<int>  m_calibIndices,     // indices of swaptions used for calibration
		         m_corrSwapIndices;  

	// parameters of the functional form of
	// the factor loadings functions. All matrices have dimensions
	// nbLibors x nbFactors
    vector<double> 
		m_diffuseTime,
		m_libor,
		m_df,           // engine discount factors 
		m_liborVols,    // caplet implied vols. associated with the model's LIBORS
		m_liborSkew,
		m_resetDates,
		m_payDates,
		m_paramDates,
		m_expReset, 
		m_initialLibor,
		m_accruals,     // array of accruals associated to libors
		m_invAccruals,  // reciprocal of time between engine ZCB maturities -
		m_previousLibor,
		m_v,               // arrays used in iterative Libor update -- see
		m_drift,     // array of Libror's drifts for predictor corrector
		m_driftNew,  // same
		m_shock, // array of Libors's diffusion term for pred. cor.
		m_aInit, m_bInit, m_dInit, 
		m_a, m_b, m_d, 
		m_fudgeA, m_fudgeD,
		m_cosTheta, m_sinTheta,
		m_targetCorr,
		m_humpParams,
		m_marketExpiries, // expiries in years offset from today     		  
		m_marketTenors,   // tenors in years offset from expiries
		m_marketCalibVol,
		CholMatrix_11, CholMatrix_21, CholMatrix_31,
		CholMatrix_22, CholMatrix_32, CholMatrix_33;

	vector< vector<double> >
		m_factorVol,  // array to store effective covariances in predictor-corrector method
		m_A, m_B, m_D, 
		m_facLdArray,   // nbLibors x 3nbFactors 
		m_marketSwaptionVols,
	    m_modelSwaptionVols, 
		m_modelCorrelations, // nbLibors x nbLibors dim. array with LIBOR correlations	 
		m_instLiborVol;    // nbFactor vectors of libors instantaneous volatilities at time t

    DateTimeArray
	    m_marketExpiriesDate; // expiries, DateTime

	simulationType 
		m_simulationType;

    RanNumClassSP ranNumLiteSP;     // random number generator
    int             diffYCIdx;
    int             discYCIdx;

	const double invSqrtTwelfth; // 1/12
};

//declare smartPtr versions	
DECLARE(SRMRatesLiborDiffuse);

DRLIB_END_NAMESPACE
#endif // SRMRATESHJMDIFF_HPP
