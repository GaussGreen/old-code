

#ifndef MLEQSHORTMODELS
#define MLEQSHORTMODELS

#include "cMatrix.h"

#include "MlEqHandles.h"
#include "solve.h"
#include "mleqobjects.h"

//#include "MonteCarlo.h"



using namespace std;


double impliedLogContractVol(MlEqAsset& asset,long nMaturity,int ngauss=50,double upper = 5.0,double lower = -5.0);

double expectedVol(MlEqAssetHandle& asset,long nMaturity,int ngauss,double upper,double lower);
double PatHaagenImpliedVol(double strike,double mat,double fwd,double correl,double beta,double volvol,double shortvol);
double PatHaagenImpliedVol(double strike,double mat,double fwd,double correl,double volvol,double shortvol);
double BlackScholesSwaptionPricer(long valuationDate, GVector<long>& SwapFixingDates,double swaptionStrike,double vol,DayCountConventionEnum & dayCount ,MlEqZeroCurve&	Z, int callPut);
double BlackScholesCapletPricer(long valuationDate, long RateSettingDate,long RateMaturityDate,double strike,double vol,DayCountConventionEnum & dayCount,MlEqZeroCurve&	discountCurve,int callput);
void AmericanDownInPut(double& price,CMatrix& AmPrices,int method,double DownBarrier,
					   double spot,double mat,double r,double q,
					   double SpotBarrierVols,double SpotBarrierVolSlope,
					   double FwdBarrierVols,double FwdBarrierVolsSlope,
					   bool PayAtTheEnd,
					   int npoints,int nplotpoints,bool bplotCummulativeHitProb=true);

void BootstrapStoppingTimes(double & price,CMatrix& HitProb,					 
					 double DownBarrier,double spot,double mat,
					 double SpotBarrierVols,double drift,double r,
					 double SpotBarrierVolSlope,
					 double FwdBarrierVols,
					 double FwdBarrierVolsSlope,bool payAtTheEnd,
					 int npoints,bool plotCummulativeHitProb);

void BootstrapStoppingTimes2(double & price,CMatrix& HitProb,					 
							 double DownBarrier,double spot,double mat,
							 double SpotBarrierVols,double drift,double r,
							 double SpotBarrierVolSlope,
							 double FwdBarrierVols,
							 double FwdBarrierVolsSlope,
							 double BarrierStrikeVols,
							 double strike,
							 double cp,
							 double rebate,
							 bool payAtTheEnd,
							 int npoints,bool plotCummulativeHitProb);

void HermiteOptionsNew(CVector& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,int ntimesteps,int npaths,int nspacialPoints,int nstdev);
void HermiteOptionsNew(CVector& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,int ntimesteps,int npaths,int nspacialPoints,int nstdev,randomGeneratorHandle& RandomGenerator);
void VolOptions(CMatrix& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,int ntimesteps,int npaths,int nspacialPoints,int nstdev);




class calibrateHermiteOptions : public LSGRGSolver
{

	protected:


	void  ConstraintFcn(double* contraintVals, double* xVals);
	void  ObjectiveFcn(double* gVals,double* xVals);

	randomGeneratorHandle		m_randomGenerator;
	CVector						m_fixedStrikes;//[nstrikes]
	GVector< MlEqStrikeHandle > m_calibStrikes;//[nstrikes]
	CVector						m_calibVols;//[nstrikes]
	CVector						m_calibPrices;//[nstrikes]
	CVector						m_calibWeights;//[nstrikes]
	CVector						m_varSwapVol;//size: 0 or 1

	MlEqDateHandle				m_Today;
	long						m_maturityDate;
	int							m_ntimesteps;
	int							m_npaths;
	int							m_nspacialPoints;
	int							m_nstdev;
	CVector						m_hermitCoeffGuess;
	CVector						m_currentHermitCoeff;
	double						m_maturity;
	double						m_forward;
	

	CVector						m_gaussPoint;
	CVector						m_gaussWeight;
	CVector						m_spacialGaussPoint;
	CVector						m_spacialGaussWeight;

//	GMatrix<CMatrix>			m_cashedVar;
	double****					m_cashedVar;// don't do this at home

	void						HermiteOptionsNew(CVector& result);

public:

	void						calibrate(CVector& x);

	void						initialize(	GVector< MlEqStrikeHandle >& calibStrikes,
											CVector& calibVols,double	forward,
											MlEqDateHandle	Today,
											long				maturityDate,
											CVector&			varSwapVol,
											int				ntimesteps,
											int				npaths,
											int				nspacialPoints,
											int				nstdev,
											CVector&		hermitCoeffGuess);

	~calibrateHermiteOptions();

};


class blackScholesBarrier
{

	double I(double n,double x);
	double J(double n,double x);

	double m_d;
	double m_lp;
	double m_lm;
	double m_la;
	double m_lap;
	double m_k;
	double m_mu;

public:

	double m_strike;
	double m_callPut;
	double m_upperBarrier;
	double m_lowerBarrier;
	double m_maturity;
	double m_vol;
	double m_forward;
	double m_spot;
	double m_rate;

	void initialize(
					double strike,
					double callPut,
					double upperBarrier,
					double lowerBarrier,
					double maturity,
					double vol,
					double forward,
					double spot,
					double rate
					);

	double price();
};





/***************************************************************
**	Class   : fitLocalVol
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/


class calibrateEffectiveLocalVol;

class fitLocalVol: public LSGRGSolver
{
	
	calibrateEffectiveLocalVol*	m_cLV;

	CVector m_currentVals;
	CVector m_currentVols;
	CVector m_weight;

	public:


	CVector m_coeff;


	void  initialize(calibrateEffectiveLocalVol& calibLV,CVector& initialGuess,	CMatrix& data,
						   CMatrix& xValBounds,CMatrix& ObjectiveBounds,
						   double InitialTolerance,double FinalTolerance,
						   double StoppingTolerance,
						   CMatrix& NonZeroPartialDerivativesSpecification,int outputFlag);

//	void  PartialDerivatives(CVector& xVal,CVector& partialDerivatives){};
	void  ObjectiveFcn(double* gVals,double* xVals);
};



/***************************************************************
**	Class   : calibrateEffectiveLocalVol
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class calibrateEffectiveLocalVol
{

	friend class fitLocalVol;

protected:

	MlEqAssetHandle  m_asset;

protected:

	fitLocalVol m_fitLV;

protected:

	MlEqAnalyticCurveWithTanhWingEdgeHandle m_quadraticLV;

protected:

	double	m_shortFwd;
	double	m_longFwd;
	double	m_deltaT;
	double  m_sqdT;
	double  m_shortMat;
	double  m_sqShortMat;
	double  m_longMat;
	CVector m_fixedCalibStrikes;//[icalib]

	CVector m_calibSpots;
	CVector m_guess;
	CVector m_targetPrices;
	CVector m_targetvols;
	CVector m_vols;//m_nGauss

	MlEqInterpolatorHandle m_skewMultipliers;

	int		m_ngauss;
	CVector m_gaussPoints;
	CVector m_gaussWeights;
	double  m_volatm;
	double  m_fwdvol;
	CVector m_currentVols;
	CMatrix m_results;

	virtual void calcOptions(CVector& prices, const CVector& vols,CVector& fixedStrikes);

	double getLocalVol(double xspot);

public:


	virtual void calcOptions(CVector& prices, const CVector& vols,GVector< MlEqStrikeHandle >& strikes);



	void init(MlEqAsset& asset,MlEqDateHandle startDate,long endDate,
							   GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,CVector& targetvols,
							   CVector& guess,int ngauss,CVector& tanhWingInfo,double finalTol,double stoppingTol,bool initSkewMultsOpt=true);


	void init(MlEqAsset& asset,MlEqAnalyticCurveWithTanhWingEdgeHandle& quadraticInterp,MlEqDateHandle startDate,long endDate,
							   GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,CVector& targetvols,
							   CVector& guess,int ngauss,bool initSkewMultsOpt=true);

	void calibrateLV(CMatrix& res,CVector& initialGuess);

	MlEqDateHandle m_startDate;
	long		   m_endDate;

	GVector< MlEqStrikeHandle > m_calibStrikes;//[icalib]

};


void calibrateEffectiveLocalVolFull(CMatrix& result,MlEqAsset& asset, GVector <MlEqAnalyticCurveWithTanhWingEdgeHandle >& quadraticInterpolator,
									MlEqDateHandle startDate,GVector< long > endDates,CMatrix& guess,int ngauss,bool usePrevRes,
									GVector< GVector< MlEqStrikeHandle > >& calibStrikes,GVector< GVector< MlEqStrikeHandle > >& calibSpots);


/***************************************************************
**	Class   : calibrateEffectiveLocalVol
**	Function: 
**	Returns : 
**	Comment : 
****************************************************************/

class calibrateEffectiveLocalVolGrid : public calibrateEffectiveLocalVol
{

protected:

	pdeLocalVolEffectiveHandle m_pdeDriver;

	virtual void calcOptions(CVector& prices, const CVector& vols,CVector& fixedStrikes);


	long m_idate;

public:

	virtual void calcOptions(CVector& prices, const CVector& vols,GVector< MlEqStrikeHandle >& strikes);

	void init(int idate,pdeLocalVolEffectiveHandle pdeDriver,
			  MlEqAsset& asset,MlEqDateHandle startDate,long endDate,
			  GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,
			  CVector& targetvols,CVector& tanhWingInfo,double finalTol,double stoppingTol);


	void init(int nt,pdeLocalVolEffectiveHandle pdeDriver,MlEqAsset& asset,MlEqAnalyticCurveWithTanhWingEdgeHandle& quadraticInterp,MlEqDateHandle startDate,long endDate,
							   GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,CVector& targetvols);


};




void calibrateEffectiveLocalVolGridFull(CMatrix& result,MlEqAsset& asset, GVector <MlEqAnalyticCurveWithTanhWingEdgeHandle >& quadraticInterpolator,
									MlEqDateHandle startDate,GVector< long > endDates,CMatrix& guess,bool usePrevRes,
									GVector< GVector< MlEqStrikeHandle > >& calibStrikes,GVector< GVector< MlEqStrikeHandle > >& calibSpots,
									double lowSpot,double highSpot,int nx,int nt);

#endif

