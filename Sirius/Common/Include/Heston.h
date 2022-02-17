
#pragma once

#include "mleqobjects.h"
#include "StaticReplication.h"
#include "solve.h"
#include "complex.h"



struct HestonParameters
{	double	v0;
	double	vas;
	double	kappa;
	double	correl;
	double	vovol;
};


class heston_pricer	:	public RCObject
{
public:
	heston_pricer(HestonParameters params);
	virtual ~heston_pricer(){}

	virtual double call_price(double mat, double forward, double strike);
	virtual double implied_vol(double fwd, double mat, double strike);

	friend class heston_fitter ;

protected:	
	double m_kappa;
	double m_vas;
	double m_v0;
	double m_cov;	// rho * vovol
	double m_var;	// vovol^2
	double m_mat;

protected:
	CVector		m_gaussWeights;
	CVector		m_gaussPoints;
	int			m_nPoints;
	double		searchUpperBound();

	virtual	complex log_fft(complex c);	
};


class heston_jump_pricer		:	public heston_pricer
{
public:
	heston_jump_pricer(HestonParameters params, double lambda, double jump);
	virtual ~heston_jump_pricer(){}

	friend class heston_jump_fitter ;

protected:	
	double m_lambda;
	double m_jump;

	virtual	complex log_fft(complex c);
};


class heston_fitter	: public CStaticReplication, public LSGRGSolver, public RCObject
{
public:
	heston_fitter( MlEqAssetHandle hUdly, long dateMaturity, int nPoints = 30);
	heston_fitter( MlEqAssetHandle hudly, long dateShort, long dateMaturity, long dateLong, int nPoints = 30 );
	virtual ~heston_fitter(){}

	virtual void	calibrate( GVector<MlEqStrikeHandle> vStrikes );	
	HestonParameters* getParams(){ return &m_params;}	


//protected:
	virtual void	calibrateToVarianceSwap();
	virtual bool	calibrateToVanilla( GVector<MlEqStrikeHandle> vStrikes );

	void	checkMaturities(bool&);

protected:	
	virtual void  ObjectiveFcn(double* gVals,double* xVals);
	bool	m_calibVanilla;


	virtual void  ObjectiveFcnVS(double* gVals,double* xVals);
	CVector m_t;
	CVector m_vs;

	double f(double x){ return 0.0;}
	double df(double x){ return 0.0;}
	double ddf(double x){ return 2./(x*x);}


	virtual void  ObjectiveFcnVanilla(double* gVals,double* xVals);
	CVector m_normStrikes;
	CVector m_strikes;
	CVector m_prices;
	CVector m_bsVega;
	virtual void	initCalibrationInputs(GVector<MlEqStrikeHandle> vStrikes);
	int		m_dataSize;

	virtual double	checkFcn();

protected:
	RCPtr<heston_pricer>	m_pricer;

	MlEqAssetHandle	m_hUnderlying ;

	long	m_nToday;
	long	m_dateShort;
	long	m_dateMiddle;
	long	m_dateLong;

	HestonParameters m_params;
};


class heston_jump_fitter	:	public heston_fitter
{
public:
	heston_jump_fitter( MlEqAssetHandle hUdly, long dateMaturity, int nPoints = 30);
	heston_jump_fitter( MlEqAssetHandle hudly, long dateShort, long dateMaturity, long dateLong, int nPoints = 30 );
	virtual ~heston_jump_fitter(){}

	void	calibrate( GVector<MlEqStrikeHandle> vStrikes );	

	virtual void	calibrateJumpsToSkew(GVector<MlEqStrikeHandle> vStrikes);

	virtual void	calibrateToVarianceSwap();
	virtual bool	calibrateToVanilla( GVector<MlEqStrikeHandle> vStrikes );

protected:	
	virtual void  ObjectiveFcnVS(double* gVals,double* xVals);
	virtual void  ObjectiveFcnVanilla(double* gVals,double* xVals);
	virtual void  initCalibrationInputs(GVector<MlEqStrikeHandle> vStrikes);

	double	checkFcn();

protected:
	RCPtr<heston_jump_pricer>	m_jump_pricer;
public:
	bool	m_calibJumps;
	double	m_lambda;
	double	m_jump;
	double	m_bsVol;
	double  MertonPrice(double mat, double forward, double strike);
	CVector m_vols;
};



struct HestonParametersTermStructure
{
	GVector<HestonParameters>	hParams;
	CVector		t;
};


class heston_term_pricer	:	public heston_pricer
{
public:
	heston_term_pricer(HestonParametersTermStructure params);
	~heston_term_pricer(){}

	double call_price(double mat, double forward, double strike);

	friend class heston_term_fitter;

protected:
	complex log_fft(complex c);

	void	updateTimeSlices(double mat);
	CVector m_dt;
	CVector m_t;

	CVector m_vKappa;
	CVector m_vCov;	
	CVector m_vVar;	
	int	m_nSlice;
	int	m_nLastSlice;

	complex C(complex c );
	complex int_d(complex c );
	complex int_d(complex c, int start, int end);
};





class heston_term_fitter	:	public heston_fitter
{

public:
	heston_term_fitter( MlEqAssetHandle udly, GVector<long> dateMaturities );
	virtual ~heston_term_fitter(){}

	void	calibrate( GVector<MlEqStrikeHandle> vStrikes );
	HestonParameters* getParams(int slice){ return &((m_params.hParams)[slice]);}
	


	void	calibrateToVarianceSwap();
	bool	calibrateToVanilla( int slice, GVector<MlEqStrikeHandle> vStrikes, bool calibCorrel = true, double correl = 0. );
protected:
	virtual void  ObjectiveFcnVanilla(double* gVals,double* xVals);
	virtual void  ObjectiveFcnVS(double* gVals,double* xVals);

protected:
	RCPtr<heston_term_pricer>		 m_pricer;
	HestonParametersTermStructure	 m_params;;
	GVector<long>					 m_dateMaturities;

	CVector		m_dt;
	int			m_nSlice;
	int			m_nCalibSlice;
	bool		m_calibCorrel;
	double		m_correl;
	double		m_avgKappa;
};




class sabr_fitter	:  public LSGRGSolver, public RCObject
{
public:
	sabr_fitter( MlEqAssetHandle hUdly, long dateMaturity);
	virtual ~sabr_fitter(){}

	virtual void	calibrate( GVector<MlEqStrikeHandle> vStrikes );	
	void getParams(CVector& param);

protected:	
	virtual void  ObjectiveFcn(double* gVals,double* xVals);
	CVector m_strikes;
	CVector m_bsVols;	
	int		m_dataSize;
	void	initCalibrationInputs(GVector<MlEqStrikeHandle> vStrikes);

protected:
	MlEqAssetHandle	m_hUnderlying ;
	long			m_dateMaturity;

	double m_forward;
	double m_mat;
//	double m_beta;	// beta is 1...

	double m_vovol;
	double m_correl;
	double m_v0;
};





#include "MonteCarlo.h"

class CStochasticVolMC	:	public CSequentialBlackMC
{
public:
	CStochasticVolMC(MlEqConstDateHandle hDate):CSequentialBlackMC(hDate){};
	virtual ~CStochasticVolMC(){};

	virtual void initialize(product& deriv, long dateMaturity, int npath,int randomNumberFlag);
	virtual void simulate(CMatrix& result,product& deriv);

protected:
	RCPtr<CMersenneTwister>		m_mt;
};

class CHestonMC	:	public CStochasticVolMC
{
public:
	CHestonMC(MlEqConstDateHandle hDate, HestonParameters params):
	CStochasticVolMC(hDate)
	{
		m_v0	=	params.v0;
		m_vas	=	params.vas;
		m_kappa	=	params.kappa;
		m_vol	=	params.vovol;
	}
	~CHestonMC(){};

	virtual CMatrix&	GetPathArray(int ipath);

protected:
	double m_v0;
	double m_vol;
	double m_vas;
	double m_kappa;	// no need for correlation here...
};


class CHestonJumpMC	:	public CHestonMC
{
public:
	CHestonJumpMC(MlEqConstDateHandle hDate, HestonParameters params, double jump, double lambda):
	  CHestonMC(hDate, params),m_jump(jump), m_lambda(lambda){}

	~CHestonJumpMC(){};

	void initialize(product& deriv, long dateMaturity, int npath,int randomNumberFlag)
	{
		CStochasticVolMC::initialize( deriv, dateMaturity, npath, randomNumberFlag );
		initRandomJump( randomNumberFlag );
	}

	CMatrix&	GetPathArray(int ipath);

protected:
	double m_jump;
	double m_lambda;

	void initRandomJump( int rngFlag);

	randomGeneratorHandle m_randomGeneratorJump;
	GVector<CVector>	  m_randomJumps;
};


class CSABR_MC	:	public CStochasticVolMC
{
public:
	CSABR_MC(MlEqConstDateHandle hDate, double v0, double vol):CStochasticVolMC(hDate),m_v0(v0),m_vol(vol){}
	~CSABR_MC(){};

	CMatrix&	GetPathArray(int ipath);

protected:
	double m_v0;
	double m_vol;// no need for correlation here...
};


class VarianceCall	:	public product
{
public:
	void init( CVector strikes, double mat, long matDate);

	virtual void payout(CMatrix& value,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& m);
	virtual void setUp(CMatrix& value,MlEqMonteCarlo& mc);


protected:
	CVector m_strikes;
	double	m_mat;
};


class PayoffTest	:	public product
{
public:
	void init( const CVector& strikes, const CVector& cp, const GVector<long>& matDates);

	virtual void payout(CMatrix& value,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& m);
	virtual void setUp(CMatrix& value,MlEqMonteCarlo& mc);


protected:
	CVector m_strikes;
	CVector m_cp;
	int m_nDates;
};











//// @@@@@@@@@@@@@@@@@@@@@@@@


#include "localvol.h"
#include "MlEqPde.h"
//#include "pdeproducts.h"


class straddleHelper : public MlEqPdeHelper
{
public:
		straddleHelper(){}
		straddleHelper(double strike,double spot);
		straddleHelper(MlEqInterpolatorHandle initialValue, double spot);
		~straddleHelper(){}

		MlEqInterpolatorHandle m_initialValue;
		double m_strike;

		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);
		void initialize( double strike, double spot );
		void initialize( MlEqInterpolatorHandle initialValue, double spot);

		void finalize( MlEqPdeDriver* pde );
};





class FittedLocalVol	:	public DupireLocalVol
{
public:
	FittedLocalVol(){}
	~FittedLocalVol(){}

	void	fitSurface();
/*
	virtual void initialize(MlEqAssetHandle pAsset,
								double lowSpot,
								double highSpot,
								int nx, 
								long maturityDate,
								double vb = 0.01,
								int gridTypeFlag = 2 );
*/
	virtual void initialize(MlEqAssetHandle pAsset,
								double lowSpot,
								double highSpot,
								int nx, 
								const std::vector<long>& payoffDates,
								double vb = 0.01,
								int gridTypeFlag = 2 );

protected:
	virtual void	initialize( int slice );
	virtual void	FitSlice( int slice );
	virtual double	fitSlice( int slice, double tol );
	virtual void	updateLocalVolSlice( int slice );

	GVector<long>	m_sliceDates;

protected:
	CVector				m_nStrikes;	// normalize strikes then interpolate
	CVector				m_mktPrices;
	CVector				m_cStrikes;
	CVector				m_mktVols;
	CVector				m_bsVega;
	int					m_nFit;

protected:
	GVector<DupireLocalVolHandle>	m_sliceLocalVol;

	double	solvePde( int slice, double strike );
	RCPtr<pdeLocalVol>		m_pde;
	RCPtr<straddleHelper>	m_straddle;	// use straddle...

	

//	void	buildTimeGrid(long startDate, long endDate, int nt, GVector<long>& timeGrid);	
	void	buildTimeGrid(long startDate, long endDate, int nt, const std::vector<long>& payoffDates, GVector<long>& timeGrid);

	int m_nSlice;
	int m_firstSlice;

protected:
//	void	createTimeGrid( long maturityDate );
	void	createTimeGrid( const std::vector<long>& payoffDates );
	virtual void	initInterpolationGrid( double stdDev = 0., int nx = 0);

protected:
	virtual void	initMarketData();
	virtual void	initInterpolator();
	virtual void	updateInterpolator( double bump );

	MlEqInterpolatorHandle m_curveBump;
};


class FittedLocalVol2	:	public FittedLocalVol
{
public:
	FittedLocalVol2(){}
	~FittedLocalVol2(){}

protected:
	void	initialize( int slice );
	double	fitSlice( int slice, double tol );
	void	updateLocalVolSlice( int slice );	
	void	initInterpolator();
	void	updateInterpolator(  const CVector& dval  );
	void	initMarketData();

protected:
//	CVector m_xdata;

	MlEqInterpolatorHandle m_leftBump;
	MlEqInterpolatorHandle m_rightBump;
	CVector m_lxdata;
	CVector m_rxdata;

	int m_nFitLeft;
	int m_nFitRight;

	void testFitter(int slice);
};










