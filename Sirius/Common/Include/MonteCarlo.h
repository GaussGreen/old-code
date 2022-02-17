


#ifndef MONTECARLOH
#define MONTECARLOH

#include <vector>

#include "CMatrix.h"
#include "utility.h"
#include "mleqdate.h"
#include "mleqshortmodels.h"
#include "edcubicspline.h"
#include "Rnd.h"

class MlEqZeroCurve;
class MlEqAsset;

#include "MlEqMaths.h"
#include "mersennetwister.h"



//	 Alex - this is the object that interfaces with handle_MonteCarlo
class CForwardSkewMC;
class CMultiAssetForwardSkewMC;
class CStochIRForwardSkewMC;
class CQuasiMC;
class CHermiteMC;
class CLocalVolMC;
class CMultiLocalVolMC;


class CMonteCarloHelper
{

public:
	const CMonteCarloHelper&			operator=(const CMonteCarloHelper& mch);


	void initialize(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday);
	void initialize(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday, const std::vector<long>& payoffDates);

	CMonteCarloHelper(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday);
	CMonteCarloHelper(){};

	const std::vector<long>&			GetMCDates(void) ;
	long								GetMCDate(int idate) ;
	const std::vector<long>&			GetPayoffDates(void) const;
	bool								CreatePayoffPath(const CMatrix& MCSpots, CMatrix& PayoffSpots);
	bool								CreatePayoffPath(const CMatrix& MCSpots, CMatrix& PayoffSpots, CVector& payoutPathDiscounts,CVector& simulatedDiscounts);	

	void postInitialize(std::vector<long>&	PayoffDates);
	void setPayoffDates(std::vector<long>&	PayoffDates);

	int mapFromPayoff(int idate);

protected:

	void								mapFromPayoff(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday,const  std::vector<long>& payoffDates);
	std::vector<int>					m_mapFromPayoff;
	std::vector<long>					m_PayoffDates;
	bool								m_bUseToday;					// true if nToday is a function of the payoff but we need to simulate	it
	CMatrix								m_vectorFixedSpots;				// [iasset][ifixingdate]	
	std::vector<long>					m_vectorMCDates;				// these are the dates that need to be simulated by the Monte-Carlo process
};





void getCorrelationMatrix( const std::vector<MlEqAssetHandle>& assets, CMatrix& correl );


////////////////////////////////////////////////////////////////////////////////////
//
//
//				Generic Monte Carlo
//
//
////////////////////////////////////////////////////////////////////////////////////



class MlEqMonteCarlo : public RCObject
{
public:
	explicit MlEqMonteCarlo(MlEqConstDateHandle hDate) : m_hDate(hDate), m_randomsAreExternal(0), m_centerFwds(false) {}	
	
protected:

	CMonteCarloHelper m_mMCHelper;
	virtual void createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray);
	virtual void createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray,CVector& payoutPathDiscounts,CVector& simulatedDiscounts);

	//void initMC(MlEqConstDateHandle hDate);
	void initMC(const vector<long>& simDates);

	void					initializeRandomGenerator(int dimensionalityFactor);
	virtual void			GenerateRandomNumbers(int idate);
	virtual void			GenerateRandomNumbers(Cholesky& cholesky);

	randomGeneratorHandle	m_randomGenerator;
	int						m_randomNumberFlag;
	GMatrix < CVector >		m_allMultiRandoms;//[ipath][iasset][m_numberOfFactors*idate]
	CMatrix					m_allRandoms;//[ipath]m_numberOfFactors*idate]
	CVector					m_oneStepRandoms;//[m_numberOfFactors*ipath]
	int						m_offset;
	GVector < int >			m_numberOfFactorsPerAsset;//[m_nAssets]
	int						m_seed;
	int						m_iasset;
	GMatrix < CVector >*	m_pallMultiRandoms;//[ipath][iasset][m_numberOfFactors*idate]
	int						m_randomsAreExternal;
	cTensor*				m_pPath_array;//[ipath][iasset][idate]  
	cTensor					m_path_array;//[ipath][iasset][idate]

	CVector m_discount;//[idate]
	CVector m_discountPrepended;//[idate+1] this includes a discount factor of 1 for valuation

	CVector m_t;//[idate]
	CVector m_dt;//[idate]
	CVector m_sqrtdt;//[idate]

	CMatrix m_termFwds;// [idate][iasset]

	bool	m_centerFwds;
	void	centerForwards(int idate,int iasset,double fwd);
	int		m_currentPath;

	vector<MlEqAssetHandle>			m_assets;	// [iasset]	
	int		findAsset(MlEqAssetHandle & asset);

protected:

	double	GetNextRandom(int& iRandom,int& noffset);
	CMatrix m_fwds;// [idate][iasset]
	virtual int GetNumberOfVolStates(int ipath,int idate,int iasset=0){return 1;};


	void basicInitializeMcHelper(product& deriv,CMatrix& fixings);



public:
	std::string			GetName(void) const;
	
	virtual double AmericanDigital(CVector& prices,int ipath,int idate,double DownBarrier,double finalSpot,
									 int iasset,int method,double stressForwardBarrierVol,
									 double stressForwardBarrierVolSlope,int npoints,double blend);

	virtual double		GetBrownianBridgeHitProb(int ipath,int idate,double U,double L);

	const vector<MlEqAssetHandle>&	getAssets() const{return m_assets;};

	void				setExternalRandoms(int iasset,GMatrix < CVector >&	m_pallMultiRandoms,cTensor& PathArray);
	void				initializeMcHelper(product& deriv,CMatrix& fixings);
	void				initializeMcHelper(product& deriv,CMatrix& fixings,const vector<long >& additionalSimDates);


	const				MlEqConstDateHandle	m_hDate;
	int					m_nPaths;
	int					m_nDates;// number of timesteps that need to be simulated in the monte carlo
	int					m_nAssets;
	double				m_numberOfFactors;

	virtual long		GetNumberOfFutureDates();


	virtual void		averageMonteCarloResults(CVector& results);
	virtual void		averageMonteCarloResults(double& result);
	void				completePayoff(double& payoff,int ipath);

	virtual void		simulate(CMatrix& result,product& deriv);
	virtual CMatrix&	GetPathArray(int ipath);
	virtual void		getPath(int ipath,CMatrix& path_array){};

	virtual double		GetDiscount(int ipath,int idate);
	virtual double		GetDiscount(int ipath,int idate,int idateAsOf);
	virtual double		GetDiscountToDate(int ipath,long To,int idateAsOf);

	virtual double		GetDiscount(int ipath,int idate,int jdate,int idateAsOf);
	virtual double		GetStochasticDiscount(int ipath,int idate);

	virtual const		CVector& GetDiscounts(int ipath);
	virtual void		generatePaths(){};

	double				GetPathValue(int ipath,int idate,int iasset=0);
	double				Getdt(int idate);
	double				GetTime(int idate);
	double				GetTimeAtPayoffDate(int ipayoffDate);


	int					GetCurrentPath(){return m_currentPath;}

	virtual double		GetBridgeVolatility(int ipath,int idate,double finalSpot,int iasset = 0);
	virtual double		GetForwardValue(int idate,int iasset=0){return m_termFwds[idate][iasset];};

	virtual double		GetImpliedLogContractVol(int ipath,int idate,int iasset,int ngauss,double lowerStdev,double upperStdev);
	
	virtual double		GetForwardForward(int idate,int iasset){return m_fwds[idate][iasset];};
	virtual void		calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1.0,int iasset=0);
	virtual void		calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1.0);

	double				GetImpliedVol(int iStart,int iEnd,double optionPrice,double strike,int iasset = 0);


	void				PutPathValue(int ipath,int idate, int iasset, double fValue);

	const std::vector<long>&			GetMCDates(void) ;
	long								GetMCDate(int idate) ;
	
private:

	CVector								m_vectorDummy;
};




class CSequentialBlackMC :public MlEqMonteCarlo
{
protected:
	GVector < int > m_numberOfFactors;//[iasset]
	GVector < CVector > m_randoms;//[iasset][istep*numberOfFactors[iasset]+jfactor]
	CVector  m_rndTemp; //[iasset]
	CMatrix  m_vols;//[iasset][idate]

	Cholesky m_cholesky;
	CMatrix  m_currentPathValues;//[iasset][idate]

	CVector m_spots;//[iasset]
	CMatrix m_fwdfwds;//[iasset][idate]

public:
	
	virtual CMatrix&	GetPathArray(int ipath);

	CSequentialBlackMC(MlEqConstDateHandle hDate) : MlEqMonteCarlo(hDate) {}

	void initialize(const vector<MlEqAssetHandle>& assets, product& deriv, CMatrix& fixings, CMatrix& correl,int npath,int randomNumberFlag);
};

// Quasi-Monte Carlo. This gives you one path with all
// the values defined explicitly. This is useful if you
// want to keep a set of defined paths that can be used
// in exactly the same manner as a genuine Monte Carlo
// object.
class CQuasiMC : public MlEqMonteCarlo
{
public:
	CQuasiMC(MlEqConstDateHandle hDate) : MlEqMonteCarlo(hDate) {}
	void Initialize(std::vector<MlEqSpotScheduleHandle>& ahSpotSchedules, const std::vector<long>& anDates);
};


class HermiteCoeff;

class CHermiteMC : public CSequentialBlackMC
{

protected:

	GVector < HermiteCoeff >  m_HermiteCoeff; //[idate]
	GMatrix < CVector >  m_b;//[idate][ihermite][jhermite]

	void createTransitionCoeff();
	int m_hermite;
	int m_ngauss;
	CVector m_gaussPoints;//[m_ngauss]
	CVector m_gaussWeights;//[m_ngauss]
	double  m_sq2pi;

	double normalizeWiener(int idate,CVector& prevHermite);
	double Wiener(CVector& currentHermite,int idate,double omega,CVector& prevHermite);


public:

	void initialize(GVector < HermiteCoeff > & HermiteCoeff,const vector<MlEqAssetHandle>& assets, product& deriv, CMatrix& fixings, CMatrix& correl,int npath,int randomNumberFlag,int ngauss,int nstdev);
	void initialize(const vector<MlEqAssetHandle>& assets, product& deriv, CMatrix& fixings, CMatrix& correl,int npath,int randomNumberFlag,int ngauss,int nstdev);

	CMatrix&	GetPathArray(int ipath);


};


////////////////////////////////////////////////////////////////////////////////////
//
//
//				Forward Skew Monte Carlo
//
//
////////////////////////////////////////////////////////////////////////////////////


class CForwardSkewMC: public MlEqMonteCarlo
{

	friend class CMultiAssetForwardSkewMC;	

	protected:


	CVector					m_currentVols;//[ipath]

	int getSlope(double& slope,double & smile,double& maxDiff,double& minDiff,CVector& impliedVols,CVector& termVols,int atmIndex,CVector& accuracy)	;
	int updateSlope(CVector& impliedVols,CVector& termVols,int atmIndex,CVector& accuracy)	;

	void GlobalCalibrate(CMatrix& TermSkews,vector<vector<MlEqStrikeHandle> >& calibStrikes,
					 long& idum,CVector& accuracy);// calibrates paths from idate->idate+1

	void calibrateTermsSkewsGlobal(CMatrix& TermSkews,	vector<vector<MlEqStrikeHandle> > & Strikes,
								   iVector& calibDateIndex,long& idum,int shortIndex,int longIndex,CVector& accuracy);// calibrates paths from idate->idate+1

	int calibrateTimeSlice(int idateCalibIndex,CVector& TermSkews,	vector< MlEqStrikeHandle > & Strikes,
					 long& idum,CVector& accuracy,int calibGlobalFlag,CMatrix& calibrationInfo,int maxiter=12,int succeedType=1,int istartDate=-1);// calibrates paths from idate->idate+1

	CMatrix								m_calibVols; // [idate][istrike]
	vector<vector<MlEqStrikeHandle> >	m_pCalibStrikes; // [idate][istrike]
	iVector								m_calibTimeIndex;//[idate]
	CMatrix								m_calibrationInfo;//[m_calibTimeIndex][icalibVariable]

	CVector								m_volvolMonitor;//[2]

	int		m_globalCalib;
	double	m_rev0;	
	double	m_rev1;
	double	m_xslope0;
	double 	m_xslope1;
	int		m_mrUpdateOnly;

	CVector* m_pcurrentVols;
	CVector* m_pnewcurrentVols;
	CVector  m_newcurrentVols;//[ipath]
	CVector* m_psavecurrentVols;

	CMatrix m_vols; // [idate][istrike]


	CMatrix m_normalizedStrikes;//[idate][istrike]
	CVector m_beta;//[idate]
	CVector m_meanReversionRate;//[idate]
	CVector m_meanReversionLevel;//[idate]
	CVector m_volvol;//[idate]

	vector<EDCubicSpline> m_ForwardVolcurve;//[idate]
	vector<EDCubicSpline> m_bareForwardVolcurve;//[idate]

	void	multiplyForwardVols(int idate,double factor);


	int		m_numberVolStates;
	double	m_lowestNd2;
	double	m_highestNd2;
	double  m_bareVolAtm;
	int		m_localVolFlag;
	int		m_saveVolofVolInfo;
	double	m_initialSpot;
	double  m_currentVariance;
	double  m_currentAverage;
	double  m_volmultiplyer;
	int		m_numberGridpoints;
	double  m_initialShortVol;
	int		m_calibflag;
	CVector m_accuracy;
	CMatrix m_guessGrid;//[iGridpoint][ivolstate]
	CMatrix m_omegaGridPoints;//[m_numberVolStates][numberGridpoints];
	CVector	m_discreteVol;//[iGridpoint]
	CVector m_discreteVolProb;//[iGridpoint]
	CVector m_VolAtm;//[idate]
//	cTensor m_path_array;//[ipath][iasset][idate]sos
	CVector m_prevVol;//[ipath]
	CMatrix m_saveDiscreteVols;//[ndates][iGridpoint] 
	CVector m_spotvolvol;//[ndates]
	CVector m_averageBeta;//[ndates]
	CVector m_discreteVolArray;//[ndates]
	CVector m_termvolvol;//[ndates]

	GVector<MlEqInterpolatorHandle> m_excessVolVolInterp;//[idate]

	CMatrix m_spotGrid;//[numberGridpoints][m_numberVolStates]
	CVector m_logFwdFactors;//[ndates]
	iMatrix m_volMap;//[ndates][npaths]  saves map into m_discreteVol
	CMatrix	m_discreteVolMatrix;//[ndates][m_numberVolStates] saves m_discreteVol for all times slices
	CMatrix m_continuousVolMatrix;//[ndates][npaths]
	int		m_saveVolGrid;


	double	getNewVolSkew(EDCubicSpline& baseVolcurve,double xStrike,double newVolAtm,double baseVolAtm,double excessVolVol);
	double	getExcessVolOfVol(int idate,double normStrike);
	void	getLocalVolFactor(double& localVolFactor,double& beta,int idate,int ipath,double volold);		
	void	generatePath(int idate,int calibflag,long& idum);// generates paths from idate->idate+1

	double	getCummulativeProb(double logK,double currentAtmVol,int idate);
	void	initTable(int numberGridpoints,int idate);
	double	getNewLogSpotIncrement(double omega,int volindex,int idate,int ipath);
	void	saveVolofVolInformation(int idate,int ipath,double vol,double& termvolvol,double&z,double& spotvolvol,double&z1);
	void	calculateVolofVolInformation(int idate,double termvolvol,double z,double spotvolvol,double z1);
	double	getGuess(double omega,int volindex);
	void	calibrate(int idate,CVector& currentAtmVols,double factor,int calibflag);
	void	calibrateHelper(CVector& values,CVector& simulatedVols,int idate,double volmultiplyer,CVector& strikes,int calibflag);
	double	calibrateFcn(double strike,int idate);



	virtual void initializeNumberOfFactors();
	void	GeneratePath(int calibflag);
	void	GeneratePath(int calibflag,long idum,CVector& accuracy);

	virtual double getlogFwdFactors(int idate,int ipath);
	virtual void setUpNextTimeStep(int idateCurrent){};
	virtual void initIRStateVars(int idate){};
	
	void generateControlVariatePath(int idate,int calibflag,long& idum);// generates paths from idate->idate+1
	int									m_controlVariate;
	int									m_controlVariateCalcFlag;
	CVector								m_controlVariateCurrentSpot;//[ipath]
	CVector								m_controlVariateVols;//[m_ndates]
	CMatrix								m_controlVariatePrices;
	
protected:

	GVector< hermiteSkew >				m_hermiteSkew;//[idate]
	CVector								m_hermiteElasticity;//[idate]
	GVector< hermiteSkew >				m_currentHermiteSkews;//[m_numberVolStates]
	long								m_useHermites;

	void		scaleHermiteSkew(CVector& hermiteCoeff, double currentVol,double elasticity);
	void		initHermiteTable(int numberGridpoints,int idate);

	void		setupHermites(CVector& volatm,CVector& targetSlopes);
	void		setupHermites(MlEqAsset& asset,const std::vector<long>& mcDates,CVector& strikeVol );
	

protected:

	CVector m_gaussPoint;
	CVector m_gaussWeight;// integration points for logcontracxt calibration
	GVector < CVector >	m_logContractVols;


public:
	
	explicit CForwardSkewMC(MlEqConstDateHandle hDate) : MlEqMonteCarlo(hDate) {};
	explicit CForwardSkewMC(MlEqConstDateHandle hDate, MlEqAsset& asset,product& deriv,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,const std::vector<long>& mcDates,MlEqStochBetaVolatilityStructure& betaVolStruct,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,const std::vector<long>& mcDates,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);

	void		generatePaths();	

	virtual void		calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1.0,int iasset=0);
	virtual void		calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1);


	
	void Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CVector& beta,const CVector& volvol,const CVector& meanRevRate,const CVector& meanRevLevel,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEq2DInterpolatorHandle pAccessVolInterp=NULL);
	void Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CMatrix& modelMarketParameters,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEq2DInterpolatorHandle pAccessVolInterp=NULL);
	void Initialize(MlEqAsset& asset,product& deriv,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(MlEqAsset& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);

	void Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);


	void Initialize(MlEqAsset& asset,const std::vector<long>& mcDates,MlEqStochBetaVolatilityStructure& betaVolStruct,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(MlEqAsset& asset,const std::vector<long>& mcDates,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);


	
	double	GetBrownianBridgeHitProb(int ipath,int idate,double U,double L);
	double AmericanDigital(CVector& prices,int ipath,int idate,double DownBarrier,double finalSpot,
									 int iasset,int method,double stressForwardBarrierVol,
									 double stressForwardBarrierVolSlope,int npoints,double blend);

//	virtual const CVector& GetDiscounts(int ipath);
	
	double	GetSimulatedVols(int ipath,int idate);
	double	GetSimulatedVolSquareDt(int ipath,int iStart,int iEnd,int normalize=1);
	double  GetBridgeVolatility(int ipath,int idate,double strike,int iasset=0);
	double	GetBridgeVolatility(double newVolAtm,int idate,double finalSpot,double initialSpot,int iasset);

	double	GetImpliedLogContractVol(int ipath,int idate,int iasset,int ngauss,double lowerStdev,double upperStdev);

	double	GetDt(int idate);
	int GetNumberOfVolStates(int ipath,int idate,int iasset=0);

	void	CalculateForwardVols(CVector& impliedVols,const CVector& Strikes,int iStart,int iEnd);
	
	CMatrix m_calibResults;
	
	
	~CForwardSkewMC();	
	
	const CVector&						GetBetas(void)					{return m_beta;}
	const CVector&						GetVolVols(void)				{return m_volvol;}
	const CVector&						GetMeanReversionRates(void)		{return m_meanReversionRate;}
	const CVector&						GetMeanReversionLevels(void)	{return m_meanReversionLevel;}
	int									GetSeed(void)					{return m_seed;}
	int									GetNumberVolStates(void)		{return m_numberVolStates;}
	int									GetLocalVolFlag(void)			{return m_localVolFlag;}
	int									GetSaveVolofVolInfo(void)		{return m_saveVolofVolInfo;}
	int									GetNumberGridPoints(void)		{return m_numberGridpoints;}
	int									GetSaveVolGrid(void)			{return m_saveVolGrid;}
	int									GetRandomNumberFlag(void)		{return m_randomNumberFlag;}
	int									GetControlVariate(void)			{return m_controlVariate;}
	int									GetGlobalCalib(void)			{return m_globalCalib;}
	const iVector&						GetCalibTimeIndex(void)			{return m_calibTimeIndex;}
	const CMatrix&						GetCalibVols(void)				{return m_calibVols;}
	const CMatrix&						GetControlVariatePrices(void)	{return m_controlVariatePrices;}

};	
	
	


////////////////////////////////////////////////////////////////////////////////////
//
//
//				Forward Skew Monte Carlo with stochastic interest rates
//
//
////////////////////////////////////////////////////////////////////////////////////



static bool  swaptionRoot(double x,void* vp,double* f);
static int swaptionCalibRoot (const CVector& x, void* vp, int calc_for_grad, CVector& f);

class CStochIRForwardSkewMC;
	
class MlEqHullAndWhite : public RCObject
{	
	friend class CStochIRForwardSkewMC;	
	friend bool  swaptionRoot(double x,void* vp,double* f);
	friend int	 swaptionCalibRoot (const CVector& x, void* vp, int calc_for_grad, CVector& f);
	friend bool	 capletRoot(double x,void* vp,double* f);	

public:
	explicit MlEqHullAndWhite(MlEqConstDateHandle hDate) : m_hDate(hDate) {}
		
protected:
		
	GVector<long>			m_calibrationDates;//[nCalibDates]
	CVector					m_times;//[nCalibDates]
	MlEqConstDateHandle			m_hDate;
	CVector					m_gamma_ti;//[nCalibDates]
	CVector					m_beta_ti;//[nCalibDates]
	double					m_lambda;
	MlEqZeroCurveHandle		m_discountCurve;
	
	double B(double t,double T,int indexStart,int indexEnd);
	double B(int& indexStart,int& indexEnd,long t,long T);
	double integral_g(double t,double T,int indexStart,int indexEnd);
	double integral_g(int& indexStart,int& indexEnd,long t,long T);
	double integral_g2(double t,double T,int indexStart,int indexEnd);
	double integral_g2(int& indexStart,int& indexEnd,long t,long T);
	void Locate(int& indexStart,int& indexEnd,long tDate,long TDate);


//	double P(long t, long T,double stateVar,int indexStart,int indexEnd);// this function uses indexStart and indexEnd
//	double P(int& indexStart,int& indexEnd,long t, long T,double stateVar);// this function also returns indexStart and IndexEnd

	double P(int& indexStart,int& indexEnd,double& fwdDisc,double& B,double& var,
						   long t, long T,double stateVar);// this function also returns indexStart and indexEnd B,var, fwdDisc

	double P(long t, long T,double stateVar,int indexStart,int indexEnd,double fwdDisc,double B,double var);
								// this function uses returns indexStart and indexEnd B,var, fwdDisc

	double CapletPrice(int icalibDate,double strike);
	double SwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,DayCountConventionEnum & dayCount,double integ_g2 ,int callPut=1);

	void calibrateToCaplets(CVector& capletVols,CVector& capletStrikes,DayCountConventionEnum & dayCount);
	void calibrateToSwaptions(CVector& swaptionVols,CVector& swaptionStrikes,GVector < GVector< long > > & swaptionDates,CVector& capletVols,CVector& capletStrikes,DayCountConventionEnum & dayCount);

	double CapletPrice(double B, long RateSettingDate,long RateMaturityDate,double strike,DayCountConventionEnum & dayCount );
	
public:
	
	const MlEqZeroCurveHandle getYieldCurve()const{return m_discountCurve;};

	void   initialize(GVector<long>& Dates,CVector& beta,CVector& gamma,double lambda,MlEqZeroCurveHandle discountCurve, CVector& capletVols,CVector& capletStrikes,CVector& swaptionVols,CVector& swaptionStrikes,CMatrix& swaptionDates,DayCountConventionEnum & dayCount);
	
	double Z(long t,long T);
	double P(long t, long T,double stateVar);// zero bond under Q_t
	double B(long t,long T);
	double integral_g(long t,long T);
	double integral_g2(long t,long T);

	double CapletPrice(long RateSettingDate,long RateMaturityDate,double strike,DayCountConventionEnum & dayCount);
	double SwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,DayCountConventionEnum & dayCount,int callPut = 1 );
	double BlackScholesSwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,double vol,DayCountConventionEnum & dayCount ,int callPut);
	double BlackScholesCapletPrice(long RateSettingDate,long RateMaturityDate,double strike,double vol,DayCountConventionEnum & dayCount,int callput);
	
	void   initialize(CVector& calibrationDates, long rateSetting,long rateMaturity);
	const  CVector& getBetas()const{return m_beta_ti;};
	const  CVector& getGammas()const{return m_gamma_ti;} ;
	
};
	
class swaptionRootDataContainer
{	
public:
	
	MlEqHullAndWhiteHandle pHullWhite;  
	double swaptionStrike;
	GVector<long> SwapFixingDates;
	DayCountConventionEnum  dayCount ;
	
	CVector swaptionVals;
	CVector swaptionStrikes;
	GVector< GVector< long > > AllSwapFixingDates;
	CVector capletVols;
	CVector capletStrikes;
	
	long RateSettingDate;	
	long RateMaturityDate;
	double capletStrike;
	double price;	
};	
	
	
class CStochIRForwardSkewMC: public CForwardSkewMC
{	
protected:

	bool	m_saveStochasticDiscounts;
	CMatrix m_stochasticDiscounts;//[ipath][idate]

protected:
	
	MlEqHullAndWhiteHandle	pHullWhite;  
	MlEqHullAndWhiteHandle	m_pPayHullWhite;  // do this later
		
	CVector			m_stateVarVols;//[m_ndates] vol of statevar between t[i] and t[i+1]
	CVector			m_stochasticDiscountsAlongPath;//[m_ndates]
	CVector			m_oneStepDiscounts;//[m_nPaths]
	CMatrix			m_StateVar;//[ipaths][idate] holds S_t{idate+1}
	double			m_corr;
	double			m_xcorr;
	
	void			initializeNumberOfFactors();
	void			initIR(double fCorrelation,bool isMinimalHybrid);
	void			createStateVariables(int idate);
	void			createStateVariables();
	void			calculateOneStepDiscounts(int idate);
	virtual double	getlogFwdFactors(int idate,int ipath);
	virtual void	setUpNextTimeStep(int idateCurrent);
	double			getStateVar(int ipath,int idate);
	virtual void	initIRStateVars(int idate);

	double			getDiscountFactor_5(int& indexStart,int& indexEnd,double& fwdDisc,double& B,double& var,
										int idateAsOf,int ipath,long To);

	double			getDiscountFactor_5(int idateAsOf,int ipath,int indexStart,int indexEnd,double fwdDisc,double B,double var);

	virtual void    GenerateRandomNumbers(int idate);


	CVector			m_convexityAdj;//[idate]
	CVector			m_bareConvexityAdj;//[idate]

	double			calculateConvextityAdjustment(int idate,long T);
	double			calculateConvextityAdjustment(int iDate,int jDate);


public:
					
	void			Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CVector& beta,const CVector& volvol,const CVector& meanRevRate,const CVector& meanRevLevel,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEqZeroCurveHandle	discountCurve,MlEq2DInterpolatorHandle pAccessVolInterp=NULL);
	void			Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CMatrix& modelMarketParameters,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEqZeroCurveHandle	discountCurve,MlEq2DInterpolatorHandle pAccessVolInterp=NULL);
	void			Initialize(MlEqAsset& asset,const vector<long>& simdates,MlEqStochBetaVolatilityStructure& betaVolStruct,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle&  phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite);
	void			Initialize(MlEqAsset& asset,product& deriv,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle&  phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite);
	void			Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle& phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite);
	void			Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void			Initialize(MlEqAsset& asset,const vector<long>& simdates,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle&  phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite, double fCorrelation,bool isMinimalHybrid);
	void			Initialize(long nToday,MlEqAsset& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle&  phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite);

	virtual const	CVector& GetDiscounts(int ipath);// gets discounts along one monteCarlo path
	double			getDiscountFactor(int idateAsOf,int ipath,long To);
	double			GetStochasticDiscount(int ipath,int idate);
	double			GetDiscount(int ipath,int idate,int jdate,int idateAsOf);

	double			GetDiscount(int ipath,int idate,int idateAsOf);
	double			GetDiscountToDate(int ipath,long To,int idateAsOf);

	bool			m_minimalHybrid;

	virtual void	calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1.0,int iasset=0);
	virtual void	calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1);

	explicit CStochIRForwardSkewMC(MlEqConstDateHandle hDate, double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CVector& beta,const CVector& volvol,const CVector& meanRevRate,const CVector& meanRevLevel,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEqZeroCurveHandle	discountCurve,MlEq2DInterpolatorHandle pAccessVolInterp=NULL);
	explicit CStochIRForwardSkewMC(MlEqConstDateHandle hDate) : CForwardSkewMC(hDate) {};
	virtual ~CStochIRForwardSkewMC();	
};	
	
		
		
	
	
class CMultiAssetForwardSkewMC : public MlEqMonteCarlo
{	
	
protected:
	
	vector < ForwardSkewMCHandle >	m_fwdSkewMC; // [iasset]
	void GeneratePath(int calibflag,long idum,CVector& accuracy);
	
	
	CVector m_accuracy;
	int     m_calibflag;
	void	Initialize(long nToday,vector < MlEqAssetHandle >& asset,GVector<int>& rateIsStochastic, const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);

protected:

	bool				m_calibToBasket;
	void				initializeCalibrationToBasket();
	void				CalibrateToBasket();
	CVector				m_calibBasketValue;
	MlEqStrikeHandle	m_pBasketCalibrationStrike;
	CVector				m_basketWeights;
	double				CalculateBasket();
	CVector				m_correlationShift;
	CVector				m_correlationFactor;
	double				GetCorrelation(MlEqAssetHandle& asset1,MlEqAssetHandle& asset2);
	
	int					GetNumberOfVolStates(int ipath,int idate,int iasset);


public:
	
	void generatePaths();

	void Initialize(vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,GVector<int>& rateIsStochastic, const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(vector < MlEqAssetHandle >& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(vector < MlEqAssetHandle >& asset,product& deriv,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);

	void Initialize(vector < MlEqAssetHandle >& asset,product& deriv,const std::vector<long>& additionalSimDates, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	void Initialize(vector < MlEqAssetHandle >& asset,product& deriv,const std::vector<long>& additionalSimDates,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);

//	virtual void  calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1.0,int iasset=0);
//	virtual void  calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor=1);


	explicit CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,GVector<int>& rateIsStochastic, const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,product& deriv,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes);
	explicit CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate) : MlEqMonteCarlo(hDate){};

	virtual const	CVector& GetDiscounts(int ipath);// gets discounts along one monteCarlo path
	double			GetDiscount(int ipath,int idate);
	double			GetDiscount(int ipath,int idate,int idateAsOf);

	double  GetBridgeVolatility(int ipath,int idate,double strike,int iasset);

	double	GetForwardValue(int idate,int iasset=0);

	double	AmericanDigital(CVector& prices,int ipath,int idate,double DownBarrier,double finalSpot,
									 int iasset,int method,double stressForwardBarrierVol,
									 double stressForwardBarrierVolSlope,int npoints,double blend);

};		
	
	
	
// This is only a beta version of a 
// forward skew MonteCarlo with stochastic interest rates !!!
	
		
	
	
	


class QuickMultiMc
{	

	void choleskySetup();
	double m_correl;
	cTensor m_choleskyVols;

	void getNewSpots(CVector& newSpots,int idate,CVector& currentSpots,CVector& randoms);

	CVector m_dt;//[m_mDates]
	CVector m_times;
	CVector m_fwdRates;//[idate][iasset]
	CVector m_vols;//[iasset]
	CVector m_blends;

	double getLevFactor(int idate);

	int m_nPayoffDates;
	int m_numberStepsBetweenSettings;

	randomGeneratorHandle m_randomGenerator;
	
public:

	int m_seed;
	double getCorrelation(int iasset,int jasset){return m_correl;};
	double payoff(CMatrix& pathArray,CVector& params);

	void simulate(CVector& results,CVector& params);
	void initialize(CVector& times,CMatrix& inputs,CVector& params);

	int m_nPaths;
	int m_nassets;
	int m_nDates;
	int m_mDates;


};




class LocalVolGrid
{	
	
protected :
		
	MlEqVolatilityStructureHandle m_mktVols;
	MlEqVolatilityStructureHandle m_localVolConstructor;

	GVector < MlEqInterpolatorHandle > m_localVolInterp;
		
	GVector< GVector < MlEqStrikeHandle > >	m_calibStrikes;//[idate][istrike]
	
	GVector<long>				  m_calibDates;//[icalibdate]
	CVector						  m_calibTimes;//[idate]
	long					      m_horizonDate;
	
	CVector						  m_spotGrid;//[nspots]
	CVector						  m_timeGrid;//[ntimes] // as double
	GVector<long>				  m_dateGrid;//[times]
	
	CMatrix						  m_localVolGrid;//[ntimes][nspots]
	CVector						  m_fwdValues;//[ntimes]
	CVector						  m_spotLogValues;//[nspots]
	double						  m_spot;
	int							  m_ngauss;
	int							  m_npoints;//number of integration points for variance integration
	int							  m_nToday;
	
	double localVolToImpliedVol(double K,double T);
	
	void constructGrid();
	void mostLikelyPath(CVector& spotPath,double t, double T, double finalSpot);
	
	GVector < long > m_map;//[idatecalib] : idatecalib->idate
	GVector < long > m_invmap;//[idatecalib] : idate->idatecalib
	
	int m_stepsPerDayRate;// 5 means one step every five days
	
	CMatrix		m_gaussWeight;//[icalibDate]]
	CMatrix		m_gaussPoint;//[icalibDate]
	
	
public :
	

	void initialise(MlEqAsset asset,GVector< GVector < MlEqStrikeHandle > >& calibStrikes,GVector<long> calibDates,GVector < MlEqInterpolatorHandle >& localVolGridInterp,long horizonDate,long stepsPerDayRate,int numberSpacialPoints,int ngauss);
	
	void iterateMostLikelyPath(CMatrix&	mostlikelyPath,CMatrix&	mostlikelyVariance,CVector& impliedVol,
							   GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,
							   CMatrix&  previous_mostlikelyPath,CMatrix&  previous_mostlikelyVariance,bool use_prev_most_likely_variance,bool use_prev_most_likely_path);


	void iterateMostLikelyPath(CMatrix&	mostlikelyPath,CMatrix&	mostlikelyVariance,CVector& impliedVol,
							   GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,
							   CMatrix&  previous_mostlikelyPath,CMatrix&  previous_mostlikelyVariance,bool use_prev_most_likely_variance,bool use_prev_most_likely_path,GVector < MlEqInterpolatorHandle >& localVolInterp);



	void localToImplied(CVector& impliedVols,CMatrix& mostlikelyPaths, GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,double accuracy,bool returnMostlikelypath);
	void localToImplied(CVector& impliedVols,CMatrix& mostlikelyPaths, int idateCalib,int idateCalibStart,double accuracy,bool returnMostlikelypath);
	void localToImplied(CVector& impliedVols,CMatrix& mostlikelyPaths, GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,double accuracy,bool returnMostlikelypath,GVector < MlEqInterpolatorHandle >& localVolInterp);

	void createLocalVol(CVector& impliedVols,double tolerance,double accuracy);


	LocalVolGrid(){};
};	





/**********************************************
***											***
***			Local Vol Monte Carlo			***
***											***
***********************************************/



class CLocalVolMC	:	public CSequentialBlackMC
{
public:
	CLocalVolMC( MlEqConstDateHandle hDate, MlEqAssetHandle hUdly );
	~CLocalVolMC(){}


	virtual void			initialize( product& deriv, const CMatrix& Fixings, int npath, int rngFlag, int stepAYear);	
	virtual void			initialize( const std::vector<long>& payoffDates, int npath, int rngFlag, int stepAYear);	
	
//	virtual void			initialize( product& deriv, const CMatrix& Fixings, int npath, int rngFlag);
//	virtual void			initialize( const std::vector<long>& payoffDates, int npath, int rngFlag);

	virtual void			initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag, int stepAYear);	// use MT
//	virtual void			initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag);

 	virtual void			createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray);
	virtual void			createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray,CVector& payoutPathDiscounts,CVector& simulatedDiscounts);
	virtual CMatrix&		GetPathArray(int ipath);
	virtual CMatrix&		MakePathArray(CVector& randoms);
	virtual const CVector&	GetDiscounts(int ipath);
	void					generatePaths();
	long					GetNumberOfFutureDates();

	double					GetDiscount(int ipath,int idate);

	void					computeGreeks() { m_computeGreeks = true; }

protected:
	CVector				m_forwards;		// may differ from the usual stuff
	CVector				m_quantoDrift; 
	
	int					m_nMcDates;	
	GVector<int>		m_mapToFixings ;
	GVector<long>		m_dates; // to build lv grid...

protected:
	double				computeLocalVol(int time, double rspot, const double* lv_0 );

	CVector				m_reducedSpotGrid;	
	GVector<CVector>	m_localvolgrid;		
	CVector				m_arraylocalvolgrid;		
	double				m_min_yGrid;
	double				m_dy;
	double				m_max_yIndex;
	int					m_space_size;

	bool				m_storedPaths;
	int					m_nStart;

protected:
	void			initTimeGrids( const std::vector<long>& payoffDates, int stepAYear);
	void			initTimeGrids( const std::vector<long>& payoffDates );
	void			initDiscounts(const std::vector<long>& payoffDates, const CMatrix& Fixings);
	virtual void	initLocalVolGrid();
	virtual void	initLocalVolGrid( const std::vector<long>& payoffDates );
	void			initRandomGenerator(int npath, int rngFlag);

protected:
	bool			m_computeGreeks;
	CVector			m_arraybumpedlocalvolgrid;

	CMatrix&		MakeGreeksPathArray(CVector& randoms);

protected:
	randomGeneratorHandle	m_quasiGenerator;
	GVector < CVector >		m_quasis;

	CVector					m_bbVar;
	CVector					m_bbPtr;
	GVector<int>			m_mapToBBfixings;
	int						m_nBBDates ;

	void			initBBvar();
	void			BrownianBridgePath(CVector& randoms, const CVector& quasi);
};


class CMultiLocalVolMC	:	public CLocalVolMC
{
public:
	CMultiLocalVolMC( MlEqConstDateHandle hDate, MlEqAssetHandle hUdly );
	CMultiLocalVolMC( MlEqConstDateHandle hDate, const std::vector<MlEqAssetHandle>& vhUdly);
	~CMultiLocalVolMC(){}

	virtual void			initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag, int stepAYear);	// use MT
 	virtual CMatrix&		GetPathArray(int ipath);
	virtual double			GetForwardValue(int idate,int iasset);

protected:
	GVector< RCPtr<CLocalVolMC> >		m_lvmcComponent;

	CVector			m_qRndTemp;		// for quasi generator...
	Cholesky		m_qCholesky;
};


#endif
