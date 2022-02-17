//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CIDParameters.hpp
//
//   Description : Parameters needed for sCID model 
//                           (for all names, not per name)
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_sCID_PARAMETERS_HPP
#define QLIB_sCID_PARAMETERS_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/SCIDbaseWorld.hpp"
#include "edginc/SCIDconvolution.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/Array.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(SCIDparameters)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads);
typedef array<ICDSParSpreadsWrapper, ICDSParSpreadsWrapper> ICDSParSpreadsWrapperArray;
typedef smartPtr<ICDSParSpreadsWrapperArray>                ICDSParSpreadsWrapperArraySP;
typedef smartPtr<SCIDparameters>							SCIDparametersSP;
typedef smartConstPtr<SCIDparameters>						SCIDparametersConstSP;

class MARKET_DLL SCIDparameters: public MarketObject
{
public:

    static CClassConstSP const TYPE;
    virtual ~SCIDparameters();

    virtual string getName() const;

	  /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);
    virtual void validatePop2Object();

// get some parameters
	int getNbWorlds() { return m_nbWorlds; };
	int getNbNames() { return m_nbNames; }
	double getInitialWeights(int indexWorld) { return m_initialWeights[indexWorld]; };
	void   getInitialWeights(DoubleArray &weights) { weights = m_initialWeights; };
	void   changeInitialWeights(DoubleArray &newWeights) { if (newWeights.size()==newWeights.size()) m_initialWeights = newWeights; }
	void   getMultShifts(array<double> &multShifts) { multShifts = m_linearShift; };
	void   getParallelShifts(array<double> &parallelShifts) { parallelShifts = m_parallelMultShift; };
	double getRecovery(int name) { return m_recovery[name]; };
	double getNotinoal(int name) { return m_notional[name]; };
	void   getSingleNameCalibrationDates(DateTimeArray &dates);

	ICDSParSpreadsSP getCDSparSpread(int name);
	double getMarketPortfolioEL(DateTime date, bool notionalEL=false);
	double getMarketPortfolioEL(double date, bool notionalEL=false);
	DoubleArraySP getMarketPortfolioEL(DateTimeArray dates, bool notionalEL=false);

	void MarketSurvProb(DoubleArray &time, DoubleMatrix &survProbability); 
	void MarketSurvProb(DateTimeArray &dates, DoubleMatrix &survProbability); 	

// Help to create TimeLine for products
	void push_backTimeLine(DateTimeArray &timeLine,
		  				   DateTime firstDate, 
						   DateTime lastDate, 
						   MaturityPeriodSP period, 
						   bool includeFirstDate = true);
	DoubleArray DateAsDouble(DateTimeArray &dates);
    
    /** compute risky annuities and default legs from a start date to a bunch of maturities. It uses flags to compute Risky annuities, and default legs **/
    void computeRiskyAnnuityLeg(DateTime& valueDate, YieldCurveSP discount, DateTimeArray& simDates,
                                    DoubleArray& expectedLoss, double trancheDiv, DateTime& startDate, DateTimeArray& endDates,
                                    int coupon, const DayCountConvention & dcc, DoubleArray & RAvalues);  
    void computeDefaultLeg(DateTime& valueDate, YieldCurveSP discount, DateTimeArray& simDates,
                                    DoubleArray& expectedLoss, double trancheDiv, DateTime& startDate, DateTimeArray& endDates,
                                    DoubleArray & DLvalues); 


// Closed Form
	void NamesSurvProbinGivenWorld(int world, DoubleArray &times, DoubleArray &survProbability); 	// survProbabilility(name,time) = survProbabilility(name*time.size()+index time)
	void NamesSurvProbinGivenWorld(int world, DateTimeArray &dates, DoubleArray &survProbability); 	// survProbabilility(name,time) = survProbabilility(name*time.size()+index time)

	void NamesSurvProb(DoubleArray &time, DoubleMatrix &survProbability); 	
	void NamesSurvProb(DateTimeArray &dates, DoubleMatrix &survProbability); 	

	void PortfolioELinGivenWorld(int world, DoubleArray &times, double *pel);  // pel must be of size times.size()...
	void PortfolioELinGivenWorld(int world, DateTimeArray &dates, double *pel);  // pel must be of size dates.size()... 

	void PortfolioEL(DoubleArray &times, double *pel);  // pel must be of of size times.size()
	void PortfolioEL(DateTimeArray &dates, double *pel);  // pel must be of of size dates.size()

	void reInitializeSurvProba(DoubleMatrix &survProbaAtCalibrationSingleNameDates);
// Fast MC
	void setFastMC(int seed, double timeSteps, int nbSteps, int noJumpPaths, int jumpsPaths);
	void resetIpathFastMC(int iPathBM = 0, int iPathPoisson = 0, int iPathJumpTime = 0);
	void getIpathFastMC(int &iPathBM, int &iPathPoisson, int &iPathJumpTime);
	void setConvolution(DoubleArray &kmin, DoubleArray &kmax);
	void ComputeTELinAllWorlds(DoubleArray &times,  // times at which we compute these TEL,
							   vector< DoubleMatrix > &TEL,  // size nbWorlds, DoubleMatrix of size times.size(), kmin.size().
							   int convolutionNoJump, int convolutionJumps); 

	void ComputeSurvProbinAllWorld_fastMCtest(DoubleArray &times,                  // compute the survival probability of each names, using the fast MC methodology
											  vector< DoubleMatrix > &survProba);  // only used for testing
																				

// Full MC
	void setFullMC(int seed, DoubleArray &discTimes);
	void setFullMC(int seed, DateTimeArray &discDates);
	void resetIpathFullMC(int iPathBM, int iPathPoisson, int iPathJumpTime, int iPathJumpSize, int iPathDeafult);
	void resetSpreads();
	void FullMCSimulation() { SimulateNbJumps(false); SimulateSpread(); SimulateDefaults(); }
	void SimulateNbJumps(bool conditional);
		// simulate the NbJumps of the Poisson processes
		// if conditional==true, conditionally on at least one jump by the last date
	void SimulateSpread();  // always call SimulateNbJumps first, as it assumes it knows how many jumps there are
							// simulate the spreads in all worlds
	void SimulateDefaults(); // given the spreads, simulate defaults and hence the weights
	void getPortfolioLoss(DoubleArray &loss) { loss = m_loss; };
	void getNotionalLoss(DoubleArray &loss) { loss = m_notionalLoss; };
	void getNbDefaultedNames(IntArray &nbDefNames) { nbDefNames= m_nbDefaultedNames; };

	void getDefaultTime(DoubleArray &defaultTime) { defaultTime = m_defaultTime; };

    //methods to store path-dependent data in full MC
    void InitializeDataStorage(int timeSlices, int nrPaths, int nrStoredValues);
    void StorePathValues(int timeSlice, int nrPath, DoubleArray & values);
    DoubleMatrixArray getStoredValues();
    //get intensities for a name
    DoubleArray getNameIntensity(int nameIndex);
    //get intensities for all names
    DoubleArrayArray getIntensities();

// HYBRID METHODS
	
	void getSurvProba(int indexLittle_t,
				  vector<int> &indexBig_t,                  
				  DoubleMatrix &survProba,
				  bool add = false);	  // return the survival probabilty between discTimes[indexLittle_t] and discTimes[indexBig_t],
							  // as computed in the given simulation (no average is being taken
							  // survProba must be of size indexBig_t.size() * m_nbNames;
							  // if add, add to survProba the answer // only for testing purposes

	// given the values of the simulation, compute at m_discTimes[index] the portfolio EL at times
	double getFuturePortfolioEL(int indexFutureTime,
		  					    DoubleArray &times,
								DoubleArray &EL);

	double getFuturePortfolioNotionalEL(int indexFutureTime,
		  					    DoubleArray &times,
								DoubleArray &notionalEL);

	/** given the values of the simulation, compute at m_discTimes[index] the survival probability at times */
	void getFutureSurvProba(int indexFutureTime, DoubleArray &times, DoubleMatrix &survProba);

	/** given the values of the simulation, compute at m_discTimes[index] the 
	 TEL at times.. uses the loss before m_discTimes[index] if countPastLosses
	 survProba must be of size times.size() * m_nbNames
	 if add, add to the previous values */
	double getFutureETL(int indexFutureTime, 
						DoubleArray &times, int begin,
						DoubleMatrix &ETL,  // ETL of size begin+times.size() * kmin.size()
						bool countPastLosses,
						int convolutionNoJump, int convolutionJumps);  

protected:  
    SCIDparameters(const CClassConstSP& clazz);
	SCIDparameters(const SCIDparameters& v);
	SCIDparameters();

    static void load(CClassSP& clazz);
	static IObject* defaultSCIDparameters();

    string               name;
private: 
	void ComputePELcurve(); // Compute the implied PEL curve from the CDS market
							 // Compute PEL at singleNameExpiries

	
	// Used Model parameters
//    YieldCurveWrapper               discount;        
    ICDSParSpreadsWrapperArraySP parSpreads;  // Porfolio of names within the model
	int m_nbNames, m_nbWorlds, m_nbJumpProcesses;
	DoubleArray m_linearShift, m_parallelMultShift, m_initialWeights;    // size m_nbWorlds
	BaseWorld m_world;          // contains the parameters for one CID world
	DoubleArray m_notional, m_recovery;
	DoubleArrayConstSP diffParameters;
	CDoubleMatrixConstSP jumpParameters;
	ExpiryConstSP cfExpiry;
	DateTime today;
	ExpiryArrayConstSP singleNameExpiries;
	DoubleArray m_singleNameTimesWithToday;
	DoubleMatrix baseWorldCalibSurvProba; // optional, read from CDS curve if not given

	
// storing of the portfolio expected loss curve
//	EffectiveCurveSP pelCurveMarket;
// variable to avoid creating in every call of a function new arrays or doubles
	DoubleArray survProb; // of size m_nbNames.. used in a few different places, in particular as in the input to convolution
	DoubleMatrix telWork1, telWork2;

//  Convolution parameters
	SCIDconvolution convol;

// Monte Carlo parameters
	// fast
	PoissonProcess m_PoissonFast;         // simulates the number of jump of the Poisson process 
	double probaZeroJump;
	BrownianMotion m_BMFast;			  // the Brownian motion leading to the common factor
	Uniform m_jumpTimeFast;				  // simulates the jumptime of the Poisson process
	double m_lastTimeFast;
	DoubleArray m_BMincrementsFast;		  // save increments of the BM (driving the cf)
	IntArray m_nbJumpsFast;				  // NbMarkets, save the (random) nb of jump of the different poisson process
	int nbPathsNoJump, nbPathsJumps;     // m_nbPaths[0] nb Paths conditionally on no jumps, m_nbPaths[1] conditionally on at least one jump in the Poisson processes
	// full
	int m_nbPathsFull;
	DoubleArray m_discTimes, m_expDecayCF;     // its size defines nbTimes
	vector<Spread> m_intensity;  // store spreads
	DoubleMatrix m_jumpIntensity;  // size nbMarkets * (nbNames * nbTimes); // store jump spreads
	IntArray m_nbJumpsFull;
	DoubleMatrix m_weightsDynamic;	// weights for the worlds at all time
	DoubleArray m_defaultTime; // size nbNames // time of defaults of all name // stupidely high is no default by maturity
	DoubleArray m_loss, m_notionalLoss; // size nbTimes // portfolio loss at m_discTimes, portfolio notional loss
	IntArray m_nbDefaultedNames;
	DoubleArray m_bmIncrementsFull, m_intIntensityTemp, m_forUpdate;// m_intensityTemp;
	DoubleMatrix m_lambdaFutureTimeJump;
	DoubleArray m_lambdaFutureTimeIdio;
	double m_lastDateFull;
	PoissonProcess m_poissonFull;
	BrownianMotion m_bmFull;
	Uniform m_defaultFull, m_jumpTimeFull, m_jumpSizeFull;

    //stored path-dependent data for full MC simulations
    DoubleMatrixArray dataStorage; 
};



DRLIB_END_NAMESPACE

#endif
