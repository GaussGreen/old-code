//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CID Tree Parameters.hpp
//
//   Description : Parameters needed for sCID model 
//                           (for all names, not per name)
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SCIDtreeParameters_HPP
#define QLIB_SCIDtreeParameters_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/SCIDTreeBaseWorld.hpp"
#include "edginc/SCIDTreeWorld.hpp"
#include "edginc/SCIDconvolution.hpp"
#include "edginc/EffectiveCurve.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(SCIDtreeParameters)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads);
typedef array<ICDSParSpreadsWrapper, ICDSParSpreadsWrapper> ICDSParSpreadsWrapperArray;
typedef smartPtr<ICDSParSpreadsWrapperArray>                ICDSParSpreadsWrapperArraySP;


typedef smartPtr<SCIDtreeParameters> SCIDtreeParametersSP;
typedef smartConstPtr<SCIDtreeParameters> SCIDtreeParametersConstSP;

class MARKET_DLL SCIDtreeParameters: public MarketObject
{
public:

    static CClassConstSP const TYPE;
    virtual ~SCIDtreeParameters();

    virtual string getName() const;

	  /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);
    virtual void validatePop2Object();

// get some parameters
	int			getNbNames() { return m_nbNames; }
	double		getRecovery(int name) { return m_recovery[name]; };
	void		getSingleNameCalibrationDates(DateTimeArray &dates);
	double		getFastLastTime() { return fast.lastTime; };

	ICDSParSpreadsSP getCDSparSpread(int name);
	double		getMarketPortfolioEL(DateTime date, bool notionalEL=false);
	double		getMarketPortfolioEL(double date, bool notionalEL=false);
	DoubleArraySP getMarketPortfolioEL(DateTimeArray dates, bool notionalEL=false);

	void 		MarketSurvProb(DoubleArray &time, DoubleMatrix &survProbability); 
	void 		MarketSurvProb(DateTimeArray &dates, DoubleMatrix &survProbability); 	

// Help to create TimeLine for products
	void		push_backTimeLine ( DateTimeArray &timeLine,
		  									DateTime firstDate, 
											DateTime lastDate, 
											MaturityPeriodSP period, 
											bool includeFirstDate = true );
	DoubleArray	DateAsDouble(const DateTimeArray &dates);

// Calibration of the Tree
	void		setTELtimesInTree(DoubleArray &times) { m_tree.setTELtimes(times); };
	void		setTELtimesInTree(DateTimeArray &dates) { DoubleArray times = DateAsDouble(dates); m_tree.setTELtimes(times); };

	void 		CalibrateBaseWorld(DoubleMatrix &survProba);

	void		readTree(DoubleMatrix &valuesAtTreeTime) { m_tree.readTree(valuesAtTreeTime); };

// Closed Form
	void		TESTSurvProbinGivenWorld ( LinearInterpolantSP world, DoubleArray &T, DoubleMatrix &survProbability);

// Fast MC
	void 		setFastMC ( int seed,                 // seed used by the random number	
							double timeSteps,		  // time steps to discretise the CIR process	
							double lastTime,		  // until when do we simulate
							int noJumpPaths,		  // how many paths conditional on no jumps
							int jumpsPaths,			  // how many paths conditional on at least one jumps
							int _convolutionNoJump,   // convolution method conditional on no jumps
							int _convolutionJumps,    // convolution method conditional on at least one jumps
							DoubleArray &kmin,        // lower attachment points
							DoubleArray &kmax,        // upper attachment points 
							int		  coupon,	      // coupon paid every ? months
							DateTimeArray maturities, // tree dates
							YieldCurveSP &discount,   // discount factor use
							DayCountConventionSP &dcc);  // day count convention used

	void		computeAllTELeasy(DoubleMatrixArray *outputTEL);  // for testing
	void		InitializeTree(DoubleMatrix *pos = 0, DoubleArray *w = 0) { if ((pos==0)||(w==0)) m_tree.globalFeed(positions,weights); 
																							 else m_tree.globalFeed(*pos,*w); }

	void 		computeAllTELandPEL(bool feedTheTree);
	void		computeAllLegsFromTEL();

	double		getStressPEL() { return stressPEL; }
	void		CalibrateTreeToTranches(const DoubleMatrix &parSpread,
										const DoubleArray &PEL,
										double smoothing,
										size_t lowGen,
										size_t highGen);												
	//void		ImprovePosition(const DoubleMatrix &parSpread,
	//							const DoubleArray &PEL,
	//							int lowGen,
	//							int highGen,
	//							bool reconcile,
	//							double epsilon,
	//							int maxNbLoops,
	//							double ftol);  
	//void		CalibratorHelp (vector<SimpleTreeSP> &leaves, 
	//							const DoubleArray &ratios, 
	//							DoubleArray &MTMandPEL, // (O)
	//							const DoubleMatrix &parSpreadMarket,
	//							const DoubleArray PELmarket,
	//							bool changeWeights,
	//							double smoothing,
	//							bool deleteZeroWorlds = false);
	void		guessPositionAtAllMaturities(const DoubleMatrix &parSpread,
								const DoubleArray &PEL,
								double Dx,
								double smoothing,
								DoubleArrayArray &selectedPos);


	void	    getLeaves(vector< SimpleTreeSP > &s, size_t gen = -1) { m_tree.getLeaves(s,gen); }

	void		readAllLegs(vector < SimpleTreeSP > &leaves, DoubleMatrixArray &DL, DoubleMatrixArray &RA);
	void		readAllPEL(vector < SimpleTreeSP > &leaves, DoubleArrayArray &PEL);
	void		readAllTEL(vector < SimpleTreeSP > &leaves, DoubleMatrixArray &TEL);

	void		readAverageLegs(vector < SimpleTreeSP > &leaves, DoubleMatrix &DL, DoubleMatrix &RA);
	void		readAveragePEL(vector < SimpleTreeSP > &leaves, DoubleArray &PEL);
	void		readAverageTEL(vector < SimpleTreeSP > &leaves, DoubleMatrix &TEL);
	vector<LinearInterpolantSP>  getWorlds(vector < SimpleTreeSP > &leaves) { vector<LinearInterpolantSP> li; m_tree.getWorldsFromTree(li,leaves); return li;};
	void	     getWorlds(DoubleArray &x, vector<DoubleArray> &y, vector<SimpleTreeSP> &leaves) { m_tree.getWorldsFromTree(x,y,leaves); }; 



	void		TESTSurvProbinGivenWorldFastMC (DoubleArray &times,
												LinearInterpolantSP world,
												DoubleMatrix &survProbability);

/*	void ComputeSurvProbinAllWorld_fastMCtest(DoubleArray &times,                  // compute the survival probability of each names, using the fast MC methodology
											  vector< DoubleMatrix > &survProba);  // only used for testing
*/																				

// Full MC
	void setFullMC(long seed, DoubleArray &discTimes);
	void setFullMC(long seed, DateTimeArray &discDates);
	void resetIpathFullMC(long _iPath) { m_world.setIpathFull(_iPath); }
	void FullMCSimulation() { SimulateNbJumps(false); SimulateSpread(); SimulateDefaults(); }

	/** simulate the NbJumps of the Poisson processes
	    if conditional==true, conditionally on at least one jump by the last simulation date */
	void SimulateNbJumps(bool conditional) { m_world.SimulateJumpsFull(full.spread.discTimesSpread.back(),conditional); }

	/** simulate the spreads in all worlds
	    always call SimulateNbJumps first, as it assumes it knows how many jumps there are */
	void SimulateSpread();  

	/** simulate default times */
	void SimulateDefaults(); 


	void getPortfolioLoss(DoubleArray &loss) { loss = full.loss; };
	void getNotionalLoss(DoubleArray &loss) { loss = full.notionalLoss; };
	void getNbDefaultedNames(IntArray &nbDefNames) { nbDefNames= full.nbDefaultedNames; };
	void getDefaultTime(DoubleArray &defaultTime) { defaultTime = full.defaultTime; };



protected:  
    SCIDtreeParameters(const CClassConstSP& clazz);
	SCIDtreeParameters(const SCIDtreeParameters& v);
	SCIDtreeParameters();

    static void load(CClassSP& clazz);
	static IObject* defaultSCIDtreeParameters();

    string               name;
private: 
	void 		computeTELinGivenWorlds(DoubleArray &times,
										vector< LinearInterpolantSP > &f,
										DoubleMatrixArray & TEL);  // must be initialized
	void 		computePELinGivenWorlds(double time,
										vector< LinearInterpolantSP > &f,
										vector<double> & PEL) 
	{ 
		vector<double *> pel2(PEL.size()); 
		for (size_t i=0; i<PEL.size(); ++i) pel2[i]=&PEL[i];
		computePELinGivenWorlds(time, f, pel2); 
	}

	void 		computePELinGivenWorlds(double time,
										vector< LinearInterpolantSP > &f,
										vector<double *> & PEL); // must be initialized
	void		computeTELandPELNGeneration(vector<SimpleTreeSP> &leaves, int lastGen);
	void 		computeTELandPELOneGeneration (vector<SimpleTreeSP> &leaves);  



	// Used Model parameters
	int m_nbNames, m_nbWorlds, m_nbJumpProcesses;

/* worlds parameters **/
	public:
		DoubleArrayArrayArray a_t;  // to MOVE BACK TO PRIVATE 
	private:
	DoubleArray bounds;
	double stressPEL;
	DoubleMatrix positions;
	DoubleArray  weights;


public: // temporary, for testing!!!
	SCIDtreeWorld m_tree;
private:

/* base CID parameters **/
	DoubleArray m_notional, m_recovery;
	double volDiff, decayDiff, idioRatio, cfRatio;
	DoubleArray jumpDecay, jumpImpact, jumpFreq, jumpRatios;
	DoubleArrayArrayArray jumpWeights, jumpMeans;
/* those parameters will feed the base world : **/
	TreeBaseWorld m_world;          

/*  pure market data **/
	DateTime today;
	ExpiryArrayConstSP singleNameExpiries;
	ICDSParSpreadsWrapperArraySP parSpreads;  // Portfolio of names within the model

/*  Convolution parameters **/
	SCIDconvolution convol;

// Monte Carlo parameters
	// fast
	struct FastMCparameters
	{
		int nbPathsNoJump, nbPathsJumps;     // Number of Paths conditionally on no jumps, m_nbPaths[1] conditionally on at least one jump in the Poisson processes
		int convolutionNoJump, convolutionJumps;
		double lastTime;
	};
	struct FullMCparameters
	{
		int nbPaths;
		DoubleArray loss, notionalLoss, defaultTime;
		IntArray nbDefaultedNames;
		baseWorldSpread spread;
		IntArray discTimesDefaults;
		DoubleMatrix weightsDynamic;
	};

	FastMCparameters fast;
	FullMCparameters full;



};



DRLIB_END_NAMESPACE

#endif
