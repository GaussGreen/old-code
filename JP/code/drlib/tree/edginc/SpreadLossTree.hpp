//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : SpreadLossTree.hpp
//
//   Description : Tree for bivariate spread and portfolio loss process
//
//   Author      : Matthias Arnsdorf
//
//----------------------------------------------------------------------------

#ifndef QLIB_SPREAD_LOSS_TREE_HPP
#define QLIB_SPREAD_LOSS_TREE_HPP

#include "edginc/FDModel.hpp"
//#include "edginc/DateTime.hpp"
#include "edginc/SkewSurface.hpp"
//#include "edginc/CreditTree.hpp"
#include "edginc/LossTreeRecovery.hpp"
#include "edginc/IAuxSpreadProcess.hpp"
#include "edginc/Optimizer.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/IForwardRatePricer.hpp"


DRLIB_BEGIN_NAMESPACE



class SpreadLossTree;

typedef smartPtr<SpreadLossTree> SpreadLossTreeSP;
typedef smartConstPtr<SpreadLossTree> SpreadLossTreeConstSP;

class TrancheContingentLeg;

typedef smartPtr<TrancheContingentLeg> TrancheContingentLegSP;
typedef smartConstPtr<TrancheContingentLeg> TrancheContingentLegConstSP;

class TrancheFeeLeg;

typedef smartPtr<TrancheFeeLeg> TrancheFeeLegSP;
typedef smartConstPtr<TrancheFeeLeg> TrancheFeeLegConstSP;

class TREE_DLL SpreadLossTree:	public FDModel, // base tree model
								virtual public FDModel::IFDSolver
{
public:


    static CClassConstSP const TYPE;
    static void load(CClassSP& );

	
    virtual ~SpreadLossTree();

	virtual void validatePop2Object();

    /** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);    

	///////////////////// SOLVER INTERFACE //////////////////////////
	/** spread loss tree specific roll function (contains the time loop) */
	// should be implemented at credit tree level
	virtual void roll();
  


	//////////// SPREAD LOSS TREE METHODS /////////////////////////////////////////////////////////////

	/**Get a reference to a specific contingentleg's price at a given time step.
	To be used for effective curve pricing only */
	virtual void getCDOContingentLegValue(
		TreeSlice& valueSlice,		// (O) slice to fill
		const TrancheContingentLeg& protLeg, int step);

	

	/** Get a reference to a specific fee leg's price at a given time step.
		This is only for pricing using effective curves */	
	virtual void getCreditFeeLegValue(
		TreeSlice& valueSlice,	// (O) slice to fill
		const TrancheFeeLeg& feeLeg, 
		int step,
		IForwardRatePricerSP model);

	/** return slice of total losses for each node for maximum loss range*/
	TreeSliceSP getLossSlice();

	/** return slice of total recovery for each node for maximum loss range*/
	TreeSliceSP getRecoverySlice();

    /** return default levels */
    IntArraySP getDefaultLevels()
    {
        return IntArraySP(new IntArray(defaults));
    }

	/** adds strikes to either el..Strikes or elWithRR...Strikes if not there already.
	returns the index in the array at which the strikes are located 
	This is needed if want to price using effective curves.
	Needs to be called after initModel */
	int addStrikes(
		double lowStrike,
		double highStrike,
		bool recoverNotional
		);

	/** Can limit the number of defaults used in the tree. To do this need to add the high strike of the tranche
	to be priced prior to finaliseModel() */
	void addMaxStrike(double strike);

    
	
friend class TrancheContingentLegProd;
friend class TrancheFeeLegProd;
	
protected:
	

	// FDModel INTERFACE //////////////////////////////////////////////
	
	/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
	virtual double getPrice0( const TreeSlice & price ) const;

	/** collect model initialization data, set up timeline */
	// This is called after product->init() but before product->initProd()
	virtual void initModel();

	/** finalize model initialization */
	// This is called after product->initProd()
	virtual void finaliseModel(CControl*    control);

	/** retrieve market data into model */
	virtual void retrieveFactor();

public :

	/** creates new slice */
	// This should ultimately be done at CreditTree level
	virtual TreeSliceSP createSlice(
		const string & curveToDEV = "",
		const string & factorName1 = "",
		const string & factorName2 = "",
		const string & factorName3 = "" ) const;


	/** as above but return a TreeSliceGeneral */
	virtual TreeSliceGeneralSP createGeneralSlice() const;

	/** override FDModel::makeProduct to retrieve IndexSpec product */
	// Use FDModel::makeProduct for now
	virtual FDProductSP makeProduct(const IProdCreatorSP & creator)
	{
		return FDModel::makeProduct(creator);
	}

	/** create the tree specific solver */
	virtual IFDSolverSP createSolver() { 
		return IFDSolverSP(this, NullDeleter()); 
	}

	/** register with the model required zero/discounters */
	// In this case nothing to do since rates are treated as deterministic
	virtual void registerZero(
		DateTime useDate, 
		DateTime matDate, 
		string curveName) {}

	/** Get slice representing zero/discounter between useDate and matDate 
	*  If "slice" is null it creates a new slice (no ownership kept), 
	*  otherwise it uses the slices passed to return the zero */
	virtual void getZero(
		DateTime useDate, 
		DateTime matDate, 
		string curveName, 
		TreeSliceSP &slice);

	/** rates tree products distinguish between the today date and the value date -
	they may or may not be the same date.  The getToday returns the 
	reference date in the market that defines today.  If there is any possibility
	that valueDate != today, derived FDModel should override CModel::getValueDate
	to return the appropriate value date */
	virtual DateTime getToday() const {
		return valueDate;
	}


    // GET CONDITIONAL QUANTITIES --------------------------------
    /** get the state probabilities ((spread , default) state)for a given time */
    TreeSliceGeneralSP getStateProbSlice(int step);

    /** get probability of defaults at given step */
    DoubleArraySP getMarginalDefaultProbs(int step);

    /** get values of slice conditional on loss states */
    DoubleArraySP getValConditionalOnLoss(const TreeSlice & valueSlice, int step);


	/**
	* Whether to enable RiskMapping when computing sensitivities for
	* instruments priced using this model
	*
	* Returns riskMappingIrrelevant --- not sure if this is always right in
	* the long run.  See IModel::wantsRiskMapping().
	*/
	virtual IModel::WantsRiskMapping wantsRiskMapping() const
	{
		return riskMappingIrrelevant;
	}

protected :

	/** Update tree variables at the "top" of the time loop */
	void update(int t, FDProduct::UpdateType type);

	/** Implement any updating logic that must occur at the "bottom" of the time loop after
	all product updates */
	void updateFinalize(int t, FDProduct::UpdateType type);

	/** Calculates discounted expected value of slice.
    * Can be used in fwd or bkwd direction. In Bckwd direction we obtain the standard DEV from
    * t+1 to t. In fwd case we obtain the values of the slice at t+1 given the slice values at t, ie the 
    * the value at node (sIdx, dIdx) at t+1 is given by the weighted average of all values at t with 
    * branches that end at (sIDx, dIdx). The weighting is 1/df * (prob of branche).
    */
	void sliceDev(
        TreeSlice& slice, 
        const string& curveName, 
        FDProduct::UpdateType type = FDProduct::BWD) const;

	/** This is same as sliceDEV with df = 1 */
	void sliceEv(TreeSlice& slice, FDProduct::UpdateType type = FDProduct::BWD) const;

	/** Expand the dimension of the slice so that it can be DEV'd. */
	void sliceExpand(TreeSlice& slice, const string& curveName) const;

	

private:
	
	/** default value for minimum condDefTransProb */
	static const double DEFAULT_MIN_TRANS_PROB;

	/** default max expansion order */
	static const int DEFAULT_MAX_ORDER;

	/** default pade order */
	static const int DEFAULT_PADE_ORDER;

	/** accuracy for key in generator caching */
	static const double DT_ACCURACY;

	/** maximum number of elements for the cache */
	static const int MAX_CACHE_SIZE;

	/** default constructor */
	static IObject* defaultConstructor(void);

	/** private constructor */
	SpreadLossTree(); // for derived classes

	// copy constructor only - only defined
	SpreadLossTree(const SpreadLossTree& rhs);

	// operator = (definition only)
	SpreadLossTree& operator=(const SpreadLossTree& rhs);
	
	
	//////////// METHODS ///////////////////////////////

	/** method to initialise model parameters that are independent of the products */
	// called first before any other product or model initialisation
	// parameters defined here will not be recalculated unless the model is changed, i.e.
	// they will be cached between different instrument pricings
	void setupModel();

    /** updates the conditional transition probabilities and returns the key */
    int updateCondProbs(int timeIdx, FDProduct::UpdateType type);

	/** update the transition probabilities in  condDefTransProbs. This vector contains all conditional default
		transition probabilities conditional on the spread given by spreadIdx and the 
		default level given by defaultIdx */
	void updateCondDefTransProb_expansion(
        int tIdx, // timeStep
        DoubleArrayArrayArraySP probsToUpdate,
		IntArrayArraySP defBranchesToUpdate);

	void updateCondDefTransProb_Pade(
        int tIdx, //Time step
        DoubleArrayArrayArraySP probsToUpdate,
		IntArrayArraySP defBranchesToUpdate);

	/** update the expected loss curves */
	void updateExpLoss();


	/** update a contingent leg slice*/
	void updateContLeg(TreeSlice& valueSlice,
		double lowStrike,
		double highStrike);


	/** factor adjustment to next-to-default intensity that depends only on the default levels*/
	// This function needs to be calibrated
	double localFactor(
		int numDefaults,    // number of defaults
		int timeIndex
		) const;


	/** searches for pair of values in pair of arrays
		returns index of paris in the arrays if found otherwise returns -1
		*/
	int findPair(double val1, double val2, const DoubleArray & array1, const DoubleArray & array2) const;

	/** returns index in spreads array (0 offset) given tree slice index (centred around 0) */
	int sOffset(int sIdx) const
	{
		return sIdx - minSpreadRangeBot;
	}


#if 0 // not needed for now
	/** returns index for 2-d spread/default slice */
	// (spread, default)
	 int sdOffset( int s, int d ) const
	 {
		 return (numNames + 1)*spreadProcess->sOffset(s) + d;
	 }
#endif

	 /** create range for 2d spread - default slice */
	 TreeSliceGeneral::RangeSP createRange2D(
		 int bot1, int top1,
		 int bot2, int top2
		 );

	// INPUT FIELDS /////////////////////////////////////

	/** number of names in portfolio */
	int numNames;

	/** number of tree steps per year*/
	int stepsPerYear;

	/** Maximum number of defaults per tree time period */
	int maxNumDefs;


	/** use the pade approximation (or expansion) */
	bool usePade;

	/** order for pade approx */
	int padeOrder;

	/** max loss level for which want to calculate generator matrix */
	// need this for calibration until have CDO collection
	double maxLossForGenerator;

	/** maximum order in expansion of conditional transition probs */
	int maxExpansionOrder;

	/** Minimum conditional default transition probability */
	double minTransProb;

	/** recovery model */
	LossTreeRecoverySP recoveryModel;

	/** Spread process which can depend on default event in time period */
	IAuxSpreadProcessSP	   spreadProcess;

	/** curve used for discounting */
	YieldCurveWrapper curveToDiscount;

	/** calibrated local contagion factor */
	SkewSurfaceSP localContagion;

    /** loss level at which we start bucketing the defaults.
    * [optional (default = 1)]
    */
    double bucketStartLoss;

    /** number of defaults corresponding to bucketStartLoss */
    int bucketStart;

    /** default bucket size in numbers of defaults 
    [optional default = 1] 
    */
    int bucketSize;


	/** use effective curve pricing */
	bool useEffCurves;
	
	// TRANSIENT FIELDS ///////////////////////////////
	
	/** current index in timeline */
	int currStep;

	/** value date */
	DateTime valueDate;

	/** discount curve */
	YieldCurveConstSP discYC;


	/** underlying value of defaults for given default idx */
    // not that defaults[i] != i in general because of bucketing
	IntArray  defaults;

    /**  index such that defaults[maxDefIdxForCache] >= maxDefsForCache */
    int maxDefIdxForCache;

	/** Bottom index for spreads array at currStep*/
	IntArray spreadRangeBot;
	/** Top index for spreads array at currSrep */
	IntArray spreadRangeTop;

	/** maximum number of spread nodes in spread process tree */
	int maxSpreadSize;

	/** number of spread process branches */
	int numSpreadBranches;

	/** minimum value of bottom spread range = range->bot at final slice */
	int minSpreadRangeBot;
	/** Bottom index for defaults array per timeLine point*/
	IntArray defaultRangeBot;
	/** Top index for defaults array per timeLine points */
	IntArray defaultRangeTop;

	/** maximum size of default slice */
	int maxDefRangeSize;

	/** maximum strike up to which we want to calc defaults */
	double maxStrike;

	/** number of default branches (= size of default range) in outer loop at state (sIdx, dIdx) */
	IntArrayArraySP numDefaultBranches;

	/** cache for numDefaultBranches corresponding to condDefTransProbs */
	map< int, IntArrayArraySP > numDefBrancheCache;

	/** Range for tree slices */
	TreeSliceGeneral::RangeSP range;

	/** slice that can be used to store temporary (copied values) */
	TreeSliceGeneralSP tempSlice;

	/** holder for coefficients needed in updateCondDefTransProb */
	DoubleArray coeffs;

    /** holder for generator matrix */
    CDoubleMatrix generator;

    /** holder for exp(generator) */
    CDoubleMatrix expGen;

    /** holder for the auxiliary variables needed to compute exp(generator) */
    CDoubleMatrix aux4exp1;
    CDoubleMatrix aux4exp11;
    CDoubleMatrix aux4exp2;
    CDoubleMatrix aux4exp21;
    CDoubleMatrix aux4exp3;

	/** conditional default transition probabilities */
	// condDeftransProbs[sIdx][numDefs][dIdx] = Prob (def(t') = dIdx+numDefs | spreads(t) = s, defs(t) = dIdx)
	DoubleArrayArrayArraySP condDefTransProbs;
	

	/** holder for expected loss curves for tranches. Write down according to loss only */
	// effectiveCurves[i][j] is slice of effective curve values at time currStep + j for strikes if index i
	vector< vector < TreeSliceGeneralSP > > trancheEL;

	/** holder for expected loss curves for tranches. 
	Write down according to loss and recovery (index convention for fee leg) */
	vector< vector < TreeSliceGeneralSP > > trancheELwithRR;

	/** set of low strikes corresponding to trancheEL */
	DoubleArray elLowStrikes;
	/** set of high strikes for tranchEL */
	DoubleArray elHighStrikes;

	/** set of low strikes for tranchELwithRR*/
	DoubleArray elWithRRLowStrikes;
	/** set of high strikes for trancheELwithRR */
	DoubleArray elWithRRHighStrikes;

	/** total losses per default */
	DoubleArrayConstSP totalLosses;

	/** total recoveries per default */
	DoubleArrayConstSP totalRecRates;


	/** cache for contagion factors corresponding to expGenerator */
	map< int, DoubleArraySP > contagionCache;

	/** cache for condDefTransProbs */
	map< int, DoubleArrayArrayArraySP > condDefProbCache;

	/** key for generator at last time step. Used to check if key has changed*/
	int lastKey;

	
	/** dates for localContagion surface (assume piecewise constant contagion in between dates) */
	DateTimeArray contagionDates;
	
	/** dates at which spread slice values change */
	DateTimeArray spreadDates;

	/** ceiling projection (see DateTime) of contagionDates in timeLine */
	vector<int> contagionDateProjection;

	/** ceiling projection of spreadDates (dates where spread values change )in timeLine */
	vector<int> spreadDateProjection;

	/** maximum number of defaults corresponding to maxLossForGenerator */
	// is -1 if max loss not defined 
	int maxDefsForCache;


    /// conditional output fields
    

    /** cache for state probs at given timeline steps */
    map< int, TreeSliceGeneralSP > stateProbCache;

    /** cache for marginal default probs */
    map< int , DoubleArraySP > marginalDefProbCache;

    
};


//====================================================================================
//
// Elementary Instruments / Products for SpreadLossTree
//
//====================================================================================
/************************************************************************/
/* Contingent leg of CDO with pair of strikes                           */
/************************************************************************/
class TREE_DLL TrancheContingentLeg :	public CObject, 
										virtual public FDModel::IIntoProduct
{
public:
	static CClassConstSP const TYPE;

	/** public constructor */
	TrancheContingentLeg(
		ICreditContingentLegSP contingentLeg,			// underlying contingent leg
		double lowStrike,							// as fraction of 1
		double highStrike,							// as fraction of 1
		const IBadDayAdjusterConstSP  bda
		);

	/** returns critical dates */
	virtual DateTimeArraySP getCritDates() const;

	const string& getName() const { return name; }

	friend class SpreadLossTree;
	friend class TrancheContingentLegProd;

	

protected :
	/// FIELDS
	/** unique name */
	string name; 
	/** underlying contingent leg */
	ICreditContingentLegSP  contingentLeg;
	/** low strike of tranche as fraction of 1 */
	double lowStrike;
	/** high strike of tranche as fraction of 1 */
	double highStrike;



	/** bad day adjuster */
	const IBadDayAdjusterConstSP  bda;

private:
	static void load(CClassSP& clazz);

	/** private constructor */
	TrancheContingentLeg();

	static IObject* defaultConstructor(void) { return new TrancheContingentLeg(); }

	virtual FDProductSP createProduct( FDModel * model ) const;

	

};



/************************************************************************/
/* product for TrancheContingentLeg                                     */
/************************************************************************/
class TREE_DLL TrancheContingentLegProd : public FDProduct 
{
public:
	TrancheContingentLegProd(const TrancheContingentLegConstSP &inst, FDModel* model);

	virtual void init(Control*) const;	// model init
	virtual void initProd();			 // product init
	virtual const TreeSlice & getValue(int step) const;
	virtual void     update(int & step, FDProduct::UpdateType);
	virtual DateTime getStartDate(void) const { return startDate;}
	virtual bool     isElementary(void) const {return true;}
	virtual void     recordOutput(Control*, YieldCurveConstSP, Results*) {}

private:
	TreeSliceSP mainSlice;
	TrancheContingentLegConstSP inst;
	SpreadLossTree*  tree;

    DateTime startDate;

};
typedef refCountPtr<TrancheContingentLegProd> TrancheContingentLegProdSP;



/************************************************************************/
/* Fee leg of CDO with pair of strikes                           */
/************************************************************************/
class TREE_DLL TrancheFeeLeg :	public CObject, 
	virtual public FDModel::IIntoProduct
{
public:
	static CClassConstSP const TYPE;

	/** public constructor */
	TrancheFeeLeg(
		ICreditFeeLegSP feeLeg,				// underlying fee leg
		double lowStrike,					// as fraction of 1
		double highStrike,					// as fraction of 1
		bool recoverNotional,
		const IBadDayAdjusterConstSP bda
		);

	/** returns critical dates */
	virtual DateTimeArraySP getCritDates() const;

	const string& getName() const { return name; }

	friend class SpreadLossTree;
	friend class TrancheFeeLegProd;

	

protected:

	// FIELDS
	/** unique name */
	string name; 
	/** underlying contingent leg */
	ICreditFeeLegSP feeLeg;
	/** low strike of tranche as fraction of 1 */
	double lowStrike;
	/** high strike of tranche as fraction of 1*/
	double highStrike;

	/** Do we write down the notional by the recovery amount (as well as loss) */
	// e.g. is true for index tranches
	bool recoverNotional;

	/** bad day adjuster */
	const IBadDayAdjusterConstSP bda;

private:
	static void load(CClassSP& clazz);

	/** private constructor */
	TrancheFeeLeg();

	static IObject* defaultConstructor(void) { return new TrancheFeeLeg(); }

	virtual FDProductSP createProduct( FDModel * model ) const;
	

};




/************************************************************************/
/* product for TrancheFeeLeg                                     */
/************************************************************************/
class TREE_DLL TrancheFeeLegProd : public FDProduct 
{
public:
	TrancheFeeLegProd(const TrancheFeeLegConstSP &inst, FDModel* model);

	virtual void init(Control*) const;	// model init
	virtual void initProd();			 // product init
	virtual const TreeSlice & getValue(int step) const;
	virtual void     update(int & step, FDProduct::UpdateType);
	virtual DateTime getStartDate(void) const 
    { 
        return startDate;
    }
	virtual bool     isElementary(void) const {return true;}
	virtual void     recordOutput(Control*, YieldCurveConstSP, Results*) {}

private:
	/** slice containing fee leg values */
	TreeSliceSP mainSlice;

	TrancheFeeLegConstSP inst;
	
	SpreadLossTree* tree;

	/** start date for trade */
    DateTime startDate;


	/** projection of risky dates within  model timeline*/
	// i.e  returns an array of size timLine.size() giving the index of timeLine[i]
	// inside riskyDates[] or -1 if this element is out of riskyDates().
	vector<int>  riskyDateIdx;

	/** riskless date projection within timeline*/
	// see above for details
	vector<int>  risklessDateIdx;

	/** risky cashflow amounts for the fee leg */
	DoubleArray  riskyFlowAmounts; 

	/** riskless cashflow amounts for the fee leg */
	DoubleArray  risklessFlowAmounts; 

	/** slice containing outstanding fraction at each tree node (= 1 - outstandingNotional/trancheSize) */
	TreeSliceSP outstandingFraction;

};
typedef refCountPtr<TrancheFeeLegProd> TrancheFeeLegProdSP;


DRLIB_END_NAMESPACE
#endif
