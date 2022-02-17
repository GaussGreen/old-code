//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : OUAuxSpreadProcess.hpp
//
//   Description : OU auxiliary spread process. Instance of IAuxSpreadProcess
//
//   Author      : Matthias Arnsdorf
//
//   Date        : September 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_OU_AUX_SPREAD_PROCESS_HPP
#define QLIB_OU_AUX_SPREAD_PROCESS_HPP

#include "edginc/IAuxSpreadProcess.hpp"
#include "edginc/ICDSParSpreads.hpp"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL OUAuxSpreadProcess: public CObject,
	public virtual IAuxSpreadProcess 
{

public:
	static CClassConstSP const TYPE;

	/** Destructor */
	virtual ~OUAuxSpreadProcess();

	/** Called immediately after object constructed */
	virtual void validatePop2Object();    

	/// INTERFACE METHODS /////////////////////////////////////
	
	/** chance to pass market data to the spread process */
	virtual void getMarketData(const MarketData* market,
		int stepsPerYear // number of timeline steps per year
		);   

	/** initialise the model. Called prior to any requests. */
	virtual void setupModel(
		TimeLineSP timeLine		// timeline to use
		);

	/** updates tree variables at each time step */
	virtual void update(int timeStep, FDProduct::UpdateType type);

    /** updates the allowed spread values at the time step */
    virtual void updateSpread(int timeStep, FDProduct::UpdateType type);

	/** get the top range limits for all time points */
	virtual IntArray getTopLimits() const
	{
		return topLimits;
	}

	/** get the bottom range limits for all time points */
	virtual IntArray getBotLimits() const
	{
		return botLimits;
	}
	
	/** get range of current spread slice */
	virtual TreeSliceGeneral::RangeSP getRange() const 
	{
#ifdef DEBUG
		if(!timeline)
		{
			throw ModelException("setupModel has not been called","OUAuxSpreadProcess::getRange" );
		}
#endif
		return range;
	}

	/** get dates at which spread values change 
	* can be called after getMarketData
	*/
	virtual DateTimeArray getSpreadDates() const
	{
#ifdef DEBUG
		if(fwdDates.size() == 0)
		{
			throw ModelException("fwdDates have size 0" );
		}
#endif
		return fwdDates;
	}

	/** returns transition probability between two spread states at time step T1 (assumed to be currStep)
		and T2 (assumed to be currStep+1). The possible spread states at T2 are labeled by the branche Index
		which ranges from 0 to numBranches-1. 
		The Probability is conditional on the number of defaults
		over this period (transient contagion) */
	virtual double transProb(
		int sIdxT1,			// spread index at time index T1
		int branchIdx,		// index of tree branch (0 offset)
		int numDefaults		// number of defaults that have occurred between T1 and T2
		) const
	{
#ifdef DEBUG
		if(!timeline)
		{
			throw ModelException("setupModel has not been called","OUAuxSpreadProcess::transProb" );
		}
#endif

		return (*transitionProbs)[sOffset(sIdxT1)][Maths::min(numDefaults, maxDefaults)][branchIdx];
	}

	/** returns entire transitionProbs array */
	virtual DoubleArrayArrayArraySP getTransitionProbs()
	{
		return transitionProbs;
	}

	/** return the spread value at current time step and for given index.
		Can be called after update */
	virtual double spread(int index) const 
	{
#ifdef DEBUG
		if(!timeline)
		{
			throw ModelException("setupModel has not been called","OUAuxSpreadProcess::spread" );
		}
#endif
		return spreads[sOffset(index)];
	}


	/** returns index in spreads array (0 offset) given tree slice index (centred around 0) */
	virtual int sOffset(int sIdx) const
	{
		return sIdx - minRangeBot;
	}

	/** returns the node at T+1 reached by top branche emanating from node sIdx at T */ 
	virtual int outerRangeTop(int sIdx, int numDefs = 0) const
	{
#ifdef DEBUG
		if(!timeline)
		{
			throw ModelException("setupModel has not been called","OUAuxSpreadProcess::outerRangeTop" );
		}
#endif
		return sIdx + shifts[sOffset(sIdx)][Maths::min(numDefs, maxDefaults)] + (DEFAULT_NUM_BRANCHES - 1)/2;
	}

	/** returns the node at T+1 reached by bottom branche emanating from node sIdx at T */ 
	virtual int outerRangeBot(int sIdx, int numDefs = 0) const
	{
#ifdef DEBUG
		if(!timeline)
		{
			throw ModelException("setupModel has not been called","OUAuxSpreadProcess::outerRangeBot" );
		}
#endif
		return sIdx + shifts[sOffset(sIdx)][Maths::min(numDefs, maxDefaults)] - (DEFAULT_NUM_BRANCHES - 1)/2;
	}

	/** returns the entire outer range array */
	virtual IntArrayArraySP getOuterRangeBot()
	{
		return outerBot;
	}

	/** returns the number of branches (assumed fixed) */
	virtual int getNumBranches()
	{
		return DEFAULT_NUM_BRANCHES;
	}

	/** returns maximum number of defaults used for contagion = 2nd dim of transitionPtobs array */
	virtual int getMaxDefaults()
	{
		return maxDefaults;
	}

private:

	// METHODS /////////////////////////////////////////////

	

	/** calculates tree limits in fwd loop */
	void calcLimits();

	/** calculate shifts array */
	void calcShifts(int timeStep, double dt);

	/** calculate transition probabilites */
	void calcTransProbs(
		double scaledMean,			// mean over period divided by dR
		double scaledVarVal,		// vol^2dt/dR^2
		int shift,					// shift at current node
		DoubleArray & transProbs	// array of transition probs for each branch	
		);

	/** default value for maximum number of defaults */
	static const int DEFAULT_MAX_DEFAULTS;

	/** default value for number of branches for tree (assumed fixed) */
	static const int DEFAULT_NUM_BRANCHES;

	/** default val for minSpreadVal */
	static const int DEFAULT_MIN_SPREAD_VAL;

	/** default val for q */
	static const double DEFAULT_Q;


	/** if |prevDT - dt| < DT_ACCURACY they will be assumed same */
	static const double DT_ACCURACY;

	/** scale factor for dR: dR^2 = SCALE * avegrageDT * vol^2  */
	static const double SCALE;

	// For reflection
	static void load (CClassSP& clazz);

	static IObject* defaultOUAuxSpreadProcess();

	/** private constructor */
	OUAuxSpreadProcess();

	/// FIELDS ////////////////////////////

	/** inital index fwd curve */
	ICDSParSpreadsWrapper indexCurve;

	

	/** volatility */
	double vol;

	/** mean reversion */
	double meanReversion;

    /** transient contagion jump parameter as % of level [default = 0]*/
    double contagionJump;

	/** q-mapping parameter */
	double q;

	/** minimum value of allowed spreads */
	double minSpreadVal;

	// TRANSIENT FIELDS /////////////////	

	/** dates for the fwd curve */
	DateTimeArray fwdDates;

	/** clean spreads associated to fwd curve */
	DoubleArraySP cleanSpreads;

	/** model timeline */
	TimeLineSP timeline;

	/** maximum number of defaults which we consider for contagion.
		A higher of number of defaults will no longer change the transitition probabilities */
	int maxDefaults;

	/** transition probabilities 
		(*transitionProbs)[sIdxT1][numDefs][branchIdx] */
	DoubleArrayArrayArraySP transitionProbs;

	/** range for spread tree */
	TreeSliceGeneral::RangeSP range;

	/** lowest value for bottom of ranges (= rangeBot at final slice) */
	int minRangeBot;

	/** if are at node sIdx at T then center node at T+1 = sIdx + shifts[sIdx][numDefs] */
	vector<IntArray> shifts;

	/** outer range bottom limits outerBot[sIdx][numDefs] */
	IntArrayArraySP outerBot;

	/** allowed spread values at any time step */
	DoubleArray spreads;

	/** q mapping of OU process values */
	DoubleArray qMap;
	

	/** fwds intensities corresponding to index spreads at timeline points*/
	DoubleArray fwds;

	/** average time step in timeline */
	double avgDT;
	/* spacing between points Gaussian auxilliary process R */
	double dR;

    /* contagionJump/dR */
    double scaledJump;

	/** top tree limits for each timeline point */
	IntArray topLimits;

	/** bottom tree limits for each timeline point */
	IntArray botLimits;

	/** variance over time step scaled by dR^2
    *   scaledVar[i] = vol^2(t_{i+1} - t_i ) / dR^2 */
	DoubleArray scaledVar;

	/** previous time step. Used to cehck if need to recalc probabilities */
	double prevDT;



};

DECLARE(OUAuxSpreadProcess);

DRLIB_END_NAMESPACE

#endif


