//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ConstAuxSpreadProcess.hpp
//
//   Description : Constant (trivial) auxiliary spread process. Instance of IAuxSpreadProcess
//
//   Author      : Matthias Arnsdorf
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONST_AUX_SPREAD_PROCESS_HPP
#define QLIB_CONST_AUX_SPREAD_PROCESS_HPP

#include "edginc/IAuxSpreadProcess.hpp"
#include "edginc/ICDSParSpreads.hpp"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL ConstAuxSpreadProcess: public CObject,
	public virtual IAuxSpreadProcess 
{

public:
	static CClassConstSP const TYPE;

	/** public constructor */
	ConstAuxSpreadProcess(double intensity);

	/** Destructor */
	virtual ~ConstAuxSpreadProcess();

	/** Called immediately after object constructed */
	virtual void validatePop2Object();    

	/// INTERFACE METHODS /////////////////////////////////////
	
	/** chance to pass market data to the spread process */
	virtual void getMarketData(const MarketData* market, int stepsPerYear); 
	

	/** initialise the model. Called prior to any requests. */
	virtual void setupModel(
		TimeLineSP timeLine		// timeline to use
		);

	/** updates tree variables at each time step */
    virtual void update(int timeStep, FDProduct::UpdateType type) 
    {
        //update the intensity
        intensity = fwds[timeStep];
    };

    /** updates the allowed spread values at the time step */
    virtual void updateSpread(int timeStep, FDProduct::UpdateType type)
    {
        //update the intensity
        intensity = fwds[timeStep];
    };

	/** get range of current spread slice */
	virtual TreeSliceGeneral::RangeSP getRange() const 
	{
		if(!range)
		{
			throw ModelException("range not initalised", "ConstAuxSpreadProcess::getRange");
		}
		return range;
	}

	/** get dates at which spread values change 
	* can be called after getMarketData
	*/
	virtual DateTimeArray getSpreadDates() const;

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
		return 1.0;
	}

	/** returns entire transitionProbs array */
	virtual DoubleArrayArrayArraySP getTransitionProbs()
	{
		return transProbs;
	}

	/** return the spread value at current time step and for given index.
	Can be called after update */
	virtual double spread(int index) const;

    

	/** returns index in spreads array (0 offset) given tree slice index (centred around 0) */
	virtual int sOffset(int sIdx) const
	{
		return 0;
	}

	/** returns the node at T+1 reached by top branche emanating from node sIdx at T */ 
	virtual int outerRangeTop(int sIdx, int numDefs = 0) const
	{
		return sIdx;
	}

	/** returns the node at T+1 reached by bottom branche emanating from node sIdx at T */ 
	virtual int outerRangeBot(int sIdx, int numDefs = 0) const
	{
		return sIdx;
	}

	/** returns the entire outer range array */
	virtual IntArrayArraySP getOuterRangeBot()
	{
		return outerBot;
	}

	/** get the top range limits for all time points */
	virtual IntArray getTopLimits() const
	{
		IntArray topLimits;
		topLimits.resize(numDates,0);
		return topLimits;
	}

	/** get the bottom range limits for all time points */
	virtual IntArray getBotLimits() const
	{
		IntArray botLimits;
		botLimits.resize(numDates,0);
		return botLimits;
	}

	/** returns the number of branches (assumed fixed) */
	virtual int getNumBranches()
	{
		return 1;
	}
	/** returns maximumum number of defaults used for contagion = 2nd dim of transitionPtobs array */
	virtual int getMaxDefaults()
	{
		return 0;
	}
private:

	// For reflection
	static void load (CClassSP& clazz);

	static IObject* defaultConstAuxSpreadProcess();

	/** private constructor */
	ConstAuxSpreadProcess();

	///  INPUT FIELDS ////////////////////////////

	/** initial default intensity value */
	double intensity;

    /** inital index fwd curve */
    ICDSParSpreadsWrapper indexCurve;

	// OTHER FIELDS ////////////////

	DateTime valueDate;

    /** dates for the fwd curve */
    DateTimeArray fwdDates;

    /** clean spreads associated to fwd curve */
    DoubleArraySP cleanSpreads;

    /** fwds intensities corresponding to index spreads at timeline points*/
    DoubleArray fwds;

	/** number of dates in timeline */
	int numDates;

	/** spreads slice range */
	TreeSliceGeneral::RangeSP range;

	DoubleArrayArrayArraySP transProbs;

	IntArrayArraySP outerBot;
};

DECLARE(ConstAuxSpreadProcess);

DRLIB_END_NAMESPACE

#endif


