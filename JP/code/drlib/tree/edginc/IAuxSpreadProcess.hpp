//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IAuxSpreadProcess.hpp
//
//   Description : Interface for auxiliary spread porcess needed by spreadLossTree
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IAUX_SPREAD_PROCESS_HPP
#define QLIB_IAUX_SPREAD_PROCESS_HPP

//#include "edginc/Object.hpp"
//#include "edginc/AtomicArray.hpp"

#include "edginc/TreeSlice.hpp"
#include "edginc/FDProduct.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class IAuxSpreadProcess: public virtual IObject {

public:
	static CClassConstSP const TYPE;

	virtual ~IAuxSpreadProcess();

	/** chance to pass market data to the spread process */
	virtual void getMarketData(const MarketData* market, int stepsPerYear) = 0;   

	/** initialise the model. Called prior to any requests. */
	virtual void setupModel(
		TimeLineSP timeLine		// timeline to use
		) = 0;

	/** updates tree variables at each time step */
	virtual void update(int timeStep, FDProduct::UpdateType type) = 0;

    /** updates the allowed spread values at the time step */
    virtual void updateSpread(int timeStep, FDProduct::UpdateType type) = 0;

	/** get range of current spread slice */
	virtual TreeSliceGeneral::RangeSP getRange() const = 0;

	/** get dates at which spread values change 
	* can be called after getMarketData
	*/
	virtual DateTimeArray getSpreadDates() const = 0;


	/** returns transition probability between two spread states at time step T1 (assumed to be currStep)
	and T2 (assumed to be currStep+1). The possible spread states at T2 are labeled by the branche Index
	which ranges from 0 to numBranches-1. 
	The Probability is conditional on the number of defaults
	over this period (transient contagion) */
	virtual double transProb(
		int sIdxT1,			// spread index at time index T1
		int branchIdx,		// index of tree branch (0 offset)
		int numDefaults		// number of defaults that have occurred between T1 and T2
		) const = 0;

	/** returns entire transitionProbs array */
	virtual DoubleArrayArrayArraySP getTransitionProbs() = 0;

	/** return the spread value at current time step and for given index.
	Can be called after update */
	virtual double spread(int index) const =0;

	/** returns index in spreads array (0 offset) given tree slice index (centred around 0) */
	virtual int sOffset(int sIdx) const = 0;

	/** returns the node at T+1 reached by top branche emanating from node sIdx at T */ 
	virtual int outerRangeTop(int sIdx, int numDefs = 0) const = 0;

	/** returns the node at T+1 reached by bottom branche emanating from node sIdx at T */ 
	virtual int outerRangeBot(int sIdx, int numDefs = 0) const = 0;

	/** get the top range limits for all time points */
	virtual IntArray getTopLimits() const = 0;
	/** get the bottom range limits for all time points */
	virtual IntArray getBotLimits() const = 0;

	/** returns the entire outer range array */
	virtual IntArrayArraySP getOuterRangeBot() = 0;

	/** returns the number of branches (assumed fixed) */
	virtual int getNumBranches() = 0;
	/** returns maximumum number of defaults used for contagion = 2nd dim of transitionPtobs array */
	virtual int getMaxDefaults() =0;

private:
	// For reflection
	static void load (CClassSP& clazz);
};

DECLARE(IAuxSpreadProcess);

DRLIB_END_NAMESPACE

#endif
