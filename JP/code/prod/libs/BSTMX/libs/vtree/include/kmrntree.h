/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kmrntree.h
 * Function:	
 * Author:	
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kmrntree_H
#define	_kmrntree_H

#define FLOORIDX    0   //spot vol interpolation flag
#define CEILIDX     1

#include "kvtree.h"
#include "kcbank.h"

/**
 * A class to be used as base class for any model based
 * on an N-factor mean reverting tree.
 */
class	KMrNTree : public KVTree  {
public:
	//------------------------------------------------------
	// KVTree base class methods
	//------------------------------------------------------

	/**  Default constructor. */
			KMrNTree();

	/** Destructor. */
virtual			~KMrNTree();

	/** Inserts a critical date. */
virtual void            Insert(TDate date);

	/** Inserts a KZeroReset (dummy at this level). */
virtual void		Insert(const KZeroReset&, bool isCrit=TRUE);

	/** Inserts a KRateReset, and return the modelReset date  (dummy).*/
virtual TDate		Insert(const KRateReset&, bool isCrit=TRUE);

	/** Inserts a KRateReset between reset date and end date
	 *  and return the modelReset date  (dummy).*/
virtual TDate		Insert(const KRateReset&, TDate, bool isCrit=TRUE);

	/** Gets a KZeroReset (dummy). */
virtual void		Get(KTSlice& ts, const KZeroReset&);

	/** Gets a KRateReset (dummy). */
virtual void		Get(KTSlice& ts, const KRateReset&);

	/** Get the value date of specified curve (dummy) */
virtual TDate            GetValueDate(int idx) {return TPToday();}

	/** Return tree dates within specified range */
virtual KVector(TDate)	GetTDatesRange(TDate stDate, TDate endDate);

	/** Return the last tree date */
virtual TDate		GetLastDate() {return mLastDate;}

	/**
	 * Setup tree timeline.
	 */
virtual	void		SetUpTimeline();

	/**
	 * Calibrate should be called after all the product cirtical
	 * dates are inserted. It performs simply the TreeSetUp
	 */
virtual	void		Calibrate();

	/** Updating MR tree sets time point and
	 *  computes the probabilities 
	 */
virtual	void		Update(int tpIdx);

	/** Virtual function of KVTree */
virtual	int		TPNum()
			{ return NbTP;}

	/** Virtual function of KVTree */
virtual	TDate		TPToday()
			{ return mTodayDate;}

	/** Virtual function of KVTree */
virtual	int		TPIdxCurrent()
			{ return tpIdxCurrent;}

	/** Virtual function of KVTree */
virtual	TDate		TPDateCurrent()
			{ return TPDates[tpIdxCurrent];}

	/** Returns the year fraction to the current timepoint. */
virtual double          TPTimeCurrent();

	/** Virtual function of KVTree */
virtual int             TPIdx(TDate date);

	/** Virtual function of KVTree */
virtual	TDate		TPDate(int tpIdx)
			{ return TPDates[tpIdx];}

	/** Virtual function of KVTree */
virtual	int		TPIsCriticalDate(int tpIdx)
			{ return TRUE; }

	/** Virtual function of KVTree */
virtual KTSlice&	TSliceDev(KTSlice&, const String&);

	/** Virtual function of KVTree */
virtual KTSlice&	TSliceEv(KTSlice& ts);


	/** Virtual function of KVTree. */
virtual KTSlice&	TSliceFwd(KTSlice&);

	/** Allocate a slice of dimension mNbFactor.
	 *  curve index and slice name are not defined.
	 */
virtual	KTSlice&	TSliceCreate(KTSlice&); 

	/** Virtual function of KVTree. */
virtual	void		TSliceDestroy(KTSlice&);

	/** Virtual function of KVTree. */
virtual	double		TSliceGetCenter(KTSlice&);

	/** Virtual function of KVTree. */
virtual	void		TSliceSetCenter(KTSlice &ts, double value);


	/** Virtual function of KVTree. */
virtual bool            TSliceCompare(KTSlice&, 
				      KTSlice&, 
				      KTSComp);

	/** Virtual function of KVTree. */
virtual	KTSlice&	TSliceScalarOper(KTSlice&, double, KOper);

	/** Virtual function of KVTree. */
virtual	KTSlice&	TSliceUnaryOper(KTSlice&, const KTSlice&, KOper);

	/** Virtual function of KVTree. */
virtual	void		TSlicePut(KTSlice&, ostream& os, int minMaxOnly);

	/** Virtual function of KVTree. */
virtual	void		TSliceSpecialOper(KTSlice&, char *what, ...);

	/**
	 * Create and set the begining of slice node.
	 */
virtual KTSliceNode&	TSliceNodeBegin(KTSlice&, KTSliceNode&);
 
	/**
	 * Set the end of slice node.
	 */
virtual bool		TSliceNodeEnd(KTSliceNode &tsNode)
			{ return ((int *)tsNode.mIndex)[0];}

	/**
	 * Move to next slice node.
	 */
virtual KTSliceNode&	TSliceNodeNext(KTSliceNode&);

	/** Virtual function of KVTree. */
virtual double&         TSliceAccessNodeValue(KTSlice&, KTSliceNode&);
 
	/** Virtual function of KVTree. */
virtual double          TSliceGetNodeStepMax(KTSlice&, KTSliceNode&);

	/** Virtual function of KVTree. */
virtual	void		TSliceNodePut(KTSliceNode&, ostream& os);


	//------------------------------------------------------
	// KMrNTree specific public methods to initialize
	//------------------------------------------------------

	/**
	 * Initialize the tree parameters. Mainly set up model parameters.
	 * The computation of timeline, jump size, orthogonalize factors, etc.
	 * will be done in Calibrate AFTER all the product
	 * related critical dates are inserted. 
	 */
	void	Initialize(
		TDate todayDt,		// (I) today's date
		int ppy,		// (I) period per year
		double smoothFact,	// (I) smoothing factor
		char EoI,		// (I) Equal or increasing time steps
		int numStdev,		// (I) number stdev to cut tree
		int numFact,		// (I) number of factors
		double *factMr,		// (I) factor mean reversions
		double *factWeight,	// (I) factor weights
		double *factCorr,	// (I) factor correlation [nbFact]
		KVector(TDate) &volDates,// (I) volatility dates
		KVector(KVector(double)) &factVol);// (I) factor spot vol 

	
	/**
	 * Initialize factor volatilities, simply assign
	 * the factor vols.  Used by CET only.
	 */
	void	InitializeFactorVol(
		KVector(TDate) &volDates,	    // (I) volatility dates
		KVector(KVector(double)) &factVol); // (I) factor spot vol 

	/**
	 * Update factor volatilities, including memory allocation for 
	 * tree limits and probabilities that depends on factor vols. 
	 * Used by CET only.
	 */
	void	UpdateFactorVol(
		KVector(TDate) &volDates,	    // (I) volatility dates
		KVector(KVector(double)) &factVol); // (I) factor spot vol 


	/**
	 * Free tree memories (such as tree limits, and transitional
	 * probabilities, etc.) that depends on factor vols.
	 * This is used only by CET routine for effective
	 * tree memory management.
	 */
	void	ClearTreeVolMem();	

	/**
	 * Return the vector of spot volatility dates.
	 */
	KVector(TDate) SpotVolDates() const
		{ return mVolDates; }

    /**
     * Return the factor volatilities.
     */
    KVector(KVector(double)) SpotVolatilities() const {
        return mFactVol;
    }


	/**
	 * Returns the matrix of spot volatility (numFact * numVolDates).
	 (/
	KVector(KVector(double)) SpotVols() const
		{ return mFactVol; }


protected:
	//------------------------------------------------------
	// Protected methods used internally
	//------------------------------------------------------

	/** Allocate a slice of specified dimension.
	 *  curve index and slice name are not defined.
	 */
virtual	KTSlice&	TSliceDimCreate(KTSlice &ts, int nDim); 

	/** Performs a slice dimension expansion */
virtual	KTSlice&	TSliceExpand(
				int tpIdx,	    // (I) time point index
			     	const KTSlice& tsLo,// (I) slice in lower dim
				KTSlice& tsHi);	    // (I/0) slice in higher dim

	/** Performs a slice dimension degeneration */
virtual	KTSlice&	TSliceDegenerate(
				int tpIdx,	    // (I) time point index
			     	const KTSlice& tsHi,// (I) slice in hiher dim
			     	KTSlice& tsLo);	    // (I/0) slice in lower dim

	/** Calculate array offsets. There are three components:<BR>
	 *  1)  node indices run into [-halfwidth, halfwidth] in each dimension 
	 *  whereas memory is allocated as [0, 2 * halfwidth]. This gives an 
	 *  offset of halfwidth.<BR>
	 *  2)  memory is allocated linearly in dimension 2 and 3. This gives an
	 *  offset of i * width[1] in dimension 2 and i * width[1] * width[2]
	 *  + j * width[1] in dimension 2.<BR>
	 *  3)  the axis of the ellipse has a slope. The coordinate of the axis
	 *  is given by (bottom2[i]+top2[i])/2 at index [i] in dimension 2 and
	 *  (bottom3[i][j]+top3[i][j])/2 at index [i][j] in dimension 3.<BR>
	 */
	int		NodeOffset(int Dim, int i, int j, int t);

	/**
	 *
	 */
	void		NodeCheck(int Dim, int j0, int j1, int j2,
				int t, int inOut);



	/**
	 * Computes the the local max diff. Used for smoothing.
	 */
	double		GetIndexStep(double *Index, int Dim,
				int i, int j, int k, int t);


	/** Set up tree perform the following:
	 *  1. Interpolate the spot vols on each time step and Compute the
	 *     jump step size.
	 *  2. Compute the orthogonal factor covariances 
	 *     (once the timeline is setup)
	 *  3. Compute the treelimits.
	 */
virtual	void		TreeSetUp();

	/** Free all the memory allocated during tree initialization
	 *  and calibration.
	 *  It is used only by CET tool to free the tree memory
	 *  during each iteration, so to allow re-initialization
	 *  of tree with new volatilities.
	 *  It is identical to destructor ~KPirTree(), but keep
	 *  the object intact.
	 */
virtual	void		DeleteMemory();

	/**
	 * Allocate a new slice as a double* to be acessed using NodeOffset.
	 */
	double*		sliceNew(int sliceDim);

	/** Deletes a pointer allocated by newSlice */
	void		sliceDelete(double *ptr);

	/**
	 * Prints a pointer allocated by newSlice.
	 * @param outLim If TRUE, prints the outer limits,
	 * otherwise Prints the inner limits.
	 */
	void		slicePrint(double *ptr, int sliceDim,
				int tpIdx, int outLim, ostream& os);


	/** Performs a scalar operation on a slice . */
	void		sliceScalarOper(
			double *data,		// (B) slice data 
			const int sliceDim,	// (I) slice dimension
			double argument,	// (I) argument of operation
			KOper oper);		// (I) operation type

	/** Performs a unary operation. */
	void		sliceUnaryOper(
			double *data,		// (B) slice data
			int sliceDim,		// (I) slice dimension
			const double *data1,	// (I) oper argument (same dim)
			KOper oper);		// (I) operation type

	/** Performs a special slice operation. */
	void		sliceSpecialOper(
			double *data,		// (I) slice data
			int sliceDim,		// (I) slice dimension
			const char *what,	// (I) operaton type
			...);			// (I) arguments ...

	/** Performs a slice dimension degeneration */
	void		sliceDegenerate(
			int	tpIdx,	  // (I) time point index
			int	sHiDim,	  // (I) higher slice dimensions
			int	sLoDim,	  // (I) lower degenerated slice dims
			double	*sHiVal,  // (I) slice values in sHiDim dims.
			double	*sLoVal); // (O) slice values in sLoDim dims.
	
	/** Performs a slice dimension expansion */
	void		sliceExpand(
			int	tpIdx,	  // (I) time point index
			int	sHiDim,	  // (I) higher slice dimensions
			int	sLoDim,	  // (I) lower degenerated slice dims
			double	*sLoVal,  // (I) slice values in sLoDim dims.
			double	*sHiVal); // (O) slice values in sHiDim dims.
	
	/**
	 * Takes the expected (possibly discounted) value of a slice.
	 */
	void		sliceEv(
			double *sliceVal,	// (B) values to be discounted 
			int sliceDim,		// (I) slice dimension
			double *discVal,	// (I) discounting ts
			int tpIdx);		// (I) current time period 

	/**
	 * Forwards a time slice (adjoint operation of sliceFw)
	 */
	void		sliceFw(
			double *sliceVal,	// (B) values to be forwarded
			int sliceDim,		// (I) slice dimension
			double *discVal,	// (I) discounting ts
			int tpIdx);		// (I) current time period 

protected:
	//------------------------------------------------------
	// Internal tree data
	//------------------------------------------------------
					// --- Time line
	TDate	mTodayDate;		// todatDate
	int	NbTP;			// Nb of time points in the tree
	TDate	*TPDates;		// Date of each time point
	double	*TPTimes;		// Time offset of each time point
	double	*Length;		// Length of time steps (ACT/365)
	double	*LengthJ;		// Time step for jump size [nbTp]

	TDate	mLastDate;		// Last tree date

					// --- Model parameters
	int	mNbFactor;		// Number of factors
	double	*mAlpha;		// Factors weight [nbFact]
	double	*mBeta;			// Mean reversions [nbFact]
	double	**mRho;			// Correlation [nbFact][nbFact]
	double	**mSigma;		// Spot volatilities [nbTp][nbFact]

					// --- Model and volatility 
	double	**mAweight;		// Weights at time point [nbTp][idx]

					// --- Tree geometry 
	int	mNbSigmaMax;		// Nb of std dev to cut the tree
	int	*mWidth;		// Ellipsoid width [nbFact]
	int	*mHalfWidth;		// Ellipsoid half width [nbFact]

	int	*mTop1,
		*mBottom1;		// Limits of 1D tree
	int	**mTop2,
		**mBottom2;		// Limits of 2D tree
	int	***mTop3,
		***mBottom3;		// Limits of 3D tree
	int	*mOutTop1,
		*mOutBottom1;		// Outer limits
	int	**mOutTop2,
		**mOutBottom2;   
	int	***mOutTop3,
		***mOutBottom3; 


					// --- DEV DATA STRUCTURE 
	int	tpIdxCurrent;
	int     *Shift1;
	int     *Shift2;
	int     *Shift3;

					// Probabilities
	double  *pu, *p0, *pd;		// 1D
	double	*qu, *q0, *qd;		// 2D
	double	*ru, *r0, *rd;		// 3D

					// Working memory
	double	*NewPrice;
	int	NbTPMax;		// max number of time points

	KVector(TDate) 	mVolDates;      	// benchmark volatility dates
	KVector(KVector(double)) mFactVol;	// Benchmark spot vols

					/** Smoothing factor.*/
	double	mSmoothFact;

                    // To solve the inconsistency between fix3 and tmx3
                    // Add additional flags in build MRTree
    int     spotVolInterp;  // spot vol interoplation flag: FLOORIDX or CEILIDX
    int     outStdDevs;     // width of the outer ellipse in std devs.

private:
	//------------------------------------------------------
	// Private methods and data 
	//------------------------------------------------------

	/**
	 * Sorts all the critical dates in ascending order 
	 * and sets up the time line.
	 * @param EoI Specifies if the times steps are equal ('E')
	 * or increasing ('I').
	 * @param ppy Number of periods per year.
	 */
	void	InitializeTimeline(
		char EoI,		// (I) Equal or increasing time steps
		int ppy);		// (I) period per year

	/**
	 * Checks time line is valid.
	 */
	void	CheckTimeline(
		int ppy);		// (I) period per year


	/**
	 * Setup factor volatilities and tree geometry. 
	 * Used by TreeSetUp and CET.
	 */
	void	SetUpFactorVol();

	/** 
	 * Initialize tree memory allocation that depends only
	 * on timeline, but independent on factor vols.
	 * The break-down of TreeSetUp routine is for efficient
	 * memory management in the CET.
	 */
	void	MrNTreeLimitsNFactMemAlloc();

	/**
	 * Compute factor vol and tree geometry only, assuming memories
	 * that only depends on timeline have been setup. 
	 * Memory for tree limits and probabilities that depends
	 * on factor vol will be allocated on the run.
	 * This routine will be used by CET as well as Initialize()
	 * for effeciency.
	 */
	void	ComputeFactorVolAndTreeLimit();
	
	/**
	 * Compute orthogonal variance
	 */
	void	OrthogonalizeVol();

	/**
	 * Checks limits are valid.
	 */
	void	LimitsCheck();


protected:
					// --- temp stored before all critical
					//     dates are inserted
	int  	mPPY;                   // period per year
	char 	mEoI;                   // Equal or increasing time steps

private:

friend	int	MrNTreeComputeTransitionProba(KMrNTree &,
			int, int, int, double*, double**, double*, double*,
			int, int, int*, int*, int**, int**,
			int, int, int*, int*, int**, int**,
			double*, double*, double*, double*,
			double*, double*, double*, double*,
			double*, int*, int*, int*);


};


typedef KMap(TDate, KDateInterval)::value_type TDate_KDateInterval;
typedef KMap(TDate, KRateReset)::value_type TDate_KRateReset;
typedef KMap(TDate, KTSlice*)::value_type TDate_KTSlice;

//
// Macros used in source
//


#define	DEBUG_LEVEL_CET		05
#define	DEBUG_LEVEL_TIMELINE	12
#define	DEBUG_LEVEL_VOLATILITY	13
#define	DEBUG_LEVEL_DRIFT	15
#define	DEBUG_LEVEL_GEOMETRY	25
#define	DEBUG_LEVEL_PROBA	30
#define	DEBUG_LEVEL_CBANK	40
#define	DEBUG_LEVEL_GET		40
#define	DEBUG_LEVEL_DEV		45
#define	DEBUG_LEVEL_KIO		50












//--------------------------------------------------------------
// Principal defs used in source
//

#ifdef	_kmrntree_SRC

// Constants
#define	JUMPCOEFF	3.0	// Jump size coefficient
#define	NFMAX		8
#define MIN_SPOT_VOL	(1e-5)
#define	NUMSTDEV	10
#define	ERROR		(1e-7)  // avoid to be divided by 0


//--------------------------------------------------------------
// Calculates the vector index for the A coefficients vector
// from the matrix coefficients (i,j) using the enumeration
//	0 1 3 .
//	1 2 4 .
//	3 4 5 .
//	. . . .
//

inline int AIDX(int i, int j)
{
	return (i <= j ? \
		j*(j+1)/2 + i : 
		i*(i+1)/2 + j );
}


//--------------------------------------------------------------
// Calculates the vector index for the coreeleation coefficients vector
// from the matrix coefficients (i,j) using the enumeration
//	- 0 1 3 .
//	0 - 2 4 .
//	1 2 - 5 .
//	3 4 5 - .
//	. . . . - 

inline int RIDX(int n, int i, int j)
{
	return ((i) <= (j) ? \
	((n)*((n)-1)/2 - ((n)-(i))*((n)-(i)-1)/2 + ((j)-(i)-1)) : 
	((n)*((n)-1)/2 - ((n)-(j))*((n)-(j)-1)/2 + ((i)-(j)-1)));
}


//--------------------------------------------------------------
// Other useful macros

/**  Inline version of max. */
inline	double	Max(double a, double b)
	{ return (a <= b ? b : a); }

/**  Inline version of max. */
inline	int	Max(int a, int b)
	{ return (a <= b ? b : a); }

/**  Inline version of min. */
inline	double	Min(double a, double b)
	{ return (a <= b ? a : b); }

/**  Inline version of min. */
inline	int	Min(int a, int b)
	{ return (a <= b ? a : b); }

/**  Square root of absolute value. */
inline	double	SqrtM(double x)
	{ return (x <= DBL_EPSILON ? DBL_EPSILON : sqrt(x)); }

/**  Rounds to the nearest integer. */
inline	int	Ceil(double x)
	{ return (int) ceil(x); }

/**  Rounds to the nearest integer. */
inline	int	NearInt(double x)
	{ return (int) (x + (x >= 0. ? 0.5 : -0.5)); }

/**  Exponential decay function. */
inline double  ExpDecay(double  a, double  t)
	{ return (fabs (a*t) < DBL_EPSILON ? 1e0 :
		(1e0 - exp (- a * t)) / (a*t)); }

/**  Day count fraction ACT/365 */
inline	double	YrsDiff(TDate date1, TDate date2)
	{ return ((double) date2 - (double) date1) / 365e0; }

	/*{ double dcf; IF_FAILED_THROW( GtoDayCountFraction(date1, date2,
		GTO_ACT_365F, &dcf)); return(dcf); }*/


extern "C" {
#include "drltime.h"
};

/**  Prints a date or N/A if the date number is less than one */
inline	const char*	printDate(TDate date)
{
	return (date > 1L ? DrlTDatePrint(NULL, date) : "NA");
}


/** Convert KOper to string for printing */
inline	const char* printOper(KOper oper)
{
	switch (oper){ \
	case COPY: \
		return "COPY"; \
		break; \
	case ADD: \
		return "ADD"; \
		break; \
	case SUB: \
		return "SUB"; \
		break; \
	case MULT: \
		return "MULT"; \
		break; \
	case POW: \
		return "POW"; \
		break; \
	case DIV: \
		return "DIV"; \
		break; \
	case MAX: \
		return "MAX"; \
		break; \
	case MIN: \
		return "MIN"; \
		break; \
	case GEQ: \
		return "GEQ"; \
		break; \
	case LEQ: \
		return "LEQ"; \
		break; \
	case STVAR: \
		return "STVAR"; \
		break; \
	default: \
		throw KFailure("unknown operation (%d).\n", oper); \
	} \
	return "COPY";
}


#endif /*_kmrntree_SRC*/


#endif /* _kmrntree_H */


