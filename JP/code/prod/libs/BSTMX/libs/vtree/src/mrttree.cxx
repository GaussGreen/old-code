/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Th original model has been developped by Lionnel Pradier & London DR
 * 		This version has been adapted/modified by NY DR.
 ***************************************************************/
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"    /* Standard definitions & error hadling */

#define	_kmrntree_SRC
#include "kmrntree.h"

extern "C" {
#include "cgeneral.h"                   /* Has stdlib.h */
#include "bastypes.h"                   /* TDateList */
#include "cerror.h"                     /* GtoErrMsg */
#include "cmemory.h"                    /* MALLOC, FREE */
#include "extract.h"                    /* GtoExtractArray */
#include "macros.h"                     /* MAX */

#include "datelist.h"                   /* TDateList routines */
#include "ldate.h"                      /* GtoDtFwdAny */
#include "convert.h"                    /* GtoFormatDate */
#include "zerodate.h"                   /* TZeroDates */
#include "termtype.h"                   /* TFloatRateArray */

#include "drlmem.h"			/* DrlDoubleVectAlloc/Free */
#include "drlio.h"			/* */
#include "drltime.h"			/* DrlTDatePrint */
#include "drlvtype.h"			/* DrlVTypeVectAdd */
#include "drlsort.h"			/* DrlTDateArrayFloorIdx() */

};

int      debugLevel;

static	void	MrNTreeLimitsFree(
	int *Top1,			// (I) Upper limits of 1D tree
	int *Bottom1,			// (I) Lower limits           
	int **Top2,			// (I) Upper limits of 2D tree
	int **Bottom2,			// (I) Lower limits           
	int ***Top3,			// (I) Upper limits of 3D tree
	int ***Bottom3,			// (I) Lower limits           
	int nbFactor,			// (I) Number of factors      
	int NbTP);			// (I) Nb of time points      

static	void	MrNTreeLimitsNFactFree(
	int *Top1,			// (I) Upper limits of 1D tree
	int *Bottom1,			// (I) Lower limits           
	int **Top2,			// (I) Upper limits of 2D tree
	int **Bottom2,			// (I) Lower limits           
	int ***Top3,			// (I) Upper limits of 3D tree
	int ***Bottom3,			// (I) Lower limits           
	int nbFactor,			// (I) Number of factors      
	int NbTP);			// (I) Nb of time points      

static	void	MrNTreeLimitsNFactPrint(
	int *Width,			// (I) Ellipsoid widt
	int *HalfWidth,			// (I) Ellipsoid half width 

	int *Top1,			// (I) Upper limits of 1D tree
	int *Bottom1,			// (I) Lower limits           
	int **Top2,			// (I) Upper limits of 2D tree
	int **Bottom2,			// (I) Lower limits           
	int ***Top3,			// (I) Upper limits of 3D tree
	int ***Bottom3,			// (I) Lower limits           
	int nbFactor,			// (I) Number of factors      
	int NbTP);			// (I) Nb of time points      

//--------------------------------------------------------------
// Compute elemntary probabilities.

inline	void	ComputeProba(
	double JumpCoeff,		// (I) 
	double e,			// (I) expected value
	double *pUp,			// (O) probability up
	double *pMi,			// (O) probability middle
	double *pLo)			// (O) probability down
{
	*pUp = (JumpCoeff + e + e*e) * 0.5e0;
	*pLo = *pUp - e;
	*pMi = 1e0 - *pUp - *pLo;
}




//--------------------------------------------------------------
// Default consstructor.


KMrNTree::KMrNTree()
{

					// Timeline
	mTodayDate = -1L;
	NbTP = 0;
	TPDates = NULL;
	TPTimes = NULL;
	Length  = NULL;
	LengthJ = NULL;

	mLastDate = -1L;

					// Volatility info
	mNbFactor = 0;
	mAlpha = NULL;
	mBeta = NULL;
	mRho = NULL;
	mSigma = NULL;
	mAweight = NULL;

					// Tree geometry 
	mNbSigmaMax = -1;
	mWidth = NULL;
	mHalfWidth = NULL;

	mTop1 = NULL;
	mBottom1 = NULL;
	mTop2 = NULL;
	mBottom2 = NULL;
	mTop3 = NULL;
	mBottom3 = NULL;
	mOutTop1 = NULL;
	mOutBottom1 = NULL;
	mOutTop2 = NULL;
	mOutBottom2 = NULL;
	mOutTop3 = NULL;
	mOutBottom3 = NULL;
	
	tpIdxCurrent = -1;
					// Transition & probabilities
	Shift1 = NULL;
	Shift2 = NULL;
	Shift3 = NULL;
	pu = NULL;
	p0 = NULL;
	pd = NULL;
	qu = NULL;
	q0 = NULL;
	qd = NULL;
	ru = NULL;
	r0 = NULL;
	rd = NULL;

	NewPrice = NULL;
	NbTPMax = 0;

	mPPY = 0;
	mEoI = 'E';

	// Smoothing off by default
	mSmoothFact = 0e0;

    // spot vol interpolation flag
    spotVolInterp = FLOORIDX;

    outStdDevs = 1;
}


//--------------------------------------------------------------
//

KMrNTree::~KMrNTree()
{
static	char	routine[] = "KMrNTree::~KMrNTree";

    // Free probailities and transitions
    // Skip if mNbFactor == 0
    if (mNbFactor > 0)
    {
	switch (mNbFactor) {
	case 3:
		sliceDelete(ru);
		sliceDelete(r0);
		sliceDelete(rd);
		delete [] Shift3;
	case 2:
		sliceDelete(qu);
		sliceDelete(q0);
		sliceDelete(qd);
		delete [] Shift2;
	case 1:
		sliceDelete(pu);
		sliceDelete(p0);
		sliceDelete(pd);
		delete [] Shift1;
		break;
	default:
		throw KFailure("%s: not implemented for mNbFactor=%d.\n",
			routine, mNbFactor);
	}

	sliceDelete(NewPrice);


	// Free limits
	//
	MrNTreeLimitsNFactFree(
		mTop1,
		mBottom1,
		mTop2,
		mBottom2,
		mTop3,
		mBottom3,
		mNbFactor,
		NbTP);

	MrNTreeLimitsNFactFree(
		mOutTop1,
		mOutBottom1,
		mOutTop2,
		mOutBottom2,
		mOutTop3,
		mOutBottom3,
		mNbFactor,
		NbTP);

	DrlIntVectFree(mWidth,     0, mNbFactor-1);
	DrlIntVectFree(mHalfWidth, 0, mNbFactor-1);

	//
	// Use "malloc" rather than "new" to be consistent with 
	// Insert(TDate) function, resolving potential inconsistency
	// between memory allocation and free if failed in the middle.  
	//delete [] TPDates;
	//delete [] TPTimes;
	DrlTDateVectFree(TPDates, 0, NbTP+1);
	DrlDoubleVectFree(TPTimes,  0, NbTP+1);
	DrlDoubleVectFree(Length,  -1, NbTP+1);
	DrlDoubleVectFree(LengthJ, -1, NbTP+1);

	// Free vol info
	//
	DrlDoubleVectFree(mAlpha, 0, mNbFactor-1);
	DrlDoubleVectFree(mBeta,  0, mNbFactor-1);
	DrlDoubleMatrFree(mRho,   0, mNbFactor-1, 0, mNbFactor-1);
	DrlDoubleMatrFree(mSigma,    0, NbTP, 0, mNbFactor-1);
	DrlDoubleMatrFree(mAweight, -1, NbTP, 0, mNbFactor*(mNbFactor+1)/2);

    } // end if mNbFactor > 0
}



//--------------------------------------------------------------
//
void
KMrNTree::DeleteMemory()
{
static	char	routine[] = "KMrNTree::DeleteMemory";

	// Free probailities and transitions
	//
	switch (mNbFactor) {
	case 3:
		sliceDelete(ru);
		sliceDelete(r0);
		sliceDelete(rd);
		delete [] Shift3;
	case 2:
		sliceDelete(qu);
		sliceDelete(q0);
		sliceDelete(qd);
		delete [] Shift2;
	case 1:
		sliceDelete(pu);
		sliceDelete(p0);
		sliceDelete(pd);
		delete [] Shift1;
		break;
	case 0:
		// Nothing to do
		return;
	default:
		throw KFailure("%s: not implemented for mNbFactor=%d.\n",
			routine, mNbFactor);
	}

	sliceDelete(NewPrice);


	// Free limits
	//
	MrNTreeLimitsNFactFree(
		mTop1,
		mBottom1,
		mTop2,
		mBottom2,
		mTop3,
		mBottom3,
		mNbFactor,
		NbTP);

	MrNTreeLimitsNFactFree(
		mOutTop1,
		mOutBottom1,
		mOutTop2,
		mOutBottom2,
		mOutTop3,
		mOutBottom3,
		mNbFactor,
		NbTP);


	DrlIntVectFree(mWidth,     0, mNbFactor-1);
	DrlIntVectFree(mHalfWidth, 0, mNbFactor-1);


	//
	// Use "malloc" rather than "new" to be consistent with 
	// Insert(TDate) function, resolving potential inconsistency
	// between memory allocation and free if failed in the middle.  
	//delete [] TPDates;
	//delete [] TPTimes;
	DrlTDateVectFree(TPDates, 0, NbTP+1);
	DrlDoubleVectFree(TPTimes,  0, NbTP+1);

	DrlDoubleVectFree(Length,  -1, NbTP+1);
	DrlDoubleVectFree(LengthJ, -1, NbTP+1);

	// Free vol info
	//
	DrlDoubleVectFree(mAlpha, 0, mNbFactor-1);
	DrlDoubleVectFree(mBeta,  0, mNbFactor-1);
	DrlDoubleMatrFree(mRho,   0, mNbFactor-1, 0, mNbFactor-1);
	DrlDoubleMatrFree(mSigma,    0, NbTP, 0, mNbFactor-1);
	DrlDoubleMatrFree(mAweight, -1, NbTP, 0, mNbFactor*(mNbFactor+1)/2);




	// Re-initialize the tree parameters.
	//
					// Timeline
	mTodayDate = -1L;
	NbTP = 0;
	TPDates = NULL;
	TPTimes = NULL;
	Length  = NULL;
	LengthJ = NULL;

					// Volatility info
	mNbFactor = 0;
	mAlpha = NULL;
	mBeta = NULL;
	mRho = NULL;
	mSigma = NULL;
	mAweight = NULL;

					// Tree geometry 
	mNbSigmaMax = -1;
	mWidth = NULL;
	mHalfWidth = NULL;

	mTop1 = NULL;
	mBottom1 = NULL;
	mTop2 = NULL;
	mBottom2 = NULL;
	mTop3 = NULL;
	mBottom3 = NULL;
	mOutTop1 = NULL;
	mOutBottom1 = NULL;
	mOutTop2 = NULL;
	mOutBottom2 = NULL;
	mOutTop3 = NULL;
	mOutBottom3 = NULL;
	
	tpIdxCurrent = -1;
					// Transition & probabilities
	Shift1 = NULL;
	Shift2 = NULL;
	Shift3 = NULL;
	pu = NULL;
	p0 = NULL;
	pd = NULL;
	qu = NULL;
	q0 = NULL;
	qd = NULL;
	ru = NULL;
	r0 = NULL;
	rd = NULL;

	NewPrice = NULL;
	NbTPMax = 0;

	mPPY = 0;
	mEoI = 'E';

}





//--------------------------------------------------------------
// Free tree memories that depends on the factor vols, including
// tree limits, transitional probabilities. etc.
//
void
KMrNTree::ClearTreeVolMem()
{
static	char	routine[] = "KMrNTree::ClearTreeVolMem";

	// Free probailities and transitions
	//
	switch (mNbFactor) {
	case 3:
		sliceDelete(ru);
		sliceDelete(r0);
		sliceDelete(rd);
		delete [] Shift3;
	case 2:
		sliceDelete(qu);
		sliceDelete(q0);
		sliceDelete(qd);
		delete [] Shift2;
	case 1:
		sliceDelete(pu);
		sliceDelete(p0);
		sliceDelete(pd);
		delete [] Shift1;
		break;
	case 0:
		// Nothing to do
		return;
	default:
		throw KFailure("%s: not implemented for mNbFactor=%d.\n",
			routine, mNbFactor);
	}

	sliceDelete(NewPrice);


	// Free limits
	//
	MrNTreeLimitsFree(
		mTop1,
		mBottom1,
		mTop2,
		mBottom2,
		mTop3,
		mBottom3,
		mNbFactor,
		NbTP);

	MrNTreeLimitsFree(
		mOutTop1,
		mOutBottom1,
		mOutTop2,
		mOutBottom2,
		mOutTop3,
		mOutBottom3,
		mNbFactor,
		NbTP);

	Shift1 = NULL;
	Shift2 = NULL;
	Shift3 = NULL;
	pu = NULL;
	p0 = NULL;
	pd = NULL;
	qu = NULL;
	q0 = NULL;
	qd = NULL;
	ru = NULL;
	r0 = NULL;
	rd = NULL;

	// Required!!!  To properly initialize memory for probability.
	NewPrice = NULL;

}




//--------------------------------------------------------------
// Returns the year fraction to the current timepoint.
//
double
KMrNTree::TPTimeCurrent()
{
static  char    routine[] = "KMrNTree::TPTimeCurrent";

	double	yearFract;

	IF_FAILED_THROW(GtoDayCountFraction(TPToday(),
					    TPDateCurrent(),
					    GTO_ACT_365F,
					    &yearFract));
	return yearFract;

}



//--------------------------------------------------------------
// Return the time point index given a date.
//
int
KMrNTree::TPIdx(TDate date)
{
static  char    routine[] = "KMrNTree::TPIdx";

	int	idxT;
	
	for (idxT=0; idxT<=NbTP+1; idxT++) {
		if (TPDates[idxT] == date)
			break;
	}

	if (idxT == NbTP+1 && TPDates[idxT] != date )
		throw KFailure("%s: input date (%s) is not a critical date.\n",
				routine, GtoFormatDate(date));
	
	return idxT;
	
}



//--------------------------------------------------------------
// Adds a date to the list of critical dates.
// Sort the dates in InitializeTimeline after 
// all critical dates are inserted. 

void
KMrNTree::Insert(TDate date)
{
static	char	routine[] = "KMrNTree::Insert(Date)";
	// Check valid date
	if ((date <= 1L) || (date >= 500000L)) 
		throw KFailure("%s: invalid date %ld.\n", routine, (long)date);

	// Only add date if > mTodayDate
	if (date >= mTodayDate)
	{
		// Add the date to array of critical dates
		// Do not add it if it is already there.
		IF_FAILED_THROW( DrlVTypeVectAdd(
					(void**) &TPDates,
					&NbTP,
					&NbTPMax,
					512,
					(void*) &date,
					TRUE, 	// TRUE=remove double items
					DRL_TDATE_T));
	}

	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s: inserting %10s.\n",
			routine, DrlTDatePrint(NULL, date));
	}
}



//--------------------------------------------------------------
// Dummy: these are not defined for an abstract tree,
// but they are here so you don'y need to implement them for testing.
//
 
void KMrNTree::Insert(const KZeroReset &zeroReset, bool isCrit)
{
	throw KFailure("KMrNTree::Insert(KZeroReset&,bool): N/A.\n");
}
 
TDate KMrNTree::Insert(const KRateReset &rtReset, bool isCrit)
{
	{throw KFailure("KMrNTree::Insert(KRateReset&,bool): N/A.\n");}
	return((TDate)-1L);
}
 

TDate KMrNTree::Insert(const KRateReset &rtReset, TDate endDate, bool isCrit)
{
	{throw KFailure("KMrNTree::Insert(KRateReset&,TDate, bool): N/A.\n");}
	return((TDate)-1L);
}
 
 
void KMrNTree::Get(KTSlice& ts, const KZeroReset&)
{
	throw KFailure("KMrNTree::Get(KTSlice&,KZeroReset&): N/A.\n");
}
 
void KMrNTree::Get(KTSlice& ts, const KRateReset&)
{
	throw KFailure("KMrNTree::Get(KTSlice&,KRateReset&): N/A.\n");
}



// Return tree dates in specified range [stDate, endDate]
//
KVector(TDate)
KMrNTree::GetTDatesRange(TDate stDate, TDate endDate)
{
static	char	routine[] = "KMrNTree::GetTDatesRange";

	int	i;
	KVector(TDate)	dateList;

try {
	// Check
	ASSERT_OR_THROW(endDate >= stDate);

	if (NbTP != 0)
	{
	    // Loop through each date in the list
	    //
	    for (i=0; i<=NbTP-1; i++)	
	    {
		if (TPDates[i] <= endDate &&
		    TPDates[i] >= stDate)
			dateList.push_back(TPDates[i]);		
	    }

	    // Sort and merge dateList in ascending order
	    //
	    sort(dateList.begin(), dateList.end());
	    KVector(TDate)::iterator iterBadDateSt
		= unique(dateList.begin(), dateList.end());
	    dateList.erase(iterBadDateSt, dateList.end());
	}

	return dateList;

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
// Sort the critical dates in ascending order and remove duplicate
// dates.  Set up time line according to ppy 

void
KMrNTree::InitializeTimeline(
	char EoI,			// (I) Equal or increasing time steps
	int ppy)			// (I) period per year
{
static	char	routine[] = "KMrNTree::InitializeTimeline";
	TDate	*pDates = NULL;		// temporary dates
	double	*pTimes = NULL;		// temporary times (offset)
	int	nPDates,		// num dates
		nPDatesMax;

	int	idxT, i, n, k;
	long	nbDaysSoFar;
	TDate	lastDate;
	double	dt, dtMax;
	int	dayStep, dayStepMax;

    try {


	// Sort and merge the critical dates in ascending order
	//
	IF_FAILED_THROW(DrlVTypeVectSort(TPDates, &NbTP, DRL_TDATE_T, TRUE));

	if (ppy <= 365) {
		// ppy allows integer dates timepoints
		// Maximum step size from ppy
		dtMax = 1e0 / (double) ppy;
		dayStepMax = (int) (365e0 / (double)ppy);
		if (dayStepMax < 1) {
			throw KFailure("%s: maximum day step %d < 1.\n",
				routine, dayStepMax);
		}


		// Allocate enough space (at most one node per day)
		//
		nPDatesMax = (int) (TPDates[NbTP-1] - TPDates[0]);
		ASSERT_OR_THROW((pDates = new TDate  [nPDatesMax+2]) != NULL);
		ASSERT_OR_THROW((pTimes = new double [nPDatesMax+2]) != NULL);
		nPDates = 0;

		if (EoI == 'I') {
		    // Set up increasing time step (from Lionnel)
		    idxT = 0;
		    dayStep = 1;
		    for (i=0; i<NbTP-1; i++) {
			pDates[idxT] = TPDates[i];
			pTimes[idxT] = YrsDiff(pDates[0], pDates[idxT]);
			idxT++;

			// If the next critical date is already
			// in the proper increasing interval
			// (TPDates[i+1]-TPDates[i] = dayStep+1),
			// then use it as the next time point,
			// rathan than insert an extra date 
			// in between.  This would happen in CET
			// when inserting product dates < 2y.
			//
		    	if ((TPDates[i+1]-TPDates[i]) == dayStep+1 &&
				dayStep < dayStepMax && i > 0)
				dayStep++;
			else {
			    while (pDates[idxT-1] + dayStep < TPDates[i+1])
		  	    {
			    	//
			    	// The last step between TPDates[i] and 
				// TPDates[i+1] is the residual left.  
				// To avoid it being too small, we divide the 
				// last TWO steps equally, with both 
				// steps <= dayStep and the (n-1)th step is 
				// either equal to or 1 day shorter than the 
				// nth step.
			    	//
			    	if(pDates[idxT-1] + 2*dayStep < TPDates[i+1])
				{
				    pDates[idxT] = pDates[idxT-1] + dayStep;
			    	    if (dayStep < dayStepMax) dayStep++;
				}
			    	else
				{
				    pDates[idxT] = (TPDates[i+1]
						  + pDates[idxT-1]) / 2;
				}

			    	pTimes[idxT] = YrsDiff(pDates[0], pDates[idxT]);
			    	idxT++;
			    } // end while
			} // end if

		    } // end for

		    pDates[idxT] = TPDates[i];	// don't forget last point
		    pTimes[idxT] = YrsDiff(pDates[0], pDates[idxT]);
		    idxT++;
		    nPDates = idxT;
		} else {
		    // Set up equal time step 
		    idxT = 0;
		    for (i=0; i<NbTP-1; i++) {
			pDates[idxT] = TPDates[i];
			pTimes[idxT] = YrsDiff(pDates[0], pDates[idxT]);
			idxT++;

			// Compute # of days between steps
			lastDate = pDates[idxT-1];
			
			n = (int) floor(YrsDiff(TPDates[i], TPDates[i+1]) / dtMax);
			dayStep = (int) (TPDates[i+1] - TPDates[i]) / (n+1);
			for(k=1; k<=n; k++) {
				//pDates[idxT] = pDates[idxT-1] + dayStep;
				pDates[idxT] = lastDate + 
					(long) (((double)k)
					* ((double)(TPDates[i+1] - TPDates[i]))
					/ ((double)(n+1)));
				pTimes[idxT] = YrsDiff(pDates[0], pDates[idxT]);
				idxT++;
			}
		    }

		    pTimes[idxT] = YrsDiff(pDates[0], TPDates[i]);
		    pDates[idxT] = TPDates[i];	// don't forget last point
		    idxT++;

		    nPDates = idxT;

		}
	} else {

		// ppy > 365, fraction dates, is not done properly here.
		// The problem has to do setting pDates = -1L in between.
		// These invalid dates are later used in the program
		// such as zerobank, etc. 
		throw KFailure("%s: ppy (%d) > 365 is not supported.\n", 
				routine,
				ppy);

		// ppy is larger than 365, allow fractional time points
		// Equal time steps form the moment
		dtMax = 1e0 / (double) ppy;

		// Allocate enough space (at most one node per day)
		//
		nPDatesMax = (int) ((TPDates[NbTP-1] - TPDates[0])/365e0/dtMax);
		ASSERT_OR_THROW((pDates=new TDate  [nPDatesMax+NbTP]) != NULL);
		ASSERT_OR_THROW((pTimes=new double [nPDatesMax+NbTP]) != NULL);
		nPDates = 0;
	
		// Set up equal time step 
		idxT = 0;
		for (i=0; i<NbTP-1; i++) {
			pDates[idxT] = TPDates[i];
			pTimes[idxT] = YrsDiff(pDates[0], pDates[idxT]);
			idxT++;
			n = (int) floor(YrsDiff(TPDates[i], TPDates[i+1]) / dtMax);
			dt = YrsDiff(TPDates[i], TPDates[i+1]) / (double) (n+1);
			for(k=1; k<=n; k++) {
				pTimes[idxT] = pTimes[idxT-1] + dt;
				pDates[idxT] = -1L; 	// not a date
				idxT++;
			}
		}

		pTimes[idxT] = YrsDiff(pDates[0], TPDates[i]);
		pDates[idxT] = TPDates[i];	// don't forget last point
		idxT++;

		nPDates = idxT;
	}

	// free initial critical dates array
	// so that we can replace it.
	DrlVTypeVectFree((void*) TPDates, NbTP, DRL_TDATE_T);
	TPDates = NULL;

	// Debug
	if (debugLevel >= DEBUG_LEVEL_TIMELINE) {
	    dppLog << format("%s: Preprocessed time nodes.\n", routine);
	    for (idxT=0; idxT<nPDates; idxT++) 
		dppLog << format(" %3d %10s %12.8f\n",
				idxT,
				printDate(pDates[idxT]),
				pTimes[idxT]);
	}






	// We now know the timepoints, so we can allocate the
	// memory in the tree structure and copy over.

	NbTP = nPDates-1;
	NbTPMax = MAX(NbTP+2, NbTPMax);

	// Use "malloc" rather than "new" to be consistent with 
	// Insert(TDate) function, resolving potential inconsistency
	// between memory allocation and free if failed in the middle.  
	ASSERT_OR_THROW((TPDates = DrlTDateVectAlloc(0, NbTP+1)) != NULL);
	ASSERT_OR_THROW((TPTimes = DrlDoubleVectAlloc(0, NbTP+1)) != NULL);

	ASSERT_OR_THROW((Length  = DrlDoubleVectAlloc(-1, NbTP+1)) != NULL);
	ASSERT_OR_THROW((LengthJ = DrlDoubleVectAlloc(-1, NbTP+1)) != NULL);
	for (idxT=0; idxT<=NbTP; idxT++) {
		TPDates[idxT] = pDates[idxT];
		TPTimes[idxT] = pTimes[idxT];
	}

	// We know compute the true time step size (Length) and
	// time step size used for jump computation (LengthJ)

	nbDaysSoFar = 0;
	for (idxT = 0; idxT < NbTP; idxT++)
	{
		// Compute actual time step ACT/365F
		if (idxT < NbTP-1)
			Length[idxT] = (double) (TPTimes[idxT+1] - TPTimes[idxT]);
		else
			Length[idxT] = Length[idxT-1];


		// If  Length is OK, then assign delta t for jumpsize. If on 
		// 'I' it is the step size corresponding to the current time
		// and if 'E', then it is a constant. 
		switch (EoI) {
		case 'I':
			LengthJ[idxT] = Min(0.5*(1.0+sqrt(1.0+8.0*nbDaysSoFar)) / 365e0,
						dtMax);
			break;
		case 'E':
			LengthJ[idxT] = 1e0 / ppy;
			break;
		}

        	nbDaysSoFar += (long)((TPTimes[idxT+1] - TPTimes[idxT])*365e0);
	}

	// We need 2 extra time points with arbitrary length for the calibration 
	Length[NbTP]  = Length[NbTP - 1];
	Length[-1]    = Length[0];

	LengthJ[NbTP] = LengthJ[NbTP - 1];
	LengthJ[-1]   = LengthJ[0];

	TPDates[NbTP+1] = TPDates[NbTP] + (long) (Length[NbTP] * 365e0 + 0.5e0);



	/*// Print debug info
	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s: TIMELINE\n", routine);
		dppLog << format("TIDX DATE       TIME         "
			"LENGTH       LENGTH-J    \n");
		for (idxT=0; idxT<=NbTP; idxT++) {
			dppLog << format(" %3d %10s %12.8f %12.8f %12.8f\n",
				idxT,
				DrlTDatePrint(NULL, TPDates[idxT]),
				TPTimes[idxT],
				Length[idxT],
				LengthJ[idxT]);
		}
	}*/


	// free memory
	delete [] pDates;
	delete [] pTimes;
    }
    catch (KFailure) {
	// free memory
	delete [] pDates;
	delete [] pTimes;
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
//

void
KMrNTree::CheckTimeline(
	int ppy)			// (I) period per year
{
static	char	routine[] = "KMrNTree::CheckTimeline";
	int	idxT;
	double	dt,
		dtMax = 1e0 / (double) ppy;

	int	dateStep;
	int	dayStepMax = (int) (365e0 / (double)ppy);

    try {

	// Print debug info
	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s: TIMELINE\n", routine);
		dppLog << format("TIDX DATE          TIME         "
			"LENGTH       LENGTH-J    \n");
		for (idxT=0; idxT<=NbTP; idxT++) {
			dppLog << format(" %3d %10s %12.8f %12.8f %12.8f\n",
				idxT,
				printDate(TPDates[idxT]),
				TPTimes[idxT],
				Length[idxT],
				LengthJ[idxT]);
		}
	}


	// Check no critical dates before today
	if (mTodayDate <= 0L)
		throw KFailure("%s: today date not set.\n", routine);

	if (mTodayDate > TPDates[0])
		throw KFailure("%s: today date (%10s) > first critical date (%10s).\n",
			routine,
			DrlTDatePrint(NULL, mTodayDate),
			DrlTDatePrint(NULL, TPDates[0]));

	// Test timeline is correct (increasing with steps
	// no larger than those specified by ppy.
	for (idxT=0; idxT<=NbTP; idxT++) {
		// If CoB date, check consistent offset time
		if (TPDates[idxT] > 0L) {
			dt = YrsDiff(TPDates[0], TPDates[idxT]);
			if (!IS_ALMOST_ZERO(dt - TPTimes[idxT])) {
				throw KFailure("%s: inconsistent offset %12.8f "
					" at date %d (%10s) from date 0 (%10s).\n",
					routine, dt, idxT,
					DrlTDatePrint(NULL, TPDates[idxT]),
					DrlTDatePrint(NULL, TPDates[0]));
			}
		}

		if (idxT < NbTP) {
                	dt       = TPTimes[idxT+1] - TPTimes[idxT];
                	dateStep = TPDates[idxT+1] - TPDates[idxT];

			// There might be cases (high ppy~70) where 
			// dateStepMax = 5, while two dates calculated are 
			// occasionally off by 6 days, and 6/5=1.2 > 1.1.
			// So we added one more condition for check
			// dateStep <= dayStepMax+1.
                	if (dt > dtMax*1.1e0 &&
			    dateStep > dayStepMax+1) {
                        	throw KFailure("%s: too large offset between "
					"date %d (%12.8f, %10s) >= "
					"date %d (%12.8f, %10s).\n",
                                	routine,
					idxT, TPTimes[idxT],
					printDate(TPDates[idxT]),
					idxT+1, TPTimes[idxT+1],
					printDate(TPDates[idxT+1]));
			}
                }
        }
    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}





//--------------------------------------------------------------
// MrNTreeLimitsNFact:
// Calculates tree limits: we calculate the volatility of each process
// and cut the tree according to a parabola in each dimension.
// Also allocate on the go memory for the tree limits.

static	void
MrNTreeLimitsNFact(
	int *Width,			// (O) Ellipsoid widt
	int *HalfWidth,			// (O) Ellipsoid half width 

	int *Top1,			// (O) Upper limits of 1D tree
	int *Bottom1,			// (O) Lower limits           
	int **Top2,			// (O) Upper limits of 2D tree
	int **Bottom2,			// (O) Lower limits           
	int ***Top3,			// (O) Upper limits of 3D tree
	int ***Bottom3,			// (O) Lower limits           

	int nbFactor,			// (I) Number of factors      
	int Ns,				// (I) Max nb of std devs     
	int NbTP,			// (I) Nb of time points      
	double *Beta,			// (I) Mean reversions [nbFactor]
	double **Aweight,		// (I) Orthogonal weights     
	double *Length,			// (I) Length of time steps   
	double *LengthJ)		// (I) Time step for jump size 
{
static	char	routine[] = "MrNTreeLimitsNFact";

	int	nFact = nbFactor;	// Used in some macros.
	double	vol[NFMAX],		// cumulative vol and corrs.
		rho[NFMAX][NFMAX],
		cov[NFMAX][NFMAX];
	double	x[NFMAX],		// reference stochastic variable
		z[NFMAX],		// normalized stochastic variable
		jump[NFMAX][NFMAX],	// Jump coeff [0<i<=j<nbFactor]
		a[NFMAX][NFMAX],	// ellipsoid axis [i<j]
		c[NFMAX][NFMAX],	// c coeffs [i<j]
		d[NFMAX][NFMAX],	// d coeffs [i<j]
		e[NFMAX][NFMAX];	// e coeffs [i<j]

	int	i1, i2, i3;

	int	t, i, j, k;		// Time and node indices
	int	h;			/* Local value of half ellipse width */
	int	h2, h3;			/* Minimum width of ellipsoid        */
	double  du;			/* Jump coefficient                  */
	int     Center;			/* Center of the ellipse             */
	double  Delta;			/* Intermediate results              */


    try {
	if (debugLevel > DEBUG_LEVEL_DRIFT) {
		dppLog << format("%s:\n", routine);
	}

	// Inialize cumulative vols and covaraiances.
	//
	for (i1=0; i1<nbFactor; i1++)
	for (i2=0; i2<=i1     ; i2++) {
		cov[i1][i2] = cov[i2][i1] = 0e0;
	}


	// Treat the case of the first node.
	// Start with 5 nodes in each dimension.
	for (i1=0; i1<nbFactor; i1++) {
		Width[i1]     = 5;
		HalfWidth[i1] = 2;
	}


	Bottom1[0] = -2;
	Top1[0]    =  2;

	if (nbFactor >= 2)
	{  
	    Bottom2[0] = DrlIntVectAlloc(-2, 2);
	    Top2[0]    = DrlIntVectAlloc(-2, 2);

	    ASSERT_OR_THROW(Bottom2[0] != NULL);
	    ASSERT_OR_THROW(Top2[0]    != NULL);

	    if (nbFactor >= 3)
	    {    
		Bottom3[0] = DrlIntPVectAlloc(-2, 2);
		Top3[0]    = DrlIntPVectAlloc(-2, 2);

		ASSERT_OR_THROW(Bottom3[0] != NULL);
		ASSERT_OR_THROW(Top3[0]    != NULL);

	    }

	    for (i = -2; i <= 2; i++)
	    {
		Bottom2[0][i] = -2;
		Top2[0][i]    =  2;
                
		if (nbFactor >= 3)
		{    
		    Bottom3[0][i] = DrlIntVectAlloc(-2, 2);
		    Top3[0][i]    = DrlIntVectAlloc(-2, 2);

		    ASSERT_OR_THROW(Bottom3[0][i] != NULL);
		    ASSERT_OR_THROW(Top3[0][i]    != NULL);

		    for (j = -2; j <= 2; j++)
		    {
			Bottom3[0][i][j] = -2;
			Top3[0][i][j]    =  2;
		    }
		}
	    }
	}


	// Forward in time
	//
	for (t = 0; t <= NbTP; t++)                  
	{   

	  if (t > 0)
	  {

	    // Compute all coefficients needed
	    // for computations below.

	    //---  Jump size
	    du = sqrt (JUMPCOEFF * LengthJ[t-1]);
	    for (i1=0; i1<nbFactor; i1++)
	    for (i2=0; i2<=i1     ; i2++) {
		jump[i1][i2] = Aweight[t-1][AIDX(i1, i2)] * du;
	    }

	    //--- c are the (Z1,..,Zk)  to Xk coefficients
	    c[0][0] = vol[0];
	    if (nbFactor >=2)
	    {
		c[1][0] = vol[1] * rho[0][1];
		c[1][1] = vol[1] * sqrt(1e0 - rho[0][1]*rho[0][1]);

		if (nbFactor >= 3)
		{
		    c[2][0] = vol[2] * rho[0][1];
		    c[2][1] = vol[2] * (rho[1][2] - rho[0][1]*rho[0][2])
			/ sqrt(1e0 - rho[0][1]*rho[0][1]);
		    c[2][2] = vol[2] * (1e0 - 
			- rho[0][1]*rho[0][1]
			- rho[0][2]*rho[0][2]
			- rho[1][2]*rho[1][2]
			+ 2e0 * rho[0][1]*rho[0][2] * rho[1][2])
			/ sqrt(1e0 - rho[0][1]*rho[0][1]);
		}
	    }

	    //---  e are the (X1,..,Xk)  to Zk coefficients (e=c^-1)
	    //---  We need them only up to n-1
	    if (nbFactor >=2)
	    {
	        e[0][0] = 1e0 / vol[0];
		if (nbFactor >= 3)
		{
	            e[1][0] = - rho[0][1] 
		        / (vol[0] * (1e0 - rho[0][1]*rho[0][1]));
	            e[1][1] = 1e0
		        / (vol[1] * (1e0 - rho[0][1]*rho[0][1]));
		}
	    }


	    //---  d are the (X1,..,Xk-1,Zk) to Xk coefficients
	    if (nbFactor >=2)
	    {
	        d[1][0] = vol[1] * rho[0][1]
		        / vol[0];

		if (nbFactor >= 3)
		{
	    	    d[2][0] = vol[2] * (rho[0][2] - rho[0][1]*rho[1][2])
		        / (vol[0] * (1e0 - rho[0][1]*rho[0][1]));
	    	    d[2][1] = vol[2] * (rho[1][2] - rho[0][1]*rho[0][2])
		        / (vol[1] * (1e0 - rho[0][1]*rho[0][1]));
		}
	    }

	    //---  a are the Z to X coefficients
	    if (nbFactor >=2)
	    {
	        a[1][0] = (d[1][0]*jump[0][0] - jump[1][0])
			/ fabs(jump[1][1]);

		if (nbFactor >= 3)
		{
	    	    a[2][0] = (d[2][0]*jump[0][0] + d[2][1]*jump[1][0] - jump[2][0])
			/ fabs(jump[2][2]);
	    	    a[2][1] = (d[2][1]*jump[1][0]                      - jump[2][1])
			/ fabs(jump[2][2]);
		}
	    }


	    // Limits of 1D Tree
	    //

	    // fabs () is needed as Jumps can be negative
	    // Also make sure there is at least 2 nodes

	    h = Max(Ceil(Ns * vol[0] / fabs(jump[0][0])), 2);

	    // No shift of axis in first dimension 
	    Center = 0;

	    Bottom1[t] = Center - h;
	    Top1[t]    = Center + h;
            
	    // Record the maximum width over time 
	    HalfWidth[0] = MAX (HalfWidth[0], h);
	    Width[0]     = 2 * HalfWidth[0] + 1;


	    // Limits of 2D tree. They depend on first dimension.
	    //

	    if (nbFactor >= 2) 
	    {  
		// Allocate memory for the 2D limits
		// and the 3D matrices of limits.

		Bottom2[t] = DrlIntVectAlloc(Bottom1[t], Top1[t]);
		Top2[t]    = DrlIntVectAlloc(Bottom1[t], Top1[t]);

		ASSERT_OR_THROW(Bottom2[t] != NULL);
		ASSERT_OR_THROW(Top2[t]    != NULL);

		if (nbFactor >= 3)
		{   
		    Bottom3[t] = DrlIntPVectAlloc(Bottom1[t], Top1[t]);
		    Top3[t]    = DrlIntPVectAlloc(Bottom1[t], Top1[t]);

		    ASSERT_OR_THROW(Bottom3[t] != NULL);
		    ASSERT_OR_THROW(Top3[t]    != NULL);
		}


		// Minimum half size of the ellipse:
		// slope + 1 node. This is to avoid case(1)
		// where there is no square of 9 nodes to branch to.
		// *   (case(2) is correct):
		// *                                       x
		// *                                   x   x
		// *               x               O   O   O
		// *           x   x               O   O   O
		// *       x   x   x   =======>    O   O   O
		// *       x   x                   x   x
		// *       x                       x
		// *
		// *           (1)                     (2)
		//

                h2 = Ceil(fabs(a[1][0])) + 2;

		for (i = Bottom1[t]; i <= Top1[t]; i++)
		{

		    x[0] = i*jump[0][0];
		    z[0] = e[0][0]*x[0];

		    Delta = c[1][1] / fabs(jump[1][1]) *
					SqrtM(Ns*Ns - z[0]*z[0]);

		    h = Max(Ceil(Delta), h2);

		    Center = NearInt(i * a[1][0]);

		    Bottom2[t][i] = Center - h;
		    Top2[t][i]    = Center + h;
                
		    HalfWidth[1] = MAX (HalfWidth[1], h);
		    Width[1]     = 2 * HalfWidth[1] + 1;


		    // Limits of 3D tree. They depend on first two dimensions.
		    //
		    if (nbFactor >= 3)
		    {    
			// Allocate memory for limits
			Bottom3[t][i] = DrlIntVectAlloc(
						Bottom2[t][i], Top2[t][i]);
			Top3[t][i]    = DrlIntVectAlloc(
						Bottom2[t][i], Top2[t][i]);
			ASSERT_OR_THROW(Bottom3[t][i] != NULL);
			ASSERT_OR_THROW(Top3[t][i]    != NULL);

			//  Minimum half size of the ellipsoid
			h3 = Ceil(fabs(a[2][0]) + fabs(a[2][1])) + 2;


			for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
			{
			    x[1] = i*jump[1][0] + j*jump[1][1];
			    z[1] = e[1][0]*x[0] + e[1][1]*x[1];

			    Delta = c[2][2] / fabs(jump[2][2]) *
					SqrtM(Ns*Ns - z[0]*z[0] - z[1]*z[1]);

			    h = Max(Ceil(Delta), h3);

			    Center = NearInt(a[2][0] * i + a[2][1] * j);

			    Bottom3[t][i][j] = Center - h;
			    Top3[t][i][j]    = Center + h;
                
			    HalfWidth[2] = Max(HalfWidth[2], h);
			    Width[2]     = 2 * HalfWidth[2] + 1;
			}
		    }
		} 
	    }

	  } // if (t > 0)


	    //  Compute correlation matrix at next time step.
	    //
	    for (i1=0; i1<nbFactor; i1++) 
	    for (i2=0; i2<=i1     ; i2++) {
	        cov[i1][i2] *= exp(-(Beta[i1]+Beta[i2])*Length[t]);

		i3 = Min(i1, i2);
	        for (k=0; k<=i3; k++) {
		    cov[i1][i2] += 
			Aweight[t][AIDX(i1,k)] * Aweight[t][AIDX(i2,k)]
			* Length[t] * ExpDecay(Beta[i1] + Beta[i2], Length[t]);
	        }

	        cov[i2][i1] = cov[i1][i2];
	    }
	    for (i1=0; i1<nbFactor; i1++)  {
		vol[i1] = sqrt(cov[i1][i1]);
		for (i2=0; i2< i1     ; i2++) {
			rho[i1][i2] = rho[i2][i1] = cov[i1][i2] / (vol[i1] * vol[i2]);
		}
	    }

	    // Debug info
	    if (debugLevel > DEBUG_LEVEL_DRIFT) {
		dppLog << format(" %3d/%3d |", t, NbTP);
		for (i1=0;  i1<nbFactor; i1++) 
		for (i2=i1; i2<nbFactor; i2++) 
			dppLog << format(" c[%d%d]=%12.6f",
						i1, i2, cov[i1][i2]);
		dppLog << format(" |");
		for (i1=0;  i1<nbFactor; i1++) 
			dppLog << format(" v[%d]=%12.6f", i1, vol[i1]);
		dppLog << format(" |");
		for (i1=0;    i1<nbFactor; i1++) 
		for (i2=i1+1; i2<nbFactor; i2++) 
			dppLog << format(" r[%d%d]=%12.6f",
				i1, i2, rho[i1][i2]);
		dppLog << format("\n");
	    }
	}


    }
    catch (KFailure) {
        MrNTreeLimitsNFactFree(
            Top1,
            Bottom1,
            Top2,
            Bottom2,
            Top3,
            Bottom3,
            nbFactor,
            NbTP);

        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
// Frees the tree limits.
// Separate the memories that depends on factor vols
// from that independent of factor vols, to be used
// by general Initialize() and CET() respectively. 

static	void
MrNTreeLimitsNFactFree(
	int *Top1,			// (I) Upper limits of 1D tree
	int *Bottom1,			// (I) Lower limits           
	int **Top2,			// (I) Upper limits of 2D tree
	int **Bottom2,			// (I) Lower limits           
	int ***Top3,			// (I) Upper limits of 3D tree
	int ***Bottom3,			// (I) Lower limits           

	int nbFactor,			// (I) Number of factors      
	int NbTP)			// (I) Nb of time points      
{
static	char	routine[] = "MrNTreeLimitsNFactFree";

	// Memory depending on factor vols
	//
	MrNTreeLimitsFree(
		Top1,
		Bottom1, 
		Top2,   
		Bottom2,
		Top3,  
		Bottom3,
 
		nbFactor, 
		NbTP);

	// Memeory independent of factor vols.
	DrlIntVectFree(Top1,    -1, NbTP+1);
	DrlIntVectFree(Bottom1, -1, NbTP+1);
	DrlVoidPVectFree((void**) Top2,    -1, NbTP+1);
	DrlVoidPVectFree((void**) Bottom2, -1, NbTP+1);
	DrlVoidPVectFree((void**) Top3,    -1, NbTP+1);
	DrlVoidPVectFree((void**) Bottom3, -1, NbTP+1);

}



//--------------------------------------------------------------
// Frees the tree limits depending on vols, not the global pointers.

static	void
MrNTreeLimitsFree(
	int *Top1,			// (O) Upper limits of 1D tree
	int *Bottom1,			// (O) Lower limits           
	int **Top2,			// (O) Upper limits of 2D tree
	int **Bottom2,			// (O) Lower limits           
	int ***Top3,			// (O) Upper limits of 3D tree
	int ***Bottom3,			// (O) Lower limits           

	int nbFactor,			// (I) Number of factors      
	int NbTP)			// (I) Nb of time points      
{
static	char	routine[] = "MrNTreeLimitsNFactFree";
	int	t, i;


	// Free limits
	//

	if ((nbFactor >= 3) 
            && (Top1 != NULL) && (Bottom1 != NULL)
            && (Top2 != NULL) && (Bottom2 != NULL)
            && (Top3 != NULL) && (Bottom3 != NULL))
	{
	    for (t = 0; t <= NbTP; t++)
	    {    
		if ((Top2[t]    != NULL) 
		    && (Bottom2[t] != NULL)
		    && (Top3[t]    != NULL)
		    && (Bottom3[t] != NULL))
		{
		    for (i = Bottom1[t]; i <= Top1[t]; i++)
		    {
			DrlIntVectFree(Top3[t][i],
				Bottom2[t][i],
				Top2[t][i]);

			DrlIntVectFree( Bottom3[t][i],
				Bottom2[t][i],
				Top2[t][i]);
		    }
	    	}  /* if */

		DrlIntPVectFree(Top3[t],
			Bottom1[t],
			Top1[t]);

		DrlIntPVectFree(Bottom3[t],
			Bottom1[t],
			Top1[t]);
	    }  /* for t */
	}  /* if */


	if ((nbFactor >= 2) 
	    && (Top1 != NULL) && (Bottom1 != NULL)
	    && (Top2 != NULL) && (Bottom2 != NULL))
	{
	    for (t = 0; t <= NbTP; t++)
	    {    
		DrlIntVectFree(Top2[t],
			Bottom1[t],
			Top1[t]);

		DrlIntVectFree(Bottom2[t],
			Bottom1[t],
			Top1[t]);
	    }  /* for t */
	}  /* if */



}


//--------------------------------------------------------------
// Frees the tree limits.

static	void
MrNTreeLimitsNFactPrint(
	int *Width,			// (I) Ellipsoid widt
	int *HalfWidth,			// (I) Ellipsoid half width 

	int *Top1,			// (I) Upper limits of 1D tree
	int *Bottom1,			// (I) Lower limits           
	int **Top2,			// (I) Upper limits of 2D tree
	int **Bottom2,			// (I) Lower limits           
	int ***Top3,			// (I) Upper limits of 3D tree
	int ***Bottom3,			// (I) Lower limits           
	int nbFactor,			// (I) Number of factors      
	int NbTP)			// (I) Nb of time points      
{
static	char	routine[] = "MrNTreeLimitsNFactPrint";
	//int	t, k, i1, i2;
	int	t, k, j0, j1;


	dppLog << format("%s: \n", routine);


	dppLog << format(" HALF-W   WIDTH \n");
	for (k=0; k<nbFactor; k++) {
		dppLog << format("    %4d    %4d\n", HalfWidth[k], Width[k]);
	}


	for (t = 0; t <= NbTP; t++) {
	    dppLog << format("TP: %4d/%4d\n", t, NbTP);
	    dppLog << format(" [%4d,%4d]",
			Bottom1[t], Top1[t]);
	    dppLog << format("\n");

	    if (nbFactor >=2) {
		for (j0 = Bottom1[t]; j0 <= Top1[t]; j0++) {
		    dppLog << format("            ");
		    dppLog << format("%4d [%4d,%4d]", 
				j0, 
				Bottom2[t][j0],
				Top2[t][j0]);
		    dppLog << format("\n");

		    if (nbFactor >= 3) {
			for (j1 = Bottom2[t][j0]; j1 <= Top2[t][j0]; j1++) {
		    	    dppLog << format("            ");
		    	    dppLog << format("                ");
			    dppLog << format("%4d [%4d,%4d]", 
				j1, 
				Bottom3[t][j0][j1],
				Top3[t][j0][j1]);
		    	    dppLog << format("\n");
			}
		    }
		}
	    }
	}
}

//--------------------------------------------------------------
// Checks the tree limits.

void
KMrNTree::LimitsCheck()
{
static	char	routine[] = "KMrNTree::LimitsCheck";
	int	t, j0, j1;

	int	nFact = mNbFactor;		// Used in some macros.
	int	*inTop0,			// Inner limits at t+1
		*inBot0,
		**inTop1,
		**inBot1,
		***inTop2,
		***inBot2;
	int	*ouTop0,			// Outer limits at t+1
		*ouBot0,
		**ouTop1,
		**ouBot1,
		***ouTop2,
		***ouBot2;

	int	ob,  ib,  it,  ot;

    try {

	inTop0 = this->mTop1;
	inTop1 = this->mTop2;
	inTop2 = this->mTop3;
	inBot0 = this->mBottom1;
	inBot1 = this->mBottom2;
	inBot2 = this->mBottom3;

	ouTop0 = this->mOutTop1;
	ouTop1 = this->mOutTop2;
	ouTop2 = this->mOutTop3;
	ouBot0 = this->mOutBottom1;
	ouBot1 = this->mOutBottom2;
	ouBot2 = this->mOutBottom3;


	for (t = 0; t <= NbTP; t++) {

	    // Check In <= Out
	    //
	    if (nFact >= 1) {
		ob = ouBot0[t];
		ib = inBot0[t];
		it = inTop0[t];
		ot = ouTop0[t];
		if ((ob > ib) || (ib >= it) || (it > ot)) {
			throw KFailure("%s: tpIdx=%d DIM1 :"
			"ouBot (%d) <= inBot (%d) < inTop (%d) <= ouTop (%d) failed.\n",
			routine, t, ob, ib, it, ot);
		}
		
	    if (nFact >= 2) {
		for (j0 = inBot0[t]; j0 <= inTop0[t]; j0++) {
		    ob = ouBot1[t][j0];
	    	    ib = inBot1[t][j0];
	    	    it = inTop1[t][j0];
	    	    ot = ouTop1[t][j0];
	    	    if ((ob > ib) || (ib >= it) || (it > ot)) {
		        throw KFailure("%s: tpIdx=%d DIM2 (j0=%d) :"
			"ouBot (%d) <= inBot (%d) < inTop (%d) <= ouTop (%d) failed.\n",
			routine, t, j0, ob, ib, it, ot);
	    
	    	    }

	    if (nFact >= 3) {
		for (j1 = inBot1[t][j0]; j1 <= inTop1[t][j0]; j1++) {
	    	    ob = ouBot2[t][j0][j1];
		    ib = inBot2[t][j0][j1];
		    it = inTop2[t][j0][j1];
		    ot = ouTop2[t][j0][j1];

		    if ((ob > ib) || (ib >= it) || (it > ot)) {
		        throw KFailure("%s: tpIdx=%d DIM3 (j0=%d j1=%d) :"
			"ouBot (%d) <= inBot (%d) < "
				"inTop (%d) <= ouTop (%d) failed.\n",
			 routine, t, j0, j1, ob, ib, it, ot);
	   	     
	    	    }
		}
	    }
	        }
	    }
	    }
	}

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}





//--------------------------------------------------------------
// Initialize the tree parameters. Mainly set up model parameters.
// The computation of timeline, jump size, orthogonalize factors, etc. 
// will be done in Calibrate AFTER all the product elated critical 
// dates are inserted.  

void
KMrNTree::Initialize(
	TDate	todayDt,		// (I) today's date
	int ppy,			// (I) period per year
	double smoothFact,		// (I) smoothing factor
	char EoI,			// (I) Equal or increasing time steps
	int numStdev,			// (I) number stdev to cut tree

	int numFact,			// (I) number of factors
	double *factMr,			// (I) factor mean reversions
	double *factWeight,		// (I) factor weights
	double *factCorr,		// (I) factor correlation [nbFact]

	KVector(TDate) &volDates,	// (I) volatility dates
	KVector(KVector(double)) &factVol)// (I) spot vols 

{
static	char	routine[] = "KMrNTree::Initialize";

	int	i1, i2, numDates;

 try {

	// Initialize the tree parameters
	this->mNbSigmaMax = numStdev;
/*
	if(numFact == 3 && numStdev < NUMSTDEV)
		this->mNbSigmaMax = NUMSTDEV;
	else
		this->mNbSigmaMax = numStdev;
*/

	this->mPPY = ppy;
	this->mEoI = EoI;
	this->mSmoothFact = smoothFact;
	

	// Check there are dates and consistency
	if (volDates.empty())
		throw KFailure("%s: no input volatility dates.\n", routine);

	if (factVol.size() != numFact) {
		throw KFailure("%s: factVol.size() (%d) != numFact (%d).\n",
			routine, factVol.size(), numFact);
	}
	numDates = volDates.size();
#ifdef	_SKIP
	for (i1=0; i1<numFact; i1++) {
	    if (volDates.size() != factVol[i1].size()) {
		throw KFailure("%s: volDates.size() (%d) != "
			"factVol[%d].size() (%d).\n", routine, 
	    		volDates.size(), i1, factVol[i1].size());
	    }
	}
#endif


	// today date
	mTodayDate = todayDt;
	Insert(mTodayDate);

	// Allocate the memory to store the vol information
	// in the tree and copy it over by doing "step" interpolation.
	// Because we have included all the vol dates in the critical dates
	// we have no problem.
	//
	if (numFact <1 || numFact >3)
		throw KFailure("%s: invalid tree dimension (%d).\n",
			       routine, numFact);
	mNbFactor = numFact;
	ASSERT_OR_THROW((mAlpha = DrlDoubleVectAlloc(0, mNbFactor-1)) != NULL);
	ASSERT_OR_THROW((mBeta  = DrlDoubleVectAlloc(0, mNbFactor-1)) != NULL);
	ASSERT_OR_THROW((mRho   = DrlDoubleMatrAlloc(0, mNbFactor-1,
						    0, mNbFactor-1)) != NULL);
	mVolDates = volDates;
	mFactVol = factVol;

	// Correlations
	for (i1=0; i1<mNbFactor; i1++) {
		mBeta[i1]  = factMr[i1];
		mAlpha[i1] = factWeight[i1];
		mRho[i1][i1] = 1e0;
		for (i2=i1+1; i2<mNbFactor; i2++) {
		    mRho[i1][i2] = mRho[i2][i1] 
				 = factCorr[RIDX(mNbFactor,i1,i2)];
		}
	}


	//
	// Debug info
	//
	if (debugLevel > DEBUG_LEVEL_VOLATILITY) {
		dppLog << format("%s:\n", routine);

		dppLog << format("Today:\n    %s",
			DrlTDatePrint(NULL, todayDt)) << endl;

		dppLog << "Ppy:\n    " << ppy << endl;
		dppLog << "EoI:\n    " << EoI << endl;
		dppLog << "Num Stdev:\n    " << numStdev<< endl;

		dppLog << "NumFact:\n    " << numFact << endl;

		dppLog << "Beta:" << endl;
		for (i1=0; i1<mNbFactor; i1++) 
			dppLog << format(" %12.8f", mBeta[i1]);
		dppLog << endl;

		dppLog << "Alpha:" << endl;
		for (i1=0; i1<mNbFactor; i1++) 
			dppLog << format(" %12.8f", mAlpha[i1]);
		dppLog << endl;

		dppLog << "Rho:" << endl;
		for (i1=0; i1<mNbFactor; i1++) {
		    for (i2=i1+1; i2<mNbFactor; i2++) 
				dppLog << format(" %12.8f", mRho[i1][i2]);
		    dppLog << endl;
		}

		dppLog << "Factor Volatilities:" << endl;
		for (i1=0; i1<numDates; i1++) {
		    dppLog << format(" %3d (%s) ", i1,
			DrlTDatePrint(NULL, volDates[i1]));

		    for (i2=0; i2<mNbFactor; i2++) 
			dppLog << format(" %12.8f", factVol[i2][i1]);
		    dppLog << endl;
		}



	}

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Initialize tree memory allocation that only depends on timeline,
// but indepedent of factor vols.  Timeline has been setup.
//
void
KMrNTree::MrNTreeLimitsNFactMemAlloc()
{

	ASSERT_OR_THROW((mSigma = DrlDoubleMatrAlloc(0, NbTP,
						     0, mNbFactor-1)) != NULL);
	ASSERT_OR_THROW((mAweight= DrlDoubleMatrAlloc(
			-1, NbTP, 0, mNbFactor*(mNbFactor+1)/2)) != NULL);

	//
	// Allocate memory for tree limits
	//
	mTop1    =          DrlIntVectAlloc(-1, NbTP+1);
	mBottom1 =          DrlIntVectAlloc(-1, NbTP+1);
	mTop2    = (int**)  DrlVoidPVectAlloc(-1, NbTP+1);
	mBottom2 = (int**)  DrlVoidPVectAlloc(-1, NbTP+1);
	mTop3    = (int***) DrlVoidPVectAlloc(-1, NbTP+1);
	mBottom3 = (int***) DrlVoidPVectAlloc(-1, NbTP+1);


	mOutTop1    =          DrlIntVectAlloc(-1, NbTP+1);
	mOutBottom1 =          DrlIntVectAlloc(-1, NbTP+1);
	mOutTop2    = (int**)  DrlVoidPVectAlloc(-1, NbTP+1);
	mOutBottom2 = (int**)  DrlVoidPVectAlloc(-1, NbTP+1);
	mOutTop3    = (int***) DrlVoidPVectAlloc(-1, NbTP+1);
	mOutBottom3 = (int***) DrlVoidPVectAlloc(-1, NbTP+1);

	mWidth     = DrlIntVectAlloc(0, mNbFactor-1);
	mHalfWidth = DrlIntVectAlloc(0, mNbFactor-1);


}



//--------------------------------------------------------------
// Compute factor vol and tree geometry only, assuming
// memory that only depends on timeline and indepedent of 
// factor vol have been setup.  
// This routine will be used by CET as well as Initialize() 
// for effeciency.
//
void
KMrNTree::ComputeFactorVolAndTreeLimit()
{
static	char	routine[] = "KMrNTree::ComputeFactorVolAndTreeLimit";

	int	size;
	int	i1, i2, idxT, idxV;

    try {


	// Set up spot vol at each time interval using step interpolation.
	for (i1=0; i1<mNbFactor; i1++) {
		for (idxT=0; idxT<=NbTP; idxT++) {
			// Perform Step interpolation
            if ( spotVolInterp == FLOORIDX)
            {
			    DrlTDateArrayFloorIdx(
				&(mVolDates[0]),
				mVolDates.size(),
				TPDates[idxT],
				&idxV);
            }
            else
            {
			    DrlTDateArrayCeilIdx_Exclude(
				&(mVolDates[0]),
				mVolDates.size(),
				TPDates[idxT],
				&idxV);
            }

			mSigma[idxT][i1] = (mFactVol[i1])[idxV];

		}
	}


	// Compute the ortthogonalized factor local covariances.
	//
	OrthogonalizeVol();


	// COMPUTE TREE LIMITS

	// Inner limits
	MrNTreeLimitsNFact(
		mWidth,
		mHalfWidth,
		mTop1,
		mBottom1,
		mTop2,
		mBottom2,
		mTop3,
		mBottom3,
		mNbFactor,
		mNbSigmaMax+1,
		NbTP,
		mBeta,
		mAweight,
		Length,
		LengthJ);


	if (debugLevel > DEBUG_LEVEL_GEOMETRY) {
	    dppLog << format("%s: Inner limits:\n", routine);
	    MrNTreeLimitsNFactPrint(
		mWidth,
		mHalfWidth,
		mTop1,
		mBottom1,
		mTop2,
		mBottom2,
		mTop3,
		mBottom3,
		mNbFactor,
		NbTP);
	}

    // Outer limits
	MrNTreeLimitsNFact(
		mWidth,
		mHalfWidth,
		mOutTop1,
		mOutBottom1,
		mOutTop2,
		mOutBottom2,
		mOutTop3,
		mOutBottom3,
		mNbFactor,
		mNbSigmaMax+1+outStdDevs,
		NbTP,
		mBeta,
		mAweight,
		Length,
		LengthJ);

	if (debugLevel > DEBUG_LEVEL_GEOMETRY) {
	    dppLog << format("%s: Outer limits:\n", routine);
	    MrNTreeLimitsNFactPrint(
		mWidth,
		mHalfWidth,
		mOutTop1,
		mOutBottom1,
		mOutTop2,
		mOutBottom2,
		mOutTop3,
		mOutBottom3,
		mNbFactor,
		NbTP);
	}

	// Check limits consistent
	//
	LimitsCheck();

	// Allocate storage for probabilities & transitions
	//
	if (NewPrice == NULL) {

		switch (mNbFactor) {
		case 3:
			ASSERT_OR_THROW((ru     = sliceNew(3)) != NULL);
			ASSERT_OR_THROW((r0     = sliceNew(3)) != NULL);
			ASSERT_OR_THROW((rd     = sliceNew(3)) != NULL);

			size = mWidth[0] * mWidth[1] * mWidth[2];
			ASSERT_OR_THROW((Shift3 = new int [size]) != NULL);
		case 2:
			ASSERT_OR_THROW((qu     = sliceNew(2)) != NULL);
			ASSERT_OR_THROW((q0     = sliceNew(2)) != NULL);
			ASSERT_OR_THROW((qd     = sliceNew(2)) != NULL);

			size = mWidth[0] * mWidth[1];
			ASSERT_OR_THROW((Shift2 = new int [size]) != NULL);
		case 1:
			ASSERT_OR_THROW((pu     = sliceNew(1)) != NULL);
			ASSERT_OR_THROW((p0     = sliceNew(1)) != NULL);
			ASSERT_OR_THROW((pd     = sliceNew(1)) != NULL);

			size = mWidth[0];
			ASSERT_OR_THROW((Shift1 = new int [size]) != NULL);
			break;
		default:
			throw KFailure("%s: not implemented.\n", routine);
		}

		ASSERT_OR_THROW((NewPrice = sliceNew(mNbFactor)) != NULL);
	}

	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << routine << ":Tree Critical Dates" << endl;
		dppLog << "    Index   Date" << endl;
		for (idxT=0; idxT<=NbTP; idxT++) {
			dppLog << format("    %5d   %s\n",
				idxT, DrlTDatePrint(NULL, TPDates[idxT]));
		}
	}

	// Debug info
	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s: \n", routine);

		dppLog << format(" NF  BETA         ALPHA       \n");
		for (i1=0; i1<mNbFactor; i1++) {
			dppLog << format(" %d    %12.8f  %12.6f\n",
				i1+1, mBeta[i1], mAlpha[i1]);
		}

		dppLog << format(" CORR\n");
		for (i1=0; i1<mNbFactor; i1++) {
		    for (i2=0; i2<mNbFactor; i2++) {
			dppLog << format("  %8.4f", mRho[i1][i2]);
		    }
		    dppLog << endl;
		}

		dppLog << format(" SPOT VOL\n");
		for (idxT=0; idxT<=NbTP; idxT++) {
			dppLog << format(" %3d/%3d %10.6f",
				idxT, NbTP,
				TPTimes[idxT]);
			for (i1=0; i1<mNbFactor; i1++) {
				dppLog << format("%12.6f", mSigma[idxT][i1]);
			}
			dppLog << endl;
		}
	}

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
// Called by CET ONLY.  Initialize the first nDim factors of 
// pure interest tree vols.
//
void
KMrNTree::InitializeFactorVol(
	KVector(TDate) &volDates,	// (I) volatility dates
	KVector(KVector(double)) &factVol)// (I) spot vols 

{
	mVolDates = volDates;
	
	//
	// Only the pure interest factor vols
	// in the first nIRDim factors
	//
	for(int i=0; i<factVol.size(); i++) 
		mFactVol[i] = factVol[i];
}


//--------------------------------------------------------------
// Update factor vol.  Used by CET only. 
//
void
KMrNTree::UpdateFactorVol(
	KVector(TDate) &volDates,	// (I) volatility dates
	KVector(KVector(double)) &factVol)// (I) spot vols 

{

	// Initialize factor vol
	//
	InitializeFactorVol(volDates,
			    factVol);

	// Compute tree vol and limits
	//
	ComputeFactorVolAndTreeLimit();

}



//--------------------------------------------------------------
// Setup factor vol.  
// Timeline has been setup.  Initialize memories for tree factor parameters 
// and tree limits that independent of factor vols, and compute
// tree vol and limits.  This should be called once in the
// tree Initialize().
//
void
KMrNTree::SetUpFactorVol()
{

	//
	// Initialize memories for tree parameters and 
	// tree limits that independent of factor vols.
	//
	MrNTreeLimitsNFactMemAlloc();

	// Compute tree vol and limits
	//
	ComputeFactorVolAndTreeLimit();

}


 

//--------------------------------------------------------------
// Setup tree timeline.  
//
void
KMrNTree::SetUpTimeline()
{
static	char	routine[] = "KMrNTree::SetUpTimeline";

    try {


	// After all critical dates are inserted,
	// sort and merge the critical dates in ascending order
	//
	IF_FAILED_THROW(DrlVTypeVectSort(TPDates, &NbTP, DRL_TDATE_T, TRUE));

	// Last critical date
	mLastDate = TPDate(NbTP-1);


	// Add all the benchmark vol calibration dates which are less than
	// the last critical dates.
	for(KVector(TDate)::iterator iterVolDate=mVolDates.begin();
				iterVolDate != mVolDates.end(); ++iterVolDate)
	{
		if (*iterVolDate < mLastDate)
			Insert(*iterVolDate);
	}

	// Set up timeline based on all the critical dates and PPY.
	InitializeTimeline(mEoI, mPPY);

	// Check the timeline
	CheckTimeline(mPPY);

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Calibrate should be called AFTER all the product critical dates
// are inserted.  It performs simply the TreeSetUp.

void
KMrNTree::Calibrate()
{
static	char	routine[] = "KMrNTree::Calibrate";

 try{

	TreeSetUp();

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
// TreeSetUp should be called AFTER all the product critical dates
// are inserted.  It performs:
// 1. Set up the timeline according to the ppy rule.
// 2. Interpolate the spot vols on each time step and Compute the 
//    jump step size.
// 3. Compute the orthogonal factor covariances (once the timeline is setup)
// 4. Compute the treelimits.

void
KMrNTree::TreeSetUp()
{

	//SetUpTimeline();

	SetUpFactorVol();

}




//--------------------------------------------------------------
// Compute the orthogonal factor covariances. 
// Called after all the timeline is set up.

void
KMrNTree::OrthogonalizeVol()
{
static	char	routine[] = "KMrNTree::OrthogonalizeVol";
	int	i1, i2, idxT;
	int	nFact = mNbFactor;	// Used in some macros.
	double	rhodet;

  try {

	for (idxT=0; idxT<=NbTP; idxT++) {
	    switch (mNbFactor) {
	    case 3:
		rhodet = 1e0 
				- mRho[0][1]*mRho[0][1]
				- mRho[0][2]*mRho[0][2]
				- mRho[1][2]*mRho[1][2]
				+ 2e0 * mRho[0][1]*mRho[0][2]*mRho[1][2];
		if (rhodet <= 0e0) {
			throw KFailure("%s: TP %3d/%3d (%10.6f) "
				" non positive 3D correlation matrix (det=%lf).\n",
				routine, idxT, NbTP, TPTimes[idxT],
				rhodet);
		}

			mAweight[idxT][AIDX(0, 2)] = mSigma[idxT][2] * mAlpha[2] *
				mRho[0][2];
			mAweight[idxT][AIDX(1, 2)] = mSigma[idxT][2] * mAlpha[2] *
				(mRho[1][2] - mRho[0][1]* mRho[0][2])
				/ sqrt(1e0 - mRho[0][1]*mRho[0][1]);

			mAweight[idxT][AIDX(2, 2)] = mSigma[idxT][2] * mAlpha[2] *
				sqrt(1e0 
					- mRho[0][1]*mRho[0][1]
					- mRho[0][2]*mRho[0][2]
					- mRho[1][2]*mRho[1][2]
					+ 2e0 * mRho[0][1]*mRho[0][2]*mRho[1][2]) 
				/ sqrt(1e0 - mRho[0][1]*mRho[0][1]);
		    case 2:
			rhodet = 1e0 - mRho[0][1]*mRho[0][1];
			if (rhodet <= 0e0) {
				throw KFailure("%s: TP %3d/%3d (%10.6f) "
					" non positive 2D correlation matrix (det=%lf).\n",
					routine, idxT, NbTP, TPTimes[idxT],
					rhodet);
			}

			mAweight[idxT][AIDX(0, 1)] = mSigma[idxT][1] * mAlpha[1] *
				mRho[0][1];
			mAweight[idxT][AIDX(1, 1)] = mSigma[idxT][1] * mAlpha[1] *
				sqrt(1e0 - mRho[0][1]*mRho[0][1]);
		    case 1:
			mAweight[idxT][AIDX(0, 0)] = mSigma[idxT][0] * mAlpha[0];
			break;
		    default:
			throw KFailure("%s: not implemented.\n", routine);
		    }
		}

		// Need Aweight at time -1 to compute grid at time 0.
		//
		for (i1=0;  i1<mNbFactor; i1++) 
		for (i2=i1; i2<mNbFactor; i2++) {
			mAweight[-1][AIDX(i1, i2)] = mAweight[0][AIDX(i1, i2)];
	}




	// Print debug information
	//
	if (debugLevel > DEBUG_LEVEL_DRIFT) {
		dppLog << format("%s:      ", routine);
		for (i1=0;  i1<mNbFactor; i1++) 
		for (i2=i1; i2<mNbFactor; i2++) 
			dppLog << format("     AWt%1d%2d[%d]", 
					i1, i2, AIDX(i1, i2));
		dppLog << endl;

		for (idxT=0; idxT<=NbTP; idxT++) {
			dppLog << format(" %3d/%3d %10s %5d (%10.6f)",
				idxT, NbTP,
				DrlTDatePrint(NULL, TPDates[idxT]),
				(int) (TPTimes[idxT] *365e0),
				TPTimes[idxT]);
			for (i1=0;  i1<mNbFactor; i1++) 
			for (i2=i1; i2<mNbFactor; i2++) {
				dppLog << format(" %12.6f",
					mAweight[idxT][AIDX(i1, i2)]);
			}
			dppLog << endl;
		}
	}
    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
// Lattice


int
MrNTreeComputeTransitionProba(
	KMrNTree &mrTree,		// (I) 
	int t,				// (I) Current time point index

	int nbFactor,			// (I) Number of factors      
	int NbTP,			// (I) Nb of time points      
	double *Beta,			// (I) Mean reversions        
	double **Aweight,		// (I) Orthogonal weights     
	double *Length,			// (I) Length of time steps   
	double *LengthJ,		// (I) Time step for jump size 

	int Top1,			// (I) Inner tree limits at t
	int Bottom1,			// (I) 
	int *Top2,			// (I) 
	int *Bottom2,			// (I) 
	int **Top3,			// (I) 
	int **Bottom3,			// (I) 
	int OutTop1,			// (I) Outer tree limits ate t+1
	int OutBottom1,			// (I) 
	int *OutTop2,			// (I) 
	int *OutBottom2,		// (I) 
	int **OutTop3,			// (I) 
	int **OutBottom3,		// (I) 

	double *puO,			// (O) output probabilities
	double *p0O,			// (O)
	double *pdO,			// (O)
	double *quO,			// (O)
	double *q0O,			// (O)
	double *qdO,			// (O)
	double *ruO,			// (O)
	double *r0O,			// (O)
	double *rdO,			// (O)

	int *Shift1O,			// (O) node transition
	int *Shift2O,			// (O)
	int *Shift3O)			// (O)
{
static	char	routine[] = "MrNTreeComputeTransitionProba";

	int	nFact = nbFactor;	// Used in some macros.

	double	jump[NFMAX][NFMAX],	// jump coeff [0<i<=j<nbFactor]
		pjump[NFMAX][NFMAX],	// old jump coeff [0<i<=j<nbFactor]
		zeta[NFMAX][NFMAX],
		d[NFMAX][NFMAX],	// used in transition proba calc
		b[NFMAX][NFMAX];	// used in transition proba calc
	int	shift[NFMAX];

	double	JumpCoeff,  du,		// Jump coefficients 
		dt,			// step size
		dtJ;			// step size for lattice 

	int	i, k;

	int	j0, j1, j2, ji,		// current time point index
		l0, l1, l2;		// center transition index

	int	outBotMin1, outTopMax1,
		outBotMin2, outTopMax2,
		outBotMin3, outTopMax3;

	int	offset;
	int	*Shift1L, *Shift2L, *Shift3L;	// intermediate pointers
	double	*puL, *p0L, *pdL,
		*quL, *q0L, *qdL,
		*ruL, *r0L, *rdL;


    try {

	// Nothing to do at the back of the tree
	//
	if (t == NbTP) {
		return (SUCCESS);
	}


	// Previous period jump size.
	//
	du = sqrt(JUMPCOEFF * LengthJ[t-1]);
	for (i=0; i<nbFactor; i++)
	for (k=0; k<=i      ; k++) {
		pjump[i][k] = Aweight[t-1][AIDX(i, k)] * du;
	}

	// Current period coefficients.
	//
	dtJ = LengthJ[t];                    
	dt  = Length[t];
	JumpCoeff = dt / (dtJ * JUMPCOEFF);

	du = sqrt (JUMPCOEFF * dtJ);
	for (i=0; i<nbFactor; i++)
	for (k=0; k<=i      ; k++) {
		jump[i][k] = Aweight[t][AIDX(i, k)] * du;
	}


	for (i=0; i<nbFactor; i++)
	for (k=0; k<=i      ; k++) {
		d[i][k] = (pjump[i][k] - jump[i][k]) / jump[i][i];
		b[i][k] = Beta[i] * dt * pjump[i][k] / jump[i][i];
	}


	//
	// Compute the transition probabilities and branching.
	//


        Shift1L = Shift1O + mrTree.NodeOffset(1, 0, 0, t);

	outBotMin1 = OutBottom1 + 1;
	outTopMax1 = OutTop1    - 1;

	if (outBotMin1 > outTopMax1)
		throw KFailure("%s: lMin > lMax.\n");


	// We do it for i=0 (ji = j0)
	//

	for (j0 = Bottom1; j0 <= Top1; j0++)
	{ 
	    i  = 0;
	    ji = j0;

	    zeta[i][i] = (d[i][i] - b[i][i]) * ji;
	    for (k=i+1; k<=nFact-1; k++) {
		zeta[k][i] = (d[k][i] - b[k][i]) * ji
				- zeta[i][i] * (jump[k][i] / jump[k][k]);
	    }

            shift[i] = NearInt(zeta[0][0]);
            shift[i] = Min(Max(outBotMin1 - j0 , shift[i]), outTopMax1 - ji);
	    l0 = j0 + shift[0];

	    offset = mrTree.NodeOffset(1, 0, 0, t);
	    puL     = puO     + offset;
	    p0L     = p0O     + offset;
	    pdL     = pdO     + offset;
	    Shift1L = Shift1O + offset;



            Shift1L[ji] = shift[i];
            zeta[i][i] -= shift[i];

	    ComputeProba(
		JumpCoeff,
		zeta[i][i],
		&puL[ji],
		&p0L[ji],
		&pdL[ji]);

	    if (nbFactor >= 2) {

	    outTopMax2 =                    OutTop2[j0+shift[0] -1];
	    outTopMax2 = Min(outTopMax2,    OutTop2[j0+shift[0]   ]);
	    outTopMax2 = Min(outTopMax2,    OutTop2[j0+shift[0] +1]) - 1;

	    outBotMin2 =                 OutBottom2[j0+shift[0] -1];
	    outBotMin2 = Max(outBotMin2, OutBottom2[j0+shift[0]   ]);
	    outBotMin2 = Max(outBotMin2, OutBottom2[j0+shift[0] +1]) + 1;

	    if (outBotMin2 > outTopMax2)
		throw KFailure("%s: outBotMin2 > outTopMax2\n");


	    // We do it for i=1 (ji = j1)
	    //
            for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++)
            {
		i = 1;
		ji = j1;

		zeta[i][i] = zeta[i][i-1] + (d[i][i] - b[i][i]) * ji;
		for (k=i+1; k<=nFact-1; k++) {
			zeta[k][i] = zeta[k][i-1] 
				+ (d[k][i] - b[k][i]) * ji
				- zeta[i][i] * (jump[k][i] / jump[k][k]);
		}

		shift[i] = NearInt(zeta[i][i]);
                shift[i] = Min(Max(shift[i], outBotMin2 - ji), outTopMax2 - ji);

		l1 = j1 + shift[1];

		offset = mrTree.NodeOffset(2, j0, 0, t);
		quL     = quO     + offset;
		q0L     = q0O     + offset;
		qdL     = qdO     + offset;
		Shift2L = Shift2O + offset;


                Shift2L[ji] = shift[i];
		zeta[i][i] -= shift[i];

		ComputeProba(JumpCoeff,
			zeta[i][i],
			&quL[ji],
			&q0L[ji],
			&qdL[ji]);

		if (nbFactor >= 3) {

                outTopMax3 =                    OutTop3[l0-1][l1-1];
                outTopMax3 = Min(outTopMax3,    OutTop3[l0-1][l1  ]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0-1][l1+1]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0  ][l1-1]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0  ][l1  ]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0  ][l1+1]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0+1][l1-1]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0+1][l1  ]);
                outTopMax3 = Min(outTopMax3,    OutTop3[l0+1][l1+1]) - 1;

                outBotMin3 =                 OutBottom3[l0-1][l1-1];
                outBotMin3 = Max(outBotMin3, OutBottom3[l0-1][l1  ]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0-1][l1+1]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0  ][l1-1]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0  ][l1  ]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0  ][l1+1]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0+1][l1-1]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0+1][l1  ]);
                outBotMin3 = Max(outBotMin3, OutBottom3[l0+1][l1+1]) + 1;


                if (outBotMin3 > outTopMax3)
		    throw KFailure("%s: outBotMin3 > outTopMax3.\n", routine);


                offset = mrTree.NodeOffset(3, j0, j1, t);
                ruL     = ruO     + offset;
                r0L     = r0O     + offset;
                rdL     = rdO     + offset;
                Shift3L = Shift3O + offset;


		// We do it for i=2 (ji = j2)
		//
		for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++)
		{
		    i = 2;
		    ji = j2;


		    zeta[i][i] = zeta[i][i-1] + (d[i][i] - b[i][i]) * ji;
		    for (k=i+1; k<=nFact-1; k++) {
			zeta[k][i] = zeta[k][i-1] 
				+ (d[k][i] - b[k][i]) * ji
				- zeta[i][i] * (jump[k][i] / jump[k][k]);
		    }

		    shift[i] = NearInt(zeta[i][i]);
                    shift[i] = Min(Max(shift[i], outBotMin3 - ji), outTopMax3 - ji);

		    l2 = j2 + shift[2];

                    Shift3L[ji] = shift[i];
		    zeta[i][i] -= shift[i];

		    ComputeProba(JumpCoeff,
			    zeta[i][i],
			    &ruL[ji],
			    &r0L[ji],
			    &rdL[ji]);


                }				// for j2
		}				// if (nbFactor >= 3)
            } 					// for j1
	    }					// if (nbFactor >= 2)
        }					// for j0


	// Print debug info
	//
	if (debugLevel >= DEBUG_LEVEL_PROBA) {


	    dppLog << format("%s:TP:%d", routine, t);
	    dppLog << format("\n");
        dppLog << format("Aweight = %f\n", Aweight[t][AIDX(0,0)]);


	    for (j0 = Bottom1; j0 <= Top1; j0++) {

		dppLog << format("(%4d,----,----)", j0);

		offset = mrTree.NodeOffset(1, 0, 0, t);
		puL     = puO     + offset;
		p0L     = p0O     + offset;
		pdL     = pdO     + offset;
		Shift1L = Shift1O + offset;
		l0 = j0 + Shift1L[j0];

		dppLog << format("<-(%4d,----,----)", l0);

		dppLog << format(" DIM1 SH=%4d  PU=%12.6f P0=%12.6f PD=%12.6f"
			" (BOT=%3d  TOP=%3d)\n",
			Shift1L[j0], puL[j0], p0L[j0], pdL[j0],
			OutBottom1, OutTop1);

	    if (nbFactor >=2) {
            for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {

		dppLog << format("(%4d,%4d,----)", j0, j1);

		offset = mrTree.NodeOffset(2, j0, 0, t);
		quL     = quO     + offset;
		q0L     = q0O     + offset;
		qdL     = qdO     + offset;
		Shift2L = Shift2O + offset;
		l1 = j1 + Shift2L[j1];

		dppLog << format("<-(%4d,%4d,----)", l0, l1);

		dppLog << format(" DIM2 SH=%4d  PU=%12.6f P0=%12.6f PD=%12.6f"
			" (BOT=%3d  TOP=%3d)\n",
			Shift2L[j1], quL[j1], q0L[j1], qdL[j1],
			OutBottom2[l0], OutTop2[l0]);


	    if (nbFactor >= 3) {
	    for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {

		dppLog << format("(%4d,%4d,%4d)", j0, j1, j2);

                offset = mrTree.NodeOffset(3, j0, j1, t);
                ruL     = ruO     + offset;
                r0L     = r0O     + offset;
                rdL     = rdO     + offset;
                Shift3L = Shift3O + offset;
		l2 = j2 + Shift3L[j2];

		dppLog << format("<-(%4d,%4d,%4d)", l0, l1, l2);

		dppLog << format(" DIM3 SH=%4d  PU=%12.6f P0=%12.6f PD=%12.6f"
			" (BOT=%3d  TOP=%3d)\n",
			Shift3L[j2], ruL[j2], r0L[j2], rdL[j2],
			OutBottom3[l0][j1], OutTop3[l0][j1]);
	    }}
	    }}
	    }
	}

	return(SUCCESS);

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}


//--------------------------------------------------------------
//

void
KMrNTree::Update(int tpIdx)
{
static	char	routine[] = "KMrNTree::Update";
	//int	size;

    try {

	// Set the tree time point
	tpIdxCurrent = tpIdx;


	// Compute probabilities
	//

	MrNTreeComputeTransitionProba(
		*this,
		tpIdx,
		mNbFactor,		// (I) Number of factors      
		NbTP,			// (I) Nb of time points      
		mBeta,			// (I) Mean reversions        
		mAweight,		// (I) Orthogonal weights     
		Length,			// (I) Length of time steps   
		LengthJ,		// (I) Time step for jump size 

		mTop1   [tpIdx],
		mBottom1[tpIdx],
		mTop2   [tpIdx],
		mBottom2[tpIdx],
		mTop3   [tpIdx],
		mBottom3[tpIdx],

		mOutTop1   [tpIdx+1],
		mOutBottom1[tpIdx+1],
		mOutTop2   [tpIdx+1],
		mOutBottom2[tpIdx+1],
		mOutTop3   [tpIdx+1],
		mOutBottom3[tpIdx+1],

		pu,
		p0,
		pd,
		qu,
		q0,
		qd,
		ru,
		r0,
		rd,

		Shift1,
		Shift2,
		Shift3);



    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
// Compute the expected value of a slice at the current
// timepoint tpIdx.



void
KMrNTree::sliceEv(
	double *sliceVal,		// (B) values to be discounted 
	int sliceDim,			// (I) slice dimension
	double *discVal,		// (I) discounting ts
	int tpIdx)			// (I) current time period 
{
static	char	routine[] = "KMrNTree::sliceEv";
	int	i0, i1, i2,
		l0,
		l1, l1u, l1d,
		l2, l2u, l2d;

	int	*shift1T,		// tmp with offset
		*shift2T,
		*shift3T,
		offset1,		// offset for 1-D slices
		offset2,		// etc.
		offset3;

	int	Top0,			// Inner limits at t
		Bot0,
		*Top1,
		*Bot1,
		**Top2,
		**Bot2;
	int	InTop0,			// Inner limits at t+1
		InBot0,
		*InTop1,
		*InBot1,
		**InTop2,
		**InBot2;
	int	OuTop0,			// Outer limits at t+1
		OuBot0,
		*OuTop1,
		*OuBot1,
		**OuTop2,
		**OuBot2;

	double	*puT, *p0T, *pdT,	// offset probabilities arrays
		*quT, *q0T, *qdT,
		*ruT, *r0T, *rdT;

					// tmp storage of 1-D 3 probabilities
	double	puI,  p0I,  pdI,
		quuI, qu0I, qudI,
					// tmp storage of 2-D 9 probabilities
		q0uI, q00I, q0dI,
		qduI, qd0I, qddI;
					// tmp storage of 3-D 27 probabilities
		/*ruuuI, ruu0I, ruudI,
		ru0uI, ru00I, ru0dI,
		ruduI, rud0I, ruddI,

		r0uuI, r0u0I, r0udI,
		r00uI, r000I, r00dI,
		r0duI, r0d0I, r0ddI,

		rduuI, rdu0I, rdudI,
		rd0uI, rd00I, rd0dI,
		rdduI, rdd0I, rdddI;*/


	double	PFlat,
		*discValT,		// local slice pointers
		*sliceValT,
		*newValT,
		*svuT, *sv0T, *svdT,
		*svuuT, *sv0uT, *svduT,
		*svu0T, *sv00T, *svd0T,
		*svudT, *sv0dT, *svddT;

 try{

	// Nothing to do at the back of the tree 
	if (tpIdx == TPNum()) return;


	// Store bounds for efficiency
	Top0   = this->mTop1[tpIdx];
	Top1   = this->mTop2[tpIdx];
	Top2   = this->mTop3[tpIdx];
	Bot0   = this->mBottom1[tpIdx];
	Bot1   = this->mBottom2[tpIdx];
	Bot2   = this->mBottom3[tpIdx];

	InTop0 = this->mTop1[tpIdx+1];
	InTop1 = this->mTop2[tpIdx+1];
	InTop2 = this->mTop3[tpIdx+1];
	InBot0 = this->mBottom1[tpIdx+1];
	InBot1 = this->mBottom2[tpIdx+1];
	InBot2 = this->mBottom3[tpIdx+1];

	OuTop0 = this->mOutTop1[tpIdx+1];
	OuTop1 = this->mOutTop2[tpIdx+1];
	OuTop2 = this->mOutTop3[tpIdx+1];
	OuBot0 = this->mOutBottom1[tpIdx+1];
	OuBot1 = this->mOutBottom2[tpIdx+1];
	OuBot2 = this->mOutBottom3[tpIdx+1];



	if(debugLevel > DEBUG_LEVEL_DEV) {
		dppLog << endl;
		dppLog << format("============ %s ============", routine)<<endl;
		dppLog << "Slice BEFORE Dev:" << endl;
 
        	slicePrint(
                   	  sliceVal,
                     	  sliceDim,
                   	  tpIdx+1,
                   	  FALSE,
                   	  dppLog);
	}




	// Do it according the slice dimension


	if (sliceDim == 1)
	{
		// Flat values for nodes between the inner ellipsoid 
		// where values have been calculated and the outer 
		// ellipsoid defined by OuBot and OuTop.

		sliceValT = sliceVal + NodeOffset(1, 0, 0, tpIdx+1);

		for (i0 = OuBot0; i0 < InBot0; i0++)
		{
		    sliceValT[i0] = sliceValT[InBot0];
		}

		for (i0 = InTop0+1;  i0 <= OuTop0; i0++)
		{
		    sliceValT[i0] = sliceValT[InTop0];
		}

		// Discounted expected value
		//
		offset1 = NodeOffset(1, 0, 0, tpIdx);
		shift1T = this->Shift1 + offset1;

		puT = this->pu + offset1;
		p0T = this->p0 + offset1;
		pdT = this->pd + offset1;

		discValT = discVal        + offset1;
		newValT  = this->NewPrice + offset1;

		if (discVal) {
			for (i0 = Bot0; i0 <= Top0; i0 ++)                  
			{
				l0 = i0 + shift1T[i0];

				newValT[i0] =
					puT[i0] * sliceValT[l0+1] +
					p0T[i0] * sliceValT[l0  ] +
					pdT[i0] * sliceValT[l0-1];

				newValT[i0] *= discValT[i0];
			}
		} else {
			for (i0 = Bot0; i0 <= Top0; i0 ++)                  
			{
				l0 = i0 + shift1T[i0];

				newValT[i0] =
					puT[i0] * sliceValT[l0+1] +
					p0T[i0] * sliceValT[l0  ] +
					pdT[i0] * sliceValT[l0-1];

			}
		}


		// Copy values back to original slice
		sliceValT = sliceVal + offset1;
		for (i0 = Bot0; i0 <= Top0; i0 ++)
		{
			sliceValT[i0] = newValT[i0];
		}


	}
	else if (sliceDim == 2)
	{
		//
		// Copy values between inner and outer limits.
		//


		sliceValT = sliceVal + NodeOffset(2, InBot0, 0, tpIdx+1);
		PFlat = sliceValT[InTop1[InBot0]];
		for (i0 = OuBot0; i0 < InBot0; i0++)
		{
		    sliceValT = sliceVal + NodeOffset(2, i0, 0, tpIdx+1);
		    for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++)
		    {
			sliceValT[i1] = PFlat;
		    }
		} 


		for (i0 = InBot0; i0 <= InTop0; i0++)
		{
		    sliceValT = sliceVal + NodeOffset(2, i0, 0, tpIdx+1);

		    PFlat = sliceValT[InBot1[i0]];
		    for (i1 = OuBot1[i0]; i1 < InBot1[i0]; i1++)
		    {
			sliceValT[i1] = PFlat;
		    }

		    PFlat = sliceValT[InTop1[i0]];
		    for (i1 = InTop1[i0]+1; i1 <= OuTop1[i0];    i1++)
		    {
			sliceValT[i1] = PFlat;
		    }
		}


		sliceValT = sliceVal + NodeOffset(2, InTop0, 0, tpIdx+1);

		PFlat = sliceValT[InTop1[InTop0]];
		for (i0 = InTop0+1; i0 <= OuTop0; i0++)
		{
		    sliceValT = sliceVal + NodeOffset(2, i0, 0, tpIdx+1);
		    for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++)
		    {
			sliceValT[i1] = PFlat;
		    }
		}


		//
		// Perform actual expected value
		//

		offset1 = NodeOffset(1, 0, 0, tpIdx);
		shift1T = this->Shift1 + offset1;
		puT     = this->pu     + offset1;
		p0T     = this->p0     + offset1;
		pdT     = this->pd     + offset1;


		for (i0 = Bot0; i0 <= Top0; i0++)
		{
		    // compute 1-D 3 probabilities and transitions
		    puI = puT[i0];
		    p0I = p0T[i0];
		    pdI = pdT[i0];
		    l0  = i0 + shift1T[i0];

		    // offset arrays for next loop 
		    offset2 = NodeOffset(2, i0, 0, tpIdx);
		    shift2T = this->Shift2 + offset2;
		    quT     = this->qu     + offset2;
		    q0T     = this->q0     + offset2;
		    qdT     = this->qd     + offset2;


		    // Inner loop optimized: offset all needed slices
		    // and store to tmp pointers
		    discValT = discVal         + offset2;
		    newValT  = this->NewPrice  + offset2;

		    svuT = sliceVal + NodeOffset(2, l0+1, 0, tpIdx+1);
		    sv0T = sliceVal + NodeOffset(2, l0  , 0, tpIdx+1);
		    svdT = sliceVal + NodeOffset(2, l0-1, 0, tpIdx+1);


		    for (i1 = Bot1[i0]; i1 <= Top1[i0]; i1++)
		    {
			l1  = i1 + shift2T[i1];
			l1u = l1 + 1;
			l1d = l1 - 1;

			newValT[i1] = 
				puI * (quT[i1] * svuT[l1u] +
				       q0T[i1] * svuT[l1 ] +
				       qdT[i1] * svuT[l1d]) +
				p0I * (quT[i1] * sv0T[l1u] +
				       q0T[i1] * sv0T[l1 ] +
				       qdT[i1] * sv0T[l1d]) +
				pdI * (quT[i1] * svdT[l1u] +
				       q0T[i1] * svdT[l1 ] +
				       qdT[i1] * svdT[l1d]);

			if (discVal)
				newValT[i1] *= discValT[i1];
		    }
		}


		for (i0 = Bot0; i0 <= Top0; i0++)
		{        
		    offset2    = NodeOffset(2, i0, 0, tpIdx);
		    sliceValT = sliceVal       + offset2;
		    newValT   = this->NewPrice + offset2;

		    for (i1 = Bot1[i0]; i1 <= Top1[i0]; i1++)
		    {	
			sliceValT[i1] = newValT[i1];
		    }
		}

	}
	else if (sliceDim == 3)
	{
		//
		// Copy values between inner and outer limits.
		//

		sliceValT = sliceVal + NodeOffset(3, InBot0, InTop1[InBot0], tpIdx+1);
		PFlat = sliceValT[InTop2[InBot0][InTop1[InBot0]]];

		for (i0 = OuBot0;     i0 <  InBot0;     i0++)
		for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++)
		{
		    sliceValT = sliceVal + NodeOffset(3, i0, i1, tpIdx+1);

		    for (i2 = OuBot2[i0][i1]; i2 <= OuTop2[i0][i1]; i2++)
		    {
			sliceValT[i2] = PFlat;
		    }
		}


		for (i0 = InBot0; i0 <= InTop0; i0++)
		{
		    sliceValT = sliceVal + NodeOffset(3, i0, InBot1[i0], tpIdx+1);

		    PFlat = sliceValT[InTop2[i0][InBot1[i0]]];

		    for (i1 = OuBot1[i0]; i1 < InBot1[i0]; i1++)
		    {
			sliceValT = sliceVal + NodeOffset(3, i0, i1, tpIdx+1);

			for (i2 = OuBot2[i0][i1]; i2 <= OuTop2[i0][i1]; i2++)
			{
			    sliceValT[i2] = PFlat;
			}
		    }

		    for (i1 = InBot1[i0]; i1 <= InTop1[i0]; i1++)
		    {
			sliceValT = sliceVal + NodeOffset(3, i0, i1, tpIdx+1);

			PFlat = sliceValT[InBot2[i0][i1]];
			for (i2 = OuBot2[i0][i1]; i2 < InBot2[i0][i1]; i2++)
			{
			    sliceValT[i2] = PFlat;
			}


			PFlat = sliceValT[InTop2[i0][i1]];
			for (i2 = InTop2[i0][i1]+1; i2 <= OuTop2[i0][i1]; i2++)
			{
			    sliceValT[i2] = PFlat;
			}
		    }


		    sliceValT = sliceVal + NodeOffset (3, i0, InTop1[i0], tpIdx+1);
		    PFlat = sliceValT[InTop2[i0][InTop1[i0]]];

		    for (i1 = InTop1[i0]+1; i1 <= OuTop1[i0]; i1++)
		    {
			sliceValT = sliceVal + NodeOffset(3, i0, i1, tpIdx+1);
			for (i2 = OuBot2[i0][i1]; i2 <= OuTop2[i0][i1]; i2++)
			{
			   sliceValT[i2] = PFlat;
			}
		    }
		}


		sliceValT = sliceVal + NodeOffset(3, InTop0, InTop1[InTop0], tpIdx+1);
		PFlat = sliceValT[InTop2[InTop0][InTop1[InTop0]]];

		for (i0 = InTop0+1;   i0 <= OuTop0;     i0++)
		for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++)
		{
                    sliceValT = sliceVal + NodeOffset (3, i0, i1, tpIdx+1);

                    for (i2 = OuBot2[i0][i1]; i2 <= OuTop2[i0][i1]; i2++)
                    {
			sliceValT[i2] = PFlat;
                    }
		}

		//
		// Perform actual expected value
		//

		offset1 = NodeOffset(1, 0, 0, tpIdx);
		shift1T = this->Shift1 + offset1;
		puT     = this->pu     + offset1;
		p0T     = this->p0     + offset1;
		pdT     = this->pd     + offset1;


		for (i0 = Bot0; i0 <= Top0; i0++)
		{
		    // compute 1-D 3 probabilities and transitions
		    puI = puT[i0];
		    p0I = p0T[i0];
		    pdI = pdT[i0];
		    l0 = i0 + shift1T[i0];

		    // offset arrays for next loop 
		    offset2 = NodeOffset(2, i0, 0, tpIdx);
		    shift2T = this->Shift2 + offset2;
		    quT     = this->qu     + offset2;
		    q0T     = this->q0     + offset2;
		    qdT     = this->qd     + offset2;

		    for (i1 = Bot1[i0]; i1 <= Top1[i0]; i1++)
		    {

			// compute 2-D 9 probabilities and transitions
			quuI = puI * quT[i1];
			qu0I = puI * q0T[i1];
			qudI = puI * qdT[i1];
			q0uI = p0I * quT[i1];
			q00I = p0I * q0T[i1];
			q0dI = p0I * qdT[i1];
			qduI = pdI * quT[i1];
			qd0I = pdI * q0T[i1];
			qddI = pdI * qdT[i1];
			l1 = i1 + shift2T[i1];

			// offset arrays for next loop 
			offset3 = NodeOffset(3, i0, i1, tpIdx);
			shift3T = this->Shift3 + offset3;
			ruT     = this->ru     + offset3;
			r0T     = this->r0     + offset3;
			rdT     = this->rd     + offset3;

			// Inner loop optimized: offset all needed slices
			// and store to tmp pointers
			discValT = discVal        + offset3;
			newValT  = this->NewPrice + offset3;

			svuuT = sliceVal + NodeOffset(3, l0+1, l1+1, tpIdx+1);
			svu0T = sliceVal + NodeOffset(3, l0+1, l1  , tpIdx+1);
			svudT = sliceVal + NodeOffset(3, l0+1, l1-1, tpIdx+1);
			sv0uT = sliceVal + NodeOffset(3, l0  , l1+1, tpIdx+1);
			sv00T = sliceVal + NodeOffset(3, l0  , l1  , tpIdx+1);
			sv0dT = sliceVal + NodeOffset(3, l0  , l1-1, tpIdx+1);
			svduT = sliceVal + NodeOffset(3, l0-1, l1+1, tpIdx+1);
			svd0T = sliceVal + NodeOffset(3, l0-1, l1  , tpIdx+1);
			svddT = sliceVal + NodeOffset(3, l0-1, l1-1, tpIdx+1);


			for (i2 = Bot2[i0][i1]; i2 <= Top2[i0][i1]; i2++) {
			    l2  = i2 + shift3T[i2];
			    l2u = l2 + 1;
			    l2d = l2 - 1;

			    if ((l2u > OuTop2[l0][l1]) || (l2d < OuBot2[l0][l1])) {
			    	printf("ERROR EV:%3d (%3d %3d %3d)<-(%3d %3d %3d)\n",
					tpIdx, i0, i1, i2, l0, l1, l2);

				printf("CURRENT IN LIMITS: "
					"(%3d %3d %3d) (%3d %3d %3d)\n",
					Top0, Top1[i0], Top2[i0][i1],
					Bot0, Bot1[i0], Bot2[i0][i1]);

				printf("NEXT    IN LIMITS: "
					"(%3d %3d %3d) (%3d %3d %3d)\n",
					InTop0, InTop1[i0], InTop2[i0][i1],
					InBot0, InBot1[i0], InBot2[i0][i1]);

				printf("NEXT    OU LIMITS: "
					"(%3d %3d %3d) (%3d %3d %3d)\n",
					OuTop0, OuTop1[i0], OuTop2[i0][i1],
					OuBot0, OuBot1[i0], OuBot2[i0][i1]);
			    }


			    ASSERT_OR_THROW (l2u <= OuTop2[l0][l1]);
			    ASSERT_OR_THROW (l2d >= OuBot2[l0][l1]);

/*
#define	PRT(z)	printf("%s = %lf\n", #z, z);
			    printf("EV:%3d (%3d %3d %3d)<-(%3d %3d %3d)\n",
					tpIdx, i0, i1, i2, l0, l1, l2);

				PRT(quuI);
				PRT(ruT[i2]);
				PRT(svuuT[l2u]);
				PRT(r0T[i2]);
				PRT(svuuT[l2 ]);
				PRT(rdT[i2]);
				PRT(svuuT[l2d]);
				PRT(qu0I);
				PRT(ruT[i2]);
				PRT(svu0T[l2u]);
				PRT(r0T[i2]);
				PRT(svu0T[l2 ]);
				PRT(rdT[i2]);
				PRT(svu0T[l2d]);
				PRT(qudI);
				PRT(ruT[i2]	);
				PRT(svudT[l2u]);
				PRT(r0T[i2]);
				PRT(svudT[l2 ]);
				PRT(rdT[i2]);
				PRT(svudT[l2d]);

				PRT(q0uI);
				PRT(ruT[i2]);
				PRT(sv0uT[l2u]);
				PRT(r0T[i2]);
				PRT(sv0uT[l2 ]);
				PRT(rdT[i2]);
				PRT(sv0uT[l2d]);
				PRT(q00I);
				PRT(ruT[i2]);
				PRT(sv00T[l2u]);
				PRT(r0T[i2]);
				PRT(sv00T[l2 ]);
				PRT(rdT[i2]);
				PRT(sv00T[l2d]);
				PRT(q0dI);
				PRT(ruT[i2]);
				PRT(sv0dT[l2u]);
				PRT(r0T[i2]);
				PRT(sv0dT[l2 ]);
				PRT(rdT[i2]);
				PRT(sv0dT[l2d]);

				PRT(qduI);
				PRT(ruT[i2]);
				PRT(svduT[l2u]);
				PRT(r0T[i2]);
				PRT(svduT[l2 ]);
				PRT(rdT[i2]);
				PRT(svduT[l2d]);
				PRT(qd0I);
				PRT(ruT[i2]);
				PRT(svd0T[l2u]);
				PRT(r0T[i2]);
				PRT(svd0T[l2 ]);
				PRT(rdT[i2]);
				PRT(svd0T[l2d]);
				PRT(qddI);
				PRT(ruT[i2]);
				PRT(svddT[l2u]);
				PRT(r0T[i2]);
				PRT(svddT[l2 ]);
				PRT(rdT[i2]);
				PRT(svddT[l2d]);
//xxxxxxxxxxxxxxxxxxxxxxx
*/






			    newValT[i2] =
				quuI * (ruT[i2] * svuuT[l2u] +
					r0T[i2] * svuuT[l2 ] +
					rdT[i2] * svuuT[l2d]) +
				qu0I * (ruT[i2] * svu0T[l2u] +
					r0T[i2] * svu0T[l2 ] +
					rdT[i2] * svu0T[l2d]) +
				qudI * (ruT[i2] * svudT[l2u] +
					r0T[i2] * svudT[l2 ] +
					rdT[i2] * svudT[l2d]) +

				q0uI * (ruT[i2] * sv0uT[l2u] +
					r0T[i2] * sv0uT[l2 ] +
					rdT[i2] * sv0uT[l2d]) +
				q00I * (ruT[i2] * sv00T[l2u] +
					r0T[i2] * sv00T[l2 ] +
					rdT[i2] * sv00T[l2d]) +
				q0dI * (ruT[i2] * sv0dT[l2u] +
					r0T[i2] * sv0dT[l2 ] +
					rdT[i2] * sv0dT[l2d]) +

				qduI * (ruT[i2] * svduT[l2u] +
					r0T[i2] * svduT[l2 ] +
					rdT[i2] * svduT[l2d]) +
				qd0I * (ruT[i2] * svd0T[l2u] +
					r0T[i2] * svd0T[l2 ] +
					rdT[i2] * svd0T[l2d]) +
				qddI * (ruT[i2] * svddT[l2u] +
					r0T[i2] * svddT[l2 ] +
					rdT[i2] * svddT[l2d]);

			    if (discVal)
				newValT[i2] *= discValT[i2];
			} // for (i2
		    } // for (i1
		} // for (i0

		// Copy values back to slice
		for (i0 = Bot0; i0 <= Top0; i0++)
		{
		    for (i1 = Bot1[i0]; i1 <= Top1[i0]; i1++)
		    {	
			offset3 = NodeOffset(3, i0, i1, tpIdx);
			sliceValT = sliceVal       + offset3;
			newValT   = this->NewPrice + offset3;

			for (i2 = Bot2[i0][i1]; i2 <= Top2[i0][i1]; i2++)
			{
				sliceValT[i2] = newValT[i2];
			} // for (i2
		    } // for (i1
		} // for (i0

	}

	if(debugLevel > DEBUG_LEVEL_DEV) {
		dppLog << "Slice AFTER Dev:" << endl;
 
        	slicePrint(
                   	  sliceVal,
                     	  sliceDim,
                   	  tpIdx,
                   	  FALSE,
                   	  dppLog);
	}

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Forwards a time slice between tpIdx and tpIdx+1, i.e.
// performs the adjoint of sliceEv.
//


void
KMrNTree::sliceFw(
	double *sliceVal,		// (B) values to be forwarded
	int sliceDim,			// (I) slice dimension
	double *discVal,		// (I) discounting ts
	int tpIdx)			// (I) current time period 
{
static	char	routine[] = "KMrNTree::sliceFw";
	int	i0, i1, i2,
		l0,
		l1, l1u, l1d,
		l2, l2u, l2d;

	int	*shift1T,		// tmp with offset
		*shift2T,
		*shift3T,
		offset1,		// offset for 1-D slices
		offset2,		// etc.
		offset3;

	int	Top0,			// Inner limits at t
		Bot0,
		*Top1,
		*Bot1,
		**Top2,
		**Bot2;
	int	InTop0,			// Inner limits at t+1
		InBot0,
		*InTop1,
		*InBot1,
		**InTop2,
		**InBot2;
	int	OuTop0,			// Outer limits at t+1
		OuBot0,
		*OuTop1,
		*OuBot1,
		**OuTop2,
		**OuBot2;

	double	*puT, *p0T, *pdT,	// offset probabilities arrays
		*quT, *q0T, *qdT,
		*ruT, *r0T, *rdT;

					// tmp storage of 1-D 3 probabilities
	double	puI,  p0I,  pdI,
		quuI, qu0I, qudI,
					// tmp storage of 2-D 9 probabilities
		q0uI, q00I, q0dI,
		qduI, qd0I, qddI;
					// tmp storage of 3-D 27 probabilities
		/*ruuuI, ruu0I, ruudI,
		ru0uI, ru00I, ru0dI,
		ruduI, rud0I, ruddI,

		r0uuI, r0u0I, r0udI,
		r00uI, r000I, r00dI,
		r0duI, r0d0I, r0ddI,

		rduuI, rdu0I, rdudI,
		rd0uI, rd00I, rd0dI,
		rdduI, rdd0I, rdddI;*/


	double	*newVal = this->NewPrice;// tmp slice
	double	*discValT,		// local slice pointers
		*sliceValT,
		*newValT,
		*nvuT, *nv0T, *nvdT,
		*nvuuT, *nv0uT, *nvduT,
		*nvu0T, *nv00T, *nvd0T,
		*nvudT, *nv0dT, *nvddT;

	double	x;			// discounting applied

 try{
	// Nothing to do at the back of the tree 
	if (tpIdx == TPNum()) return;

	if(debugLevel > DEBUG_LEVEL_DEV) {
		dppLog << endl;
		dppLog << format("============ %s ============", routine)<<endl;
		dppLog << "Slice BEFORE being forwarded:" << endl;
 
        	slicePrint(
                   	  sliceVal,
                     	  sliceDim,
                   	  tpIdx,
                   	  FALSE,
                   	  dppLog);
	}


	// Store bounds for efficiency
	Top0   = this->mTop1[tpIdx];
	Top1   = this->mTop2[tpIdx];
	Top2   = this->mTop3[tpIdx];
	Bot0   = this->mBottom1[tpIdx];
	Bot1   = this->mBottom2[tpIdx];
	Bot2   = this->mBottom3[tpIdx];

	InTop0 = this->mTop1[tpIdx+1];
	InTop1 = this->mTop2[tpIdx+1];
	InTop2 = this->mTop3[tpIdx+1];
	InBot0 = this->mBottom1[tpIdx+1];
	InBot1 = this->mBottom2[tpIdx+1];
	InBot2 = this->mBottom3[tpIdx+1];

	OuTop0 = this->mOutTop1[tpIdx+1];
	OuTop1 = this->mOutTop2[tpIdx+1];
	OuTop2 = this->mOutTop3[tpIdx+1];
	OuBot0 = this->mOutBottom1[tpIdx+1];
	OuBot1 = this->mOutBottom2[tpIdx+1];
	OuBot2 = this->mOutBottom3[tpIdx+1];




	// Do it according the slice dimension


	if (sliceDim == 1)
	{
		//
		// Clear tmp slice in order to store new values
		//
		newValT = newVal + NodeOffset(1, 0, 0, tpIdx+1);
		for (i0 = OuBot0; i0 <= OuTop0; i0++)
		{
		    newValT[i0] = 0e0;
		}

		//
		// Calculate forward value.
		//

		offset1 = NodeOffset(1, 0, 0, tpIdx);
		shift1T = this->Shift1 + offset1;
		puT = this->pu + offset1;
		p0T = this->p0 + offset1;
		pdT = this->pd + offset1;

		discValT  = discVal  + offset1;
		sliceValT = sliceVal + offset1;

		for (i0 = Bot0; i0 <= Top0; i0 ++)                  
		{
			l0 = i0 + shift1T[i0];

			if (discVal)
				x = sliceValT[i0] * discValT[i0];
			else
				x = sliceValT[i0];

			newValT[l0+1] += puT[i0] * x;
			newValT[l0  ] += p0T[i0] * x;
			newValT[l0-1] += pdT[i0] * x;

		}


		//
		// Copy values back to output slice
		//
		sliceValT = sliceVal + offset1;
		for (i0 = OuBot0; i0 <= OuTop0; i0 ++)
		{
			sliceValT[i0] = newValT[i0];
		}


	}
	else if (sliceDim == 2)
	{
		//
		// Clear tmp slice in order to store new values
		//
		for (i0 = OuBot0; i0 <= OuTop0; i0++)
		{        
		    offset2    = NodeOffset(2, i0, 0, tpIdx+1);
		    newValT   = newVal   + offset2;
		    for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++)
		    {	
			newValT[i1] = 0e0;
		    }
		}

		//
		// Calculate forward value
		//

		offset1 = NodeOffset(1, 0, 0, tpIdx);
		shift1T = this->Shift1 + offset1;
		puT     = this->pu     + offset1;
		p0T     = this->p0     + offset1;
		pdT     = this->pd     + offset1;


		for (i0 = Bot0; i0 <= Top0; i0++)
		{
		    // compute 1-D 3 probabilities and transitions
		    puI = puT[i0];
		    p0I = p0T[i0];
		    pdI = pdT[i0];
		    l0  = i0 + shift1T[i0];

		    // offset arrays for next loop 
		    offset2 = NodeOffset(2, i0, 0, tpIdx);
		    shift2T = this->Shift2 + offset2;
		    quT     = this->qu     + offset2;
		    q0T     = this->q0     + offset2;
		    qdT     = this->qd     + offset2;


		    // Inner loop optimized: offset all needed slices
		    // and store to tmp pointers
		    discValT   = discVal  + offset2;
		    sliceValT  = sliceVal + offset2;

		    nvuT = newVal + NodeOffset(2, l0+1, 0, tpIdx+1);
		    nv0T = newVal + NodeOffset(2, l0  , 0, tpIdx+1);
		    nvdT = newVal + NodeOffset(2, l0-1, 0, tpIdx+1);


		    for (i1 = Bot1[i0]; i1 <= Top1[i0]; i1++)
		    {
			l1  = i1 + shift2T[i1];
			l1u = l1 + 1;
			l1d = l1 - 1;

			if (discVal)
				x = sliceValT[i1] * discValT[i1];
			else
				x = sliceValT[i1];



			nvuT[l1u] += puI * quT[i1] * x;
			nvuT[l1 ] += puI * q0T[i1] * x;
			nvuT[l1d] += puI * qdT[i1] * x;
			nv0T[l1u] += p0I * quT[i1] * x;
			nv0T[l1 ] += p0I * q0T[i1] * x;
			nv0T[l1d] += p0I * qdT[i1] * x;
			nvdT[l1u] += pdI * quT[i1] * x;
			nvdT[l1 ] += pdI * q0T[i1] * x;
			nvdT[l1d] += pdI * qdT[i1] * x;

		    }
		}


		//
		// Copy values back to output slice
		//
		for (i0 = OuBot0; i0 <= OuTop0; i0++)
		{        
		    offset2    = NodeOffset(2, i0, 0, tpIdx+1);
		    sliceValT = sliceVal       + offset2;
		    newValT   = this->NewPrice + offset2;

		    for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++)
		    {	
			sliceValT[i1] = newValT[i1];
		    }
		}

	}
	else if (sliceDim == 3)
	{
		//
		// Clear tmp slice in order to store new values
		//
		for (i0 = OuBot0; i0 <= OuTop0; i0++)
		for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++) {	
			offset3 = NodeOffset(3, i0, i1, tpIdx+1);
			newValT   = newVal + offset3;
			for (i2 = OuBot2[i0][i1]; i2 <= OuTop2[i0][i1]; i2++)
			{
				newValT[i2] = 0e0;
			}
		}


		//
		// Calculate forward value.
		//

		offset1 = NodeOffset(1, 0, 0, tpIdx);
		shift1T = this->Shift1 + offset1;
		puT     = this->pu     + offset1;
		p0T     = this->p0     + offset1;
		pdT     = this->pd     + offset1;


		for (i0 = Bot0; i0 <= Top0; i0++)
		{
		    // compute 1-D 3 probabilities and transitions
		    puI = puT[i0];
		    p0I = p0T[i0];
		    pdI = pdT[i0];
		    l0 = i0 + shift1T[i0];

		    // offset arrays for next loop 
		    offset2 = NodeOffset(2, i0, 0, tpIdx);
		    shift2T = this->Shift2 + offset2;
		    quT     = this->qu     + offset2;
		    q0T     = this->q0     + offset2;
		    qdT     = this->qd     + offset2;

		    for (i1 = Bot1[i0]; i1 <= Top1[i0]; i1++)
		    {

			// compute 2-D 9 probabilities and transitions
			quuI = puI * quT[i1];
			qu0I = puI * q0T[i1];
			qudI = puI * qdT[i1];
			q0uI = p0I * quT[i1];
			q00I = p0I * q0T[i1];
			q0dI = p0I * qdT[i1];
			qduI = pdI * quT[i1];
			qd0I = pdI * q0T[i1];
			qddI = pdI * qdT[i1];
			l1 = i1 + shift2T[i1];

			// offset arrays for next loop 
			offset3 = NodeOffset(3, i0, i1, tpIdx);
			shift3T = this->Shift3 + offset3;
			ruT     = this->ru     + offset3;
			r0T     = this->r0     + offset3;
			rdT     = this->rd     + offset3;

			// Inner loop optimized: offset all needed slices
			// and store to tmp pointers
			discValT  = discVal  + offset3;
			sliceValT = sliceVal + offset3;

			nvuuT = newVal + NodeOffset(3, l0+1, l1+1, tpIdx+1);
			nvu0T = newVal + NodeOffset(3, l0+1, l1  , tpIdx+1);
			nvudT = newVal + NodeOffset(3, l0+1, l1-1, tpIdx+1);
			nv0uT = newVal + NodeOffset(3, l0  , l1+1, tpIdx+1);
			nv00T = newVal + NodeOffset(3, l0  , l1  , tpIdx+1);
			nv0dT = newVal + NodeOffset(3, l0  , l1-1, tpIdx+1);
			nvduT = newVal + NodeOffset(3, l0-1, l1+1, tpIdx+1);
			nvd0T = newVal + NodeOffset(3, l0-1, l1  , tpIdx+1);
			nvddT = newVal + NodeOffset(3, l0-1, l1-1, tpIdx+1);


			for (i2 = Bot2[i0][i1]; i2 <= Top2[i0][i1]; i2++) {
			    l2  = i2 + shift3T[i2];
			    l2u = l2 + 1;
			    l2d = l2 - 1;

			    if (discVal)
				x = sliceValT[i2] * discValT[i2];
			    else
				x = sliceValT[i2];


			    nvuuT[l2u] += quuI * ruT[i2] * x;
			    nvuuT[l2 ] += quuI * r0T[i2] * x;
			    nvuuT[l2d] += quuI * rdT[i2] * x;
			    nvu0T[l2u] += qu0I * ruT[i2] * x;
			    nvu0T[l2 ] += qu0I * r0T[i2] * x;
			    nvu0T[l2d] += qu0I * rdT[i2] * x;
			    nvudT[l2u] += qudI * ruT[i2] * x;
			    nvudT[l2 ] += qudI * r0T[i2] * x;
			    nvudT[l2d] += qudI * rdT[i2] * x;
			    nv0uT[l2u] += q0uI * ruT[i2] * x;
			    nv0uT[l2 ] += q0uI * r0T[i2] * x;
			    nv0uT[l2d] += q0uI * rdT[i2] * x;
			    nv00T[l2u] += q00I * ruT[i2] * x;
			    nv00T[l2 ] += q00I * r0T[i2] * x;
			    nv00T[l2d] += q00I * rdT[i2] * x;
			    nv0dT[l2u] += q0dI * ruT[i2] * x;
			    nv0dT[l2 ] += q0dI * r0T[i2] * x;
			    nv0dT[l2d] += q0dI * rdT[i2] * x;
			    nvduT[l2u] += qduI * ruT[i2] * x;
			    nvduT[l2 ] += qduI * r0T[i2] * x;
			    nvduT[l2d] += qduI * rdT[i2] * x;
			    nvd0T[l2u] += qd0I * ruT[i2] * x;
			    nvd0T[l2 ] += qd0I * r0T[i2] * x;
			    nvd0T[l2d] += qd0I * rdT[i2] * x;
			    nvddT[l2u] += qddI * ruT[i2] * x;
			    nvddT[l2 ] += qddI * r0T[i2] * x;
			    nvddT[l2d] += qddI * rdT[i2] * x;


			} // for (i2
		    } // for (i1
		} // for (i0

		//
		// Copy values back to output slice
		//
		for (i0 = OuBot0; i0 <= OuTop0; i0++)
		for (i1 = OuBot1[i0]; i1 <= OuTop1[i0]; i1++) {	
			offset3 = NodeOffset(3, i0, i1, tpIdx+1);
			sliceValT = sliceVal + offset3;
			newValT   = newVal   + offset3;

			for (i2 = OuBot2[i0][i1]; i2 <= OuTop2[i0][i1]; i2++)
			{
				sliceValT[i2] = newValT[i2];
			}
		}

	}

	if(debugLevel > DEBUG_LEVEL_DEV) {
		dppLog << "Slice AFTER being forwarded:" << endl;
        	slicePrint(
                   	  sliceVal,
                     	  sliceDim,
                   	  tpIdx+1,
                   	  FALSE,
                   	  dppLog);
	}

    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}


