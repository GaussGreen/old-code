/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	E. Ben-Artzi - C. Daher (from D. Gallager)
 * Revision:	$Header$
 ***************************************************************/
#include <stdarg.h>
#include "kvtalib.h"	// Prototype Consistency 
#include "kutilios.h"
#include "kstlutil.h"
#include "kmemutil.h"

extern	"C" {
#include "datelist.h"
#include "date_sup.h"

#include "modltype.h"
#include "ir1fact.h"
#include "ir1fini.h"                    /* Gto1FIRModelInit */
#ifdef UNIX
#include "ir2fact.h"                    /* Gto2FIRCalibInfoNewWrap */
#include "ir2fini.h"                    /* Gto2FIRModelInit */
#endif
#include "modlcons.h"                   /* GTO_MODEL_XXX constants */
#include "calibir.h"                    /* TCalibInfoIR  */
#include "parfixed.h"                    /* TCalibInfoIR  */



/*$$$
 * Using internals of ALIB 1F to access reference to node
 */
#include "f1timept.h"                   /* TTimePoint1F */
#include "f1vector.h"			/* To access timepoints */

typedef struct
{
    TTimePoint1F *timePoints1f;             /* 1-factor timepoints */
    TIndex        idx;                      /* Node index */
    long          tpIdx;                    /* TimePoint index */
} TTimeSlice1FIndex;                           

#include "f1.h"
static	void GtoTimeSlice1FPrintOS
    (OTreeContext treeContext,          /* (I) Tree context */
     long         tpIdx,                /* (I) TimePoint index */
     OTimeSlice   timeSlice,            /* (I) Timeslice to be printed out */
     char        *name,                 /* (I) String to be printed out  */
     TBoolean     minMaxOnly,           /* (I) Only print min,max */
     ostream& os);



#include "drlio.h"		/* FScanStruct */
#include "drltime.h"
#include "drlts.h"
#include "vnfmanly.h"


};


//--------------------------------------------------------------

inline	static	void	
ensureTsExists(KTSlice& ts, TVirtualTree vt, int tpIdx)
{
	ts.SetTpIdx(tpIdx); 

	if (ts.mData == NULL) 
		ts.mData = CALL(vt,tsNew)(vt.treeContext, tpIdx);
	if (ts.mData == NULL)
		throw KFailure("Failed to allocate slice data.\n");
   
}


//--------------------------------------------------------------

KVTreeAL::KVTreeAL()
{

	criticalDL = NULL;
	zeroDates = NULL;
	zeroBank = NULL;

	modelInitFunc = NULL;
	calibInfo = NULL;
	calibFreeFunc = NULL;
	tlInfo = NULL;
	mToday = 0L;
	criticalDL = NULL;

	vt.timeLine= NULL;
	vt.treeContext = NULL;


	mZcCurves = NULL;

	mBeta  = NULL;
	mAlpha = NULL;
	mRho   = NULL;

	mVolDates    = NULL;
	mVolMat      = NULL;
	mVolFreq     = NULL;
	mVolRates    = NULL;

	//
	criticalDL = GtoNewEmptyDateList(0);
	ASSERT_OR_THROW(criticalDL != NULL);



}


//--------------------------------------------------------------

KVTreeAL::~KVTreeAL()
{

	GtoFreeDateList(criticalDL);
	GtoZeroDatesFree(zeroDates);
	GtoZeroPriceFree(zeroBank);

	GtoTimeLineFree(vt.timeLine);
	if (vt.treeContext)
		CALL(vt, treeFree)(vt.treeContext);
	GtoTimeLineInfoFree(tlInfo);
	if (calibFreeFunc != NULL)
		(*calibFreeFunc)(calibInfo);



	criticalDL = NULL;
	zeroDates = NULL;
	zeroBank = NULL;

	modelInitFunc = NULL;
	calibInfo = NULL;
	tlInfo = NULL;
	mToday = 0L;
	criticalDL = NULL;


	if (mZcCurves) {
	    for (int curveIdx = 0; curveIdx < mNumZcCurves; curveIdx++) {
		GtoFreeTCurve(mZcCurves[curveIdx]);
	    }
	    delete[] mZcCurves;
	}

	delete [] mBeta;
	delete [] mAlpha;
	delete [] mRho;

	delete [] mVolDates;
	delete [] mVolMat;
	delete [] mVolFreq;
	delete [] mVolRates;


}


//--------------------------------------------------------------
//

void
KVTreeAL::Insert(TDate date)
{

	TDateList	*dl1 = NULL,
			*dl2 = NULL,
			*dl3 = NULL;

	dl1 = GtoNewDateListFromDates(&date, 1);
	ASSERT_OR_THROW(dl1 != NULL);

	dl2 = GtoMergeDateLists(dl1, criticalDL);
	ASSERT_OR_THROW(dl1 != NULL);

	GtoFreeDateList(criticalDL);
	GtoFreeDateList(dl1);
	criticalDL = dl2;
}


//--------------------------------------------------------------
//

void
KVTreeAL::Insert(const KZeroReset& z, bool isCrit) // all critical for now
{
	long		curveIdx; 
	
	TZeroDates	*newZeroDates = NULL;
	TDate		matDate, earlyDate;
	
	try {
		matDate   = z.mMaturityDate;
		earlyDate = z.mEarliestDate;

		curveIdx = GetCurveIdx(z.mCurveName);
		newZeroDates = GtoZeroDatesNew(
			&matDate,
			&earlyDate,
			&curveIdx,
			1);
		ASSERT_OR_THROW(newZeroDates != NULL);
		
		IF_FAILED_THROW( GtoZeroDatesAdd(
			newZeroDates,
			&this->zeroDates));

		GtoZeroDatesFree(newZeroDates);

	}

	catch (KFailure) {
		GtoZeroDatesFree(newZeroDates);
		
		throw KFailure("KVTreeAL::Insert: failed.\n");
	}
	
	return;
}

//--------------------------------------------------------------
//


void
KVTreeAL::Get(KTSlice& ts, const KZeroReset& z)
{
static	char	routine[] = "KVTreeAL::Get ZeroReset";
	OTimeSlice	tsz;
	TDate maturityDate = z.mMaturityDate;
	long		curveIdx = GetCurveIdx(z.mCurveName); 

	tsz = GtoZeroPriceGet(
		this->zeroBank,
		maturityDate,
		curveIdx); 
	if (tsz == NULL) {
		throw KFailure("%s: failed at date %s curve %s.\n",
			routine, DrlTDatePrint(NULL, maturityDate),
			GetCurveName(curveIdx).c_str());
	}

	ensureTsExists(ts, this->vt, this->mTpIdx);

	CALL(this->vt, tsCopy)(
		this->vt.treeContext,
		mTpIdx,
   		tsz,
		(OTimeSlice) ts.mData);

	return;

}


//--------------------------------------------------------------
//

TDate
KVTreeAL::Insert(const KRateReset& r, bool isCrit) // all crit for now
{
	static	char	routine[] = "KVTreeAL::InsertRateReset";
	
	TDateList	startDL;
	TZeroDates	*newZeroDates;
	TFloatRateArray	fltArray;
	TFloatRate      rate = (TFloatRate &)r;
	long		curveIdx;

	TDate		resetDate = r.ResetDate();
	
	try {
		
		if (!IS_ALMOST_ZERO(rate.weight)) {
			curveIdx = GetCurveIdx(r.Rate().CurveName());

			//
			// TRUE floating rate
			//
			startDL.fNumItems = 1;
			startDL.fArray = &resetDate;
			
			fltArray.defs = &rate;
			fltArray.curveIndices = &curveIdx;
			fltArray.numRates = 1L;
			fltArray.rateInfo = NULL;
			fltArray.spread = 0e0;
			fltArray.stubRateType = GTO_STUB_RATE_INDEX;
			
			
			/*GtoFloatRateArrayPrint(&fltArray, routine);
			  GtoPrintDateList(&startDL, routine);
			  if (this->zeroDates)
			  GtoZeroDatesPrint(this->zeroDates, routine);*/
			
			newZeroDates = GtoZeroDatesBuild1(
				this->zeroDates,
				&startDL,
				&fltArray,
				NULL);
			ASSERT_OR_THROW(newZeroDates != NULL);
			
			
			// Replace old one
			GtoZeroDatesFree(this->zeroDates);
			this->zeroDates = newZeroDates;
			
			/*if (this->zeroDates)
			  GtoZeroDatesPrint(this->zeroDates, routine);*/
			
		} else {
			//
			// No weight = fixed rate.
			//
			Insert(r.EffDate());
		}
		
		
		
		return (r.EffDate());
	}
	catch (KFailure) {
		throw KFailure("%s: failed inserting at %s.\n",
				routine,
				GtoFormatDate(r.EffDate()));
	}
	
}




//--------------------------------------------------------------
//

TDate
KVTreeAL::Insert(const KRateReset& r, TDate endDate, bool isCrit) 
{
	static	char	routine[] = "KVTreeAL::InsertRateResetAmerican";

	TDate	lastZeroDate;

	try {
		Insert(r, isCrit);

		lastZeroDate = endDate + r.Rate().SpotOffset() 
			     + r.Rate().Maturity();

                Insert( KZeroReset(r.Rate().CurveName(),
                                   r.ResetDate(),
                                   lastZeroDate),
                        isCrit);

		return (r.EffDate());

	}
	catch (KFailure) {
		throw KFailure("%s: failed inserting between date %s and %s\n", 
				 routine,
			  	 GtoFormatDate(r.EffDate()),
			  	 GtoFormatDate(endDate));
	}
	
}




//--------------------------------------------------------------


void
KVTreeAL::Get(KTSlice& ts, const KRateReset& r)
{
static	char	routine[] = "KVTreeAL::Get RateReset";
	TFloatRate      rate = (TFloatRate &)r;
	TDate		currentDate = TPDateCurrent();
	double		resetFwd, currentFwd;

    try {

	ensureTsExists(ts, this->vt, this->mTpIdx);

	if (!IS_ALMOST_ZERO(rate.weight)) {
	    if (r.ResetDate() == currentDate) {
		//
		// We are at the reset date.
		//
		long	curveIdx = GetCurveIdx(r.Rate().CurveName());

		IF_FAILED_THROW (GtoTreeRateSingle(
			&this->vt,
			this->zeroBank, 
			mTpIdx,
			r.ResetDate(),
			&rate,
			curveIdx,
			(OTimeSlice) ts.mData));

	    } else if (r.ResetDate() < currentDate) {
		//
		// Requested reset is in the past: need to adjust by
		// ratios of forwards.
		//
		long	curveIdx = GetCurveIdx(r.Rate().CurveName());

		IF_FAILED_THROW (GtoTreeRateSingle(
			&this->vt,
			this->zeroBank, 
			mTpIdx,
			currentDate,
			&rate,
			curveIdx,
			(OTimeSlice) ts.mData));

		// Compute deterministic forwards at current and
		// reset dates.

		resetFwd = r.Forward(mZcCurves[curveIdx]);

		currentFwd = r.Forward(mZcCurves[curveIdx]);

		ts *= (resetFwd / currentFwd);


	    } else {
		throw KFailure();
	    }
	} else {
		ts = rate.spread;
	}

	return;
    } catch (KFailure) {
	throw KFailure("%s: failed at date %s on reset:\n",
		routine, DrlTDatePrint(NULL, currentDate));
    }
}



//--------------------------------------------------------------
//

void
KVTreeAL::Update(int tpIdx)
{
static	char	routine[] = "KVTreeAL::Update";

	if (mTpIdx == -1) {
		mTpIdx = tpIdx;
	} else {
		mTpIdx = tpIdx;
	}

	dppLog << format("Updating TP index %4d (%10s)",
		tpIdx,
		(TPDateCurrent() > 1L ? 
		DrlTDatePrint(NULL, TPDateCurrent()) : "-")) << endl;

	if (GtoZeroPriceUpdate(zeroBank, mTpIdx) != SUCCESS) {
		throw KFailure("%s: failed.\n", routine);
	}

}



//--------------------------------------------------------------
//

KTSlice&
KVTreeAL::TSliceDev(KTSlice& ts, const String& discCurveName)
{
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	int		tpIdx = mTpIdx;

	int	discCurveIdx;

	if (ts.mData != NULL) {

		/*int	lastTpIdx = GtoTimeLineLastTP(vt.timeLine);
		GtoErrMsg("tpIdx = %d  (lastTP %d)\n", tpIdx, lastTpIdx);
	    	CALL(vt,tsPrint)(
		vt.treeContext,
		tpIdx,
		ots,
		"before dev",
		TRUE);
		*/


	    discCurveIdx = GetCurveIdx(discCurveName);

	    if (CALL(vt,dev)(
		vt.treeContext,
		vt.timeLine, 
		tpIdx,
		discCurveIdx,
		ots) IS FAILURE)
			throw KFailure();

	    /*CALL(vt,tsPrint)(
		vt.treeContext,
		tpIdx,
		ots,
		"after  dev",
		TRUE);*/




	}

	return(ts);
}

//--------------------------------------------------------------
//

KTSlice&
KVTreeAL::TSliceEv(KTSlice& ts)
{
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	int		tpIdx = mTpIdx;

	if (ts.mData != NULL) {
	    if (CALL(vt,ev)(
		vt.treeContext,
		vt.timeLine, 
		tpIdx,
		ots) IS FAILURE)
			throw KFailure();
	}
	return(ts);

}


//--------------------------------------------------------------
// This constructor is empty since we want to allocate data
// in the tree only when needed
//

KTSlice&
KVTreeAL::TSliceCreate(KTSlice& ts)
{
	return (ts);

}

//--------------------------------------------------------------
//

void
KVTreeAL::TSliceDestroy(KTSlice& ts)
{
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	CALL(vt,tsFree)(ots, vt.treeContext);
}

//--------------------------------------------------------------
//

double
KVTreeAL::TSliceGetCenter(KTSlice& ts)
{
	OTimeSlice	ots = (OTimeSlice) ts.mData;

	ts.CheckNonEmpty("KVTreeAL::TSliceGetCenter");

	return CALL(vt,tsIndexGetValueToday)(ots);
}

//--------------------------------------------------------------
//

bool
KVTreeAL::TSliceCompare(KTSlice& ts1, KTSlice& ts2, KTSComp compType)
{
	return(TRUE);
}

//--------------------------------------------------------------
//

KTSlice&
KVTreeAL::TSliceScalarOper(KTSlice& ts, double value, KOper oper)
{
static	char	routine[] = "KVTreeAL::TSliceScalarOper";
	OTimeSlice	ots = (OTimeSlice) ts.mData;


	switch (oper) {
	case COPY:
		ensureTsExists(ts, this->vt, this->mTpIdx);
		ots = (OTimeSlice) ts.mData;
		CALL(vt,tsSet)(
			vt.treeContext,
			mTpIdx,
			value,
			ots);
		break;
	case ADD:
		ts.CheckNonEmpty("KVTreeAL::TSliceScalarOper(ADD)");
		CALL(vt,tsAddConstant)(
			vt.treeContext,
			mTpIdx,
			value,
			ots,
			ots);
		break;
	case SUB:
		ts.CheckNonEmpty("KVTreeAL::TSliceScalarOper(SUB)");
		CALL(vt,tsAddConstant)(
			vt.treeContext,
			mTpIdx,
			(-value),
			ots,
			ots);
		break;
	case MULT:
		ts.CheckNonEmpty("KVTreeAL::TSliceScalarOper(MULT)");
		CALL(vt,tsMultConstant)(
			vt.treeContext,
			mTpIdx,
			value,
			ots,
			ots);
		break;
	case DIV:
		ts.CheckNonEmpty("KVTreeAL::TSliceScalarOper(DIV)");
		CALL(vt,tsMultConstant)(
			vt.treeContext,
			mTpIdx,
			(1e0/value),
			ots,
			ots);
		break;

	case MAX:
		ts.CheckNonEmpty("KVTreeAL::TSliceScalarOper(MAX)");
		CALL(vt,tsMax)(
			vt.treeContext,
			mTpIdx,
			TRUE,
			0e0,
			ots,
			ots);
		break;


	case MIN:
	case GEQ:
	case LEQ:
	default:
		throw KFailure("%s: operation not implemented.\n", routine);
	}

	return(ts);
}

//--------------------------------------------------------------
//

KTSlice&
KVTreeAL::TSliceUnaryOper(KTSlice& ts, const KTSlice& ts1, KOper oper)
{
static	char	routine[] = "KVTreeAL::TSliceUnaryOper";
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	OTimeSlice	ots1 = (OTimeSlice) ts1.mData;

	ts1.CheckNonEmpty("KVTreeAL::TSliceUnaryOper");

	switch (oper) {
	case COPY:
		ensureTsExists(ts, this->vt, this->mTpIdx);
		ots = (OTimeSlice) ts.mData;
		CALL(vt,tsCopy)(
			vt.treeContext,
			mTpIdx,
			ots1,
			ots);
		break;
	case ADD:
		ts.CheckNonEmpty("KVTreeAL::TSliceUnaryOper(ADD)");
		CALL(vt,tsWeightedSum)(
			vt.treeContext,
			mTpIdx,
			1e0,
			ots,
			1e0,
			ots1,
			ots);
		break;
	case SUB:
		ts.CheckNonEmpty("KVTreeAL::TSliceUnaryOper(SUB)");
		CALL(vt,tsWeightedSum)(
			vt.treeContext,
			mTpIdx,
			1e0,
			ots,
			-1e0,
			ots1,
			ots);
		break;
	case MULT:
		ts.CheckNonEmpty("KVTreeAL::TSliceUnaryOper(MULT)");
		CALL(vt,tsProduct)(
			vt.treeContext,
			mTpIdx,
			ots,
			ots1,
			ots);
		break;

	case MAX:
		ts.CheckNonEmpty("KVTreeAL::TSliceUnaryOper(MAX)");
		CALL(vt,tsOption)(
			vt.treeContext,
			mTpIdx,
			TRUE,
			0e0,
			ots,
			ots1,
			ots);
		break;

	case DIV: {
		double	tsVal, ts1Val;
		OTSIndex nodeIdx;
		nodeIdx = CALL(vt, tsIndexNew)(vt.treeContext, mTpIdx);
		ASSERT_OR_THROW(nodeIdx != NULL);

		do
		{
		    tsVal  = CALL(vt,tsIndexGetValue)(nodeIdx, ots);
		    ts1Val = CALL(vt,tsIndexGetValue)(nodeIdx, ots1);
		    tsVal /= ts1Val;
		    CALL(vt,tsIndexSetValue)(nodeIdx, ots, tsVal);
		} while (CALL(vt,tsIndexNext) (nodeIdx));
		}
		break;

	case MIN:
	case GEQ:
	case LEQ:
	default:
		throw KFailure("%s: operation not implemented.\n", routine);
	}

	return(ts);
}


//--------------------------------------------------------------
//

void	
KVTreeAL::TSlicePut(KTSlice& ts, ostream& os, int minMaxOnly)
{
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	long		tpIdx = mTpIdx;

	GtoTimeSlice1FPrintOS(
		vt.treeContext,
     		tpIdx,			/* TimePoint index */
		ots,
		"",
		minMaxOnly, 			/* Only print min,max */
		os);
}

//--------------------------------------------------------------
//

void
KVTreeAL::TSliceSpecialOper(KTSlice& ts, char* what, ...)
{
static	char	routine[] = "KVTreeAL::TSliceSpecialOper";

	va_list	ap;
	va_start(ap, what);

	if (!strcmp(what, "EXTIMESET")) {
		//
		// Hack ! This should not be here !
		// Ask me for explanations !
		//

		double	undVal, optVal;

		KTSlice *undKs = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *optKs = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt0Ks = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt1Ks = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt2Ks = (KTSlice*) va_arg(ap, KTSlice*);
		double	te     = (double)   va_arg(ap, double);


		OTimeSlice undTs = (OTimeSlice) undKs->mData;
		OTimeSlice optTs = (OTimeSlice) optKs->mData;
		OTimeSlice xt0Ts = (OTimeSlice) xt0Ks->mData;
		OTimeSlice xt1Ts = (OTimeSlice) xt1Ks->mData;
		OTimeSlice xt2Ts = (OTimeSlice) xt2Ks->mData;


		OTSIndex nodeIdx;

		nodeIdx = CALL(vt, tsIndexNew)(vt.treeContext, mTpIdx);
		ASSERT_OR_THROW(nodeIdx != NULL);

		do
		{
		    undVal = CALL(vt,tsIndexGetValue)(nodeIdx, undTs);
		    optVal = CALL(vt,tsIndexGetValue)(nodeIdx, optTs);

		    if (undVal > optVal) {
			CALL(vt,tsIndexSetValue)(nodeIdx, xt0Ts, 1e0);
			CALL(vt,tsIndexSetValue)(nodeIdx, xt1Ts, te);
			CALL(vt,tsIndexSetValue)(nodeIdx, xt2Ts, te*te);
		    }
		} while (CALL(vt,tsIndexNext) (nodeIdx));

	} else if (!strcmp(what, "TSPRT_MINMAX")) {
		//
		//
		//

		double	val,
			minVal =  1e30,
			maxVal = -1e30;

		ostream	*os = (ostream*) va_arg(ap, ostream*);

		OTimeSlice	oTs = (OTimeSlice) ts.mData;
		OTSIndex	nodeIdx;

		nodeIdx = CALL(vt, tsIndexNew)(vt.treeContext, mTpIdx);
		ASSERT_OR_THROW(nodeIdx != NULL);

		do
		{
		    val = CALL(vt,tsIndexGetValue)(nodeIdx, oTs);
	 	    minVal = MIN(minVal, val);
	 	    maxVal = MAX(maxVal, val);

		} while (CALL(vt,tsIndexNext) (nodeIdx));


		(*os) << format("min: %lf  max:%lf",
			minVal, maxVal);


	} else 
		throw KFailure("%s: operation `%s' not supported.\n",
			routine, what);

	va_end(ap);

    try {
    }
    catch(...) {
	va_end(ap);
	throw KFailure("%s: failed.\n", routine);
    }
}


//**************************************************************
//
//	Node operations
//
//**************************************************************


//--------------------------------------------------------------
// Create and set the begining of slice node 
//

KTSliceNode&
KVTreeAL::TSliceNodeBegin(KTSlice& ts, KTSliceNode& tsNode)
{
	OTSIndex nodeIdx;
	nodeIdx = CALL(vt, tsIndexNew)(vt.treeContext, mTpIdx);
	ASSERT_OR_THROW(nodeIdx != NULL);

	tsNode.mIndex = (void*) nodeIdx;
	tsNode.mSlice = &ts;

	return (tsNode);
}


//--------------------------------------------------------------
// Test the end of of slice node 
//

bool
KVTreeAL::TSliceNodeEnd(KTSliceNode& tsNode)
{
	OTSIndex nodeIdx = (OTSIndex) tsNode.mIndex;
	return (nodeIdx == NULL);
}


//--------------------------------------------------------------
// Move to next slice node 
//

KTSliceNode&
KVTreeAL::TSliceNodeNext(KTSliceNode& tsNode)
{
	OTSIndex nodeIdx = (OTSIndex) tsNode.mIndex;
	if (CALL(vt,tsIndexNext) (nodeIdx) == FALSE)
		tsNode.mIndex = NULL;
	return (tsNode);
}

//--------------------------------------------------------------
// Get the value of slice node 
//

double&	
KVTreeAL::TSliceAccessNodeValue(KTSlice& ts, KTSliceNode& tsNode)
{

	OTSIndex	nodeIdx = (OTSIndex) tsNode.mIndex;
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	TTimeSlice1FIndex *idxp = (TTimeSlice1FIndex *)nodeIdx;
	double 		*tsDouble;

	tsDouble = TS_TO_REF_PTR(ots);
	return (tsDouble[idxp->idx]);
}

//--------------------------------------------------------------
//
// Get the maximum difference of node values of slice 
// in nearest neighbors of specified node 

double
KVTreeAL::TSliceGetNodeStepMax(KTSlice& ts, KTSliceNode& tsNode)
{
	OTSIndex	nodeIdx = (OTSIndex) tsNode.mIndex;
	OTimeSlice	ots = (OTimeSlice) ts.mData;
	return CALL(vt,tsIndexGetMaxDiff)(nodeIdx, ots);
}

//--------------------------------------------------------------
// Prints a slice node on a stream.
//
void
KVTreeAL::TSliceNodePut(KTSliceNode& tsNode, ostream& os)
{
	OTSIndex	nodeIdx = (OTSIndex) tsNode.mIndex;
	os << "KVTreeAL::TSliceNodePut: N/A." << endl;

}



//**************************************************************
//
//	Time line operations
//
//**************************************************************

//--------------------------------------------------------------
//

int
KVTreeAL::TPNum()
{
	return GtoTimeLineLastTP(vt.timeLine);
}

//--------------------------------------------------------------
//

int
KVTreeAL::TPIdxCurrent()
{
	return(mTpIdx);
}

//--------------------------------------------------------------
//

TDate
KVTreeAL::TPDateCurrent()
{
	return GtoTPGetDate(vt.timeLine, mTpIdx);
}


//--------------------------------------------------------------
//

double
KVTreeAL::TPTimeCurrent()
{
	return GtoTPGetYearsSinceValue(vt.timeLine, mTpIdx);
}


//--------------------------------------------------------------
//

TDate
KVTreeAL::TPDate(int tpIdx)
{
	return (int)GtoTPGetClosestDate(vt.timeLine, tpIdx);
}


//--------------------------------------------------------------
//

int
KVTreeAL::TPIdx(TDate date)
{
	long tpIdx;

	IF_FAILED_THROW(GtoTimeLineDateToTimePoint(vt.timeLine, date, &tpIdx));

	return (int)tpIdx;
}



//--------------------------------------------------------------
//

int
KVTreeAL::TPIsCriticalDate(int tpIdx)
{
	
	return(TP_IS_CRITICAL_DATE(vt.timeLine, tpIdx));
}





//--------------------------------------------------------------
//

void
KVTreeAL::Initialize(
	TDate todayDate,		// (I) 
	int   numZcCurves,	        // (I) 
	TCurve **zcCurves,		// (I) [numZcCurves] first = DISCOUNT 
	KVector(String) curveNames,     // (I) [numZcCurves]
	const String& diffuseCurve,     // (I) one of the above names
	int numFact,			// (I) number of factors
	const double *beta,		// (I) model param
	const double *alpha,		// (I) model param
	const double *rho,		// (I) model param
	int ppy,			// (I) ppy
	double smoothFact,		// (I) 
	double numStdevCut,		// (I) 
	double q1,			// (I) 
	double q2,			// (I) 
	double fwdSh,			// (I) 
					/*     VOLATILITY DATA */
	int numVolDates,		/* (I) # of volatility dates */
	const TDate *volDates,		/* (I) array of vol dates */
	const double *volMat,		/* (I) array of vol fwd maturities */
	const int *volFreq,		/* (I) array of vol frequency */
	const double *volRates)		/* (I) array of vol */
{
static	char	routine[] = "KVTreeAL::Initialize";

	mToday = todayDate;

	//
        // Init curve index map - MAX allowed num curves = 3
	//
	if(numZcCurves != curveNames.size()){
		throw KFailure("KVTreeAL:: Num of curves (%d) "
			       "not equal num of curve names (%d).\n",
			       numZcCurves, curveNames.size());
	}
	mZcCurves = new TCurve*[numZcCurves];
	ASSERT_OR_THROW(mZcCurves != NULL)
	for (int curveIdx = 0; curveIdx < numZcCurves; curveIdx++) {
		mZcCurves[curveIdx] = GtoCopyCurve(zcCurves[curveIdx]);
		ASSERT_OR_THROW(mZcCurves[curveIdx] != NULL);
	}
	mNumZcCurves = numZcCurves;



	MapCurveName(curveNames[0],GTO_CURVE_DISCOUNT);
	if(numZcCurves>1) MapCurveName(curveNames[1],GTO_CURVE_INDEX1);
	if(numZcCurves>2) MapCurveName(curveNames[2],GTO_CURVE_INDEX2);
	if(numZcCurves>3) MapCurveName(curveNames[3],GTO_CURVE_INDEX3);
	
	mDiffuseIdx = GetCurveIdx(diffuseCurve);
	if(mDiffuseIdx == K_DEFAULT_IDX) {
		throw KFailure("KVTreeAL::Initialize: Unknown diffuse "
				"curve name `%s'.\n", diffuseCurve.c_str());
	}
	MapCurveName(K_DEFAULT_NAME,K_DEFAULT_IDX);  // defaults
	
	mNumFact = numFact;
	mBeta    = DppNewArrayCopy(numFact, beta);
	mAlpha   = DppNewArrayCopy(numFact, alpha);
	mRho     = DppNewArrayCopy(numFact*(numFact-1)/2, rho);
	mPpy     = ppy;
	mSmoothFact = smoothFact;
	mNumStdevCut = numStdevCut;
	mQ1 = q1;
	mQ2 = q2;

	mNumVolDates = numVolDates;
	mVolDates    = DppNewArrayCopy(numVolDates, volDates);
	mVolMat      = DppNewArrayCopy(numVolDates, volMat);
	mVolFreq     = DppNewArrayCopy(numVolDates, volFreq);
	mVolRates    = DppNewArrayCopy(numVolDates, volRates);
}






//--------------------------------------------------------------
//

void
KVTreeAL::Calibrate()
{
static	char	routine[] = "KVTreeAL::Calibrate()";


	VnfmData	*tfData = NULL;
	int		ppy;
	TDateInterval	bvMat;
	double		vnfmDistType = -1e0;
	int		distType;	// 0=N, 1=LN
	TCurve		*bvCurve = NULL;

    try {


	//$$$ Insert at least 3M !!!
	// 
	// Otherwise, weird stuff happening:
	// The tree will  NOT RUN WITH LESS THAN TRHREE NODES
	{
		KDateInterval	itv = KDateInterval(2e0/((double) mPpy));
		TDate endDate = mZcCurves[0]->fBaseDate + itv;
		Insert(KZeroReset(GetCurveName(GTO_CURVE_DISCOUNT),
				  mZcCurves[0]->fBaseDate,endDate));
	}




	//
	// (1) Initialize tlInfo (ppy).
	//
	dppLog << "Starting " << routine << endl;
	dppLog << "  Today date: " << GtoFormatDate(mToday) << endl; \
	dppLog << "  Ppy: " << mPpy << endl;
	{ int idx; \
	    dppLog << format( "   /     volDates  volMat   volFreq  vol\n");\
	    for (idx=0; idx<=mNumVolDates-1; idx++) { \
		dppLog << format( "%3d/%3d  %10s %8.4f %7d  %8.6f", \
		    idx, mNumVolDates, GtoFormatDate(mVolDates[idx]), \
		    mVolMat[idx], mVolFreq[idx], mVolRates[idx]) << endl;}}


	if (IS_ALMOST_ZERO(mQ1-1e0) && IS_ALMOST_ZERO(mQ1-1e0)) {
		vnfmDistType = VNFM_NORMAL_DIST;
		distType = 0;
	} else
	if (IS_ALMOST_ZERO(mQ1-0e0) && IS_ALMOST_ZERO(mQ1-0e0)) {
		vnfmDistType = VNFM_LOGNORMAL_DIST;
		distType = 1;
	} else {
		throw KFailure("%s: Q1=%lf Q2=%lf not supported.\n",
			routine, mQ1, mQ2);
	}


	// Check supported.
	if ((mNumFact != 1) &&
	    (mNumFact != 2))
		throw KFailure("%s: mNumFact=%d not supported.\n",
			routine, mNumFact);


	ppy = mPpy;
	this->tlInfo = GtoTimeLineInfoNew(ppy, (TDate)0, 0, FALSE);
	ASSERT_OR_THROW(this->tlInfo != NULL);

	//
	// (2) Calibrate swaption volatility 
	//

	// implied base vol maturity 
	IF_FAILED_THROW(GtoFreq2TDateInterval((long) 4, &bvMat));

	/* create tf calib data */
	tfData = VnfmNew2FactSimple(
		mToday,
		mNumVolDates,
		mVolDates,
		NULL,
		mZcCurves[mDiffuseIdx], 
		vnfmDistType,
		mBeta[0],
		(mNumFact != 2 ? 1e0 : mBeta[1]),
		1e0,
		(mNumFact != 2 ? 0e0 : mAlpha[1] / mAlpha[0]),
		0.20e0,
		0.20e0,
		(mNumFact != 2 ? 1e0 : mRho[0]));
	ASSERT_OR_THROW(tfData != NULL)

	IF_FAILED_THROW(VnfmComputeCoeff(tfData));


	/* Call bootstraping routine.
	 * If valu date has been added, we may need
	 * to offset by 1 the input vectors
	 *
	 */
	if (tfData->fNDates == mNumVolDates) {
		IF_FAILED_THROW(VnfmVolCalib1VArbitrary(
			tfData,
			0, tfData->fNDates-1,
			mVolMat,
			mVolFreq,
			mVolRates,
			LOGVOL,
			TRUE,
			NULL));
	} else if (tfData->fNDates == (mNumVolDates+1)) {
		IF_FAILED_THROW(VnfmVolCalib1VArbitrary(
			tfData,
			0, tfData->fNDates-1,
			mVolMat-1,
			mVolFreq-1,
			mVolRates-1,
			LOGVOL,
			TRUE,
			NULL));
	} else {
		THROW_BUG;
	}




	IF_FAILED_THROW(VnfmComputeCoeff(tfData));

	// generate base vol curve
	IF_FAILED_THROW(VnfmGenerateVolTCurve(
		tfData,
		bvMat,
		0,
		0,
		mToday,
		mNumVolDates,
		mVolDates,
		&bvCurve));

#ifdef	__DEBUG__
	dppLog << format("Done Calibrating Vnfm:\n", routine);
	VnfmFpWrite(tfData, dppFpLog);
	DrlTCurveFpWrite(bvCurve, dppFpLog, 9999L);
#endif

	//
	// Initialize model type
	//
/*#define	USE_SPOT_VOLS*/
#if !defined(USE_SPOT_VOLS)


	if (distType == 1) {
	    /*
	     * Log Normal distribution
	     */
	    if (mNumFact == 1) {
		this->calibInfo = Gto1FIRCalibInfoNew(
			mBeta[0],	// Mean reversion coeff
			mSmoothFact,	// (0=none; 1=norm)
			mNumStdevCut,	// Num stdev 
			0e0,		   // max spotvol
			(long)mDiffuseIdx, // ZC diffusion index
#ifndef	CLIB7X
			bvCurve,
#endif
#ifndef	ALIB87
			1e0,		/* Process power */
#endif
			GTO_CONST_SPOT_VOL_INTERP,
			GTO_NO_STATE_PRICES);
		ASSERT_OR_THROW(this->calibInfo != NULL);

	    	this->modelInitFunc = Gto1FIRModelInit;
	    	this->calibFreeFunc = Gto1FIRCalibInfoFree;

	    } else {
	    	throw KFailure("%s: Unknown model type.\n", routine);
	    }
	} else if (distType == 0) {
#else
	if ((distType == 0) || (distType == 1)) {
#endif


	    /*
	     * Normal distribution
	     */
	    if (mNumFact == 1) {
		TCalibInfoIR	*infop = NULL;	/* To be returned */
		TVolDefIR	*volDef = NULL;	/* Volatility definition */
#if !defined(USE_SPOT_VOLS)
		double		processPower = 0e0;	/* normal */
#else
		double		processPower = (double) distType;
#endif
		double		maxSpotVols = 0e0;

		TDate		todayDate,
				spvolDates[256];
		double		spvolValue1[256];
		double		*spvolValue = &spvolValue1[0];
		int		numSpvol, idxS;

		ASSERT_OR_THROW(tfData->fNDates < 256);

		/*
		 * Convert to Gto vol array conventions
		 */
		todayDate = tfData->fDate[0];
		numSpvol = tfData->fNDates;
		for (idxS=0; idxS<=numSpvol-2; idxS++) {
			spvolDates[idxS] = tfData->fDate[idxS+1];
			spvolValue[idxS] = tfData->fSigma[0][idxS] *
						tfData->fAlpha[0];
		}
		spvolDates[numSpvol-1] = spvolDates[numSpvol-2] + 3650L;
		spvolValue[numSpvol-1] = tfData->fSigma[0][numSpvol-1] *
						tfData->fAlpha[0];

#ifdef	__DEBUG__
#endif
		dppLog << format( "%s: normal spot vol: %15s\n",
			routine, GtoFormatDate(todayDate));
		for (idxS=0; idxS<=numSpvol-1; idxS++) {
			dppLog << format( "\t%2d\t%12s\t%lf\n",
				idxS,
				GtoFormatDate(spvolDates[idxS]),
				spvolValue[idxS]);
		}

		/* Use spot vol: Set up TVolDefIR using spotVols */
		volDef = GtoVolDefIRNew(
			&mBeta[0],	// array of MR
			1L,		// num fact 
			todayDate,
			spvolDates,
			&spvolValue,
			NULL,
			numSpvol,
#ifndef	ALIB87
			processPower,	/* Process power */
#endif
			GTO_CONST_SPOT_VOL_INTERP);
		ASSERT_OR_THROW(volDef != NULL);
		GtoVolDefIRPrint(volDef, routine);



		/* Create a new Calibration structure using the TVolDefIR.  */
		infop =  GtoCalibInfoIRNewFromVolDef(
			volDef, 
#ifdef	ALIB87
			processPower,	   /* Process power */
#endif
			&maxSpotVols,
			mSmoothFact,	   /* (0=none; 1=norm)*/
			mNumStdevCut,	   /* Num stdev to cut tree */
			(long)mDiffuseIdx, /* ZC diffusion index  */
			GTO_NO_STATE_PRICES);
		ASSERT_OR_THROW(infop != NULL);
		GtoCalibInfoIRPrint(infop, routine);


		this->calibInfo = (OCalibInfo) infop;
	    	this->modelInitFunc = Gto1FIRModelInit;
	    	//this->calibFreeFunc = Gto1FIRCalibInfoFree;

		GtoVolDefIRFree(volDef);


	    } else {
	    	throw KFailure("%s: Unknown model type.\n", routine);
	    }
	} else {
	    throw KFailure("%s: unknown irDistType (%d).\n",
			routine, distType);
	}

	dppLog << "Done creating OCalibInfo" << endl;
	dppLog << "Calibrating VirtualTree" << endl;

	//
	// Initialize the model. Note that an empty zeroBank is allocated.
	//

#ifdef	__DEBUG__
	GtoErrMsg("%s: critical dates.\n", routine);
	GtoPrintDateList(this->criticalDL, routine);
	GtoErrMsg("%s: zero dates.\n", routine);
	if (this->zeroDates)
		GtoZeroDatesPrint(this->zeroDates, routine);
#endif


	if ((*(this->modelInitFunc))(
		this->criticalDL,
		this->zeroDates,
		mZcCurves[0],                     // discount curve
		mNumZcCurves==1?NULL:&mZcCurves[1], // index curves
		(long)(mNumZcCurves-1),           // num index curves
		this->tlInfo,
		this->calibInfo,
		routine,
		&this->vt,
		&this->zeroBank) IS FAILURE) {
			throw KFailure();
	}


	// Turn off
	//
	GtoZeroPriceNoWarnings(&this->zeroBank);


    
#ifdef DEBUG_TIME_LINE
	dppLog << "Generated timeline: " << endl;
	CALL(vt,tlPrint)(vt.treeContext, vt.timeLine);
#endif

	dppLog << "Done calibrating VirtualTree" << endl;

	mTpIdx = -1;

	/* Free memory */
	VnfmFree(tfData);
        GtoFreeTCurve(bvCurve);

    }
    catch (KFailure) {
	VnfmFree(tfData);
        GtoFreeTCurve(bvCurve);
	throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
//
/*
 * GtoTimeSlice1FPrint():  Prints out a 1-factor timeslice. If minMaxOnly
 * is set to TRUE, only prints min and max values.
 */
void GtoTimeSlice1FPrintOS
    (OTreeContext treeContext,          /* (I) Tree context */
     long         tpIdx,                /* (I) TimePoint index */
     OTimeSlice   timeSlice,            /* (I) Timeslice to be printed out */
     char        *name,                 /* (I) String to be printed out  */
     TBoolean     minMaxOnly,           /* (I) Only print min,max */
     ostream& os)
{
    TOneFactContext *cxt = (TOneFactContext *)treeContext;
    TTimePoint1F *timePoints1f = cxt->timePoints1f;
    double *ts;
    TIndex idx;                /* node index */

    int retFlag = 0;            /* Return flag */
    double minVal = DBL_MAX;
    double maxVal = -DBL_MAX;

    TIndex top, bottom;          /* Bounds of tree */

    if (minMaxOnly) {
	//os << "(" ;
    } else {
        os << format("\n\nTimeSlice: %20.20s\n", name);
    }

    if (timeSlice ISNT NULL)
    {
        ts = TS_TO_REF_PTR(timeSlice);
        
        top = GtoTP1FGetTop(timePoints1f, tpIdx);
        bottom = GtoTP1FGetBottom(timePoints1f, tpIdx);
        
        for (idx=bottom; idx <= top; idx++)          
        {
            
            if (minMaxOnly)
            {
                minVal = MIN(minVal, ts[idx]);
                maxVal = MAX(maxVal, ts[idx]);
            }
            else
            {
                os << format("(%2d)=%6.6f ", idx, ts[idx]);
                retFlag = ++retFlag & 3;
                if (retFlag IS 0)
                    os << endl;
            }
        }                           /* for idx */
        if (minMaxOnly)
            os << format("[%15.12f(%d), %15.12f(%d)]",
			minVal, bottom, maxVal, top);
        else
            os << endl;
    }
    else                        /* NULL pointer received. */
    {
        os << "\tNULL" << endl;
    }
}











