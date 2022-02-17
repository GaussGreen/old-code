/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	C. Daher
 * Revision:	$Header$
 ***************************************************************/
#include <stdarg.h>
#include "kvtdeter.h"	// Prototype Consistency 
#include "kutilios.h"
#include "kstlutil.h"

extern	"C" {
#include "datelist.h"
#include "date_sup.h"

#include "drlio.h"		/* FScanStruct */
#include "drltime.h"
#include "drlts.h"
#include "vnfmanly.h"
};


//--------------------------------------------------------------

static	void	
ensureTsExists(KTSlice& ts, int tpIdx)
{
	ts.SetTpIdx(tpIdx); 

	if (ts.mData == NULL) 
		ts.mData = (void*) new double;
	if (ts.mData == NULL)
		throw KFailure("Failed to allocate slice data.\n");
   
}


//--------------------------------------------------------------

KVTreeDeter::KVTreeDeter()
{

	criticalDL = NULL;


	//
	criticalDL = GtoNewEmptyDateList(0);
	ASSERT_OR_THROW(criticalDL != NULL);

}


//--------------------------------------------------------------

KVTreeDeter::~KVTreeDeter()
{

	GtoFreeDateList(criticalDL);
	criticalDL = NULL;
}


//--------------------------------------------------------------
//

void
KVTreeDeter::Insert(TDate date)
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
KVTreeDeter::Insert(const KZeroReset& z, bool isCrit) // all critical for now
{
    try {
	Insert(z.mEarliestDate);
	//Insert(z.mMaturityDate);
    }
    catch (KFailure) {
	throw KFailure();
    }
    return;
}

//--------------------------------------------------------------
//


void
KVTreeDeter::Get(KTSlice& ts, const KZeroReset& z)
{
static	char	routine[] = "KVTreeDeter::Get ZeroReset";
	TDate	currentDate = TPDateCurrent();
	TDate	maturityDate = z.mMaturityDate;
const	String&	curveName = z.mCurveName; 

	ensureTsExists(ts, mTpIdx);
	*((double*) ts.mData) = FwdDisc(
		currentDate,
		maturityDate,
		curveName);

	return;
}


//--------------------------------------------------------------
//

TDate
KVTreeDeter::Insert(const KRateReset& r, bool isCrit) // all crit for now
{
	static	char	routine[] = "KVTreeDeter::InsertRateReset";
	
	try {
		Insert(r.EffDate());
		return (r.EffDate());
	}
	catch (KFailure) {
		throw KFailure("%s: failed inserting at %s.\n", 
				routine, 
				GtoFormatDate(r.EffDate()));
	}
	
}



//--------------------------------------------------------------


void
KVTreeDeter::Get(KTSlice& ts, const KRateReset& r)
{
static	char	routine[] = "KVTreeDeter::Get RateReset";
	TFloatRate      rate = (TFloatRate &)r;
	TDate		currentDate = TPDateCurrent();
	double		resetValue;
	long		curveIdx;

    try {

	ensureTsExists(ts, this->mTpIdx);

	if (!IS_ALMOST_ZERO(rate.weight)) {
	    if (r.ResetDate() <= currentDate) {
		curveIdx = GetCurveIdx(r.Rate().CurveName());
		resetValue = r.Rate().Forward(mZcCurves[curveIdx],
						r.ResetDate());
		*((double*) ts.mData) = resetValue;

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
KVTreeDeter::Update(int tpIdx)
{
static	char	routine[] = "KVTreeDeter::Update";

	if (mTpIdx == -1) {
		mTpIdx = tpIdx;
	} else {
		mTpIdx = tpIdx;
	}


	dppLog << format("=============  Updating TP index %4d (%10s)",
		tpIdx,
		(TPDateCurrent() > 1L ? 
		DrlTDatePrint(NULL, TPDateCurrent()) : "-"))
		<< endl;

}



//--------------------------------------------------------------
//

KTSlice&
KVTreeDeter::TSliceDev(KTSlice& ts, const String& discCurveName)
{
	int		tpIdx = mTpIdx;
	double		discFact;

	// Do not dev at last point
	//
	if (mTpIdx < TPNum()) {
	    if (ts.mData != NULL) {
		discFact = FwdDisc(
				TPDate(mTpIdx),
				TPDate(mTpIdx+1),
				discCurveName);
		*((double*) ts.mData) *= discFact;

		/*dppLog << format("KVTreeDeter: Dev %10s <- %10s: %12.10f",
			DrlTDatePrint(NULL, TPDate(mTpIdx)),
			DrlTDatePrint(NULL, TPDate(mTpIdx+1)),
			discFact) << endl;*/
	    }
	}

	return(ts);
}

//--------------------------------------------------------------
//

KTSlice&
KVTreeDeter::TSliceEv(KTSlice& ts)
{
	return(ts);
}


//--------------------------------------------------------------
// This constructor is empty since we want to allocate data
// in the tree only when needed
//

KTSlice&
KVTreeDeter::TSliceCreate(KTSlice& ts)
{
	return (ts);

}

//--------------------------------------------------------------
//

void
KVTreeDeter::TSliceDestroy(KTSlice& ts)
{
	delete ((double*) ts.mData);
	ts.mData = NULL;
}

//--------------------------------------------------------------
//

double
KVTreeDeter::TSliceGetCenter(KTSlice& ts)
{
	return *((double*) ts.mData);
}

//--------------------------------------------------------------
//

bool
KVTreeDeter::TSliceCompare(KTSlice& ts1, KTSlice& ts2, KTSComp compType)
{
	return(TRUE);
}

//--------------------------------------------------------------
//

KTSlice&
KVTreeDeter::TSliceScalarOper(KTSlice& ts, double value, KOper oper)
{
static	char	routine[] = "KVTreeDeter::TSliceScalarOper";
	double&	tsVal = *((double*) ts.mData);


	switch (oper) {
	case COPY:
		ensureTsExists(ts, mTpIdx);
		*((double*) ts.mData) = value;
		break;
	case ADD:
		tsVal += value;
		break;
	case SUB:
		tsVal -= value;
		break;
	case MULT:
		tsVal *= value;
		break;
	case DIV:
		tsVal /= value;
		break;

	case MAX:
		tsVal = MAX(value, tsVal);
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
KVTreeDeter::TSliceUnaryOper(KTSlice& ts, const KTSlice& ts1, KOper oper)
{
static	char	routine[] = "KVTreeDeter::TSliceUnaryOper";
	double&		tsVal  = *((double*) ts.mData);
	double&		tsVal1 = *((double*) ts1.mData);

	ts1.CheckNonEmpty("KVTreeDeter::TSliceUnaryOper");

	switch (oper) {
	case COPY:
		ensureTsExists(ts, mTpIdx);
		*((double*) ts.mData) = tsVal1;
		break;
	case ADD:
		tsVal  += tsVal1;
		break;
	case SUB:
		tsVal  -= tsVal1;
		break;
	case MULT:
		tsVal  *= tsVal1;
		break;

	case MAX:
		tsVal  = MAX(tsVal, tsVal1);
		break;

	case DIV:
		tsVal /= MAX(tsVal, tsVal1);
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
KVTreeDeter::TSlicePut(KTSlice& ts, ostream& os, int minMaxOnly)
{
	if (ts.mData) {
		double&	tsVal  = *((double*) ts.mData);
		os << format("[%15.8f]", tsVal);
	} else {
		os << "(null)";
	}
}

//--------------------------------------------------------------
//

void
KVTreeDeter::TSliceSpecialOper(KTSlice& ts, char* what, ...)
{
static	char	routine[] = "KVTreeDeter::TSliceSpecialOper";

	va_list	ap;
	va_start(ap, what);

	if (!strcmp(what, "EXTIMESET")) {
		//
		// Hack ! This should not be here !
		// Ask me for explanations !
		//

		KTSlice *undKs = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *optKs = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt0Ks = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt1Ks = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt2Ks = (KTSlice*) va_arg(ap, KTSlice*);
		double	te     = (double)   va_arg(ap, double);


		double& undVal = *((double*) undKs->mData);
		double& optVal = *((double*) optKs->mData);
		double& xt0Val = *((double*) xt0Ks->mData);
		double& xt1Val = *((double*) xt1Ks->mData);
		double& xt2Val = *((double*) xt2Ks->mData);


		    if (undVal > optVal) {
			xt0Val = 1e0;
			xt1Val = te;
			xt2Val = te*te;
		    }

	} else {
		throw KFailure("%s: operation `%s' not supported.\n",
			routine, what);
	}

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
KVTreeDeter::TSliceNodeBegin(KTSlice& ts, KTSliceNode& tsNode)
{
	tsNode.mIndex = (void*) 0x01;
	tsNode.mSlice = &ts;

	return (tsNode);
}


//--------------------------------------------------------------
// Test the end of of slice node 
//

bool
KVTreeDeter::TSliceNodeEnd(KTSliceNode& tsNode)
{
	return (tsNode.mIndex == NULL);
}


//--------------------------------------------------------------
// Move to next slice node 
//

KTSliceNode&
KVTreeDeter::TSliceNodeNext(KTSliceNode& tsNode)
{
	if (tsNode.mIndex != NULL)
		tsNode.mIndex = NULL;
	return (tsNode);
}

//--------------------------------------------------------------
// Get the value of slice node 
//

double&	
KVTreeDeter::TSliceAccessNodeValue(KTSlice& ts, KTSliceNode& tsNode)
{
	return *((double*) ts.mData);
}

//--------------------------------------------------------------
//
// Get the maximum difference of node values of slice 
// in nearest neighbors of specified node 

double
KVTreeDeter::TSliceGetNodeStepMax(KTSlice& ts, KTSliceNode& tsNode)
{
	return 0e0;
}

//--------------------------------------------------------------
// Prints a slice node on a stream.
//
void
KVTreeDeter::TSliceNodePut(KTSliceNode& tsNode, ostream& os)
{
	os << "KVTreeDeter::TSliceNodePut: N/A." << endl;

}



//**************************************************************
//
//	Time line operations
//
//**************************************************************

//--------------------------------------------------------------
//

int
KVTreeDeter::TPNum()
{
	return criticalDL->fNumItems-1;
}

//--------------------------------------------------------------
//

int
KVTreeDeter::TPIdxCurrent()
{
	return(mTpIdx);
}

//--------------------------------------------------------------
//

TDate
KVTreeDeter::TPDateCurrent()
{
	return TPDate(mTpIdx);
}


//--------------------------------------------------------------
//

double
KVTreeDeter::TPTimeCurrent()
{
	return ((double) (TPDate(mTpIdx) - mToday)) / 365e0;
}


//--------------------------------------------------------------
//

TDate
KVTreeDeter::TPDate(int tpIdx)
{
	if (tpIdx > criticalDL->fNumItems) {
		throw KFailure("tpIdx (%d) > criticalDL->fNumItems(%d)\n",
			tpIdx, criticalDL->fNumItems);
	}
	return criticalDL->fArray[tpIdx];
}


//--------------------------------------------------------------
//

int
KVTreeDeter::TPIdx(TDate date)
{
	/*long tpIdx;
	IF_FAILED_THROW(GtoTimeLineDateToTimePoint(vt.timeLine, date, &tpIdx));
	return (int)tpIdx;*/
	THROW_NA;
	return(-1);
}



//--------------------------------------------------------------
//

int
KVTreeDeter::TPIsCriticalDate(int tpIdx)
{
	return (TRUE);
}





//--------------------------------------------------------------
//

void
KVTreeDeter::Initialize(
	TDate todayDate,		// (I) 
	int   numZcCurves,	        // (I) 
	TCurve **zcCurves,              // (I) [numZcCurves] first = DISCOUNT 
	KVector(String) curveNames,     // (I) [numZcCurves]
	const String &valueCurveName)	// (I) reference disc curve
{
	mNumZcCurves = numZcCurves;
	mToday = todayDate;
	mZcCurves = zcCurves;
	mValueCurveName = valueCurveName;

	//
        // Init curve index map - MAX allowed num curves = 3
	//
	if(numZcCurves != curveNames.size()){
		throw KFailure("KVTreeDeter:: Num of curves (%d) "
			       "not equal num of curve names (%d).\n",
			       numZcCurves, curveNames.size());
	}
	MapCurveName(curveNames[0],GTO_CURVE_DISCOUNT);
	if(numZcCurves>1) MapCurveName(curveNames[1],GTO_CURVE_INDEX1);
	if(numZcCurves>2) MapCurveName(curveNames[2],GTO_CURVE_INDEX2);
	if(numZcCurves>3) MapCurveName(curveNames[3],GTO_CURVE_INDEX3);
	
	MapCurveName(K_DEFAULT_NAME,K_DEFAULT_IDX);  // defaults

}


//--------------------------------------------------------------
//

void
KVTreeDeter::Calibrate()
{
static	char	routine[] = "KVTreeDeter::Calibrate()";
	int	tpIdx;
	int	curveIdx;

    try {
	dppLog << routine << endl;

	Insert(mToday);
	dppLog << "\ttoday " << DrlTDatePrint(NULL, mToday) << endl;

	// Insert value date of reference curve (used to
	// pv the valur to).
	//
	for (curveIdx=0; curveIdx< mNumZcCurves ; curveIdx++) {
		Insert(mZcCurves[curveIdx]->fBaseDate);
	}
	curveIdx = GetCurveIdx(mValueCurveName);
	mValueDiscZeroShift = FwdDisc(
				mToday,
				mZcCurves[curveIdx]->fBaseDate,
				mValueCurveName);

	dppLog << "\treference curve " << mValueCurveName << " "
		<< DrlTDatePrint(NULL, mZcCurves[curveIdx]->fBaseDate)
		<< " shift " << mValueDiscZeroShift << endl;

	// Log
	dppLog << "\ttimeline:" << endl;
	for (tpIdx=0; tpIdx <= TPNum(); tpIdx++) {
		dppLog << format("\t\t%3d/%3d %10s",
			tpIdx, TPNum(),
			DrlTDatePrint(NULL, TPDate(tpIdx))) << endl;
	}


	// 
	mTpIdx = -1;

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}


//--------------------------------------------------------------
//

double
KVTreeDeter::FwdDisc(TDate date1, TDate date2, const String& curveName)
{
	int		curveIdx;
	TCurve		*zcCurve;
	double		z1, z2;

	curveIdx = GetCurveIdx(curveName);
	zcCurve = mZcCurves[curveIdx];

	IF_FAILED_THROW( GtoDiscountDate(
		date1,
		zcCurve,
		GTO_LINEAR_INTERP,
		&z1));

	IF_FAILED_THROW( GtoDiscountDate(
		date2,
		zcCurve,
		GTO_LINEAR_INTERP,
		&z2));

	return (z2/z1);
}



