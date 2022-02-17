/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "vtltest.h"
#include "kvtspar.h"	// Slice parser

#include "kutilios.h"


extern "C" {
#include "drltime.h"
#include "drloptio.h"
};



const	int	nRMax=120;
const	int	nKMax=120;


//---------------------------------------------------------------

void KVTreeTestForwardRates1(
	KVTree &vt,
	KVector(TDate) &resetDates,
	KVector(TDate) &payDates,
	KVector(double) &strikes,
	KRate &floatRate,
	const String &curveName,
	const KZCurve &zcCurve,
	const String &discCurveName,
	const KZCurve &discZcCurve,
	ostream &os)
{
static	char	routine[] = "KVTreeTestForwardRates1";
	int	nR, idxR,
		nK, idxK;


	KVector(KRateReset*)	mRateResets;
	KVector(KZeroReset*)	mZeroResets;
	KVector(TDate)		mModelResetDates;
	KVector(KTSlice*)	mValue;


	TDate			currDate;
	TDate			todayDate = vt.TPToday();

	KTSlice			*mRateTs = NULL,
				*mZeroTs = NULL;
	KVector(KTSlice*)	mVM0;		// 0 moment (disc)
	KVector(KTSlice*)	mVM1;		// 1 moment
	KVector(KTSlice*)	mVM2;		// 2 moment
	KVector(KTSlice*)	mVL1;		// 1 moment of log
	KVector(KTSlice*)	mVL2;		// 2 moment of log

	KVector(KVector(KTSlice*))	mVK;

	double	mval[nRMax][6],
		kstk[nRMax][nKMax],	// strikes
		kval[nRMax][nKMax];


#define	DEV(sliceP)	{if (sliceP != NULL) {sliceP->Dev(discCurveName);}}



    try {


	ASSERT_OR_THROW(resetDates.size() == payDates.size());
	nR = resetDates.size();
	ASSERT_OR_THROW(nR < nRMax);

	nK = strikes.size();
	ASSERT_OR_THROW(nK < nKMax);


	// enough space to store the model reset dates and resets.
	mRateResets.resize(nR);
	mZeroResets.resize(nR);
	mModelResetDates.resize(nR);
	mVM0.resize(nR);
	mVM1.resize(nR);
	mVM2.resize(nR);
	mVL1.resize(nR);
	mVL2.resize(nR);
	mVK.resize(nR);
	for (idxR=0; idxR<nR; idxR++) {
		mVM0[idxR] = NULL;
		mVM1[idxR] = NULL;
		mVM2[idxR] = NULL;
		mVL1[idxR] = NULL;
		mVL2[idxR] = NULL;
		mVK[idxR].resize(nK);
		for (idxK=0; idxK<nK; idxK++) {
			mVK[idxR][idxK] = NULL;
		}
	}


	//
	// Add events in tree
	//
	for (idxR=0; idxR<nR; idxR++) {
		//
		// Add events to tree.
		// Rate reset: we collect the model reset date
		// (i.e. where the rate is known in the tree)
		// Pay date: we need the zero from pay to reset.
		// Create (and store) rate and zero resets
		// for each coupon payment
		//
		//
		if (resetDates[idxR] <  todayDate)
			throw KFailure("%s: resetDate %d %s < todayDate %s.\n",
				routine, idxR+1,
				DrlTDatePrint(NULL, resetDates[idxR]),
				DrlTDatePrint(NULL, todayDate));

		if (payDates[idxR] < resetDates[idxR])
			throw KFailure("%s: payDate %d %s < resetDate %s.\n",
				routine, idxR+1,
				DrlTDatePrint(NULL, payDates[idxR]),
				DrlTDatePrint(NULL, resetDates[idxR]));



		mRateResets[idxR] = new KRateReset(
			curveName,
			resetDates[idxR],
			resetDates[idxR], // !!! effective dates!
			floatRate);

		mModelResetDates[idxR] = vt.Insert(*mRateResets[idxR]);


		mZeroResets[idxR] = new KZeroReset(
			discCurveName,
			mModelResetDates[idxR],
			payDates[idxR]);

		vt.Insert(resetDates[idxR]);
		vt.Insert(*mZeroResets[idxR]);

	}


	//
	// Compute strikes
	//
	for (idxR=0; idxR<nR; idxR++) {
		double	ef;
		ef = floatRate.Forward(
			zcCurve,
			mRateResets[idxR]->EffDate());
		for (idxK=0; idxK<nK; idxK++)  {
			kstk[idxR][idxK] = ef + strikes[idxK];
		}
	}


	//
	// Initialize tree timeline and calibrate tree
	//
	vt.SetUpTimeline();
	vt.Calibrate();

	//
	// Rollback
	//
	for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) {
	    //
	    // Update Tree
	    //
	    vt.Update(tpIdx);


	    currDate = vt.TPDateCurrent();

	    //
	    // Dev
	    //
	    for (idxR=0; idxR<nR; idxR++) {
		DEV(mVM0[idxR]);
		DEV(mVM1[idxR]);
		DEV(mVM2[idxR]);
		DEV(mVL1[idxR]);
		DEV(mVL2[idxR]);
		for (idxK=0; idxK<nK; idxK++) {
			DEV(mVK[idxR][idxK]);
		}
	    }


	    //
	    // if Reset date: the floating rate can be computed.
	    //
	    for (idxR=0; idxR<nR; idxR++) {
	      if (currDate == mModelResetDates[idxR]) {

		//
		// Ensure slices allocated
		//
		mVM0[idxR] = new KTSlice(vt, "mM0", curveName);
		mVM1[idxR] = new KTSlice(vt, "mM1", curveName);
		mVM2[idxR] = new KTSlice(vt, "mM2", curveName);
		mVL1[idxR] = new KTSlice(vt, "mL1", curveName);
		mVL2[idxR] = new KTSlice(vt, "mL2", curveName);
		for (idxK=0; idxK<nK; idxK++) {
			mVK[idxR][idxK] = 
				new KTSlice(vt, "mVK", curveName);
		}

		if (mRateTs == NULL) {
			mRateTs = new KTSlice(vt, "mRateTs", curveName);
			mZeroTs = new KTSlice(vt, "mZeroTs", discCurveName);
		}

		// Get rate reset
		//
		vt.Get(*mRateTs, *mRateResets[idxR]);


		// if paid in future, get zero
		//
		if (payDates[idxR] != mModelResetDates[idxR]) {
			vt.Get(*mZeroTs, *mZeroResets[idxR]);
		} else {
			*mZeroTs = 1e0;
		}

		//
		// Set values
		//
		*mVM0[idxR] = *mZeroTs;
		*mVM1[idxR] = *mZeroTs;
		*mVM2[idxR] = *mZeroTs;
		*mVL1[idxR] = *mZeroTs;
		*mVL2[idxR] = *mZeroTs;

		*mVM1[idxR] *= *mRateTs;
		*mVM2[idxR] *= *mRateTs;
		*mVM2[idxR] *= *mRateTs;

		// calc option
		for (idxK=0; idxK<nK; idxK++)  {
			*mVK[idxR][idxK]  = *mRateTs;
			*mVK[idxR][idxK] -= kstk[idxR][idxK];
			(*mVK[idxR][idxK]).max(0e0);
			*mVK[idxR][idxK] *= *mZeroTs;
		}




	      }
	    }



	    //
	    // Last TP: store the result.
	    //
	    if (vt.TPIdxCurrent() == 0) {
	      for (idxR=0; idxR<nR; idxR++) {

		mval[idxR][0] = mVM0[idxR]->GetCenter();
		mval[idxR][1] = mVM1[idxR]->GetCenter();
		mval[idxR][2] = mVM2[idxR]->GetCenter();
		for (idxK=0; idxK<nK; idxK++)  {
			kval[idxR][idxK] = (*mVK[idxR][idxK]).GetCenter();
		}
	      }
	    }
	}


	//
	//
	//
	os << format("TODAY: %10s\n", DrlTDatePrint(NULL, vt.TPToday()));

	//
	// Expected value and var
	//



	os << "FORWARD RATES:" << endl;
	os << " RESET DATE  TEXP  PAY DATE     TPAY ";
	os << "| ZERO(F)     ZERO(TREE)  DIFF      ";
	os << "|   RATE(F)     RATE(TREE)    DIFF      ";
	os << "VOL";
	os << endl;
	for (idxR=0; idxR<nR; idxR++) {
		double	te,	// time to exer
			tp,	// time to pay
			zf,	// disc fact
			ef,	// fwd value
			zv,	// disc fact
			ev,	// expected value
			sv;	// stdev
		KDayCc dcc = KDayCc(GTO_ACT_365F);


		zv = mval[idxR][0];
		ev = mval[idxR][1] / zv;
		sv = mval[idxR][2] / zv;


		// calc fwd and disc fact
		ef = floatRate.Forward(
			zcCurve,
			mRateResets[idxR]->EffDate());

		zf = discZcCurve.DiscFact(payDates[idxR]);

		te = DayCountFraction(
			vt.TPToday(),
			(TDate)resetDates[idxR],
			dcc);
		tp = DayCountFraction(
			vt.TPToday(),
			(TDate)payDates[idxR],
			dcc);


		os << format(" %10s", DrlTDatePrint(NULL, resetDates[idxR]));
		os << format(" %5.2f", te);
		os << format(" %10s", DrlTDatePrint(NULL, payDates[idxR]));
		os << format(" %5.2f", tp);

		os << format(" |");
		os << format(" %10.8f", zf);
		os << format(" %10.8f", zv);
		os << format(" %11.8f", zv-zf);

		os << format(" |");
		os << format(" %12.8f", ef*1e2);
		os << format(" %12.8f", ev*1e2);
		os << format(" %11.8f", (ev-ef)*1e2);

		os << format(" |");
		os << format(" %7.4f", sqrt((sv-ev*ev)/te)/ev * 1e2);




		for (idxK=0; idxK<nK; idxK++)  {

			kval[idxR][idxK] = (*mVK[idxR][idxK]).GetCenter();
		}

		os << endl;
	}



	//
	// Options smile
	//
	os << format("IMPLIED BLACK_SCHOLES VOLATILITY SMILE:\n");
	os << format("        STRIKE OFFSET (%%)\n");
	os << format(" TEXP ");
	for (idxK=0; idxK<nK; idxK++) 
		os << format("  %5.2f  ", strikes[idxK] * 1e2);
	os << endl;
	for (idxR=0; idxR<nR; idxR++) {
		double	te,	// time to exer
			tp,	// time to pay
			zf,	// disc fact
			ef,	// fwd value
			zv,	// disc fact
			ev,	// expected value
			sv,	// stdev
			iv;	// implied bs vol

		KDayCc dcc = KDayCc(GTO_ACT_365F);

		zv = mval[idxR][0];
		ev = mval[idxR][1] / zv;
		sv = mval[idxR][2] / zv;


		// calc fwd and disc fact
		ef = floatRate.Forward(
			zcCurve,
			mRateResets[idxR]->EffDate());

		zf = discZcCurve.DiscFact(payDates[idxR]);

		te = DayCountFraction(
			vt.TPToday(),
			(TDate)resetDates[idxR],
			dcc);
		tp = DayCountFraction(
			vt.TPToday(),
			(TDate)payDates[idxR],
			dcc);


		os << format(" %5.2f", te);
		for (idxK=0; idxK<nK; idxK++)  {
			if (!IS_ALMOST_ZERO(te)) {
				//
				// Turn off error msg for impl vol.
				//
				int errStat = GtoErrMsgSet(FALSE);

				if (DrlBlackImplVol(
					te,
					ef,
					kval[idxR][idxK] / zv,
					kstk[idxR][idxK],
					"C",
					"P",
					&iv) != SUCCESS)
						iv = -0.01;

				GtoErrMsgSet(errStat);
			} else {
				iv = 0e0;
			}

			os << format("  %7.4f", iv * 1e2);

		}

		os << endl;
	}









	//
	// Free memory
	//
	delete mRateTs;
	delete mZeroTs;
	for (idxR=0; idxR<nR; idxR++) {
		delete mVM0[idxR];
		delete mVM1[idxR];
		delete mVM2[idxR];
		delete mVL1[idxR];
		delete mVL2[idxR];
		for (idxK=0; idxK<nK; idxK++) {
			delete mVK[idxR][idxK];
		}
	}
	for (idxR=0; idxR<nR; idxR++) {
		delete mRateResets[idxR];
		delete mZeroResets[idxR];
	}


    }
    catch (KFailure) {
	delete mRateTs;
	delete mZeroTs;
	for (idxR=0; idxR<nR; idxR++) {
		delete mVM0[idxR];
		delete mVM1[idxR];
		delete mVM2[idxR];
		delete mVL1[idxR];
		delete mVL2[idxR];
		for (idxK=0; idxK<nK; idxK++) {
			delete mVK[idxR][idxK];
		}
	}
	for (idxR=0; idxR<nR; idxR++) {
		delete mRateResets[idxR];
		delete mZeroResets[idxR];
	}

	throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------



void KVTreeTestCorrRates(
	KVTree &vt,			// (I) 
	KVector(TDate) &resetDates,	// (I)
	KVector(TDate) &payDates,	// (I)
	KVector(double) &strikes,	// (I)	offset 

	KRate &floatRate0,		// (I)
	const String &curveName0,	// (I)
	const KZCurve &zcCurve0,	// (I)

	KRate &floatRate1,		// (I)
	const String &curveName1,	// (I)
	const KZCurve &zcCurve1,	// (I)

	const String &discCurveName,	// (I)
	const KZCurve &discZcCurve,	// (I)
	ostream &os)			// (I)
{
static	char	routine[] = "KVTreeTestCorrRates";
	int	nR, idxR,
		nK, idxK;


	KVector(KRateReset*)	mRateResets[2];
	KVector(KZeroReset*)	mZeroResets[2];
	KVector(KZeroReset*)	mPZeroResets;
	KVector(KZeroReset*)	mEZeroResets[2];
	KVector(TDate)		mModelResetDates[2];
	KVector(TDate)		mCorrResetDates;

	KVector(KTSlice*)	mValue;


	TDate			currDate;
	TDate			todayDate = vt.TPToday();



	KVector(KTSlice*)	mRate[2];
	KTSlice			*mZero = NULL;
	KTSlice			*mZero0 = NULL;
	KTSlice			*mZero1 = NULL;

	KVector(KTSlice*)	mVM0;		// 0-moment
	KVector(KTSlice*)	mVE0;		// 1-moments
	KVector(KTSlice*)	mVE1;
	KVector(KTSlice*)	mVC00;	// 2-moments
	KVector(KTSlice*)	mVC01;
	KVector(KTSlice*)	mVC11;

	KVector(KVector(KTSlice*))	mVOpt;



	double			pVM0[nRMax],
				pVE0[nRMax],
				pVE1[nRMax],
				pVC00[nRMax],
				pVC01[nRMax],
				pVC11[nRMax],
				pVOpt[nRMax][nKMax],
				kstk[nRMax][nKMax];


#define	DEV(sliceP)	{if (sliceP != NULL) {sliceP->Dev(discCurveName);}}



    try {


	ASSERT_OR_THROW(resetDates.size() == payDates.size());
	nR = resetDates.size();
	ASSERT_OR_THROW(nR < nRMax);

	nK = strikes.size();
	ASSERT_OR_THROW(nK < nKMax);


	// enough space to store the model reset dates and resets.
	mRateResets[0].resize(nR);
	mZeroResets[0].resize(nR);
	mModelResetDates[0].resize(nR);
	mRateResets[1].resize(nR);
	mZeroResets[1].resize(nR);
	mModelResetDates[1].resize(nR);

	mPZeroResets.resize(nR);
	mEZeroResets[0].resize(nR);
	mEZeroResets[1].resize(nR);

	mCorrResetDates.resize(nR);

	//
	mRate[0].resize(nR);
	mRate[1].resize(nR);

	mVM0.resize(nR);
	mVE0.resize(nR);
	mVE1.resize(nR);
	mVC00.resize(nR);
	mVC01.resize(nR);
	mVC11.resize(nR);
	mVOpt.resize(nR);

	for (idxR=0; idxR<nR; idxR++) {
		mRate[0][idxR] = NULL;
		mRate[1][idxR] = NULL;

		mVM0[idxR] = NULL;
		mVE0[idxR] = NULL;
		mVE1[idxR] = NULL;
		mVC00[idxR] = NULL;
		mVC01[idxR] = NULL;
		mVC11[idxR] = NULL;

		mVOpt[idxR].resize(nK);
		for (idxK=0; idxK<nK; idxK++) {
			mVOpt[idxR][idxK] = NULL;
		}
	}


	//
	// Add events in tree
	//
	for (idxR=0; idxR<nR; idxR++) {
		//
		// Add events to tree.
		// Rate reset: we collect the model reset date
		// (i.e. where the rate is known in the tree)
		// Pay date: we need the zero from pay to reset.
		// Create (and store) rate and zero resets
		// for each coupon payment
		//
		//
		if (resetDates[idxR] <  todayDate)
			throw KFailure("%s: resetDate %d %s < todayDate %s.\n",
				routine, idxR+1,
				DrlTDatePrint(NULL, resetDates[idxR]),
				DrlTDatePrint(NULL, todayDate));

		if (payDates[idxR] < resetDates[idxR])
			throw KFailure("%s: payDate %d %s < resetDate %s.\n",
				routine, idxR+1,
				DrlTDatePrint(NULL, payDates[idxR]),
				DrlTDatePrint(NULL, resetDates[idxR]));

		// rate 1

		mRateResets[0][idxR] = new KRateReset(
			curveName0,
			resetDates[idxR],
			resetDates[idxR], // !!! effective dates!
			floatRate0);
		mModelResetDates[0][idxR] =
				vt.Insert(*mRateResets[0][idxR]);
		mZeroResets[0][idxR] = new KZeroReset(
			discCurveName,
			mModelResetDates[0][idxR],
			payDates[idxR]);

		vt.Insert(resetDates[idxR]);
		vt.Insert(*mZeroResets[0][idxR]);


		// rate 1

		mRateResets[1][idxR] = new KRateReset(
			curveName1,
			resetDates[idxR],
			resetDates[idxR], // !!! effective dates!
			floatRate1);
		mModelResetDates[1][idxR] =
				vt.Insert(*mRateResets[1][idxR]);
		mZeroResets[1][idxR] = new KZeroReset(
			discCurveName,
			mModelResetDates[1][idxR],
			payDates[idxR]);

		vt.Insert(resetDates[idxR]);
		vt.Insert(*mZeroResets[1][idxR]);


		mCorrResetDates[idxR] = MIN(
			mModelResetDates[0][idxR],
			mModelResetDates[1][idxR]);

		// Payment
		mEZeroResets[0][idxR] = new KZeroReset(
			discCurveName,
			mCorrResetDates[idxR],
			mModelResetDates[0][idxR]);
		vt.Insert(*mEZeroResets[0][idxR]);

		mEZeroResets[1][idxR] = new KZeroReset(
			discCurveName,
			mCorrResetDates[idxR],
			mModelResetDates[1][idxR]);
		vt.Insert(*mEZeroResets[1][idxR]);

		mPZeroResets[idxR] = new KZeroReset(
			discCurveName,
			mCorrResetDates[idxR],
			payDates[idxR]);
		vt.Insert(*mPZeroResets[idxR]);



	}


	//
	// Compute strikes
	//
	for (idxR=0; idxR<nR; idxR++) 
	for (idxK=0; idxK<nK; idxK++) {
		kstk[idxR][idxK] = 
			  mRateResets[0][idxR]->Forward(zcCurve0)
			- mRateResets[1][idxR]->Forward(zcCurve1)
			+ strikes[idxK];
	}


	//
	// Calibrate tree
	//
	vt.Calibrate();

	//
	// Rollback
	//
	for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) {
	    //
	    // Update Tree
	    //
	    vt.Update(tpIdx);


	    currDate = vt.TPDateCurrent();

	    //
	    // Dev
	    //
	    for (idxR=0; idxR<nR; idxR++) {
		DEV(mRate[0][idxR]);
		DEV(mRate[1][idxR]);

		DEV(mVM0[idxR]);
		DEV(mVE1[idxR]);
		DEV(mVE0[idxR]);
		DEV(mVC00[idxR]);
		DEV(mVC01[idxR]);
		DEV(mVC11[idxR]);

		for (idxK=0; idxK<nK; idxK++) {
			DEV(mVOpt[idxR][idxK]);
		}
	    }



	    //
	    // if Reset date: the floating rate can be computed.
	    //
	    for (idxR=0; idxR<nR; idxR++) {

	      // 
	      // Index 0
	      // 
	      if (currDate == mModelResetDates[0][idxR]) {
		if (mZero == NULL) {
			mZero  = new KTSlice(vt, "mZero", discCurveName);
			mZero0 = new KTSlice(vt, "mZero", discCurveName);
			mZero1 = new KTSlice(vt, "mZero", discCurveName);
		}
		mRate[0][idxR] = new KTSlice(vt, "mRate0", curveName0);
		vt.Get(*(mRate[0][idxR]), *mRateResets[0][idxR]);

		*mVE0[idxR]   = *mRate[0][idxR];
		*mVC00[idxR]  = *mRate[0][idxR];
		*mVC00[idxR] *= *mRate[0][idxR];

		if (payDates[idxR] != mModelResetDates[0][idxR]) {
			vt.Get(*mZero, *mZeroResets[0][idxR]);
			//*(mRate[0][idxR]) *= *mZero;
			*mVE0[idxR]  *= *mZero;
			*mVC00[idxR] *= *mZero;
		}


	      }

	      // 
	      // Index 1
	      // 
	      if (currDate == mModelResetDates[1][idxR]) {
		mRate[1][idxR] = new KTSlice(vt, "mRate1", curveName1);
		vt.Get(*(mRate[1][idxR]), *mRateResets[1][idxR]);

		*mVE1[idxR]   = *mRate[1][idxR];
		*mVC11[idxR]  = *mRate[1][idxR];
		*mVC11[idxR] *= *mRate[1][idxR];
		if (payDates[idxR] != mModelResetDates[1][idxR]) {
			vt.Get(*mZero, *mZeroResets[1][idxR]);
			//*(mRate[1][idxR]) *= *mZero;
			*mVE1[idxR]  *= *mZero;
			*mVC11[idxR] *= *mZero;
		}


	      }

	      //
	      // Payoff
	      //
	      if (currDate == mCorrResetDates[idxR]) {

		mVM0[idxR]  = new KTSlice(vt, "mVM0", curveName0);
		mVE0[idxR]  = new KTSlice(vt, "mVE0", curveName0);
		mVE1[idxR]  = new KTSlice(vt, "mVE1", curveName1);
		mVC00[idxR] = new KTSlice(vt, "mVC00", curveName0);
		mVC01[idxR] = new KTSlice(vt, "mVC01", curveName1);
		mVC11[idxR] = new KTSlice(vt, "mVC11", curveName1);
		for (idxK=0; idxK<nK; idxK++) {
			mVOpt[idxR][idxK] = new KTSlice(vt, "mVK", curveName1);
		}


		if (payDates[idxR] != mCorrResetDates[idxR]) {
			vt.Get(*mZero, *mPZeroResets[idxR]);
		} else {
			*mZero = 1e0;
		}

		if (mCorrResetDates[idxR] != mModelResetDates[0][idxR]) {
			vt.Get(*mZero0, *mEZeroResets[0][idxR]);
		} else {
			*mZero0 = 1e0;
		}


		if (mCorrResetDates[idxR] != mModelResetDates[1][idxR]) {
			vt.Get(*mZero1, *mEZeroResets[1][idxR]);
		} else {
			*mZero1 = 1e0;
		}


		// Adjust for pay
		*mRate[0][idxR] /= *mZero0;

		*mRate[1][idxR] /= *mZero1;

		*mVM0[idxR]   = *mZero;

		*mVC01[idxR]  = *mRate[0][idxR];
		*mVC01[idxR] *= *mRate[1][idxR];
		*mVC01[idxR] *= *mZero;


		// calc option
		for (idxK=0; idxK<nK; idxK++)  {
			*mVOpt[idxR][idxK]  = *mRate[0][idxR];
			*mVOpt[idxR][idxK] -= *mRate[1][idxR];
			*mVOpt[idxR][idxK] -= kstk[idxR][idxK];
			(*mVOpt[idxR][idxK]).max(0e0);
		}






	      }

	    }



	    //
	    // Last TP: store the result.
	    //
	    if (vt.TPIdxCurrent() == 0) {
	      for (idxR=0; idxR<nR; idxR++) {

		pVM0[idxR]  = mVM0[idxR]->GetCenter();
		pVE0[idxR]  = mVE0[idxR]->GetCenter();
		pVE1[idxR]  = mVE1[idxR]->GetCenter();

		pVC00[idxR] = mVC00[idxR]->GetCenter();
		pVC01[idxR] = mVC01[idxR]->GetCenter();
		pVC11[idxR] = mVC11[idxR]->GetCenter();

		for (idxK=0; idxK<nK; idxK++)  {
			pVOpt[idxR][idxK] = (*mVOpt[idxR][idxK]).GetCenter();
		}
	      }
	    }
	}


	//
	//
	//
	os << format("TODAY: %10s\n", DrlTDatePrint(NULL, vt.TPToday()));

	//
	// Expected value and var
	//



	os << "FORWARD RATES:" << endl;
	os << " RESET DATE  TEXP  PAY DATE     TPAY ";
	os << "| ZERO(F)     ZERO(TREE)  DIFF      ";
	os << "|   RATE(F)     RATE(TREE)    DIFF      ";
	os << "VOL";
	os << endl;
	for (idxR=0; idxR<nR; idxR++) {
		double	tr,		// time to reset
			te,		// time to model reset
			tp,		// time to pay
			zf,		// disc fact
			ef0,		// fwd value
			ef1,		// fwd value
			zv,		// disc fact
			ev0,		// expected value
			ev1,		// expected value
			sv0, sv1, cr,	// stdev and corr
			c00, c01, c11;	// covariance

		KDayCc dcc = KDayCc(GTO_ACT_365F);

		zv  = pVM0[idxR];
		ev0 = pVE0[idxR] / zv;
		ev1 = pVE1[idxR] / zv;

		c00 = pVC00[idxR] / zv;
		c01 = pVC01[idxR] / zv;
		c11 = pVC11[idxR] / zv;

		// calc fwd and disc fact
		ef0 = mRateResets[0][idxR]->Forward(zcCurve0);
		ef1 = mRateResets[1][idxR]->Forward(zcCurve1);

		zf = discZcCurve.DiscFact(payDates[idxR]);

		tr = DayCountFraction(
			vt.TPToday(),
			(TDate)resetDates[idxR],
			dcc);
		te = DayCountFraction(
			vt.TPToday(),
			(TDate)mCorrResetDates[idxR],
			dcc);
		tp = DayCountFraction(
			vt.TPToday(),
			(TDate)payDates[idxR],
			dcc);


		os << format(" %10s", DrlTDatePrint(NULL, resetDates[idxR]));
		os << format(" %5.2f", tr);
		//os << format(" %10s", DrlTDatePrint(NULL, payDates[idxR]));
		os << format(" %5.2f", tp);
		os << format(" %5.2f", te);

		os << format(" |");
		os << format(" %10.8f", zv);
		os << format(" %11.2e", zv-zf);

		os << format(" |");
		os << format(" %8.4f", ev0*1e2);
		os << format(" %11.2e", (ev0-ef0)*1e2);
		os << format(" %8.4f", ev1*1e2);
		os << format(" %11.2e", (ev1-ef1)*1e2);

		sv0 = sqrt((c00-ev0*ev0)/te) / ev0;
		sv1 = sqrt((c11-ev1*ev1)/te) / ev1;
		cr  = (c01 - ev0*ev1) / 
			sqrt((c00-ev0*ev0)*(c11-ev1*ev1));

		os << format(" |");
		os << format(" %8.4f", sv0 * 1e2);
		os << format(" %8.4f", sv1 * 1e2);
		os << format(" %8.4f", cr);


		os << endl;


	}


#ifdef	_SKIP
	//
	// Options smile
	//
	os << format("IMPLIED BLACK_SCHOLES VOLATILITY SMILE:\n");
	os << format("        STRIKE OFFSET (%%)\n");
	os << format(" TEXP ");
	for (idxK=0; idxK<nK; idxK++) 
		os << format("  %5.2f  ", strikes[idxK] * 1e2);
	os << endl;
	for (idxR=0; idxR<nR; idxR++) {
		double	te,	// time to exer
			tp,	// time to pay
			zf,	// disc fact
			ef,	// fwd value
			zv,	// disc fact
			ev,	// expected value
			sv,	// stdev
			iv;	// implied bs vol


		zv = mval[idxR][0];
		ev = mval[idxR][1] / zv;
		sv = mval[idxR][2] / zv;


		// calc fwd and disc fact
		ef = floatRate.Forward(
			zcCurve,
			mRateResets[idxR]->mEffDate);

		zf = discZcCurve.DiscFact(payDates[idxR]);

		te = DayCountFraction(
			vt.TPToday(),
			resetDates[idxR],
			KDayCc(GTO_ACT_365F));
		tp = DayCountFraction(
			vt.TPToday(),
			payDates[idxR],
			KDayCc(GTO_ACT_365F));


		os << format(" %5.2f", te);
		for (idxK=0; idxK<nK; idxK++)  {
			if (!IS_ALMOST_ZERO(te)) {
				//
				// Turn off error msg for impl vol.
				//
				int errStat = GtoErrMsgSet(FALSE);

				if (DrlBlackImplVol(
					te,
					ef,
					kval[idxR][idxK] / zv,
					kstk[idxR][idxK],
					"C",
					"P",
					&iv) != SUCCESS)
						iv = -0.01;

				GtoErrMsgSet(errStat);
			} else {
				iv = 0e0;
			}

			os << format("  %7.4f", iv * 1e2);

		}

		os << endl;
	}
#endif








	//
	// Free memory
	//
	delete mZero;
	for (idxR=0; idxR<nR; idxR++) {
		delete mRate[0][idxR];
		delete mRate[1][idxR];

		delete mVM0[idxR];
		delete mVE1[idxR];
		delete mVE0[idxR];
		delete mVC00[idxR];
		delete mVC01[idxR];
		delete mVC11[idxR];

		for (idxK=0; idxK<nK; idxK++) {
			delete mVOpt[idxR][idxK];
		}
	}
	for (idxR=0; idxR<nR; idxR++) {
		delete mRateResets[0][idxR];
		delete mZeroResets[0][idxR];
		delete mRateResets[1][idxR];
		delete mZeroResets[1][idxR];
		delete mPZeroResets[idxR];
	}
 


    }
    catch (KFailure) {
	delete mZero;
	for (idxR=0; idxR<nR; idxR++) {
		delete mRate[0][idxR];
		delete mRate[1][idxR];

		delete mVM0[idxR];
		delete mVE1[idxR];
		delete mVE0[idxR];
		delete mVC00[idxR];
		delete mVC01[idxR];
		delete mVC11[idxR];

		for (idxK=0; idxK<nK; idxK++) {
			delete mVOpt[idxR][idxK];
		}
	}
	for (idxR=0; idxR<nR; idxR++) {
		delete mRateResets[0][idxR];
		delete mZeroResets[0][idxR];
		delete mRateResets[1][idxR];
		delete mZeroResets[1][idxR];
		delete mPZeroResets[idxR];
	}


	throw KFailure("%s: failed.\n", routine);
    }
}


