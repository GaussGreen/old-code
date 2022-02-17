/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      David Liu
 ************************************************************************/
#include "vtbstest.h"
#include "kvtspar.h"	// Slice parser
#include "kbirtree.h"

#include "kutilios.h"


extern	KMap(String,double)	globVpConstTable;	// Constant tables


extern "C" {
#include "drltime.h"
#include "drloptio.h"
#include "vnfmanly.h"
};



const	int	nRMax=120;
const	int	nKMax=120;

int IntegJ(
	double		*beta,
	double		*alpha,
	double		*rho,
	double		s,
	KVector(TDate)	&dates,
	KVector(double)	&sigma1,
	KVector(double)	&sigma2,
	KVector(double)	&JInteg);

//---------------------------------------------------------------
//

void KVTreeTestBasisRatesVolAndCorr(
	KVTree *vt,			// (I) 
	
	KVolDiag& irVolDiag,		// (I)
	KVolDiag& bsVolDiag,		// (I)
	KMrParam&  treeModelParam,	// (I) full model parameter

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
static	char	routine[] = "KVTreeTestBasisRatesVolAndCorr";
	int	nR, idxR,
		nK, idxK;
	
	int	i;

	double	s;

	KVector(double)		JInteg;

	KBirTree *birTree = dynamic_cast<KBirTree*>(vt);

	KVector(KRateReset*)	mRateResets[2];
	KVector(KZeroReset*)	mZeroResets[2];
	KVector(KZeroReset*)	mPZeroResets;
	KVector(KZeroReset*)	mEZeroResets[2];
	KVector(TDate)		mModelResetDates[2];
	KVector(TDate)		mCorrResetDates;

	KVector(KRateReset*)	mLiborResets;	// libor at model reset dates

	KVector(KTSlice*)	mValue;


	TDate			currDate;
	TDate			todayDate = vt->TPToday();

	int			numSpotVols;
	KVector(TDate)		spotVolDates;
	KVector(double)		irSpotVols;
	KVector(double)		bsSpotVols;

	KVector(double)		spreadVol;
	KVector(double)		liborVol;
	KVector(double)		corr_LS;


	double			basisVol;

	TDate			dateL, dateR;

	double			liborL, liborR;
			

	KVector(KTSlice*)	mRate[2];

	KTSlice			*mZero = NULL;
	KTSlice			*mZero0 = NULL;
	KTSlice			*mZero1 = NULL;

	KVector(KTSlice*)	mVM0;		// 0-moment
	KVector(KTSlice*)	mVE0;		// 1-moments
	KVector(KTSlice*)	mVE1;
	KVector(KTSlice*)	mVC00;		// 2-moments
	KVector(KTSlice*)	mVCBL;
	KVector(KTSlice*)	mVC11;

	KVector(KTSlice*)	mVCLS;

	KVector(KTSlice*)	mSpread;	// spread 
	KVector(KTSlice*)	mVSpread;	// 2-moments of spread

	KVector(KTSlice*)	mLog0;		// log of rate 0 
	KVector(KTSlice*)	mLog1;		// log of rate 1
	KVector(KTSlice*)	mLogLibor;	// log of libor
	KVector(KTSlice*)	mVLogLibor;	// log of libor



	double			pVM0[nRMax],
				pVE0[nRMax],
				pVE1[nRMax],
				pSpread[nRMax],
				pVSpread[nRMax],
				pLog0[nRMax],
				pLog1[nRMax],
				pLogLibor[nRMax],
				pVLogLibor[nRMax],
				pVCLS[nRMax],
				pVC00[nRMax],
				pVCBL[nRMax],
				pVC11[nRMax],
				kstk[nRMax][nKMax];

	double			volBbq;		// backbone factor

	KDayCc 			dcc = KDayCc(GTO_ACT_365F);

#define	DEV(sliceP)	{if (sliceP != NULL) {sliceP->Dev(discCurveName);}}
#define	EV(sliceP)	{if (sliceP != NULL) {sliceP->Ev();}}



    try {


	ASSERT_OR_THROW(resetDates.size() == payDates.size());
	nR = resetDates.size();
	ASSERT_OR_THROW(nR < nRMax);

	nK = strikes.size();
	ASSERT_OR_THROW(nK < nKMax);


	// enough space to store the model reset dates and resets.
	JInteg.resize(3);

	mRateResets[0].resize(nR);
	mZeroResets[0].resize(nR);
	mModelResetDates[0].resize(nR);
	mRateResets[1].resize(nR);
	mZeroResets[1].resize(nR);
	mModelResetDates[1].resize(nR);

	mLiborResets.resize(nR);

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
	mVCBL.resize(nR);
	mVC11.resize(nR);

	mVCLS.resize(nR);

	mLogLibor.resize(nR);
	mVLogLibor.resize(nR);

	mSpread.resize(nR);
	mVSpread.resize(nR);
	mLog0.resize(nR);
	mLog1.resize(nR);

	liborVol.resize(nR);
	spreadVol.resize(nR);
	corr_LS.resize(nR);

	for (idxR=0; idxR<nR; idxR++) {
		mRate[0][idxR] = NULL;
		mRate[1][idxR] = NULL;

		mVM0[idxR] = NULL;
		mVE0[idxR] = NULL;
		mVE1[idxR] = NULL;
		mVC00[idxR] = NULL;
		mVCBL[idxR] = NULL;
		mVC11[idxR] = NULL;

		mVCLS[idxR] = NULL;
		mLogLibor[idxR] = NULL;
		mVLogLibor[idxR] = NULL;

		mSpread[idxR] = NULL;
		mVSpread[idxR] = NULL;

		mLog0[idxR] = NULL;
		mLog1[idxR] = NULL;
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

		// rate 0

		mRateResets[0][idxR] = new KRateReset(
			curveName0,
			resetDates[idxR],
			resetDates[idxR], // !!! effective dates!
			floatRate0);
		mModelResetDates[0][idxR] =
				vt->Insert(*mRateResets[0][idxR]);
		mZeroResets[0][idxR] = new KZeroReset(
			discCurveName,
			mModelResetDates[0][idxR],
			payDates[idxR]);

		vt->Insert(resetDates[idxR]);
		vt->Insert(*mZeroResets[0][idxR]);


		// rate 1

		mRateResets[1][idxR] = new KRateReset(
			curveName1,
			resetDates[idxR],
			resetDates[idxR], // !!! effective dates!
			floatRate1);
		mModelResetDates[1][idxR] =
				vt->Insert(*mRateResets[1][idxR]);
		mZeroResets[1][idxR] = new KZeroReset(
			discCurveName,
			mModelResetDates[1][idxR],
			payDates[idxR]);

		vt->Insert(resetDates[idxR]);
		vt->Insert(*mZeroResets[1][idxR]);


		// Libor reset at modelResetDate
		//
		mLiborResets[idxR] = new KRateReset(
			curveName0,
			mModelResetDates[1][idxR],
			mModelResetDates[1][idxR], // !!! effective dates!
			floatRate1);

		vt->Insert(*mLiborResets[idxR]);

		mCorrResetDates[idxR] = MIN(
			mModelResetDates[0][idxR],
			mModelResetDates[1][idxR]);

		// Payment
		mEZeroResets[0][idxR] = new KZeroReset(
			discCurveName,
			mCorrResetDates[idxR],
			mModelResetDates[0][idxR]);
		vt->Insert(*mEZeroResets[0][idxR]);

		mEZeroResets[1][idxR] = new KZeroReset(
			discCurveName,
			mCorrResetDates[idxR],
			mModelResetDates[1][idxR]);
		vt->Insert(*mEZeroResets[1][idxR]);

		mPZeroResets[idxR] = new KZeroReset(
			discCurveName,
			mCorrResetDates[idxR],
			payDates[idxR]);
		vt->Insert(*mPZeroResets[idxR]);



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
	// Inititialize tree timeline and Calibrate tree
	//
	vt->SetUpTimeline();
	vt->Calibrate();


	// Get spot vols
	//
	spotVolDates = birTree->SpotVolDates();
	irSpotVols   = birTree->IrSpotVols();
	bsSpotVols   = birTree->BsSpotVols();
	numSpotVols = spotVolDates.size();
	ASSERT_OR_THROW(irSpotVols.size() == numSpotVols);
	ASSERT_OR_THROW(bsSpotVols.size() == numSpotVols);



	//
	// Rollback
	//
	for (int tpIdx=vt->TPNum(); tpIdx >= 0; tpIdx--) {
	    //
	    // Update Tree
	    //
	    vt->Update(tpIdx);


	    currDate = vt->TPDateCurrent();

	    //
	    // Dev
	    //
	    for (idxR=0; idxR<nR; idxR++) {
		EV(mRate[0][idxR]);
		EV(mRate[1][idxR]);

		DEV(mVM0[idxR]);
		DEV(mVE1[idxR]);
		DEV(mVE0[idxR]);
		EV(mVC00[idxR]);
		EV(mVCBL[idxR]);
		EV(mVC11[idxR]);

		EV(mSpread[idxR]);
		EV(mVSpread[idxR]);

		EV(mLog0[idxR]);
		EV(mLog1[idxR]);

		EV(mVCLS[idxR]);

		EV(mLogLibor[idxR]);
		EV(mVLogLibor[idxR]);

	    }


	    if (mZero == NULL) {
		mZero  = new KTSlice(*vt, "mZero",  discCurveName);
		mZero0 = new KTSlice(*vt, "mZero0", discCurveName);
		mZero1 = new KTSlice(*vt, "mZero1", discCurveName);
	    }



	    //
	    // if Reset date: the floating rate can be computed.
	    //
	    for (idxR=0; idxR<nR; idxR++) {

	      // 
	      // Index 0
	      // 
	      if (currDate == mModelResetDates[0][idxR]) {
		mRate[0][idxR] = new KTSlice(*vt, "mRate0", curveName0);

		vt->Get(*(mRate[0][idxR]), *mRateResets[0][idxR]);

		mVE0[idxR]  = new KTSlice(*vt, "mVE0", curveName0);
		*mVE0[idxR]   = *mRate[0][idxR];
		
		mVM0[idxR]  = new KTSlice(*vt, "mVM0", curveName0);
		mVC00[idxR] = new KTSlice(*vt, "mVC00", curveName0);


		mLog0[idxR] = new KTSlice(*vt, "mLog0", curveName0);
		*mLog0[idxR]   = *mRate[0][idxR];
		mLog0[idxR]->Log();

		if (payDates[idxR] != mModelResetDates[0][idxR]) {
			vt->Get(*mZero, *mZeroResets[0][idxR]);
			*mVE0[idxR] *= *mZero;
		}


		*mVM0[idxR]   = *mZero;

		*mVC00[idxR]  = *mLog0[idxR];
		*mVC00[idxR] *= *mLog0[idxR];
	      }


	      // 
	      // Index 1
	      // 
	      if (currDate == mModelResetDates[1][idxR]) {
		mRate[1][idxR] = new KTSlice(*vt, "mRate1", curveName1);

	      	mSpread[idxR]  = new KTSlice(*vt, "mSpread", curveName1);
	      	mVSpread[idxR] = new KTSlice(*vt, "mVSpread", curveName1);

	      	mLogLibor[idxR]  = new KTSlice(*vt, "mLogLibor", curveName0);
	      	mVLogLibor[idxR] = new KTSlice(*vt, "mVLibor", curveName0);

		mVCLS[idxR] = new KTSlice(*vt, "mVCLS", curveName1);

		mLog1[idxR] = new KTSlice(*vt, "mLog1", curveName1);
		mVC11[idxR] = new KTSlice(*vt, "mVC11", curveName1);

		mVCBL[idxR] = new KTSlice(*vt, "mVCBL", curveName1);

		vt->Get(*(mRate[1][idxR]), *mRateResets[1][idxR]);
		vt->Get(*(mLogLibor[idxR]), *mLiborResets[idxR]);
		
		//if (KBirTree *birTree = dynamic_cast<KBirTree*>(vt))
			birTree->GetSpread(*mSpread[idxR], currDate);
		//else
		//	throw KFailure("%s: failed dynamic_cast to birtree.\n",
				//	routine);

		// Interpolate the basis spread vol
		//
		s = DayCountFraction(
			vt->TPToday(),
			mModelResetDates[1][idxR],
			dcc);


		// Get the spot vola information
		//


		ASSERT_OR_THROW(IntegJ(
				treeModelParam.mBeta,
				treeModelParam.mAlpha,
				treeModelParam.mRho,
				s,
				spotVolDates,	//irVolDiag.mVolDates,
				irSpotVols,	//irVolDiag.mSpotVolRates,
				bsSpotVols,	//bsVolDiag.mSpotVolRates,
				JInteg) == SUCCESS);

		volBbq = birTree->GetBasisSpreadVolBbq(vt->TPIdx(currDate));

		spreadVol[idxR] = volBbq * sqrt(JInteg[2]/s);
		corr_LS[idxR] = JInteg[1]/sqrt(JInteg[0]*JInteg[2]);

		for (i=0; i<= numSpotVols-2; i++)
		{
			dateL  = spotVolDates[i];
			dateR  = spotVolDates[i+1];

			liborL = irVolDiag.VolInterp(dateL);
			liborR = irVolDiag.VolInterp(dateR);

			if ((dateL<=currDate) && (dateR>currDate))
			{
				liborVol[idxR] = liborL+(liborR - liborL)
					   *(double)(currDate-dateL)
					   /(double)(dateR-dateL);	

				break;
			}
		}
		
		*mLog1[idxR]   = *mRate[1][idxR];
		mLog1[idxR]->Log();

		*mVC11[idxR]  = *mLog1[idxR];
		*mVC11[idxR] *= *mLog1[idxR];

		mLogLibor[idxR]->Log();

		mSpread[idxR]->Log();

		*mVSpread[idxR]  = *mSpread[idxR];
		*mVSpread[idxR] *= *mSpread[idxR];
		
		*mVLogLibor[idxR]  = *mLogLibor[idxR];
		*mVLogLibor[idxR] *= *mLogLibor[idxR];

		*mVCBL[idxR]  = *mLogLibor[idxR];
		*mVCBL[idxR] *= *mLog1[idxR];

		*mVCLS[idxR]	   = *mLogLibor[idxR];
		*mVCLS[idxR]	  *= *mSpread[idxR];

		mVE1[idxR]  = new KTSlice(*vt, "mVE1", curveName1);
		*mVE1[idxR] = *mRate[1][idxR];


		if (payDates[idxR] != mModelResetDates[1][idxR]) {
			vt->Get(*mZero, *mZeroResets[1][idxR]);
			*mVE1[idxR] *= *mZero;
		}
	      }



	    }



	    //
	    // Last TP: store the result.
	    //
	    if (vt->TPIdxCurrent() == 0) {
	      for (idxR=0; idxR<nR; idxR++) {

		pVM0[idxR]  = mVM0[idxR]->GetCenter();
		pVE0[idxR]  = mVE0[idxR]->GetCenter();
		pVE1[idxR]  = mVE1[idxR]->GetCenter();

		pSpread[idxR]   = mSpread[idxR]->GetCenter();
		pVSpread[idxR]  = mVSpread[idxR]->GetCenter();

		pLog0[idxR]  = mLog0[idxR]->GetCenter();
		pLog1[idxR]  = mLog1[idxR]->GetCenter();

		pLogLibor[idxR]  = mLogLibor[idxR]->GetCenter();
		pVLogLibor[idxR]  = mVLogLibor[idxR]->GetCenter();

		pVC00[idxR] = mVC00[idxR]->GetCenter();
		pVCBL[idxR] = mVCBL[idxR]->GetCenter();
		pVC11[idxR] = mVC11[idxR]->GetCenter();

		pVCLS[idxR] = mVCLS[idxR]->GetCenter();

	      }
	    }
	}


	//
	//
	//
	os << format("TODAY: %10s\n", DrlTDatePrint(NULL, vt->TPToday()));

	//
	// Expected value and var
	//



	os << "FORWARD RATES:" << endl;
	os << " RESET DATE   TEXP  TPAY ";
	os << "| ZERO(F)    ZERO(TREE)    DIFF     ";
	os << "|   RATE_L(F) RATE_L(T)  DIFF1  ";
	os << "|   RATE_B(F) RATE_B(T)  DIFF2  ";
	os << "|   VOL_L(T) VOL_L(ENV) DIFF_L  ";
	os << "|  VOL_S(T) VOL_S(ANL) DIFF_S   ";
	os << "|  VOL_B(T) VOL_B(ANL) DIFF_B   ";
	os << "|  CR_LS(T) CR_LS(ANL) DIFF_LS   ";
	os << "|  CR_LB(T) CR_LB(ANL) DIFF_LB";
	os << endl;
	for (idxR=0; idxR<nR; idxR++) {
		double	tr,		// time to reset
			te,		// time to model reset
			te0,		// time to model reset 0
			te1,		// time to model reset 1
			tp,		// time to pay
			zf,		// disc fact
			ef0,		// fwd value
			ef1,		// fwd value
			zv,		// disc fact
			ev0,		// expected value
			ev1,		// expected value
			elog0,		// expected log value
			elog1,		// expected log value
			espread,	// expected log spread value
			evspread,	// expected log(spread)^2 value
			eLibor,		// expected log(Libor) value
			evLibor,	// expected log(Libor)^2 value
			svS,		// vol of spread
			svB,		// vol of basis rate
			svL,		// vol of Libor rate
			sv0, 		// stdev and corr
			cr_LS_t, cr_LS_a, cr_LB_t, cr_LB_a,	// corr
			c00, cLS, cBL, c11;	// covariance


		zv  = pVM0[idxR];
		ev0 = pVE0[idxR] / zv;
		ev1 = pVE1[idxR] / zv;

		elog0 = pLog0[idxR];
		elog1 = pLog1[idxR];

		espread  = pSpread[idxR];
		evspread = pVSpread[idxR];

		eLibor  = pLogLibor[idxR];
		evLibor = pVLogLibor[idxR];

		c00 = pVC00[idxR];
		cBL = pVCBL[idxR];
		cLS = pVCLS[idxR];
		c11 = pVC11[idxR];

		// calc fwd and disc fact
		ef0 = mRateResets[0][idxR]->Forward(zcCurve0);
		ef1 = mRateResets[1][idxR]->Forward(zcCurve1);

		zf = discZcCurve.DiscFact(payDates[idxR]);

		tr = DayCountFraction(
			vt->TPToday(),
			resetDates[idxR],
			dcc);
		te = DayCountFraction(
			vt->TPToday(),
			mCorrResetDates[idxR],
			dcc);
		te0 = DayCountFraction(
			vt->TPToday(),
			mModelResetDates[0][idxR],
			dcc);
		te1 = DayCountFraction(
			vt->TPToday(),
			mModelResetDates[1][idxR],
			dcc);
		tp = DayCountFraction(
			vt->TPToday(),
			payDates[idxR],
			dcc);


		os << format(" %10s", DrlTDatePrint(NULL, resetDates[idxR]));
		os << format(" %5.2f", tr);
		os << format(" %5.2f", tp);

		os << format(" |");
		os << format(" %10.8f", zf);
		os << format(" %10.8f", zv);
		os << format(" %11.8f", zv-zf);

		os << format(" |");
		os << format(" %8.4f", ef0*1e2);
		os << format(" %8.4f", ev0*1e2);
		os << format(" %10.6f", (ev0-ef0)*1e2);

		os << format(" |");
		os << format(" %8.4f", ef1*1e2);
		os << format(" %8.4f", ev1*1e2);
		os << format(" %10.6f", (ev1-ef1)*1e2);

		sv0 = sqrt((c00-elog0*elog0)/te0);
		os << format(" |");
		os << format(" %8.4f", sv0 * 1e2);
		os << format(" %8.4f", liborVol[idxR] * 1e2);
		os << format(" %10.6f", (sv0-liborVol[idxR]) * 1e2);

		svS = sqrt((evspread - espread*espread)/te1);
		os << format(" |");
		os << format(" %8.4f", svS * 1e2);
		os << format(" %8.4f", spreadVol[idxR] * 1e2);
		os << format(" %10.6f", (svS-spreadVol[idxR]) * 1e2);

		svB = sqrt((c11-elog1*elog1)/te1);

		svL = sqrt((evLibor-eLibor*eLibor)/te1);

		cr_LS_t  = (cLS - eLibor*espread) / 
			  sqrt((evLibor-eLibor*eLibor)
				*(evspread - espread*espread));

		cr_LS_a  = corr_LS[idxR];

		basisVol = sqrt(svL*svL + svS*svS + 2e0*cr_LS_t*svL*svS);

		os << format(" |");
		os << format(" %8.4f", svB * 1e2);
		os << format(" %8.4f", basisVol * 1e2);
		os << format(" %10.6f", (svB-basisVol) * 1e2);

		os << format(" |");
		os << format(" %8.4f", cr_LS_t);
		os << format(" %8.4f", cr_LS_a);
		os << format(" %10.6f", (cr_LS_t - cr_LS_a));

		cr_LB_t  = (cBL - eLibor*elog1) / 
			  sqrt((evLibor-eLibor*eLibor)
				*(c11-elog1*elog1));

		cr_LB_a   = (cr_LS_t*svS + svL) / svB;

		os << format("  |");
		os << format(" %8.4f", cr_LB_t);
		os << format(" %8.4f", cr_LB_a);
		os << format(" %10.6f", (cr_LB_t - cr_LB_a));

		os << endl;


	}









	//
	// Free memory
	//
	delete mZero;
	delete mZero0;
	delete mZero1;

	for (idxR=0; idxR<nR; idxR++) {
		delete mRate[0][idxR];
		delete mRate[1][idxR];

		delete mLog0[idxR];
		delete mLog1[idxR];

		delete mSpread[idxR];
		delete mVSpread[idxR];

		delete mLogLibor[idxR];
		delete mVLogLibor[idxR];

		delete mVM0[idxR];
		delete mVE1[idxR];
		delete mVE0[idxR];
		delete mVC00[idxR];
		delete mVC11[idxR];
		delete mVCBL[idxR];
		delete mVCLS[idxR];

	}

	for (idxR=0; idxR<nR; idxR++) {
		delete mLiborResets[idxR];
		delete mRateResets[0][idxR];
		delete mZeroResets[0][idxR];
		delete mRateResets[1][idxR];
		delete mZeroResets[1][idxR];
		delete mPZeroResets[idxR];
		delete mEZeroResets[0][idxR];
		delete mEZeroResets[1][idxR];
	}
 


    }
    catch (KFailure) {
	delete mZero;
	delete mZero0;
	delete mZero1;
	for (idxR=0; idxR<nR; idxR++) {
		delete mRate[0][idxR];
		delete mRate[1][idxR];

		delete mLog0[idxR];
		delete mLog1[idxR];

		delete mSpread[idxR];
		delete mVSpread[idxR];

		delete mLogLibor[idxR];
		delete mVLogLibor[idxR];

		delete mVM0[idxR];
		delete mVE1[idxR];
		delete mVE0[idxR];
		delete mVC00[idxR];
		delete mVC11[idxR];
		delete mVCBL[idxR];
		delete mVCLS[idxR];

	}
	for (idxR=0; idxR<nR; idxR++) {
		delete mLiborResets[idxR];
		delete mRateResets[0][idxR];
		delete mZeroResets[0][idxR];
		delete mRateResets[1][idxR];
		delete mZeroResets[1][idxR];
		delete mPZeroResets[idxR];
		delete mEZeroResets[0][idxR];
		delete mEZeroResets[1][idxR];
	}


	throw KFailure("%s: failed.\n", routine);
    }
}



int IntegJ(
	double		*beta,
	double		*alpha,
	double		*rho,
	double		s,
	KVector(TDate)	&dates,
	KVector(double)	&sigma1,
	KVector(double)	&sigma2,
	KVector(double)	&JInteg)
{
static  char    routine[] = "IntegJ";

	int		status = FAILURE;
	
	int		i;

	FILE		*debugFp = DppLogFp();

	VnfmData        *vnfmData = NULL;
        TCurve          *zcCurve = NULL;
	double		backboneqL[2] = {1, 0e0};
	double		betaL[3] = {2, beta[0], beta[1]};
	double		alphaL[3]= {2, alpha[0], alpha[1]};
	KVector(double)	rhoL;
	KVector(TDate)	dateL;
	KVector(double)	sigmaL;

	TDate		baseDate = dates[0];
	TDate		zcDates[1] = {dates[1]};
	double		zcRates[1] = {0.1e0};


  try{

	// Create a dummy zero curve
	//
	ASSERT_OR_THROW((zcCurve = GtoMakeTCurve(baseDate, 
				    		 zcDates,	
				    		 zcRates,	
				    		 1,
				    		 1e0,
				    		 GTO_ACT_365F)) != NULL);

	// Add the size to dateL[0] 
	//
	dateL.push_back(dates.size());
	dateL.insert(dateL.end(), dates.begin(), dates.end());

	// Check size of sigams
	//
	if (dates.size()!=sigma1.size())
		throw KFailure("size of dates (%d) != size of sigma1 (%d).\n",
				dates.size(), sigma1.size());

	if (sigma1.size()!=sigma2.size())
		throw KFailure("size of sigma1 (%d) != size of sigma2 (%d).\n",
				sigma1.size(), sigma2.size());

	rhoL.push_back(dates.size());
	rhoL.insert(rhoL.end(), dates.size(), rho[0]);

	sigmaL.push_back(2*sigma1.size());

	for(i=0; i<=sigma1.size()-1;i++){
		sigmaL.insert(sigmaL.end(), sigma1[i]);
		sigmaL.insert(sigmaL.end(), sigma2[i]);
	}
	
	/* get model parameters */
        ASSERT_OR_THROW(VnfmWrapRead(&vnfmData,
			backboneqL, betaL, alphaL,
			&dateL[0], &sigmaL[0], &rhoL[0],
			zcCurve) == SUCCESS);

	//VnfmFpWrite(vnfmData, debugFp);

	/* Generate corr matrix */
        ASSERT_OR_THROW(VnfmComputeCoeff(vnfmData) == SUCCESS);


	ASSERT_OR_THROW(VnfmJ(vnfmData,
		 	      0e0,
		 	      s,
		 	      s,
		 	      &JInteg[0]) == SUCCESS);

	GtoFreeTCurve(zcCurve);
	VnfmFree(vnfmData);
	
	status = SUCCESS;

	return (status);

    }

  catch (KFailure) {
	GtoFreeTCurve(zcCurve);
	VnfmFree(vnfmData);
	throw KFailure("%s: failed.\n", routine);
    }
}
