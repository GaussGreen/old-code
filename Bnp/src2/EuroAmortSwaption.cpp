
// ------------------------------------------------------------------------------------------------------------------
// // EuroAmortSwaption.cpp
//
// Implementation of functions to price european amortizing swaptions
//
// ------------------------------------------------------------------------------------------------------------------
// //
//

#include "EuroAmortSwaption.h"
#include "LGMCalib.h"
#include "nag.h"
#include "nagc05.h"
#include "opfnctns.h"
#include "srt_h_lgmUSprotos.h"
#include "swp_h_curve_list.h"
#include "swp_h_swap_pricing.h"

#ifdef WIN32
#ifdef _DEBUG
#include "crtdbg.h"
#endif //_DEBUG
#endif // WIN32

/*
        static data for local volatility function defined to support old LGM1F
   pricer.
*/
static Err (*LGM1F_GetVolNew)(char *vol_curve_name, double start_date,
                              double end_date, double cash_strike, int zero,
                              char *ref_rate_name, double *vol, double *power);

static char *LGM1F_szYieldCurveName;
static char *LGM1F_szVolCurveName;
static char *LGM1F_szRefRate;
static LGMMarkConv LGM1F_MarkConv;
static long LGM1F_lToday;

/*
        init function for local volatility function defined to support old LGM1F
   pricer.
*/
static Err LGM1F_GetVolInit(Err (*getVolNew)(char *vol_curve_name,
                                             double start_date, double end_date,
                                             double cash_strike, int zero,
                                             char *ref_rate_name, double *vol,
                                             double *power),
                            char *szYieldCurveName, char *szVolCurveName,
                            char *szRefRate, long lToday) {
  LGM1F_GetVolNew = getVolNew;
  LGM1F_szYieldCurveName = szYieldCurveName;
  LGM1F_szVolCurveName = szVolCurveName;
  LGM1F_szRefRate = szRefRate;
  LGM1F_lToday = lToday;
  return LGMCcyDefaults(lookup_curve(szYieldCurveName), &LGM1F_MarkConv);
}

/*
        local volatility function defined to support old LGM1F pricer.
*/
static LGMErr LGM1F_GetCashLogVol(Date lStart, Date lEnd, double dStrike,
                                  SRT_Boolean boolIsLog, double *dVol) {
  /* local variables */
  SwapDP swap_dp;
  double dPower, dForward, dExpiry;
  Err err;
  long lFix;

  /* call the new vol function */
  if (err =
          (*LGM1F_GetVolNew)(LGM1F_szVolCurveName, (double)lStart, (double)lEnd,
                             dStrike, 0, LGM1F_szRefRate, dVol, &dPower)) {
    *dVol = 0.0;
    return err;
  }

  /* if the returned volatility is normal  , convert it to log */
  if (dPower == 0.0) {
    /* Set up the swap DP */
    if (err = swp_f_setSwapDP(lStart, lEnd, LGM1F_MarkConv.sfreq,
                              LGM1F_MarkConv.sbasis, &swap_dp))
      return err;
    swap_dp.spot_lag = LGM1F_MarkConv.lag;

    /* Calculate the CASH forward swap rate */
    if (err = swp_f_ForwardRate_SwapDP(&swap_dp, LGM1F_szYieldCurveName, "CASH",
                                       &dForward))
      return err;

    /* Get the log vol */
    lFix = add_unit(swap_dp.start, -swap_dp.spot_lag, SRT_BDAY, SUCCEEDING);
    dExpiry = (lFix - LGM1F_lToday) * YEARS_IN_DAY;
    if (err = srt_f_optbetavoltoblkvol(dForward, dStrike, *dVol, dExpiry, 0.0,
                                       dVol))
      return err;
  }

  return 0;
}

//
// Improved EuroAmortSwaption
//
// Author: Albert Wang
// Date: 04/04/03

Err EuropeanAmortizingSwaption2(
    long AsOfDate, long FixDate, long TheoEndDate, long NumPeriods,
    long *plAcrlStartDates, // vector(1  ,NumPeriods)  ,
    long *plAcrlEndDates,   // vector(1  ,NumPeriods)  ,
    long *plPayDates,       // vector(1  ,NumPeriods)  ,
    double *pFixCvgsAdjusted, double *pFixCvgs, double Coupon,
    double *pdNotionals,                           // vector(1  ,NumPeriods)  ,
    char *pszYieldCurveName, char *szVolCurveName, /*	vc */
    char *pszFreq, char *pszBasis, char *pszPayRec, char *pszRefRate,
    Err (*get_cash_vol)(char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*								LGMErr  (*GetVol)( Date  , Date  , double  , SRT_Boolean
       , double *)  , SrtDiffusionType srt_vol_type  ,
    */
    // Output
    double *dvPayDates, // today's zerobond values to all paydates
    double *dvReplicatingStrikes, double *dvReplicatingNotionals,
    double *dPrice, // price of
    double *dFixedPV, double *dFloatPV, double *dSwapRate,
    SigKapTS **lgmSigKapTSPtrPtr, double *dLGMVol,
    double *dvReplicatingSwaptions,

    // more input
    double *pdMargins, double *pdCoupons, double *pdFltNotionals,
    double dExerciseFee) {
  // internal function declaration
  static char *EuropeanAmortizingSwaption_(
      long, long, long, long *, long *, double *, double *, long *, double,
      double *, char *, char *, char *, char *, char *, char *,
      LGMErr (*)(long, long, double, SRT_Boolean, double *), SRT_Boolean,
      LGM_TS *, LGMCalSet *, double *, double *, double *, double *, double *,
      double *, double *, double);
  void Init_CalibInstruments(long, long, long, long, long *, double *, char *,
                             LGMCalSet *);
  void LGMFreeLGM_CalibInst(long, LGMCalSet *);
  char *pszVolType = 0, *ErrMessage = 0;

  double dLambda, dTex, dT0, dT1, vol;
  double tmpsum, tmpswaption, tmpvar;
  long ntime, tmpinc, nI;
  char *psztmpPayRec = "REC";
  SRT_Boolean boolIsLog;
  LGM_TS *pTermStruct = 0;
  LGMCalSet MyCalibInstruments;
  SrtBasisCode FloatBasis;
  SrtCompounding FloatFrequency, FixedFrequency;
  SrtReceiverType PayRec, tmpPayRec;
  SrtBasisCode basis;

  if (!pdCoupons) {
    if (Coupon < 1.e-13 || Coupon > 1.) {
      ErrMessage = "Invalid Coupon in EuropeanAmortizingSwaption2(...)!";
      return ErrMessage;
    }
  }

  else {
    for (nI = 1; nI <= NumPeriods; ++nI) {
      if (pdCoupons[nI] < 1.e-13 || pdCoupons[nI] > 1.) {
        ErrMessage = "Invalid Coupon in EuropeanAmortizingSwaption2(...)!";
        return ErrMessage;
      }
    }
  }

  // restrict floating frequency to be no larger than fixed frequency
  // note FloatFrequency returned is 12/#months
  swp_f_get_ref_rate_details(pszRefRate, &FloatBasis, &FloatFrequency);
  interp_compounding(pszFreq, &FixedFrequency);

  if (FloatFrequency < FixedFrequency) {
    ErrMessage = "Floating Frequency must not be larger than Fixed Frequency";
    return ErrMessage;
  }

  /* Set up the vol function */
  LGM1F_GetVolInit(get_cash_vol, pszYieldCurveName, szVolCurveName, pszRefRate,
                   AsOfDate);

  interp_basis(pszBasis, &basis);

  // is passed into function by user
  interp_rec_pay(pszPayRec, &PayRec);

  // always "REC"
  interp_rec_pay(psztmpPayRec, &tmpPayRec);

  // fill dvPayDates  , i.e.  , today's zerobond values to all pay dates
  dvPayDates[0] = swp_f_df(AsOfDate, plAcrlStartDates[1], pszYieldCurveName);
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    dvPayDates[ntime] =
        swp_f_df(AsOfDate, plPayDates[ntime], pszYieldCurveName);
  }

  *dFixedPV = 0;
  *dFloatPV = 0;
  tmpinc = FixedFrequency / FloatFrequency;

  if (!pdCoupons) {
    for (ntime = 1; ntime <= NumPeriods; ++ntime) {
      // = Coupon*pFixCvgs[ntime]*dvPayDates[ntime] ;
      *dFixedPV +=
          pdNotionals[ntime] * Coupon * pFixCvgs[ntime] * dvPayDates[ntime];

      // for(indflt=indfix*tmpinc+1;indflt<=(indfix+1)*tmpinc;++indflt)
      //{
      // fltinc = (dvPayDates[ntime-1]-dvPayDates[ntime]);

      // adjust for margin
      // fltinc += (pdFltSpreads[indflt] + pdMargins[indflt]) *
      // pdFltCvgs[indflt] * swp_f_df(plFltDates[indflt]  ,plFltDates[indflt+1]
      // , pszYieldCurveName);

      *dFloatPV +=
          pdNotionals[ntime] * (dvPayDates[ntime - 1] - dvPayDates[ntime]);
      //}
    }
  } else {
    for (ntime = 1; ntime <= NumPeriods; ++ntime) {
      // = Coupon*pFixCvgs[ntime]*dvPayDates[ntime] ;
      *dFixedPV += pdNotionals[ntime] * pdCoupons[ntime] * pFixCvgs[ntime] *
                   dvPayDates[ntime];

      // for(indflt=indfix*tmpinc+1;indflt<=(indfix+1)*tmpinc;++indflt)
      //{
      // fltinc = (dvPayDates[ntime-1]-dvPayDates[ntime]);

      // adjust for margin
      // fltinc += (pdFltSpreads[indflt] + pdMargins[indflt]) *
      // pdFltCvgs[indflt] * swp_f_df(plFltDates[indflt]  ,plFltDates[indflt+1]
      // , pszYieldCurveName);

      *dFloatPV +=
          pdFltNotionals[ntime] * (dvPayDates[ntime - 1] - dvPayDates[ntime]);
      //}
    }
  }

  /* Default the volatility type to lognormal
          if ( srt_vol_type == SRT_LOGNORMAL )
          {
  */
  pszVolType = "LOGNORMAL";
  boolIsLog = SRT_TRUE;
  /*	}
          else
          {
                  pszVolType = "NORMAL";
                  boolIsLog = SRT_FALSE;
          }
  */

  // initilize term structure
  // 2 zeta  , NumPeriods G(Ti)  , i = 0  ,1  ,... NumPeriods
  // not sure why LGMCreateLGM_TS disallow 1 zeta
  pTermStruct = LGMCreateLGM_TS(2, NumPeriods + 1);

  if (pTermStruct == 0) {
    ErrMessage = "Failure to instantiate term structure occurred in "
                 "EuropeanAmortizingSwaption2()!";
    return ErrMessage;
  }

  // initialize calibration instruments
  Init_CalibInstruments(NumPeriods, AsOfDate, FixDate, plAcrlStartDates[1],
                        plPayDates, pFixCvgs, pszYieldCurveName,
                        &MyCalibInstruments);

  // if swaption is of receiving fix coupon type  ,
  // do nothing and let EuropeanAmortizingSwaption_ handle dExerciseFee
  // otherwise  , flit the sign of dExerciseFee

  if (PayRec != tmpPayRec) {
    dExerciseFee *= -1.;
  }

  // main function
  // always price using "REC"

  if (ErrMessage = EuropeanAmortizingSwaption_(
          AsOfDate, FixDate, NumPeriods, plAcrlStartDates, plAcrlEndDates,
          pFixCvgsAdjusted, pFixCvgs, plPayDates, Coupon, pdNotionals,
          pszYieldCurveName, pszFreq, pszBasis, psztmpPayRec, pszRefRate,
          pszVolType, LGM1F_GetCashLogVol, boolIsLog, pTermStruct,
          &MyCalibInstruments, dvPayDates, dPrice, dvReplicatingStrikes,
          dvReplicatingNotionals, // etas
          dvReplicatingSwaptions, pdCoupons, pdFltNotionals, dExerciseFee)) {
    return ErrMessage;
  }

  // if dvReplicatingNotionals are less than tol == 1.e-10
  // then set dvReplicatingNotionals to zero

  /*for(ntime=1;ntime<=NumPeriods; ++ntime)
  {
          if(dvReplicatingNotionals[ntime-1] < 1.e-10)
          {
                  dvReplicatingNotionals[ntime-1]  = 0.;
          }
  }*/

  // adjust for "error" in dates  , etc.
  // idea here is to adjust using put-call parity
  // will come back later for a permanent fix ...
  tmpsum = 0.;
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {

    // Get Vol from Vol Cube
    (*LGM1F_GetCashLogVol)(plAcrlStartDates[1], plAcrlEndDates[ntime],
                           dvReplicatingStrikes[ntime - 1], boolIsLog, &vol);

    // Get ESO Price
    if (ErrMessage = swp_f_Swaption(
            plAcrlStartDates[1], plPayDates[ntime], pszFreq, pszBasis, vol,
            dvReplicatingStrikes[ntime - 1], "PAY", pszRefRate,
            pszYieldCurveName, "PREMIUM", pszVolType, &tmpswaption)) {
      return ErrMessage;
    }

    tmpsum += tmpswaption * dvReplicatingNotionals[ntime - 1];

    if (PayRec != tmpPayRec) {
      dvReplicatingSwaptions[ntime - 1] = tmpswaption;
    }
  }

  // restore notionals and fix/floating legs to $ units
  tmpvar = fabs(pdNotionals[1]);
  tmpsum *= tmpvar;
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    dvReplicatingNotionals[ntime - 1] *= tmpvar;
  }

// temporarily commented out
#if 0
	if(fabs(- (*dFixedPV) + (*dFloatPV))>1.e-10)
	{	
		tmpratio = (tmpsum - *dPrice)/(- (*dFixedPV) + (*dFloatPV));

		// do not adjust "legs'; 
		// adjust swaption prices instead
		//(*dFixedPV) *= tmpratio;
		//(*dFloatPV) *= tmpratio;


		if(fabs(tmpratio)>1.e-10)
		{
			tmpratio = 1./tmpratio;
			tmpsum *= 	tmpratio;
			(*dPrice) *= tmpratio;

			for(ntime=1;ntime<=NumPeriods; ++ntime)
			{
				dvReplicatingNotionals[ntime-1] *= tmpratio;
			}
		}
	}

#endif //#if 0

  // if Amort swap is of pay type  ,
  // then adjust price
  if (PayRec != tmpPayRec) {

    *dPrice = tmpsum;
    //*dPrice = *dPrice - (*dFixedPV) + (*dFloatPV);
  }

  (*dPrice) = (*dPrice) > 0. ? (*dPrice) : 0.;

  *dSwapRate = Coupon * (*dFloatPV) / (*dFixedPV);

  // get out the sigma tau term structure
  *lgmSigKapTSPtrPtr = LGMConvertZGtoSigKap(AsOfDate, pTermStruct);

  dLambda = (*lgmSigKapTSPtrPtr)->kap[0];

  // dTex = ( tEx[0] - AsOfDate ) / 365.0;
  dTex = (FixDate - AsOfDate) / 365.0;

  // dT1 = ( lvFixedEndDates[0] - today ) / 365.0;
  dT1 = (plPayDates[1] - AsOfDate) / 365.0;

  // dT0 = ( StartDate - today ) / 365.0;
  dT0 = (plAcrlStartDates[1] - AsOfDate) / 365.0;

  //*dLGMVol = (ex_G - coupon_G[1] ) * sqrt(2.0*zeta) * dLambda
  //					* sqrt( dLambda / (exp(2.0*dLambda*dTex) - 1.0)
  //) 					/ ( exp(-dLambda*dT0) - exp(-dLambda*dT1) );

  *dLGMVol = (pTermStruct->G[0] - pTermStruct->G[1]) *
             sqrt(2.0 * (pTermStruct->zeta[0])) * dLambda *
             sqrt(dLambda / (exp(2.0 * dLambda * dTex) - 1.0)) /
             (exp(-dLambda * dT0) - exp(-dLambda * dT1));

  LGMFreeLGM_CalibInst(NumPeriods, &MyCalibInstruments);
  LGMFreeLGM_TS(&pTermStruct);
  // free_vector(pdCoverages  ,1  ,NumPeriods);

  return 0;
}

static char *
UpdateCalibInstr(long NumPeriods, double *pdCoupons, char *pszFreq,
                 char *pszBasis, char *pszPayRec, char *pszRefRate,
                 char *pszVolType, char *pszYieldCurveName,
                 LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *),
                 SRT_Boolean boolIsLog, LGMCalSet *pLGMCalSet, double *pdESOs) {
  long ntime;
  double vol;
  char *ErrMessage = 0;

  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    /* Rfix[*  ,...  ,nlast] Rfix[j] is fixed rate for swaption j */
    pLGMCalSet->Rfix[ntime] = pdCoupons[ntime];

    // Get Vol from Vol Cube
    (*GetVol)(pLGMCalSet->tStart[ntime], pLGMCalSet->tPay[ntime],
              pLGMCalSet->Rfix[ntime], boolIsLog, &vol);

    // Get ESO Price
    if (ErrMessage =
            swp_f_Swaption(pLGMCalSet->tStart[ntime], pLGMCalSet->tPay[ntime],
                           pszFreq, pszBasis, vol, pLGMCalSet->Rfix[ntime],
                           pszPayRec, pszRefRate, pszYieldCurveName, "PREMIUM",
                           pszVolType, &(pLGMCalSet->Vfix[ntime]))) {
      return ErrMessage;
    }

    pdESOs[ntime] = pLGMCalSet->Vfix[ntime];
  }

  return 0;
}

static char *
SetInitCalibInstr(long NumPeriods, double Coupons, char *pszFreq,
                  char *pszBasis, char *pszPayRec, char *pszRefRate,
                  char *pszVolType, char *pszYieldCurveName,
                  LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *),
                  SRT_Boolean boolIsLog, LGMCalSet *pLGMCalSet) {
  long ntime;
  double vol;
  char *ErrMessage = 0;

  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    /* Rfix[*  ,...  ,nlast] Rfix[j] is fixed rate for swaption j */
    pLGMCalSet->Rfix[ntime] = Coupons;

    // Get Vol from Vol Cube
    (*GetVol)(pLGMCalSet->tStart[ntime], pLGMCalSet->tPay[ntime],
              pLGMCalSet->Rfix[ntime], boolIsLog, &vol);

    // Get ESO Price
    if (ErrMessage =
            swp_f_Swaption(pLGMCalSet->tStart[ntime], pLGMCalSet->tPay[ntime],
                           pszFreq, pszBasis, vol, pLGMCalSet->Rfix[ntime],
                           pszPayRec, pszRefRate, pszYieldCurveName, "PREMIUM",
                           pszVolType, &(pLGMCalSet->Vfix[ntime]))) {
      return ErrMessage;
    }
  }

  return 0;
}

void LGMFreeLGM_CalibInst(long NumPeriods, LGMCalSet *pLGMCalSet) {
  free_lngvector(pLGMCalSet->tPay, 0, NumPeriods);
  free_vector(pLGMCalSet->cvgpay, 0, NumPeriods);
  free_vector(pLGMCalSet->Dpay, 0, NumPeriods);
  free_lngvector(pLGMCalSet->tEx, 0, NumPeriods);
  free_lngvector(pLGMCalSet->tStart, 0, NumPeriods);
  free_vector(pLGMCalSet->DStart, 0, NumPeriods);
  free_vector(pLGMCalSet->Rfix, 0, NumPeriods);
  free_vector(pLGMCalSet->Vfix, 0, NumPeriods);
}

void Init_CalibInstruments(long NumPeriods, long AsOfDate, long FixDate,
                           long AcrlStartDate, long *plPayDates,
                           double *pdCoverages, char *pszYieldCurveName,
                           LGMCalSet *pLGMCalSet) {

  long ntime;
  double vol = 0., tmpDStart = 0.;
  char *ErrMessage = 0;

  // Set values that are not needed by "1 into k swaptions"
  // to 0 or garbage values

  /* 1 into k swaptions */
  pLGMCalSet->ExIndex = -1;
  pLGMCalSet->MAXCpn = -1;

  /* long swaptions */
  pLGMCalSet->longflag = 0;
  pLGMCalSet->ifirst = 0;
  pLGMCalSet->nlong = 0;
  pLGMCalSet->cvgfirst = 0;
  pLGMCalSet->FrLong = 0;
  pLGMCalSet->Rflong = 0;
  pLGMCalSet->Vlong = 0;
  pLGMCalSet->VegaLong = 0;
  pLGMCalSet->StDLong = 0;
  pLGMCalSet->BetaLong = 0;
  pLGMCalSet->CEVLong = 0;

  /* short swaptions */
  pLGMCalSet->shortflag = 0;
  pLGMCalSet->nshort = 0;
  pLGMCalSet->Rfshort = 0;
  pLGMCalSet->Vshort = 0;

  /* caplets */
  pLGMCalSet->capflag = 0;
  pLGMCalSet->tEnd = 0;
  pLGMCalSet->Dcap = 0;
  pLGMCalSet->cvgcap = 0;
  pLGMCalSet->Rfcap = 0;
  pLGMCalSet->Vcap = 0;
  pLGMCalSet->VegaCap = 0;

  // finally  , Set Calibration Instruments

  pLGMCalSet->n = NumPeriods;    /* tPay[n] is last date we need G(t) for */
  pLGMCalSet->nPay = NumPeriods; /* number of pay dates */
  pLGMCalSet->tNow = AsOfDate;   /* evaluation date */
  pLGMCalSet->nEx = NumPeriods;  /* number of exercise dates */
  pLGMCalSet->fixflag = 1;       /* 1 means 1 into k have been created */
  pLGMCalSet->ftEx = FixDate;    /* exercise dates are all ftEx */
  pLGMCalSet->ftStart = AcrlStartDate; /* start dates are all ftStart */
  pLGMCalSet->ffirst = 1; /* pay dates of swap. j are tPay[ffirst  , ...  ,j] */
  pLGMCalSet->nfirst =
      1; /* the swaptions are j = nfirst  , nfirst+1  , ...  , nlast */
  pLGMCalSet->nlast = NumPeriods;
  pLGMCalSet->cvgfix = pdCoverages[1]; /* cvgfix is the cvg for first period
                                          ftStart to tPay[ffirst] */
  pLGMCalSet->DfixStart =
      swp_f_df(AsOfDate, pLGMCalSet->ftStart,
               pszYieldCurveName); /* Dfixst is discount factor at ftStart[i] */

  // tmpDStart = swp_f_df(AsOfDate  ,AcrlStartDate  ,pszYieldCurveName);

  // create memory
  pLGMCalSet->tPay = lngvector(0, NumPeriods);
  pLGMCalSet->cvgpay = vector(0, NumPeriods);
  pLGMCalSet->Dpay = vector(0, NumPeriods);
  pLGMCalSet->tEx = lngvector(0, NumPeriods);
  pLGMCalSet->tStart = lngvector(0, NumPeriods);
  pLGMCalSet->DStart = vector(0, NumPeriods);
  pLGMCalSet->Rfix = vector(0, NumPeriods);
  pLGMCalSet->Vfix = vector(0, NumPeriods);

  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    /* tPay[*  ,1  ,...  ,nPay] dates for swaption fixed legs */
    pLGMCalSet->tPay[ntime] = plPayDates[ntime];

    /* cvgpay[*  ,1  ,...  ,nPay] cvg(tPay[i-1]  ,tPay[i]) */
    pLGMCalSet->cvgpay[ntime] = pdCoverages[ntime];

    /* Dpay[*  ,1  ,...  ,nPay] discount factor to tPay[i] */
    pLGMCalSet->Dpay[ntime] =
        swp_f_df(AsOfDate, pLGMCalSet->tPay[ntime], pszYieldCurveName);

    /* tEx[*  ,1  ,...  ,nEx] exercise (fixing) dates */
    pLGMCalSet->tEx[ntime] = FixDate;

    /* tStart[*  ,1  ,...  ,nEx] start (settlement) dates */
    pLGMCalSet->tStart[ntime] = AcrlStartDate;

    /* DStart[*  ,1  ,...  ,nEx] discount factor to tStart[j] */
    pLGMCalSet->DStart[ntime] = pLGMCalSet->DfixStart;

    /* Rfix[*  ,...  ,nlast] Rfix[j] is fixed rate for swaption j
    //pLGMCalSet->Rfix[ntime] = pdCoupons[ntime];

    // Get Vol from Vol Cube
    (*GetVol)(pLGMCalSet->tStart[ntime]  ,
                      pLGMCalSet->tPay[ntime]  ,
                      pLGMCalSet->Rfix[ntime]  ,
                      boolIsLog  ,
                      &vol);

     // Get ESO Price
    if(ErrMessage = swp_f_Swaption(	pLGMCalSet->tStart[ntime]  ,
                                                                    pLGMCalSet->tPay[ntime]
    , pszFreq  , pszBasis  , vol  , pLGMCalSet->Rfix[ntime]  , pszPayRec  ,
                                                                    pszRefRate ,
                                                                    pszYieldCurveName
    , "PREMIUM"  , pszVolType  ,
                                                                    &(pLGMCalSet->Vfix[ntime])))
    {
            return ErrMessage;
    }*/
  }

  // return 0;
}

static char *EuropeanAmortizingSwaption_(
    long AsOfDate, long FixDate, long NumPeriods,
    long *plAcrlStartDates, // vector(1  ,NumPeriods)  ,
    long *plAcrlEndDates,   // vector(1  ,NumPeriods)  ,
    double *pdCoverages,    // vector(1  ,NumPeriods)  ,
    double *pdCoveragesUnAdjusted,
    long *plPayDates, // vector(1  ,NumPeriods)  ,
    double Coupon,
    double *pdNotionals, // vector(1  ,NumPeriods)  ,
    char *pszYieldCurveName, char *pszFreq, char *pszBasis, char *pszPayRec,
    char *pszRefRate, char *pszVolType,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *),
    SRT_Boolean boolIsLog, LGM_TS *pTermStruct, LGMCalSet *pCalibInstr,
    double *pdInitZeroBonds, double *ASOPrice, double *dvReplicatingStrikes,
    double *dvReplicatingNotionals, double *dvReplicatingSwaptions,
    double *pdFixLegCoupons, double *pdFltNotionals, double dExerciseFee) {
  // internal function declaration
  void ComputeESOStrikes(long, double *, double *, double *);
  void ComputeESOWeights(long, double *, double *, double *, double *);
  void ComputeBondWeights(long, double, double *, double *, double *, double *,
                          double *);
  double AmortSwapPV(long, long, long, long *, char *, double *);
  static char *AmortSwaption(
      long, long, long *, char *, char *, char *, char *, char *, char *,
      double *, double *, SRT_Boolean,
      LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *), double *);
  static char *UpdateCalibInstr(
      long, double *, char *, char *, char *, char *, char *, char *,
      LGMErr (*)(long, long, double, SRT_Boolean, double *), SRT_Boolean,
      LGMCalSet *, double *);
  static char *SetInitCalibInstr(
      long, double, char *, char *, char *, char *, char *, char *,
      LGMErr (*)(long, long, double, SRT_Boolean, double *), SRT_Boolean,
      LGMCalSet *);
  static char *LGMFindStateStar(long, long, long, long, long *, char *,
                                double *, double, double, double *, double,
                                double *, double *, double *, double, double,
                                double, double *, double);
  void DetermineStateBounds(long, long, long, long *, char *, double, double *,
                            double, double *);
  static char *IsValidTermStructure(LGM_TS *, long);
  static char *LGM1FZeroBond(long, long, long, char *, double, double, double,
                             double, double *);

  void free_all_vectors(long, double *, double *, double *, double *, double *,
                        double *);

  // pdA[i]'s are weights that describe Amortizing Swap
  // in terms of underlying Zero-Bonds
  double *pdA = vector(0, NumPeriods), *pdEta = vector(1, NumPeriods),
         *pdCoupons = vector(1, NumPeriods),
         *pdInitCoupons = vector(1, NumPeriods),
         *pdASFZeroBond = vector(1, NumPeriods), *ptmpG = vector(0, NumPeriods),
         *pdESOs = vector(1, NumPeriods);

  double ASFLevel = 0, averageGDiffSqr;
  double DStart = swp_f_df(AsOfDate, plAcrlStartDates[1], pszYieldCurveName);
  double tmpsum, zero_t_T0 = DStart, adjfac = 0, tmpasflevel = 0;

  // search range for StateStar
  double dl = -10., dh = 10.;
  // initial guess of StateStart
  double StateStar = 0;
  // tol for solving StateStar
  double statestartol = 1.e-15;
  // tol for solving term structure of G(Ti)
  double termstructtol = 1.e-5;

  char *ErrMessage = 0;
  long ntime, counter, maxcounter = 10;

#ifdef WIN32
#ifdef _DEBUG
  double asfvtmp = 0, asfvleveltmp = 0, asfvzerotmp = 0, sumtmp;
  long indtmp, indtmp1;
#endif //_DEBUG
#endif // WIN32

  // compute A's  , i.e.  , bond weights
  ComputeBondWeights(NumPeriods, Coupon, pdA, pdNotionals,
                     pdCoveragesUnAdjusted, pdFixLegCoupons, pdFltNotionals);

  // Normalize A's
  // Normalization is needed to ensure numerical stability for solver
  adjfac = 1. / fabs(pdA[0]);
  for (ntime = 0; ntime <= NumPeriods; ++ntime) {
    pdA[ntime] *= adjfac;
  }

  // compute pdInitCoupons to be passed to UpdateCalibInstr
  // so that we can have initial ESO strikes and prices
  tmpsum = 0.;
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    tmpsum += pdCoveragesUnAdjusted[ntime] * pdInitZeroBonds[ntime];
  }

  tmpsum = 1. / tmpsum;
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    pdInitCoupons[ntime] = (zero_t_T0 - pdInitZeroBonds[ntime]) * tmpsum;
  }

  // Compute initial ESO prices based on pdInitCoupons
  if (ErrMessage =
          UpdateCalibInstr(NumPeriods, pdInitCoupons, pszFreq, pszBasis,
                           pszPayRec, pszRefRate, pszVolType, pszYieldCurveName,
                           GetVol, boolIsLog, pCalibInstr, pdESOs)) {

    free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons, pdInitCoupons,
                     pdESOs);
    return ErrMessage;
  }

  // Calibrate using ESOs with strikes equal pdInitCoupons
  // just so that we have an initial term structure to start off with
  if ((ErrMessage = FitAllGFwdSwap(pTermStruct, pCalibInstr, pTermStruct->zdate,
                                   pTermStruct->zeta)) ||
      (ErrMessage = IsValidTermStructure(pTermStruct, NumPeriods))) {

    // if calibration fails
    // return max(intrinsic value  , 0.)

    *ASOPrice = AmortSwapPV(NumPeriods, AsOfDate, FixDate,
                            plPayDates, // vector(1  ,NumPeriods)
                            pszYieldCurveName,
                            pdA); // vector(0  ,NumPeriods)

    *ASOPrice = max((*ASOPrice), 0.);

    // De-normalize price
    *ASOPrice /= adjfac;

    for (ntime = 1; ntime <= NumPeriods; ++ntime) {
      // will come back...
      dvReplicatingStrikes[ntime - 1] = Coupon;
      dvReplicatingNotionals[ntime - 1] = 0.; // pdEta[ntime]*fabs(pdA[0]);
      dvReplicatingSwaptions[ntime - 1] = 0.; // pdESOs[ntime];
    }

    free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons, pdInitCoupons,
                     pdESOs);
    return 0;
  }

  // extract term structucture from pTermStruct to ptmpG
  for (ntime = 0; ntime <= NumPeriods; ++ntime) {
    ptmpG[ntime] = pTermStruct->G[ntime];
  }

  // determine lower bound for state variable  , based
  // on if is based on zero bond value <=1.
  DetermineStateBounds(NumPeriods, AsOfDate, plAcrlStartDates[1], plPayDates,
                       pszYieldCurveName, pTermStruct->G[0], pTermStruct->G,
                       *pTermStruct->zeta, &dl);

  dl *= (dl > 0) ? 1 / 2. : 2.;

  // Solve for StateStar given a initial term structure
  if (LGMFindStateStar(NumPeriods, AsOfDate, plAcrlStartDates[1],
                       plAcrlStartDates[1], plPayDates, pszYieldCurveName, pdA,
                       pTermStruct->G[0], pTermStruct->G[0], pTermStruct->G,
                       *pTermStruct->zeta, pdCoverages, &ASFLevel,
                       pdASFZeroBond, dl, dh, statestartol, &StateStar,
                       dExerciseFee))
  // case 1. no solution
  {
    free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons, pdInitCoupons,
                     pdESOs);
    return "Numerical solver error in EuropeanAmortizingSwaption_()!";
  }

  // case 2: at least one zero solution
  else {

    averageGDiffSqr = 1.;
    counter = 0;

    while ((averageGDiffSqr > termstructtol) && (counter++ < maxcounter)) {

      // Compute new strikes using StateStar from above
      ComputeESOStrikes(NumPeriods, pdCoveragesUnAdjusted, pdASFZeroBond,
                        pdCoupons);

      // Set new strikes to calibration instruments  , and update vol and prices
      if (ErrMessage = UpdateCalibInstr(NumPeriods, pdCoupons, pszFreq,
                                        pszBasis, pszPayRec, pszRefRate,
                                        pszVolType, pszYieldCurveName, GetVol,
                                        boolIsLog, pCalibInstr, pdESOs)) {
        free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons,
                         pdInitCoupons, pdESOs);
        return ErrMessage;
      }

      // Calibrate to get new term structure
      if ((ErrMessage =
               FitAllGFwdSwap(pTermStruct, pCalibInstr, pTermStruct->zdate,
                              pTermStruct->zeta)) ||
          (ErrMessage = IsValidTermStructure(pTermStruct, NumPeriods))) {

        // if calibration fails
        // return max(intrinsic value  , 0.)

        *ASOPrice = AmortSwapPV(NumPeriods, AsOfDate, FixDate,
                                plPayDates, // vector(1  ,NumPeriods)
                                pszYieldCurveName,
                                pdA); // vector(0  ,NumPeriods)

        *ASOPrice = max((*ASOPrice), 0.);

        // De-normalize price
        *ASOPrice /= adjfac;

        for (ntime = 1; ntime <= NumPeriods; ++ntime) {
          // will come back...
          dvReplicatingStrikes[ntime - 1] = Coupon;
          dvReplicatingNotionals[ntime - 1] = 0.; // pdEta[ntime]*fabs(pdA[0]);
          dvReplicatingSwaptions[ntime - 1] = 0.; // pdESOs[ntime];
        }

        free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons,
                         pdInitCoupons, pdESOs);
        return 0;
      }

      // compute avrage squared difference between 2
      // term structures
      averageGDiffSqr = 0.;

      for (ntime = 0; ntime <= NumPeriods; ++ntime) {
        averageGDiffSqr += (ptmpG[ntime] - pTermStruct->G[ntime]) *
                           (ptmpG[ntime] - pTermStruct->G[ntime]);
      }

      averageGDiffSqr /= NumPeriods;

      // update ptmpG
      for (ntime = 0; ntime <= NumPeriods; ++ntime) {
        ptmpG[ntime] = pTermStruct->G[ntime];
      }

      // Solve for StateStar given a new term structure
      if (LGMFindStateStar(NumPeriods, AsOfDate, plAcrlStartDates[1],
                           plAcrlStartDates[1], plPayDates, pszYieldCurveName,
                           pdA, pTermStruct->G[0], pTermStruct->G[0],
                           pTermStruct->G, *pTermStruct->zeta, pdCoverages,
                           &ASFLevel, pdASFZeroBond, dl, dh, statestartol,
                           &StateStar, dExerciseFee)) {
        // changed by albert wang 05/12/03
        // for certain combination of amortizing structure and dates
        // numerical solver might fail !
        // in case this happens  , rather than throwing error  ,
        // get out of while loop and continue instead

        break;

        // previously
        // free_all_vectors(NumPeriods  ,pdA  ,ptmpG  ,pdEta  ,pdCoupons
        // ,pdInitCoupons  ,pdESOs); return "Numerical solver error in
        // EuropeanAmortizingSwaption_()!";
      }
    }

    // finally  , got out of while loop ...
    // update ESO prices one last time.
    if (ErrMessage = UpdateCalibInstr(NumPeriods, pdCoupons, pszFreq, pszBasis,
                                      pszPayRec, pszRefRate, pszVolType,
                                      pszYieldCurveName, GetVol, boolIsLog,
                                      pCalibInstr, pdESOs)) {
      free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons, pdInitCoupons,
                       pdESOs);
      return ErrMessage;
    }

    ComputeESOWeights(NumPeriods, pdEta, pdA, pdCoupons, pdCoveragesUnAdjusted);

#ifdef WIN32
#ifdef _DEBUG
    // check Etas' sum up to pdA[0]
    _RPT0(_CRT_WARN, "\n Test: Eta's sum up to pdA[0]");

    sumtmp = 0.;
    for (indtmp = 1; indtmp <= NumPeriods; ++indtmp) {
      sumtmp += pdEta[indtmp];
    }

    _RPT0(_CRT_WARN, "\n Sum of Eta's");
    _RPT1(_CRT_WARN, "%.10f", sumtmp);
    _RPT0(_CRT_WARN, "\n pdA[0]");
    _RPT1(_CRT_WARN, "%.10f", pdA[0]);

    for (indtmp = 1; indtmp <= NumPeriods; ++indtmp) {
      _RPT0(_CRT_WARN, "\n pdA");
      _RPT1(_CRT_WARN, "%.1f", indtmp);
      _RPT1(_CRT_WARN, "\t%.10f", pdA[indtmp]);

      sumtmp = 0.; // pdEta[indtmp];
      for (indtmp1 = indtmp; indtmp1 <= NumPeriods; ++indtmp1) {
        sumtmp += pdCoupons[indtmp1] * pdEta[indtmp1];
      }

      sumtmp *= pdCoveragesUnAdjusted[indtmp];
      sumtmp += pdEta[indtmp];
      _RPT0(_CRT_WARN, "\n reproduced ai");
      _RPT1(_CRT_WARN, "%.10f", sumtmp);
    }
#endif // _DEBUG
#endif // WIN32

    if (ErrMessage = AmortSwaption(NumPeriods, plAcrlStartDates[1],
                                   plAcrlEndDates, pszFreq, pszBasis, pszPayRec,
                                   pszRefRate, pszYieldCurveName, pszVolType,
                                   pdCoupons, // vector(1  ,NumPeriods)
                                   pdEta,     // vector(1  ,NumPeriods)
                                   boolIsLog, GetVol, ASOPrice)) {
      free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons, pdInitCoupons,
                       pdESOs);
      return ErrMessage;
    }
  }

#ifdef WIN32
#ifdef _DEBUG
  _RPT0(_CRT_WARN, "\n pdA's");

  for (ntime = 0; ntime <= NumPeriods; ++ntime) {
    _RPT1(_CRT_WARN, "\t%.10f", pdA[ntime]);
  }

  _RPT0(_CRT_WARN, "\n Gi's");

  for (ntime = 0; ntime <= NumPeriods; ++ntime) {
    _RPT1(_CRT_WARN, "\t%.10f", pTermStruct->G[ntime]);
  }

  _RPT0(_CRT_WARN, "\n Gi's");
  _RPT1(_CRT_WARN, "\t%.10f", *pTermStruct->zeta);
#endif // _DEBUG
#endif // WIN32

  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    dvReplicatingStrikes[ntime - 1] = pdCoupons[ntime];
    dvReplicatingNotionals[ntime - 1] = pdEta[ntime] * fabs(pdA[0]);
    dvReplicatingSwaptions[ntime - 1] = pdESOs[ntime];
  }

  // not necessary  , but check anyways
  *ASOPrice = max((*ASOPrice), 0.);

  // De-normalize price
  *ASOPrice /= adjfac;

  free_all_vectors(NumPeriods, pdA, ptmpG, pdEta, pdCoupons, pdInitCoupons,
                   pdESOs);
  return 0;
}

static char *IsValidTermStructure(LGM_TS *pTermStruct, long NumPeriods) {
  int tmpind;
  char *ErrMessage = 0;

  for (tmpind = 0; tmpind <= NumPeriods; ++tmpind) {
    if (_isnan(pTermStruct->G[tmpind]) || !_finite(pTermStruct->G[tmpind])) {
      ErrMessage = "Invalid term structure!";

      return ErrMessage;
    }
  }

  return 0;
}

void Init_MyData(EASData *_pMyData, long NumPeriods, long AsOfDate,
                 long FutureDate, long AcrlStartDate,
                 long *plPayDates, // vector(1  ,NumPeriods)
                 char *pszYieldCurveName,
                 double *pdA, // vector(0  ,NumPeriods)
                 double G_FutureDate, double G_AcrlStartDate,
                 double *pdG, // vector(1  ,NumPeriods)
                 double zeta, double *pdASFLevel, double *pdASFZeroBond,
                 double *pdCoverages, double dExerciseFee) {
  _pMyData->m_NumPeriods = NumPeriods;
  _pMyData->m_AsOfDate = AsOfDate;
  _pMyData->m_FutureDate = FutureDate;
  _pMyData->m_AcrlStartDate = AcrlStartDate;
  _pMyData->m_plPayDates = plPayDates;
  _pMyData->m_pszYieldCurveName = pszYieldCurveName;
  _pMyData->m_pdA = pdA;
  _pMyData->m_G_FutureDate = G_FutureDate;
  _pMyData->m_G_AcrlStartDate = G_AcrlStartDate;
  _pMyData->m_pdG = pdG;
  _pMyData->m_zeta = zeta;
  _pMyData->m_pdASFLevel = pdASFLevel;
  _pMyData->m_pdASFZeroBond = pdASFZeroBond;
  _pMyData->m_pdCoverages = pdCoverages;
  _pMyData->m_dExerciseFee = dExerciseFee;
}

double NAG_CALL _ObjectiveFunction(double dX, Nag_User *pComm) {
  static char *AmortSwapFV(double, EASData *, double *);
  double result;

  // unpack "other parameters"
  EASData *pData = (EASData *)(pComm->p);

  // pass evaluation point plus "other parameters"
  if (!AmortSwapFV(dX, pData, &result)) {
    return result;
  }

  // return garbage value that
  // should cause computation to fail
  return -99999999999.;
}

void DetermineStateBounds(long NumPeriods, long AsOfDate, long AcrlStartDate,
                          long *plPayDates, // vector(1  ,NumPeriods)
                          char *pszYieldCurveName, double G_AcrlStartDate,
                          double *pdG, // vector(1  ,NumPeriods)
                          double zeta,
                          // double *dh  ,
                          double *dl) {
  // func declaration
  void LGM1FZeroBond_DetermineBounds(long, long, long, char *, double, double,
                                     double, double *);

  long ntime;
  // double state;
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {

    LGM1FZeroBond_DetermineBounds(AsOfDate, AcrlStartDate, plPayDates[ntime],
                                  pszYieldCurveName, G_AcrlStartDate,
                                  pdG[ntime], zeta, dl);
  }
}

void LGM1FZeroBond_DetermineBounds(long AsOfDate, long StartDate, long EndDate,
                                   char *pszYieldCurveName, double G_StartDate,
                                   double G_EndDate, double zeta,
                                   double *state) {
  double price;
  // func declaration
  static char *LGM1FZeroBond(long, long, long, char *, double, double, double,
                             double, double *);

  if (StartDate > AsOfDate && StartDate < EndDate) {
    double tmp1 = G_StartDate - G_EndDate;
    double tmp2 = swp_f_df(StartDate, EndDate, pszYieldCurveName);
    price = tmp2 * exp(-tmp1 * (*state + 0.5 * tmp1 * zeta));

    if (price > 1.) {
      *state = (log(tmp2)) / tmp1 - 0.5 * tmp1 * zeta;

#ifdef _DEBUG
      LGM1FZeroBond(AsOfDate, StartDate, EndDate, pszYieldCurveName, *state,
                    G_StartDate, G_EndDate, zeta, &price);
#endif //_DEBUG
    }
  }
}

// Nag root-finding routine
static char *
LGMFindStateStar(long NumPeriods, long AsOfDate, long FutureDate,
                 long AcrlStartDate,
                 long *plPayDates, // vector(1  ,NumPeriods)
                 char *pszYieldCurveName,
                 double *pdA, // vector(0  ,NumPeriods)
                 double G_FutureDate, double G_AcrlStartDate,
                 double *pdG, // vector(1  ,NumPeriods)
                 double zeta, double *pdCoverages, double *pdASFLevel,
                 double *pdASFZeroBond, double dLowerBound, double dUpperBound,
                 double tol, double *StateStar, double dExerciseFee)

{
  // func declaration
  static char *AmortSwapFV(double, EASData *, double *);
  double NAG_CALL _ObjectiveFunction(double, Nag_User *);
  void Init_MyData(EASData *, long, long, long, long, long *, char *, double *,
                   double, double, double *, double, double *, double *,
                   double *, double);

  EASData data;
  Nag_User commParams;
  NagError fail;

#ifdef WIN32
#ifdef _DEBUG
  double asfvtmp = 0.;
  long tmpind;
#endif // _DEBUG
#endif // WIN32

  // initialize fail
  SET_FAIL(fail);
  fail.print = FALSE; // dont print!

  Init_MyData(&data, NumPeriods, AsOfDate, FutureDate, AcrlStartDate,
              plPayDates, pszYieldCurveName, pdA, G_FutureDate, G_AcrlStartDate,
              pdG, zeta, pdASFLevel, pdASFZeroBond, pdCoverages, dExerciseFee);

  commParams.p = (EASData *)(&data);
  nag_zero_cont_func_bd_1(dLowerBound, // lower bound on estimate of solution
                          dUpperBound, // uppper bound on estimate of solution
                          StateStar,   // locus of solution
                          _ObjectiveFunction, // the objective function
                          tol,                // relative accuracy
                          tol,
                          &commParams, // common block including parameters
                                       // passed to objective function
                          &fail        // error-handling parameters
  );

  if (fail.code != NE_NOERROR) {
    return fail.message;
  }

#ifdef WIN32
#ifdef _DEBUG

  // check StateStar does solve AmortSwapFV(state)=0.
  if (!AmortSwapFV(*StateStar, &data, &asfvtmp)) {
    _RPT0(_CRT_WARN, "\n Test: StateStar solve AmortSwapFV = 0");
    _RPT1(_CRT_WARN, "%.10f", asfvtmp);
    _RPT0(_CRT_WARN, "\n StateStar equals");
    _RPT1(_CRT_WARN, "%.10f", *StateStar);
    _RPT0(_CRT_WARN, "\n  Future level equals");
    _RPT1(_CRT_WARN, "%.10f", *(data.m_pdASFLevel));
    _RPT0(_CRT_WARN, "\n  Future bond values are:");
    for (tmpind = 1; tmpind <= NumPeriods; ++tmpind) {
      _RPT1(_CRT_WARN, "%.10f", data.m_pdASFZeroBond[tmpind]);
      // data.m_pdASFZeroBond[ntime];
    }

    _RPT0(_CRT_WARN, "\n ");
  }

#endif // _DEBUG
#endif // WIN32

  return 0;
}

/*void Init_a(long AsOfDate  ,
                        long NumPeriods  ,
                        long *plPayDates  ,
                        double*a  ,
                        double*pdA  ,
                        char *pszYieldCurveName)
{
        long ntime;
        for(ntime = 1; ntime<=NumPeriods; ++ntime)
        {
                a[ntime] =	pdA[ntime] *
                                        swp_f_df(AsOfDate  ,plPayDates[ntime]
,pszYieldCurveName)/ swp_f_df(AsOfDate  ,plPayDates[ntime]  ,pszYieldCurveName);
        }
}*/

// free memory
void free_all_vectors(long NumPeriods, double *pdA, double *ptmpG,
                      double *pdEta, double *pdCoupons, double *pdInitCoupons,
                      double *pdESOs) {
  free_vector(pdA, 0, NumPeriods);
  free_vector(ptmpG, 0, NumPeriods);

  free_vector(pdCoupons, 1, NumPeriods);
  free_vector(pdEta, 1, NumPeriods);
  free_vector(pdInitCoupons, 1, NumPeriods);
  free_vector(pdESOs, 1, NumPeriods);
}

static char *
AmortSwaption(long NumPeriods, long AcrlStart,
              long *plAcrlEndDates, // vector(1  ,NumPeriods)
              char *szFreq, char *szBasis, char *szPayRec, char *szRefRate,
              char *szYieldCurveName, char *szVolType,
              double *pdCoupons, // vector(1  ,NumPeriods)
              double *pdEta,     // vector(1  ,NumPeriods)
              SRT_Boolean boolIsLog,
              LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *),
              double *ASOPrice) {
  long ntime;
  double ESO = 0, vol = 0;
  char *ErrMessage = 0;

  *ASOPrice = 0.;
  for (ntime = 1; ntime <= NumPeriods; ++ntime) {

    if (ErrMessage = (*GetVol)(AcrlStart, plAcrlEndDates[ntime],
                               pdCoupons[ntime], boolIsLog, &vol)) {
      return ErrMessage;
    }

    if (ErrMessage =
            swp_f_Swaption(AcrlStart, plAcrlEndDates[ntime], szFreq, szBasis,
                           vol, pdCoupons[ntime], szPayRec, szRefRate,
                           szYieldCurveName, "PREMIUM", szVolType, &ESO)) {
      return ErrMessage;
    }

    *ASOPrice += ESO * pdEta[ntime];
  }

  return 0;
}

// returns price of underlying Amortizing Swap
static char *AmortSwapFV(double state, EASData *_pMyData,
                         double *pASFV) // vector(1  ,NumPeriods)
{
  // func declaration
  static char *LGM1FZeroBond(long, long, long, char *, double, double, double,
                             double, double *);

  long ntime;
  double zerobond = 0, ASFV = 0;
  char *ErrorMessage = 0;

  if (ErrorMessage = LGM1FZeroBond(
          _pMyData->m_AsOfDate, _pMyData->m_FutureDate,
          _pMyData->m_AcrlStartDate, _pMyData->m_pszYieldCurveName, state,
          _pMyData->m_G_FutureDate, _pMyData->m_G_AcrlStartDate,
          _pMyData->m_zeta, &zerobond)) {
    return ErrorMessage;
  }

  *pASFV = -(_pMyData->m_pdA[0] * zerobond + _pMyData->m_dExerciseFee);
  *_pMyData->m_pdASFLevel = 0;

  for (ntime = 1; ntime <= _pMyData->m_NumPeriods; ++ntime) {
    if (ErrorMessage = LGM1FZeroBond(
            _pMyData->m_AsOfDate, _pMyData->m_FutureDate,
            _pMyData->m_plPayDates[ntime], _pMyData->m_pszYieldCurveName, state,
            _pMyData->m_G_FutureDate, _pMyData->m_pdG[ntime], _pMyData->m_zeta,
            &zerobond)) {
      return ErrorMessage;
    }

    *pASFV += _pMyData->m_pdA[ntime] * zerobond;
    _pMyData->m_pdASFZeroBond[ntime] = zerobond;
    *_pMyData->m_pdASFLevel += _pMyData->m_pdCoverages[ntime] * zerobond;
  }

  return 0;
}

double AmortSwapPV(long NumPeriods, long AsOfDate, long AcrlStartDate,
                   long *plPayDates, // vector(1  ,NumPeriods)
                   char *pszYieldCurveName,
                   double *pdA) // vector(0  ,NumPeriods)
{
  long ntime;

  double ASPV = -pdA[0] * swp_f_df(AsOfDate, AcrlStartDate, pszYieldCurveName);

  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    ASPV +=
        pdA[ntime] * swp_f_df(AsOfDate, plPayDates[ntime], pszYieldCurveName);
  }

  return ASPV;
}

static char *LGM1FZeroBond(long AsOfDate, long StartDate, long EndDate,
                           char *pszYieldCurveName, double state,
                           double G_StartDate, double G_EndDate, double zeta,
                           double *price) {
  if (StartDate <= AsOfDate) {
    return "StartDate must not be smaller than AsOfDate!";
  }

  if (StartDate > EndDate) {
    return "EndDate must not be smaller than StartDate!";
  }

  if (StartDate == EndDate) {
    *price = 1.;
  }

  else //(StartDate<EndDate)
  {
    double tmp1 = G_StartDate - G_EndDate;
    double tmp2 = swp_f_df(StartDate, EndDate, pszYieldCurveName);
    *price = tmp2 * exp(-tmp1 * (state + 0.5 * tmp1 * zeta));
  }

  return 0;
}

void ComputeESOStrikes(long NumPeriods, double *pdCoverages,
                       double *pdASFZeroBond,
                       double *pdCoupons) // vector(1  ,NumPeriods)
{
  long ntime;
  double tmp = 0.;

  for (ntime = 1; ntime <= NumPeriods; ++ntime) {
    tmp += pdASFZeroBond[ntime] * pdCoverages[ntime];
    pdCoupons[ntime] = (1. - pdASFZeroBond[ntime]) / tmp;
  }
}

// pdEta MUST be  passed in with allocated memory
void ComputeESOWeights(
    long NumPeriods, double *pdEta,
    double *pdA,       // vector(0  ,NumPeriods)
    double *pdCoupons, // vector(1  ,NumPeriods)  , strikes of replicating ESOs
    double *pdCoverages) // vector(1  ,NumPeriods)

{
  long ntime;
  double sum;

  pdEta[NumPeriods] =
      pdA[NumPeriods] / (1. + pdCoverages[NumPeriods] * pdCoupons[NumPeriods]);
  sum = pdCoupons[NumPeriods] * pdEta[NumPeriods];

  for (ntime = NumPeriods - 1; ntime >= 1; --ntime) {
    pdEta[ntime] = (pdA[ntime] - pdCoverages[ntime] * sum) /
                   (1. + pdCoverages[ntime] * pdCoupons[ntime]);
    sum += pdCoupons[ntime] * pdEta[ntime];
  }
}

// it is assumed that Amorting Swaps covers period [T0  , T1  , ...  , Tn]  ,
// where T0 is first (Acrl) StartDate and Tn last (Acrl) End Date
// pdA MUST be passed in with allocated memory
void ComputeBondWeights(long NumPeriods,
                        double Coupon,       // Coupon of Amortizing Swap
                        double *pdA,         // vector(0  ,NumPeriods)
                        double *pdNotionals, // vector(1  ,NumPeriods)
                        double *pdCoverages, // vector(1  ,NumPeriods)
                        double *pdFixLegCoupons, double *pdFltNotionals) {
  long ntime;

  // 4 cases

  //// assumes that floating notional and fix notionals are the same
  //// shall come back later ...

  pdFltNotionals = 0; //// LEAK !!!! will come back later ....

  if (pdFltNotionals) {
    if (pdFixLegCoupons) {
      pdA[0] = pdFltNotionals[1];
      pdA[NumPeriods] = pdNotionals[NumPeriods] * pdCoverages[NumPeriods] *
                            pdFixLegCoupons[NumPeriods] +
                        pdFltNotionals[NumPeriods];

      for (ntime = 1; ntime < NumPeriods; ++ntime) {
        // pdA[ntime] = pdNotionals[ntime]*(1.+ pdCoverages[ntime]*Coupon) -
        // pdNotionals[ntime+1];
        pdA[ntime] =
            pdNotionals[ntime] * pdCoverages[ntime] * pdFixLegCoupons[ntime] +
            pdFltNotionals[ntime] - pdFltNotionals[ntime + 1];
      }
    } else {
      pdA[0] = pdFltNotionals[1];
      pdA[NumPeriods] =
          pdNotionals[NumPeriods] * pdCoverages[NumPeriods] * Coupon +
          pdFltNotionals[NumPeriods];

      for (ntime = 1; ntime < NumPeriods; ++ntime) {
        // pdA[ntime] = pdNotionals[ntime]*(1.+ pdCoverages[ntime]*Coupon) -
        // pdNotionals[ntime+1];
        pdA[ntime] = pdNotionals[ntime] * pdCoverages[ntime] * Coupon +
                     pdFltNotionals[ntime] - pdFltNotionals[ntime + 1];
      }
    }
  } else {
    if (pdFixLegCoupons) {
      pdA[0] = pdNotionals[1];
      pdA[NumPeriods] =
          pdNotionals[NumPeriods] *
          (1. + pdCoverages[NumPeriods] * pdFixLegCoupons[NumPeriods]);

      for (ntime = 1; ntime < NumPeriods; ++ntime) {
        pdA[ntime] = pdNotionals[ntime] *
                         (1. + pdCoverages[ntime] * pdFixLegCoupons[ntime]) -
                     pdNotionals[ntime + 1];
      }
    } else {
      pdA[0] = pdNotionals[1];
      pdA[NumPeriods] =
          pdNotionals[NumPeriods] * (1. + pdCoverages[NumPeriods] * Coupon);

      for (ntime = 1; ntime < NumPeriods; ++ntime) {
        pdA[ntime] = pdNotionals[ntime] * (1. + pdCoverages[ntime] * Coupon) -
                     pdNotionals[ntime + 1];
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------------------------
// // finds ystar  , the state variable at which a swaption is ATM.  Uses
// bisection because the existing Newton is unstable Note that the function is
// decreasing in y
static char *FindYstar(long n, double *a, double DStart, double zeta,
                       double *Gpay, double Gst, double *pdYstar) {
  double y, term, func, target, dErr, yLower, yUpper;
  long iter, itermax, i;

  itermax = 25;
  target = DStart;
  y = 0.;

  // Find a lower bound in y
  y = -0.5;
  dErr = -1.0;
  iter = 0;
  while (dErr < 0.0 && iter < itermax) {
    y *= 2.0;
    func = 0.0;
    for (i = 1; i <= n; i++) {
      term = a[i] * exp(-(Gst - Gpay[i]) * (y + 0.5 * (Gst - Gpay[i]) * zeta));
      func = func + term;
    }
    iter++;
    dErr = func - target;
  }
  if (iter == itermax)
    return "Could not find Ystar";
  yLower = y;

  // Find an upper bound
  y = 0.5;
  dErr = 1.0;
  iter = 0;
  while (dErr > 0.0 && iter < itermax) {
    y *= 2.0;
    func = 0.0;
    for (i = 1; i <= n; i++) {
      term = a[i] * exp(-(Gst - Gpay[i]) * (y + 0.5 * (Gst - Gpay[i]) * zeta));
      func = func + term;
    }
    iter++;
    dErr = func - target;
  }
  if (iter == itermax)
    return "Could not find Ystar";
  yUpper = y;

  // Now use bisection to find the solution
  while (yUpper - yLower > 1.0e-06) {
    y = 0.5 * (yUpper + yLower);
    func = 0.0;
    for (i = 1; i <= n; i++) {
      term = a[i] * exp(-(Gst - Gpay[i]) * (y + 0.5 * (Gst - Gpay[i]) * zeta));
      func = func + term;
    }
    dErr = func - target;
    if (dErr > 0.0)
      yLower = y;
    else
      yUpper = y;
  }
  *pdYstar = y;
  return 0;
}

// Local Functions
static const char *GetStrikes(int bIsConstNotional, double dCoupon,
                              long NumPeriod, double *c, double zeta,
                              double *coupon_G, double ex_G, double *tCvg,
                              double *tDFStartDate, double *strikevalue,
                              double *pdYstar) {
  double level, ystar;
  int i, j;
  char *err;

  /*Compute the valueof the state variable for which the amortising swap goes to
   * 0*/
  if (err = FindYstar(NumPeriod, c, 1, zeta, coupon_G, ex_G, pdYstar))
    return err;

  ystar = *pdYstar;

  /*Compute the strikes from 1 to NumPeriod*/
  if (bIsConstNotional) {
    for (j = 1; j < NumPeriod + 1;
         j++) /*we want to find the same y star for all vanilla swaptions by
                 changing the strike*/
      strikevalue[j] = dCoupon;
  } else {
    for (j = 1; j < NumPeriod + 1;
         j++) /*we want to find the same y star for all vanilla swaptions by
                 changing the strike*/
    {
      level = 0.0;
      for (i = 1; i <= j; i++) /*value of the vanilla swap at maturity date*/
      {
        level += tCvg[i - 1] * tDFStartDate[i - 1] *
                 exp(-(ex_G - coupon_G[i]) *
                     (ystar + 0.5 * (ex_G - coupon_G[i]) * zeta));
      }
      strikevalue[j] =
          (1 - tDFStartDate[j - 1] *
                   exp(-(ex_G - coupon_G[j]) *
                       (ystar + 0.5 * (ex_G - coupon_G[j]) * zeta))) /
          level;
    }
  }
  return 0;
}

static void freeMemory(long *x1, long n1) {
  if (x1)
    free_lngvector(x1, 0, n1);
}

// static char* PriceConstantNotional( long lNumPeriod  , double dCoupon  ,
// double* out_pdPrice

// We need to alter the pricer such that it takes in schedules.

Err EuropeanAmortizingSwaption(
    long today, long StartDate, long TheoEndDate, long NumPeriod,
    long *lvFixedStartDates, long *lvFixedEndDates, long NoticePeriod,
    double Coupon,
    double *dvNotionals, // 0  ,..  ,NumPeriod  ( dvNotionals[NumPeriod]=0 )
    char *szYieldCurveName, char *szFreq, char *szBasis, char *szPayRec,
    char *szRefRate,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *),
    SrtDiffusionType srt_vol_type,
    // Output
    double *dvPayDates, double *dvReplicatingStrikes,
    double *dvReplicatingNotionals, double *dPrice, double *dFixedPV,
    double *dFloatPV, double *dSwapRate, SigKapTS **lgmSigKapTSPtrPtr,
    double *dLGMVol,
    //								double* dvTauDates
    //, 								double* dvTau  , 								long* nTau  ,
    double *dvReplicatingSwaptions, double dMargin) {
  // Local variables
  double DStart; // discount factor at start date
  long k;
  double floatinglegPV, fixedlegPV;
  int i, j, numG;
  double ex_G;
  double beta, zeta, ystar, SwapRate, amorprice, price, vol;
  long zdatefix;
  double zetafix, Maturity, dLambda, dTex, dT1, dT0, dSpread, dFltCvg;
  int bIsConstNotional = 1;
  long lMonthsPerPeriod;
  int iNotFinished;

  // Sort structures
  //	SrtBusDayConv Conv;
  SrtBasisCode FloatBasis;

  SrtCompounding FloatFrequency;
  LGMErr (*GetBeta)(Date, Date, double *);
  LGMErr error;
  SRT_Boolean boolIsLog;
  SrtBasisCode basis;
  SrtReceiverType PayRec;
  SrtCompounding frequency;

  LGM_TS *LGM_TSPtr;
  LGM_TS LGM_TSobj;
  LGMCalSet *CSPtr;
  LGMCalSet CS;
  LGM_TS *tsPtr;

  // vectors
  //	long *tStart = lngvector(0  ,NumPeriod);
  //	long *tPremiumDates = lngvector(0  ,NumPeriod);
  //	long *tEnd = lngvector(0  ,NumPeriod);
  //	long *tPay = lngvector(0  ,NumPeriod);
  //	long *tStartDate = lngvector(0  ,NumPeriod);
  //	long *tEndDate = lngvector(0  ,NumPeriod);
  long *lvTheoEndDates = lngvector(0, NumPeriod);

  long *tEx = lngvector(0, NumPeriod);
  double *tCvg = dvector(0, NumPeriod);
  double *tExPremium = dvector(0, NumPeriod);
  double *tStrike = dvector(0, NumPeriod);
  double *tRedFirstPay = dvector(0, NumPeriod);
  double *tDFnow = dvector(0, NumPeriod);
  double *tDFStartDate = dvector(0, NumPeriod);
  double *coupon_G = dvector(0, NumPeriod + 2);
  double *a = dvector(0, NumPeriod + 1);
  double *c = dvector(0, NumPeriod + 1);
  double *b = dvector(0, NumPeriod + 1);
  double *strikevalue = dvector(0, NumPeriod + 1);
  double *A = dvector(0, NumPeriod + 2);
  char *szVolType;

  // Get the refrate details
  swp_f_get_ref_rate_details(szRefRate, &FloatBasis, &FloatFrequency);

  // Set up vol structure
  if (srt_vol_type == SRT_LOGNORMAL) {
    GetBeta = LGMReturnExp1; /* This will return beta=1 for all swaptions */
    beta = 1.0;
    szVolType = "LOGNORMAL";
    boolIsLog = SRT_TRUE;
  } else {
    GetBeta = LGMReturnExp0; /* This will return beta=0 for all swaptions */
    beta = 0.0;
    szVolType = "NORMAL";
    boolIsLog = SRT_FALSE;
  }

  // Check to see if it is a constant notional swap  , and if so  , return the
  // standard swaption value
  for (i = 0; i < NumPeriod; i++)
    bIsConstNotional = bIsConstNotional && (dvNotionals[0] == dvNotionals[i]);

  if (bIsConstNotional) {
    // set the output strike outputs
    for (i = 0; i < NumPeriod; i++) {
      dvPayDates[i] = lvFixedEndDates[i];
      dvReplicatingNotionals[i] = 0.0;
      dvReplicatingStrikes[i] = Coupon;
    }
    dvReplicatingNotionals[NumPeriod - 1] = dvNotionals[0];

    // calculate the swaption value
    (*GetVol)(StartDate, lvFixedEndDates[NumPeriod - 1], Coupon, boolIsLog,
              &vol);
    swp_f_Swaption(StartDate, TheoEndDate, szFreq, szBasis, vol, Coupon,
                   szPayRec, szRefRate, szYieldCurveName, "PREMIUM", szVolType,
                   dPrice);
    *dPrice *= dvNotionals[0];

    // return
    return 0;
  }

  // Compute Payment Dates  , coverages and df
  //	tStartDate[0] = StartDate;

  // get the Swap details
  interp_compounding(szFreq, &frequency);
  interp_basis(szBasis, &basis);
  interp_rec_pay(szPayRec, &PayRec);
  lMonthsPerPeriod = 12 / frequency;

  Maturity = (TheoEndDate - StartDate) / 365.0;
  NumPeriod = (long)(Maturity * frequency + 0.1);

  for (i = 0; i < NumPeriod; i++) {
    //		/*K =*/ tPay[i] = add_unit(TheoEndDate  , - (int) (( NumPeriod - i-1 )
    //* 12 / (frequency+0.0))   , SRT_MONTH  , NO_BUSDAY_CONVENTION); 		if
    //(i<NumPeriod-1) tStartDate[i+1] =  add_unit(tPay[i]  , 0   , SRT_DAY  ,
    //Conv); 		tEndDate[i] = add_unit(tPay[i]  , 0   , SRT_DAY  , Conv); 		tPay[i] =
    //add_unit(tPay[i]  , 0   , SRT_BDAY  , MODIFIED_SUCCEEDING);

    // Calculate the theoretical end dates for each of the underlying swaps
    lvTheoEndDates[i] = add_unit(StartDate, (i + 1) * lMonthsPerPeriod,
                                 SRT_MONTH, NO_BUSDAY_CONVENTION);

    // Calculate the spread.  For the moment  , we assume that the fixed and
    // floating frequencies are the same.  Use it to adjust the fixed coverage

    dSpread =
        swp_f_spread(lvFixedStartDates[i], lvFixedEndDates[i], szRefRate) +
        dMargin;
    dFltCvg = coverage(lvFixedStartDates[i], lvFixedEndDates[i], FloatBasis);

    tCvg[i] = coverage(lvFixedStartDates[i], lvFixedEndDates[i], basis) -
              dSpread * dFltCvg / Coupon;
    tDFnow[i] = swp_f_df(today, lvFixedEndDates[i], szYieldCurveName);
    tDFStartDate[i] = swp_f_df(StartDate, lvFixedEndDates[i], szYieldCurveName);
  }
  DStart = swp_f_df(today, StartDate, szYieldCurveName);

  for (i = 0; i < NumPeriod; i++) {
    //		tStart[i] = tStartDate[i];
    //		tEnd[i] = tEndDate[i];
    //		tPremiumDates[i] = add_unit(tStart[i]  , 0   , SRT_BDAY  ,
    //MODIFIED_SUCCEEDING);
    tEx[i] = add_unit(lvFixedStartDates[i], -NoticePeriod, SRT_BDAY,
                      MODIFIED_SUCCEEDING);
    tExPremium[i] = 0.0;
    tStrike[i] = 1.0;
    tExPremium[i] = 0.0;
    tRedFirstPay[i] = 0.0;
    strikevalue[i] = Coupon;
  }

  /*Get a calibrated term structure for the amortising deal*/
  //	lgmSigKapTSptr->TSArr = 0;
  //	lgmSigKapTSptr->NTS = 0;

  /*calibrated term structure*/
  LGM_TSobj.zdate = 0;
  LGM_TSobj.zeta = 0;
  LGM_TSobj.Gdate = 0;
  LGM_TSobj.G = 0;
  LGM_TSPtr = 0;
  LGM_TSPtr = &LGM_TSobj;

  // Set up the calibration pointer
  CS.tNow = today;
  CS.nEx = NumPeriod;
  CS.n = NumPeriod;
  CS.nPay = NumPeriod;
  CS.ffirst = 1;
  CS.nfirst = 1;
  CS.nlast = NumPeriod;
  CS.ftEx = tEx[0];
  CS.ftStart = StartDate;
  CS.cvgfix = tCvg[0];
  CS.DfixStart = DStart;

  CS.tEx = lvector(0, NumPeriod + 1);
  CS.tStart = lvector(0, NumPeriod + 1);
  CS.DStart = dvector(0, NumPeriod + 1);
  CS.tPay = lvector(0, NumPeriod + 1);
  CS.cvgpay = dvector(0, NumPeriod + 1);
  CS.Dpay = dvector(0, NumPeriod + 1);
  CS.ifirst = lvector(0, NumPeriod + 1);
  CS.nlong = lvector(0, NumPeriod + 1);
  CS.cvgfirst = dvector(0, NumPeriod + 1);
  CS.Rfix = dvector(0, NumPeriod + 1);
  CS.CEVLong = dvector(0, NumPeriod + 1);
  CS.Vfix = dvector(0, NumPeriod + 1);

  CSPtr = 0;
  CSPtr = &CS;

  for (i = 1; i <= NumPeriod; i++) {
    CSPtr->tEx[i] = tEx[0];
    CSPtr->tStart[i] = StartDate;
    CSPtr->DStart[i] = DStart;
    //		CSPtr->tPay[i] = tPay[i-1];
    CSPtr->tPay[i] = lvFixedEndDates[i - 1];
    CSPtr->cvgpay[i] = tCvg[i - 1];
    CSPtr->Dpay[i] = tDFnow[i - 1];
    CSPtr->ifirst[i] = 1;
    CSPtr->nlong[i] = i;
    CSPtr->cvgfirst[i] = tCvg[0];
  }

  // Calibrate
  /*compute the amortising swap at maturity date*/

  k = frequency / FloatFrequency; /* to take into account the case where
                                     frequency different from float frequency*/
  if (k < 1)
    k = 1; /* if float frequency > frequency  , nothing changes */

  floatinglegPV = dvNotionals[0] * DStart -
                  dvNotionals[NumPeriod - 1] * tDFnow[NumPeriod - 1];
  fixedlegPV =
      dvNotionals[NumPeriod - 1] * tCvg[NumPeriod - 1] * tDFnow[NumPeriod - 1];

  for (j = 1; j < NumPeriod; j++) {
    if (FloatFrequency >= frequency) {
      b[j] = dvNotionals[j - 1] * (1 + tCvg[j - 1] * Coupon) -
             dvNotionals[j]; /*cash-flows*/
      c[j] =
          (dvNotionals[j - 1] * (1 + Coupon * tCvg[j - 1]) - dvNotionals[j]) /
          dvNotionals[0] *
          tDFStartDate[j - 1]; /*discounted cash-flows from start date*/
      floatinglegPV += (dvNotionals[j] - dvNotionals[j - 1]) * tDFnow[j - 1];
      fixedlegPV += dvNotionals[j - 1] * tCvg[j - 1] * tDFnow[j - 1];
    } else {

      if (fmod(j, k) == 0) {
        b[j] = dvNotionals[j - 1] * (1 + tCvg[j - 1] * Coupon) -
               dvNotionals[j - 1 + k]; /*cash-flows*/
        c[j] = (dvNotionals[j - 1] * (1 + Coupon * tCvg[j - 1]) -
                dvNotionals[j - 1 + k]) /
               dvNotionals[k - 1] *
               tDFStartDate[j - 1]; /*discounted cash-flows from start date*/
        floatinglegPV += (dvNotionals[j - 1 + k] - dvNotionals[j - 1]) /
                         dvNotionals[k - 1] * tDFnow[j - 1];
        fixedlegPV += dvNotionals[j - 1] * tCvg[j - 1] * tDFnow[j - 1];
      } else {
        b[j] = dvNotionals[j - 1] * tCvg[j - 1] * Coupon; /*cash-flows*/
        c[j] = dvNotionals[j - 1] * Coupon * tCvg[j - 1] / dvNotionals[k - 1] *
               tDFStartDate[j - 1]; /*discounted cash-flows from start date*/
        floatinglegPV += 0;
        fixedlegPV += dvNotionals[j - 1] * tCvg[j - 1] * tDFnow[j - 1];
      }
    }
  }

  b[NumPeriod] = dvNotionals[NumPeriod - 1] *
                 (1 + tCvg[NumPeriod - 1] * Coupon); /*cash-flows*/
  c[NumPeriod] =
      dvNotionals[NumPeriod - 1] * (1 + Coupon * tCvg[NumPeriod - 1]) /
      dvNotionals[k - 1] *
      tDFStartDate[NumPeriod - 1]; /*discounted cash-flows from start date*/

  /* Create term structure */
  tsPtr = LGMCreateLGM_TS(2, NumPeriod + 1);

  /* if the expiry < 1M and the strike is more than 5 std dev away  , simply
   * return the intrinsic value */
  SwapRate = floatinglegPV / fixedlegPV;
  *dSwapRate = SwapRate;
  //	GetVol( StartDate  , TheoEndDate  , coupon  , SRT_FALSE  , &dVol );
  *dFixedPV = Coupon * fixedlegPV;
  *dFloatPV = floatinglegPV;

  // Calibrate

  // Find the G(t) values by calibrating on the 1 into k swaptions
  //	and fill numG  , the G dates  , and the G values
  iNotFinished = 1;
  for (j = 1; j < 10 && iNotFinished; j++) {
    for (i = 1; i <= NumPeriod; i++) {
      CSPtr->Rfix[i] = strikevalue[i - 1];
      GetVol(CSPtr->tStart[i], CSPtr->tPay[i], CSPtr->Rfix[i], SRT_FALSE,
             &CSPtr->CEVLong[i]);
      CSPtr->Vfix[i] = LGMCEVSwaptionPrice(
          today, CSPtr->tEx[i], SRT_RECEIVER, CSPtr->Rfix[i], DStart, 0.0, 0.0,
          i, &CSPtr->cvgpay[1], &CSPtr->Dpay[1], CSPtr->CEVLong[i], beta);
    }

    if (error = FitAllGFwdSwap(tsPtr, CSPtr, &zdatefix, &zetafix)) {
      // Check to see if the error is caused by short maturity  , and if so  ,
      // return the intrinsic
      if ((tEx[0] - today) / 365.0 <
          .083) /*&& (fabs(coupon-SwapRate) > 5*dVol ) */
      {
        *dPrice =
            PayRec == SRT_PAYER ? *dFloatPV - *dFixedPV : *dFixedPV - *dFloatPV;
        if (*dPrice < 0.0)
          *dPrice = 0.0;
        return 0;
      }
      // if we previously had a reasonable solution  , keep that
      if (j > 1) {
        iNotFinished = 0;
      } else
        return error;
    }

    // Check that the term structure is reasonable.  If not use the previous
    // value (if no previous value throw an error)
    else if (zetafix > .01 * (tEx[0] - today) / 365.0) {
      iNotFinished = 0;
    } else {
      tsPtr->zdate[0] = today;
      tsPtr->zeta[0] = 0;
      tsPtr->zdate[1] = zdatefix;
      tsPtr->zeta[1] = zetafix;
      LGM_TSPtr = tsPtr;
      numG = LGM_TSPtr->numG;
      ex_G = LGM_TSPtr->G[0];

      for (i = 1; i < numG; i++)
        coupon_G[i] = LGM_TSPtr->G[i];
      zeta = LGM_TSPtr->zeta[1];
      // get the strikes  , ystar value.  If it is -infinity  , then use the
      // previous strike values
      if (GetStrikes(bIsConstNotional, Coupon, NumPeriod, c, zeta, coupon_G,
                     ex_G, tCvg, tDFStartDate, strikevalue, &ystar))
        iNotFinished = 0;
    }
  }
  /*Compute now the notionals  */
  A[NumPeriod + 1] = 0; /*Convention*/
  b[NumPeriod + 1] = 0; /*Comvention*/
  tCvg[NumPeriod] = 1;  /*Convention*/

  for (j = NumPeriod; j > 0; j--) {
    A[j] = 1 / (1 + tCvg[j - 1] * strikevalue[j]) * tCvg[j - 1] *
           (A[j + 1] / tCvg[j] + b[j] / tCvg[j - 1] - b[j + 1] / tCvg[j]);
  }

  amorprice = 0;
  price = 0;

  for (j = 1; j < NumPeriod + 1;
       j++) /* once we get the strikes and notionals  , we compute the vanilla
               swaptions prices*/
  {

    (*GetVol)(StartDate, lvFixedEndDates[j - 1], strikevalue[j], boolIsLog,
              &vol);

    swp_f_Swaption(StartDate, lvTheoEndDates[j - 1], szFreq, szBasis, vol,
                   strikevalue[j], szPayRec, szRefRate, szYieldCurveName,
                   "PREMIUM", szVolType, &price);
    //	dvDiscPay[j] += tDFnow[j-1];
    //	dvDiscPay[j] = tCvg[j-1]*strikevalue[j]*tDFnow[j-1] + tDFnow[j-1];
    //	vecLGMPrice[j] = LGMRecVal(&vecYstar[j]  , j  , dvDiscPay  , DStart  ,
    ///* receiver swaption */ 			sqrt(zeta)  , coupon_G  , ex_G); 	dvDiscPay[j] -=
    //tDFnow[j-1];
    dvReplicatingSwaptions[j - 1] = price;

    amorprice += price * A[j];
  }
  // Check that the option price is positive
  if (amorprice < 0.0)
    amorprice = 0.0;

  // Get the output
  for (i = 0; i < NumPeriod; i++) {
    dvPayDates[i] = lvFixedEndDates[i];
    dvReplicatingNotionals[i] = A[i + 1];
    dvReplicatingStrikes[i] = strikevalue[i + 1];
  }

  *dPrice = amorprice;

  *lgmSigKapTSPtrPtr = LGMConvertZGtoSigKap(today, tsPtr);

  // get out the sigma tau term structure
  dLambda = (*lgmSigKapTSPtrPtr)->kap[0];
  dTex = (tEx[0] - today) / 365.0;
  dT1 = (lvFixedEndDates[0] - today) / 365.0;
  dT0 = (StartDate - today) / 365.0;
  *dLGMVol = (ex_G - coupon_G[1]) * sqrt(2.0 * zeta) * dLambda *
             sqrt(dLambda / (exp(2.0 * dLambda * dTex) - 1.0)) /
             (exp(-dLambda * dT0) - exp(-dLambda * dT1));

  // free the memory
  freeMemory(lvTheoEndDates, NumPeriod);

  return NULL;
}

char *ZeroCouponSwaption_ComputeCoupon(long Today, char *pcYieldCurveName,
                                       double dBondMult, long lNumPeriods,
                                       long *plFixedDates, long lNumFltPeriods,
                                       long *plFltDates, double *pdMargins,
                                       double *pdSpread, double *pdFixCov,
                                       double *pdFltCov, double *pdCoupon)
// long	 *pdFixCov_adjusted)
{

  // objective func  declaration
  double NAG_CALL _ObjectiveFunction_ZC_BondMult_(double, Nag_User *);

  // variable declaration
  EASData_ZC_BondMult data;
  Nag_User commParams;
  NagError fail;

  double dLowerBound = 1.e-13;
  double dUpperBound = 1.;
  double dtol = 1.e-15;

  // initialize data
  data.m_Today = Today;
  data.m_pYieldCurve = pcYieldCurveName;

  data.m_plFixedDates = plFixedDates;
  data.m_lNumPeriods = lNumPeriods;

  data.m_lNumFltPeriods = lNumFltPeriods;
  data.m_plFltDates = plFltDates;

  data.m_BondMult = dBondMult;
  data.m_pdMargins = pdMargins;
  data.m_pdSpread = pdSpread;
  data.m_pdFixCov = pdFixCov;
  data.m_pdFltCov = pdFltCov;

  data.m_pdFixCov_adjusted = dvector(0, lNumPeriods - 1);

  // initialize fail
  SET_FAIL(fail);
  fail.print = FALSE; // dont print!

  commParams.p = (EASData_ZC_BondMult *)(&data);
  nag_zero_cont_func_bd_1(
      dLowerBound,                     // lower bound on estimate of solution
      dUpperBound,                     // uppper bound on estimate of solution
      pdCoupon,                        // locus of solution
      _ObjectiveFunction_ZC_BondMult_, // the objective function
      dtol,                            // relative accuracy
      dtol,
      &commParams, // common block including parameters passed to objective
                   // function
      &fail);      // error-handling parameters

  if (fail.code != NE_NOERROR) {
    free_dvector(data.m_pdFixCov_adjusted, 0, lNumPeriods - 1);
    return fail.message;
  }

  free_dvector(data.m_pdFixCov_adjusted, 0, lNumPeriods - 1);
  return 0;
}

/*double NAG_CALL _ObjectiveFunction_ZC_Price_(double dX  ,Nag_User *pComm)
{
        int tmpind;
        double result = 1.;
        char *ErrMessage = 0;

        // unpack "other parameters"
        EASData_ZC_Price * pData = (EASData_ZC_Price * )(pComm->p);


        // call xl function to value

        SExcelRange xlLGMZCPrice = XLZCEuropeanAmortSwaption(
pData->m_und_name  , pData->m_mkt_name  , pData->m_notperiod  ,
                                                                                                                        pData->m_ref_name  ,
                                                                                                                        pData->m_freq_name  ,
                                                                                                                        pData->m_basis_name  ,
                                                                                                                        pData->m_start_date  ,
                                                                                                                        pData->m_end_date  ,
                                                                                                                        pData->m_notional  ,
                                                                                                                        pData->m_num_margins  ,
                                                                                                                        pData->m_should_be_1_margins  ,
                                                                                                                        pData->m_margins  ,
                                                                                                                        pData->m_lpxCoupon  ,
                                                                                                                        pData->m_lpxBondMult);

        return xlLGMZCPrice.GetDouble(0  ,0  ,0);


}*/

double NAG_CALL _ObjectiveFunction_ZC_BondMult_(double dX, Nag_User *pComm) {
  int tmpind;
  double result = 1.;

  // this is just to appease function
  // double tmpadjust;

  char *ErrMessage = 0;

  // unpack "other parameters"
  EASData_ZC_BondMult *pData = (EASData_ZC_BondMult *)(pComm->p);

  /*	// call ComputeCoverages

          ComputeCoverages( pData->m_Today  ,
                                            pData->m_pYieldCurve  ,
                                            pData->m_pdFixCov  , // 0  ,
     NumFixPeriods -1 pData->m_plFixedDates  , //0  ,NumFixPeriods
                                            pData->m_lNumPeriods  ,
                                            pData->m_pdFltCov   ,// 0  ,
     NumFltPeriods -1 pData->m_pdSpread  ,// 0  , NumFltPeriods -1
                                            pData->m_pdMargins  ,
                                            pData->m_plFltDates  ,
                                            pData->m_lNumFltPeriods  ,
                                            pData->m_pdFixCov_adjusted  ,
                                            &dX  ,
                                            &tmpadjust);*/

  for (tmpind = 0; tmpind < pData->m_lNumPeriods; ++tmpind) {
    result *= (1. + pData->m_pdFixCov[tmpind] * dX);
  }

  return (result - 1. - pData->m_BondMult);
}

// preliminary version
// will come back later ...
void ComputeCoverages(long Today, char *pcYCName,

                      // fixed leg info
                      double *pdFixedCvgs, // 0  , NumFixPeriods -1
                      long *plFixedDates,  // 0  ,NumFixPeriods
                      long NumFixPeriods,

                      // flt leg info
                      double *pdFltCvgs,    // 0  , NumFltPeriods -1
                      double *pdFltSpreads, // 0  , NumFltPeriods -1
                      double *pdMargins,    // 0  , NumFltPeriods -1
                      // double *pdFltNotionals  ,// 0  , NumFltPeriods -1
                      long *plFltDates, // 0  , NumFltPeriods
                      long NumFltPeriods,
                      // double *pdFltNotionals  ,

                      // output
                      double *pdFixedCvgs_adjusted,

                      double *pdCoupon, double *pdFltAdj) {

  int tmpinc, indfix, indflt;
  double tmp;
  double fixedleglevel = 0, fltleglevel = 0;

  *pdFltAdj = 0.;

  tmpinc = NumFltPeriods / NumFixPeriods;

  for (indfix = 0; indfix < NumFixPeriods; ++indfix) {
    // tmpsum=0.;

    for (indflt = indfix * tmpinc; indflt < (indfix + 1) * tmpinc; ++indflt) {
      // looks though plFltDates might be out of bound  , but this is not true
      tmp = (/*pdFltSpreads[indflt] + */ pdMargins[indflt]) *
            pdFltCvgs[indflt] *
            swp_f_df(plFltDates[indflt], plFltDates[indflt + 1], pcYCName);
      fltleglevel += tmp;
      //*pdFltAdj += tmp*pdFltNotionals[indflt];
    }

    if (pdFixedCvgs_adjusted) {
      pdFixedCvgs_adjusted[indfix] = pdFixedCvgs[indfix];
    }

    // pdFixedCvgs_adjusted[indfix] = pdFixedCvgs[indfix];// - tmpsum/(Coupon *
    // swp_f_df(plFixedDates[indfix]  ,plFixedDates[indfix+1]  , pcYCName));
    fixedleglevel +=
        swp_f_df(plFixedDates[indfix], plFixedDates[indfix + 1], pcYCName) *
        pdFixedCvgs[indfix];
  }

  (*pdCoupon) -= (fltleglevel / fixedleglevel);
}

// preliminary version
// will come back later ...
void ComputeCoverages2(long Today, char *pcYCName,

                       // fixed leg info
                       double *pdFixedCvgs, // 0  , NumFixPeriods -1
                       long *plFixedDates,  // 0  ,NumFixPeriods
                       long NumFixPeriods, double *pdFixedNotionals,

                       // flt leg info
                       double *pdFltCvgs,    // 0  , NumFltPeriods -1
                       double *pdFltSpreads, // 0  , NumFltPeriods -1
                       double *pdMargins,    // 0  , NumFltPeriods -1
                       // double *pdFltNotionals  ,// 0  , NumFltPeriods -1
                       long *plFltDates, // 0  , NumFltPeriods
                       long NumFltPeriods,
                       // double *pdFltNotionals  ,

                       // output
                       double *pdFixedCvgs_adjusted,

                       double *pdCoupon, double *pdFltAdj) {

  int tmpinc, indfix, indflt;
  double tmp;
  double fixedleglevel = 0, fltleglevel = 0;

  *pdFltAdj = 0.;

  tmpinc = NumFltPeriods / NumFixPeriods;

  for (indfix = 0; indfix < NumFixPeriods; ++indfix) {
    // tmpsum=0.;

    for (indflt = indfix * tmpinc; indflt < (indfix + 1) * tmpinc; ++indflt) {
      // looks though plFltDates might be out of bound  , but this is not true
      tmp = (/*pdFltSpreads[indflt] + */ pdMargins[indflt]) *
            pdFltCvgs[indflt] * pdFixedNotionals[indfix] *
            swp_f_df(plFltDates[indflt], plFltDates[indflt + 1], pcYCName);
      fltleglevel += tmp;
      //*pdFltAdj += tmp*pdFltNotionals[indflt];
    }

    if (pdFixedCvgs_adjusted) {
      pdFixedCvgs_adjusted[indfix] = pdFixedCvgs[indfix];
    }

    // pdFixedCvgs_adjusted[indfix] = pdFixedCvgs[indfix];// - tmpsum/(Coupon *
    // swp_f_df(plFixedDates[indfix]  ,plFixedDates[indfix+1]  , pcYCName));
    fixedleglevel +=
        pdFixedNotionals[indfix] *
        swp_f_df(plFixedDates[indfix], plFixedDates[indfix + 1], pcYCName) *
        pdFixedCvgs[indfix];
  }

  (*pdCoupon) -= (fltleglevel / fixedleglevel);
}

/*void EuropeanAmortSwaption_PreProcess(char *mkt_name  ,
                                   char *freq_name  ,
                                   char *basis_name  ,
                                   char *ref_name  ,
                                   char *)
{

                // get freqFixed
                SrtCompounding freqFloating = SPP::ToFrequencyCode(freq_name);

                // get basisFixed
                SrtBasisCode basisFixed = SPP::ToBasisCode(basis_name );

                // get freqFloating and basisFloating
                SrtBasisCode basisFloating;
                SrtCompounding freqFloating;
                swp_f_get_ref_rate_details(ref_name  , &basisFloating  ,
&freqFloating);


                // 12/ ???
                if(freqFloating > freqFloating)
                {
                        throw "AdjustMargins(...): Floating frequency must not
be smaller than fixed frequency!";
                }

                // get ptrSchedule for fixed and floating legs
                SPP::TPtr_Schedule ptrFixedSchedule =
SPP::XSchedule::Create(static_cast<long>(start_date)
,static_cast<long>(end_date)  ,freqFixed  , basisFixed  , ptrCurrency  , false
); SPP::TPtr_Schedule ptrFloatingSchedule =
SPP::XSchedule::Create(static_cast<long>(start_date)
,static_cast<long>(end_date)  ,freqFloating  , basisFloating  , ptrCurrency  ,
false );


                // get NumPeriod
                unsigned short int num_fix_notionals =
ptrFixedSchedule->m_vecDates_Payment.size(); unsigned short int
num_float_notionals = ptrFloatingSchedule->m_vecDates_Payment.size();


                // check to see if number of margins equal number of flt leg
payment periods if(*num_margins != num_float_notionals)
                {
                        throw XSort("Number of margins must equal number of
floating leg payment periods!");
                }



                // assign notional to both legs
                double *fix_notionals = new double[num_fix_notionals];
                double *float_notionals = new double[num_float_notionals];

                // fixed leg dates
                long	*lvFixedStartDates =
ptrFixedSchedule->m_vecDates_Calculation.data(); long	*lvFixedEndDates =
lvFixedStartDates + 1;

                // floating leg dates
                long	*lvFloatingStartDates =
ptrFloatingSchedule->m_vecDates_Calculation.data();
                long	*lvFloatingEndDates = lvFloatingStartDates+1;



        int tmpind;

        // compute speads and coverages
        for (int tmpind = 0; tmpind<num_fix_notionals; ++tmpind)
        {
        pdSpread[tmpind] = swp_f_spread(lvFixedStartDates[tmpind]  ,
lvFixedEndDates[tmpind]  , ref_name); pdFltCvg[tmpind] =
coverage(lvFixedStartDates[tmpind]  , lvFixedEndDates[tmpind]  , basisFloating);
        pdFixCvg[tmpind] = coverage(lvFixedStartDates[tmpind]  ,
lvFixedEndDates[tmpind]  ,basisFixed);
        }
}*/
